#
# Parallel-tempering check for Experiment II shock detection.
#
# Plain MCMC (mcmc_check.jl) mixed incompletely: chains started from the EM
# fit and from the truth did not always agree, and the flip move was rarely
# accepted. This script runs parallel_tempering_time_varying with nine
# likelihood temperatures, geometric from 1 to 0.2, again from the EM fit and
# from the truth, and reports the posterior (beta = 1) chain alongside the EM
# fit. It also reports the shock AUROC of some tempered chains.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/tempering_check.jl \
#       p=10 k=2 seeds=1:10 tau_active=1.5 minimum_active=6 \
#       out=results/experiment2/tempering_p10_k2.csv
#
# MCMC_BURNIN and MCMC_SAMPLES set the number of sweeps (default 1000, 2000).
#

include(joinpath(@__DIR__, "mcmc_check.jl"))

const BETAS = collect(0.2 .^ ((0:8) ./ 8))
# Tempered chains reported alongside the posterior chain: beta = 0.55, 0.30, 0.20.
const REPORTED_TEMPERED = [4, 7, 9]

function tempering_one(p, k, seed; tau_active, minimum_active)
    sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active)
    fit = fit_experiment2(sim, k, seed; tau_active=tau_active)

    settings = (
        betas=BETAS,
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_INACTIVE / sqrt(k),
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
        coefficient_prior_sd=COEFFICIENT_PRIOR_SD,
        n_burnin=N_BURNIN,
        n_samples=N_SAMPLES
    )
    runs = [
        "pt_em" => parallel_tempering_time_varying(
            sim.Y, sim.alive, k;
            Lambda_init=[copy(L) for L in fit.Lambda],
            Delta_init=copy(fit.Delta),
            coefficients_init=copy(fit.coefficient_vector),
            settings...,
            rng=MersenneTwister(50_000 + seed)
        ),
        "pt_truth" => parallel_tempering_time_varying(
            sim.Y, sim.alive, k;
            Lambda_init=[copy(L) for L in sim.Lambda],
            Delta_init=copy(sim.Delta),
            coefficients_init=[A0, A_EVENT, A_EVENT],
            settings...,
            rng=MersenneTwister(60_000 + seed)
        )
    ]

    base = (p=p, k=k, seed=seed, tau_active=tau_active, minimum_active=minimum_active)
    no_tempering = (
        flip_acceptance=NaN, swap_acceptance_min=NaN,
        auroc_beta055=NaN, auroc_beta030=NaN, auroc_beta020=NaN
    )
    rows = NamedTuple[]
    push!(rows, merge(
        base, (method="em",),
        shock_metrics(sim, fit.gamma),
        precision_metrics(sim, em_precisions(fit)),
        (a0=fit.coefficient_vector[1], a_birth=fit.coefficient_vector[2],
         a_death=fit.coefficient_vector[3]),
        no_tempering
    ))
    for (name, run) in runs
        tempered_auroc = [shock_metrics(sim, g).auroc for g in run.gamma_by_beta]
        push!(rows, merge(
            base, (method=name,),
            shock_metrics(sim, run.gamma),
            precision_metrics(sim, run.Omega_mean),
            (a0=run.coefficient_mean[1], a_birth=run.coefficient_mean[2],
             a_death=run.coefficient_mean[3]),
            (
                flip_acceptance=run.acceptance.flip,
                swap_acceptance_min=minimum(run.swap_acceptance),
                auroc_beta055=tempered_auroc[REPORTED_TEMPERED[1]],
                auroc_beta030=tempered_auroc[REPORTED_TEMPERED[2]],
                auroc_beta020=tempered_auroc[REPORTED_TEMPERED[3]]
            )
        ))
    end
    return rows
end

function tempering_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        rows = tempering_one(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        summary = join(
            [@sprintf("%s AUROC=%.3f recall=%.2f offRV=%.3f", r.method, r.auroc, r.recall, r.offdiag_rv)
             for r in rows], " | "
        )
        @printf("p=%d k=%d seed=%d  %s  min swap=%.2f  (%.0fs)\n",
            options.p, options.k, seed, summary, rows[2].swap_acceptance_min,
            time() - started)
        if io !== nothing
            if !header_written
                println(io, join(string.(keys(first(rows))), ","))
                header_written = true
            end
            for row in rows
                println(io, join(format_value.(values(row)), ","))
            end
            flush(io)
        end
    end

    io === nothing || close(io)
end

tempering_main()

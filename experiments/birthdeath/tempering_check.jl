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

# TEMPERING_TAU_INACTIVE sets the inactive (spike) scale tau_i of the EM fit and
# the samplers (default 0.075). TEMPERING_FIXED_PI, if set, fixes the shock
# probability at that value (a0 = logit(pi), a_b = a_d = 0) instead of
# estimating it. TEMPERING_FROM_TRUTH=0 skips the run started from the truth;
# TEMPERING_PLAIN=1 adds plain MCMC (temperature 1 only) from the EM fit.
const TAU_I = parse(Float64, get(ENV, "TEMPERING_TAU_INACTIVE", string(TAU_INACTIVE)))
const FIXED_PI = get(ENV, "TEMPERING_FIXED_PI", "")
const FROM_TRUTH = get(ENV, "TEMPERING_FROM_TRUTH", "1") == "1"
const PLAIN = get(ENV, "TEMPERING_PLAIN", "0") == "1"

"""Keyword overrides for fit_experiment2 implementing TAU_I and FIXED_PI."""
function calibration_overrides(k)
    overrides = (
        sigma_inactive=TAU_I / sqrt(k),
        sigma_inactive_start=TAU_I / sqrt(k)
    )
    isempty(FIXED_PI) && return overrides
    pi_fixed = parse(Float64, FIXED_PI)
    return merge(overrides, (
        estimate_shock_coefficients=false,
        a0_init=log(pi_fixed / (1 - pi_fixed)),
        a_birth_init=0.0,
        a_death_init=0.0
    ))
end

function tempering_one(p, k, seed; tau_active, minimum_active)
    sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active)
    fit = fit_experiment2(sim, k, seed; tau_active=tau_active, calibration_overrides(k)...)

    settings = (
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_I / sqrt(k),
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
        coefficient_prior_sd=COEFFICIENT_PRIOR_SD,
        n_burnin=N_BURNIN,
        n_samples=N_SAMPLES,
        estimate_coefficients=isempty(FIXED_PI)
    )
    from_em = (
        Lambda_init=[copy(L) for L in fit.Lambda],
        Delta_init=copy(fit.Delta),
        coefficients_init=copy(fit.coefficient_vector)
    )
    runs = Pair{String, Any}[
        "pt_em" => parallel_tempering_time_varying(
            sim.Y, sim.alive, k;
            from_em..., settings..., betas=BETAS,
            rng=MersenneTwister(50_000 + seed)
        )
    ]
    if PLAIN
        push!(runs, "mcmc_em" => parallel_tempering_time_varying(
            sim.Y, sim.alive, k;
            from_em..., settings..., betas=[1.0],
            rng=MersenneTwister(30_000 + seed)
        ))
    end
    if FROM_TRUTH
        push!(runs, "pt_truth" => parallel_tempering_time_varying(
            sim.Y, sim.alive, k;
            Lambda_init=[copy(L) for L in sim.Lambda],
            Delta_init=copy(sim.Delta),
            coefficients_init=isempty(FIXED_PI) ? [A0, A_EVENT, A_EVENT] :
                copy(fit.coefficient_vector),
            settings..., betas=BETAS,
            rng=MersenneTwister(60_000 + seed)
        ))
    end

    base = (p=p, k=k, seed=seed, tau_active=tau_active, minimum_active=minimum_active,
            tau_inactive=TAU_I, fixed_pi=isempty(FIXED_PI) ? NaN : parse(Float64, FIXED_PI))
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
        tempered(i) = i <= length(tempered_auroc) ? tempered_auroc[i] : NaN
        push!(rows, merge(
            base, (method=name,),
            shock_metrics(sim, run.gamma),
            precision_metrics(sim, run.Omega_mean),
            (a0=run.coefficient_mean[1], a_birth=run.coefficient_mean[2],
             a_death=run.coefficient_mean[3]),
            (
                flip_acceptance=run.acceptance.flip,
                swap_acceptance_min=isempty(run.swap_acceptance) ? NaN :
                    minimum(run.swap_acceptance),
                auroc_beta055=tempered(REPORTED_TEMPERED[1]),
                auroc_beta030=tempered(REPORTED_TEMPERED[2]),
                auroc_beta020=tempered(REPORTED_TEMPERED[3])
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
        flush(stdout)
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

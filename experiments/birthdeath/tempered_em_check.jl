#
# Parallel-tempered EM check for Experiment II shock detection.
#
# The fitting-strategy comparison showed that higher posterior modes than the
# default EM fit's exist and that they call almost no shocks. This script runs
# parallel_tempered_EM (nine EM replicas at temperatures geometric from 1 to
# 0.2, run in parallel on threads, exchanging temperatures) from the default
# EM fit and compares the two on shock detection, precision recovery
# and log posterior.
#
# Usage (from the repository root):
#
#   julia -t 4 --project=. experiments/birthdeath/tempered_em_check.jl \
#       p=10 k=2 seeds=1:20 tau_active=1.5 minimum_active=6 \
#       out=results/experiment2/tempered_em_p10_k2.csv
#
# TEMPERED_EM_EXCHANGE sets the exchange rule: metropolis (default), greedy,
# always or rank. TEMPERED_EM_SPIKE=1 also tempers the inactive (spike) scale.
#

include(joinpath(@__DIR__, "mcmc_check.jl"))

# Exchange rule between replicas: metropolis, greedy, always or rank.
const EXCHANGE = Symbol(get(ENV, "TEMPERED_EM_EXCHANGE", "metropolis"))
# Set TEMPERED_EM_SPIKE=1 to temper the inactive (spike) scale as well.
const TEMPER_SPIKE = get(ENV, "TEMPERED_EM_SPIKE", "0") == "1"
# Set TEMPERED_EM_LIKELIHOOD=0 to leave the likelihood untempered, so that the
# temperatures only set each replica's spike scale.
const TEMPER_LIKELIHOOD = get(ENV, "TEMPERED_EM_LIKELIHOOD", "1") == "1"
# Set TEMPERED_EM_LEVELS=1 to add one output row per temperature level.
const REPORT_LEVELS = get(ENV, "TEMPERED_EM_LEVELS", "0") == "1"

function tempered_em_one(p, k, seed; tau_active, minimum_active)
    sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active)
    fit = fit_experiment2(sim, k, seed; tau_active=tau_active)
    tempered = parallel_tempered_EM(
        sim.Y, sim.alive, k;
        Lambda_init=[copy(L) for L in fit.Lambda],
        Delta_init=copy(fit.Delta),
        coefficients_init=copy(fit.coefficient_vector),
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_INACTIVE / sqrt(k),
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
        coefficient_prior_sd=COEFFICIENT_PRIOR_SD,
        coefficient_damping=COEFFICIENT_DAMPING,
        exchange=EXCHANGE,
        temper_spike=TEMPER_SPIKE,
        temper_likelihood=TEMPER_LIKELIHOOD,
        rng=MersenneTwister(70_000 + seed)
    )

    base = (p=p, k=k, seed=seed, tau_active=tau_active,
            minimum_active=minimum_active, exchange=string(EXCHANGE),
            temper_spike=TEMPER_SPIKE, temper_likelihood=TEMPER_LIKELIHOOD)
    rows = [
        merge(
            base, (method="em",),
            shock_metrics(sim, fit.gamma),
            precision_metrics(sim, em_precisions(fit)),
            (logposterior=fit.final_logposterior, chosen_beta=NaN,
             swap_acceptance_min=NaN, replicas_visiting_beta_one=0)
        ),
        merge(
            base, (method="tempered_em",),
            shock_metrics(sim, tempered.gamma),
            precision_metrics(sim, [L * L' + Diagonal(tempered.Delta) for L in tempered.Lambda]),
            (logposterior=tempered.final_logposterior,
             chosen_beta=tempered.chosen_beta,
             swap_acceptance_min=minimum(tempered.swap_acceptance),
             replicas_visiting_beta_one=tempered.replicas_visiting_beta_one)
        )
    ]
    if REPORT_LEVELS
        for (level, gamma) in enumerate(tempered.level_gamma)
            push!(rows, merge(
                base,
                (method=@sprintf("level%d_tau_i=%.3f", level,
                    tempered.level_sigma_inactive[level] * sqrt(k)),),
                shock_metrics(sim, gamma),
                precision_metrics(sim, [L * L' + Diagonal(tempered.level_Delta[level])
                                        for L in tempered.level_Lambda[level]]),
                (logposterior=NaN, chosen_beta=NaN, swap_acceptance_min=NaN,
                 replicas_visiting_beta_one=0)
            ))
        end
    end
    return rows
end

function tempered_em_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        rows = tempered_em_one(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        em, tempered = rows[1], rows[2]
        @printf(
            "p=%d k=%d seed=%d  AUROC em=%.3f tempered=%.3f  recall em=%.2f tempered=%.2f  logpost gain=%.1f  chosen beta=%.2f  min swap=%.2f  (%.0fs)\n",
            options.p, options.k, seed, em.auroc, tempered.auroc, em.recall,
            tempered.recall, tempered.logposterior - em.logposterior,
            tempered.chosen_beta, tempered.swap_acceptance_min, time() - started
        )
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

tempered_em_main()

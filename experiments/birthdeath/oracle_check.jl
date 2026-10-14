#
# Oracle check for Experiment II shock detection.
#
# For each simulated trajectory, three solutions are compared:
#
#   truth    the true Lambda, Delta and event coefficients, with shock
#            probabilities from one E-step;
#   oracle   generalized EM started from the truth;
#   default  the usual fit (fit_experiment2, binary-splitting initialization).
#
# If the oracle keeps good detection and reaches a higher log posterior than
# the default fit, the default fit is stuck in a worse mode: an optimization
# problem. If the oracle drifts to poor detection, or the default fit reaches
# an equal or higher log posterior, the posterior itself does not single out
# the true shocks: an information or identifiability problem.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/oracle_check.jl \
#       p=10 k=5 seeds=1:20 tau_active=1.5 minimum_active=6 \
#       out=results/experiment2/oracle_p10_k5.csv
#

include(joinpath(@__DIR__, "run_experiment2.jl"))

function posterior_settings(k, tau_active)
    return (
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_INACTIVE / sqrt(k)
    )
end

"""E-step at the current state, then the log posterior used by the fit."""
function final_gamma_and_logposterior(state, coefficients, covariates, k, tau_active)
    sigmas = posterior_settings(k, tau_active)
    pi_time = BGM.shock_probabilities_from_coefficients(
        coefficients, covariates.births, covariates.deaths
    )
    gamma = BGM.E_step_independent_advanced(
        state, pi_time;
        sigma_active=sigmas.sigma_active,
        sigma_inactive=sigmas.sigma_inactive
    )
    logposterior, _ = BGM.full_log_posterior(
        state, gamma, pi_time, coefficients;
        sigmas...,
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
        coefficient_prior_sd=COEFFICIENT_PRIOR_SD
    )
    return gamma, logposterior
end

"""
Generalized EM started from the true parameters. The steps and settings are
those of EM_time_varying in fit_experiment2, except that the E-step also runs
at the first iteration, so the shock probabilities start from the truth.
"""
function oracle_em(sim, k, seed; tau_active, niter=500, tol=1e-6)
    rng = MersenneTwister(20_000 + seed)
    sigmas = posterior_settings(k, tau_active)
    births, deaths = BGM.birth_death_counts(sim.alive)
    covariates = (births=Float64.(births), deaths=Float64.(deaths))

    state = BGM.make_tensor_precision_model(
        [copy(L) for L in sim.Lambda], copy(sim.Delta), sim.Y, sim.alive
    )
    coefficients = [A0, A_EVENT, A_EVENT]
    pi_time = BGM.shock_probabilities_from_coefficients(
        coefficients, covariates.births, covariates.deaths
    )

    for iter in 1:niter
        gamma = BGM.E_step_independent_advanced(
            state, pi_time;
            sigma_active=sigmas.sigma_active,
            sigma_inactive=sigmas.sigma_inactive
        )
        before, _ = BGM.full_log_posterior(
            state, gamma, pi_time, coefficients;
            sigmas..., a_delta=A_DELTA, b_delta=B_DELTA,
            coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
            coefficient_prior_sd=COEFFICIENT_PRIOR_SD
        )

        BGM.update_Delta_bayes!(
            state, gamma, pi_time;
            a=A_DELTA, b=B_DELTA, step=1.0, Delta_min=1e-10, Delta_max=1e6,
            sigmas...
        )
        BGM.update_DiffLambda_hess!(
            state, gamma;
            step=1.0, sigmas..., max_jump=1.0, rng=rng
        )

        if iter > 1
            update = BGM.update_shock_coefficients!(
                copy(coefficients), gamma, sim.alive,
                covariates.births, covariates.deaths;
                prior_mean=COEFFICIENT_PRIOR_MEAN,
                prior_sd=COEFFICIENT_PRIOR_SD
            )
            if update.accepted
                coefficients .= (1 - COEFFICIENT_DAMPING) .* coefficients .+
                                COEFFICIENT_DAMPING .* update.coefficients
                pi_time = BGM.shock_probabilities_from_coefficients(
                    coefficients, covariates.births, covariates.deaths
                )
            end
        end

        after, _ = BGM.full_log_posterior(
            state, gamma, pi_time, coefficients;
            sigmas..., a_delta=A_DELTA, b_delta=B_DELTA,
            coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
            coefficient_prior_sd=COEFFICIENT_PRIOR_SD
        )
        iter > 1 && abs(after - before) <= tol * (1 + abs(before)) && break
    end

    gamma, logposterior = final_gamma_and_logposterior(
        state, coefficients, covariates, k, tau_active
    )
    return (state=state, gamma=gamma, coefficients=coefficients,
            logposterior=logposterior)
end

function oracle_one(p, k, seed; tau_active, minimum_active)
    sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active)
    births, deaths = BGM.birth_death_counts(sim.alive)
    covariates = (births=Float64.(births), deaths=Float64.(deaths))

    truth_state = BGM.make_tensor_precision_model(
        [copy(L) for L in sim.Lambda], copy(sim.Delta), sim.Y, sim.alive
    )
    truth_gamma, truth_logposterior = final_gamma_and_logposterior(
        truth_state, [A0, A_EVENT, A_EVENT], covariates, k, tau_active
    )
    truth = evaluate_experiment2(sim, sim.Lambda, sim.Delta, truth_gamma)

    oracle = oracle_em(sim, k, seed; tau_active=tau_active)
    oracle_metrics = evaluate_experiment2(
        sim, oracle.state.Lambda, oracle.state.Delta, oracle.gamma
    )

    fit = fit_experiment2(sim, k, seed; tau_active=tau_active)
    default = evaluate_experiment2(sim, fit.Lambda, fit.Delta, fit.gamma)

    return (
        p=p, k=k, seed=seed, tau_active=tau_active, minimum_active=minimum_active,
        shock_rate=truth.shock_rate,
        truth_auroc=truth.auroc,
        truth_logposterior=truth_logposterior,
        oracle_auroc=oracle_metrics.auroc,
        oracle_recall=oracle_metrics.recall,
        oracle_precision=oracle_metrics.precision,
        oracle_offdiag_rv=oracle_metrics.offdiag_rv,
        oracle_logposterior=oracle.logposterior,
        default_auroc=default.auroc,
        default_recall=default.recall,
        default_precision=default.precision,
        default_offdiag_rv=default.offdiag_rv,
        default_logposterior=fit.final_logposterior
    )
end

function oracle_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        row = oracle_one(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        @printf(
            "p=%d k=%d seed=%d  AUROC truth=%.3f oracle=%.3f default=%.3f  logpost oracle-default=%.1f  (%.0fs)\n",
            row.p, row.k, row.seed, row.truth_auroc, row.oracle_auroc,
            row.default_auroc, row.oracle_logposterior - row.default_logposterior,
            time() - started
        )
        if io !== nothing
            if !header_written
                println(io, join(string.(keys(row)), ","))
                header_written = true
            end
            println(io, join(format_value.(values(row)), ","))
            flush(io)
        end
    end

    io === nothing || close(io)
end

oracle_main()

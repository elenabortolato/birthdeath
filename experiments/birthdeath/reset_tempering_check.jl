#
# Parallel tempering under the reset prior (shock_center = :zero).
#
# For each seed, data are simulated with reset shocks (correctly specified)
# and with jump shocks, and the reset model is fitted by
#
#   em           generalized EM (fit_experiment2);
#   pt_mcmc      parallel_tempering_time_varying, nine likelihood
#                temperatures geometric from 1 to 0.2, started from the EM fit;
#   tempered_em  parallel_tempered_EM with the same temperatures (Metropolis
#                exchanges), started from the EM fit.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/reset_tempering_check.jl \
#       p=10 k=2 seeds=1:4 tau_active=1.5 minimum_active=6 \
#       out=results/experiment2/reset_tempering_p10_k2_tau150_min6.csv
#
# MCMC_BURNIN and MCMC_SAMPLES set the number of sweeps (default 1000, 2000).
#

include(joinpath(@__DIR__, "mcmc_check.jl"))

const BETAS = collect(0.2 .^ ((0:8) ./ 8))

function reset_tempering_one(p, k, seed; tau_active, minimum_active)
    settings = (
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_INACTIVE / sqrt(k),
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
        coefficient_prior_sd=COEFFICIENT_PRIOR_SD,
        shock_center=:zero
    )
    rows = NamedTuple[]
    for data in (:reset, :jump)
        sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active, shock_type=data)
        fit = fit_experiment2(sim, k, seed; tau_active=tau_active, shock_center=:zero)
        from_em = (
            Lambda_init=[copy(L) for L in fit.Lambda],
            Delta_init=copy(fit.Delta),
            coefficients_init=copy(fit.coefficient_vector)
        )
        pt = parallel_tempering_time_varying(
            sim.Y, sim.alive, k;
            from_em..., settings..., betas=BETAS,
            n_burnin=N_BURNIN, n_samples=N_SAMPLES,
            rng=MersenneTwister(50_000 + seed)
        )
        tempered = parallel_tempered_EM(
            sim.Y, sim.alive, k;
            from_em..., settings..., betas=BETAS,
            coefficient_damping=COEFFICIENT_DAMPING,
            rng=MersenneTwister(70_000 + seed)
        )

        base = (p=p, k=k, seed=seed, tau_active=tau_active,
                minimum_active=minimum_active, data=string(data), model="reset")
        none = (logposterior=NaN, flip_acceptance=NaN, swap_acceptance_min=NaN,
                chosen_beta=NaN)
        push!(rows, merge(
            base, (method="em",),
            shock_metrics(sim, fit.gamma),
            precision_metrics(sim, em_precisions(fit)),
            merge(none, (logposterior=fit.final_logposterior,))
        ))
        push!(rows, merge(
            base, (method="pt_mcmc",),
            shock_metrics(sim, pt.gamma),
            precision_metrics(sim, pt.Omega_mean),
            merge(none, (flip_acceptance=pt.acceptance.flip,
                         swap_acceptance_min=minimum(pt.swap_acceptance)))
        ))
        push!(rows, merge(
            base, (method="tempered_em",),
            shock_metrics(sim, tempered.gamma),
            precision_metrics(sim, [L * L' + Diagonal(tempered.Delta) for L in tempered.Lambda]),
            merge(none, (logposterior=tempered.final_logposterior,
                         swap_acceptance_min=minimum(tempered.swap_acceptance),
                         chosen_beta=tempered.chosen_beta))
        ))
    end
    return rows
end

function reset_tempering_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        rows = reset_tempering_one(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        @printf("p=%d k=%d seed=%d  (%.0fs)\n", options.p, options.k, seed, time() - started)
        for r in rows
            @printf("   data=%-5s %-11s AUROC=%.3f recall=%.2f precision=%.2f recall02=%.2f offRV=%.3f min swap=%.2f\n",
                r.data, r.method, r.auroc, r.recall, r.precision, r.recall_02,
                r.offdiag_rv, r.swap_acceptance_min)
        end
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

reset_tempering_main()

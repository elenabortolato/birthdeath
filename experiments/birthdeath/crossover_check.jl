#
# Shock-crossover EM check for Experiment II shock detection.
#
# Nine EM replicas on the untempered model start from different shock
# patterns (checkerboards with period 4, offsets 0-3, and period 5, offsets
# 0-4) and every 5 iterations exchange the loading path of one variable in
# random pairs (shock_crossover_EM). The script compares the refined best
# replica, and each replica at the end of the crossover phase, with the
# default EM fit.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/crossover_check.jl \
#       p=10 k=2 seeds=1:5 tau_active=1.5 minimum_active=6 \
#       out=results/experiment2/crossover_p10_k2.csv
#

include(joinpath(@__DIR__, "mcmc_check.jl"))

checkerboard(p, T, period, offset) =
    [mod(t - 1 + j + offset, period) == 0 ? 1.0 : 0.0 for j in 1:p, t in 1:T]

function crossover_one(p, k, seed; tau_active, minimum_active)
    sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active)
    fit = fit_experiment2(sim, k, seed; tau_active=tau_active)

    T_ = size(sim.alive, 1)
    patterns = vcat(
        [checkerboard(p, T_, 4, offset) for offset in 0:3],
        [checkerboard(p, T_, 5, offset) for offset in 0:4]
    )
    em_settings = (
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_INACTIVE / sqrt(k),
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
        coefficient_prior_sd=COEFFICIENT_PRIOR_SD,
        coefficient_damping=COEFFICIENT_DAMPING,
        a0_init=-2.0,
        a_birth_init=1.0,
        a_death_init=1.0,
        estimate_shock_coefficients=true,
        coefficient_warmup=1,
        Delta_init_value=0.10,
        step_delta=1.0,
        step_lambda=1.0
    )
    crossover = shock_crossover_EM(
        sim.Y, sim.alive, k;
        patterns=patterns,
        em_settings=em_settings,
        rng=MersenneTwister(80_000 + seed)
    )

    base = (p=p, k=k, seed=seed, tau_active=tau_active, minimum_active=minimum_active)
    extra(logposterior) = (logposterior=logposterior, chosen_replica=crossover.chosen_replica)
    precisions(Lambda, Delta) = [L * L' + Diagonal(Delta) for L in Lambda]
    rows = [
        merge(base, (method="em",), shock_metrics(sim, fit.gamma),
              precision_metrics(sim, em_precisions(fit)), extra(fit.final_logposterior)),
        merge(base, (method="crossover",), shock_metrics(sim, crossover.gamma),
              precision_metrics(sim, precisions(crossover.Lambda, crossover.Delta)),
              extra(crossover.final_logposterior))
    ]
    for r in eachindex(patterns)
        name = r <= 4 ? "replica$(r)_period4_offset$(r - 1)" : "replica$(r)_period5_offset$(r - 5)"
        push!(rows, merge(
            base, (method=name,),
            shock_metrics(sim, crossover.replica_gamma[r]),
            precision_metrics(sim, precisions(crossover.replica_Lambda[r], crossover.replica_Delta[r])),
            extra(crossover.replica_logposterior[r])
        ))
    end
    return rows
end

function crossover_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        rows = crossover_one(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        em, crossover = rows[1], rows[2]
        @printf(
            "p=%d k=%d seed=%d  AUROC em=%.3f crossover=%.3f  recall em=%.2f crossover=%.2f  offRV em=%.3f crossover=%.3f  logpost gain=%.1f  (%.0fs)\n",
            options.p, options.k, seed, em.auroc, crossover.auroc, em.recall,
            crossover.recall, em.offdiag_rv, crossover.offdiag_rv,
            crossover.logposterior - em.logposterior, time() - started
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

crossover_main()

#
# Fitting strategies for Experiment II shock detection.
#
# The oracle check (oracle_check.jl) showed that at k = 2 a better posterior
# mode than the default fit's usually exists. This script compares, on the
# same simulated trajectories, the default fit with:
#
#   anneal              the inactive scale starts at the active scale and
#                       shrinks geometrically to tau_i over 50 iterations;
#   inactive            no binary splitting; every continuation starts with
#                       shock probability 0.05;
#   inactive_anneal     both of the above;
#   restarts            five fits from different random starting loadings
#                       (the default fit and four more), keeping the one with
#                       the highest log posterior.
#
# One CSV row per trajectory and strategy. Log posteriors are comparable across
# strategies: all are evaluated with the target inactive scale after a final
# E-step.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/fit_strategies.jl \
#       p=10 k=2 seeds=1:20 tau_active=1.5 minimum_active=6 \
#       out=results/experiment2/strategies_p10_k2.csv
#

include(joinpath(@__DIR__, "run_experiment2.jl"))

const N_RESTARTS = 5
const ANNEAL_ITERATIONS = 50

function strategy_rows(p, k, seed; tau_active, minimum_active)
    sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active)
    anneal = (
        sigma_inactive_start=tau_active / sqrt(k),
        annealed_iterations=ANNEAL_ITERATIONS
    )
    inactive = (changepoint_init=false, gamma_inactive_init=0.05)

    fits = Pair{String, Any}[]
    restart_fits = [
        fit_experiment2(sim, k, seed; tau_active=tau_active, fit_seed=10_000 + seed + 1_000 * r)
        for r in 0:(N_RESTARTS - 1)
    ]
    push!(fits, "default" => restart_fits[1])
    push!(fits, "restarts" => restart_fits[argmax([f.final_logposterior for f in restart_fits])])
    push!(fits, "anneal" => fit_experiment2(sim, k, seed; tau_active=tau_active, anneal...))
    push!(fits, "inactive" => fit_experiment2(sim, k, seed; tau_active=tau_active, inactive...))
    push!(fits, "inactive_anneal" =>
        fit_experiment2(sim, k, seed; tau_active=tau_active, inactive..., anneal...))

    rows = NamedTuple[]
    for (name, fit) in fits
        metrics = evaluate_experiment2(sim, fit.Lambda, fit.Delta, fit.gamma)
        push!(rows, (
            p=p, k=k, seed=seed, tau_active=tau_active, minimum_active=minimum_active,
            strategy=name,
            auroc=metrics.auroc,
            recall=metrics.recall,
            precision=metrics.precision,
            offdiag_rv=metrics.offdiag_rv,
            rfe=metrics.rfe,
            logposterior=fit.final_logposterior,
            iterations=length(fit.mstep_gain_history)
        ))
    end
    return rows
end

function strategies_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        rows = strategy_rows(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        summary = join(
            [@sprintf("%s=%.3f", r.strategy, r.auroc) for r in rows], "  "
        )
        @printf("p=%d k=%d seed=%d  AUROC %s  (%.0fs)\n",
            options.p, options.k, seed, summary, time() - started)
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

strategies_main()

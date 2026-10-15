#
# Experiment II (paper, Section 7.4): dynamic birth-death factor model, with
# the ablation of the event-informed shock probability (Section 7.6) and the
# sensitivity to the fitted rank (Section 7.7).
#
# Data: simulate_experiment2 (k_true = 5, T = 50, tau_0 = tau_change = 1.5,
# Delta = 0.2, logit pi_t = -2.5 + 0.8 (b_t + d_t)), with at least
# k_true + 1 = 6 active variables at every time.
#
# Fits: tau_a = 0.75, k_fit = 5 by default. Precision recovery is reported on
# the full time-varying active sets and on the common node-time block used
# for the competitors (the longest interval with at least five variables
# active throughout). Shock recovery is reported on all continuation pairs and
# separately at event times (b_t + d_t > 0) and non-event times.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/paper/experiment2.jl p=10 seeds=1:5 \
#       k_fit=5 variant=full out=results/paper/experiment2_p10.csv
#
# variant = full (default) estimates a_b and a_d; variant = ablation fixes
# a_b = a_d = 0.
#

include(joinpath(@__DIR__, "common.jl"))

const K_TRUE = 5

function experiment2_one(p, seed; k_fit, variant)
    sim = simulate_experiment2(p, K_TRUE, seed; minimum_active=K_TRUE + 1)
    Omega_true = [L * L' + Diagonal(sim.Delta) for L in sim.Lambda]
    fits = fit_birthdeath(sim.Y, sim.alive, k_fit, seed;
                          tau_active=0.75, fixed_event_effects=variant == "ablation")
    block = longest_block(sim.alive; min_vars=5)
    births, deaths = birth_death_counts(sim.alive)
    event_times = [t for t in 2:size(sim.alive, 1) if births[t] + deaths[t] > 0]
    nonevent_times = [t for t in 2:size(sim.alive, 1) if births[t] + deaths[t] == 0]

    prefixed(prefix, nt) = NamedTuple{Tuple(Symbol(prefix, key) for key in keys(nt))}(values(nt))
    rows, gamma_rows = NamedTuple[], NamedTuple[]
    for (method, fit) in fits
        base = (experiment="II", variant=variant, p=p, k_true=K_TRUE, k_fit=k_fit, seed=seed)
        push!(rows, merge(
            base, (method=method,),
            full_metrics(Omega_true, fit.Omega, sim.alive),
            prefixed("block_", block_metrics(Omega_true, fit.Omega, sim.alive, block)),
            block_summary(block),
            shock_metrics_paper(fit.gamma, sim.s_true, sim.alive),
            prefixed("event_", shock_metrics_paper(fit.gamma, sim.s_true, sim.alive; times=event_times)),
            prefixed("nonevent_", shock_metrics_paper(fit.gamma, sim.s_true, sim.alive; times=nonevent_times)),
            (a0=fit.coefficients[1], a_birth=fit.coefficients[2],
             a_death=fit.coefficients[3], seconds=fit.seconds)
        ))
        append!(gamma_rows, [
            merge(r, (shock_at_t=t == 1 ? false : any(sim.s_true[:, t]),))
            for (t, r) in enumerate(gamma_time_rows(base, method, fit.gamma, sim.alive))
        ])
    end
    return rows, gamma_rows
end

function experiment2_main()
    options = parse_options(ARGS, Dict(
        "p" => "10", "seeds" => "1:5", "k_fit" => "5", "variant" => "full", "out" => ""
    ))
    options["variant"] in ("full", "ablation") || error("variant must be full or ablation")
    p = parse(Int, options["p"])
    k_fit = parse(Int, options["k_fit"])
    out = options["out"]
    for seed in parse_seeds(options["seeds"])
        rows, gamma_rows = experiment2_one(p, seed; k_fit=k_fit, variant=options["variant"])
        for r in rows
            @printf("II %s p=%d k_fit=%d seed=%d %-4s RV=%.3f blockRV=%.3f offRV=%.3f AUROC=%.3f prec=%.2f rec=%.2f Brier=%.3f (%.0fs)\n",
                options["variant"], p, k_fit, seed, r.method, r.rv, r.block_rv,
                r.offdiag_rv, r.auroc, r.precision, r.recall, r.brier, r.seconds)
        end
        flush(stdout)
        append_rows(out, rows)
        isempty(out) || append_rows(replace(out, ".csv" => "_gamma.csv"), gamma_rows)
    end
end

experiment2_main()

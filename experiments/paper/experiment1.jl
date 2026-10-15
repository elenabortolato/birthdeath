#
# Experiment I (paper, Section 7.3): piecewise-constant sparse precision
# matrices. The active set is fixed, and the precision matrix changes once,
# at t* = 26, between two independent sparse diagonally dominant matrices
# (edge probability 0.15, edge values N(0, 0.9^2), diagonal buffer 0.2).
# The fitted low-rank-plus-diagonal model is therefore misspecified.
#
# Fitted with k = 5 (p = 10) or k = 15 (p = 50), tau_a = 0.75, minimum
# segment length 6, Delta initialized at 2.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/paper/experiment1.jl p=10 seeds=1:5 \
#       out=results/paper/experiment1_p10.csv
#
# Writes one row per trajectory and method to `out`, and the mean shock
# probability by time to the same path with suffix _gamma.csv.
#

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "ControlledPrecisionExperiments.jl"))
using .ControlledPrecisionExperiments: sparse_precision

const T_EXP1 = 50
const CHANGEPOINT = 26

function simulate_experiment1(p, seed)
    rng = MersenneTwister(seed)
    regimes = [
        sparse_precision(p; edge_probability=0.15, edge_scale=0.9,
                         diagonal_buffer=0.2, rng=rng)
        for _ in 1:2
    ]
    Omega = [t < CHANGEPOINT ? regimes[1] : regimes[2] for t in 1:T_EXP1]
    alive = trues(T_EXP1, p)
    Y = zeros(T_EXP1, p)
    for t in 1:T_EXP1
        Y[t, :] .= cholesky(Symmetric(Omega[t])).U \ randn(rng, p)
    end
    return (Y=Y, alive=BitMatrix(alive), Omega=Omega)
end

function experiment1_one(p, seed)
    k = p == 10 ? 5 : 15
    sim = simulate_experiment1(p, seed)
    fits = fit_birthdeath(sim.Y, sim.alive, k, seed;
                          tau_active=0.75, min_segment_length=6, Delta_init_value=2.0)
    rows, gamma_rows = NamedTuple[], NamedTuple[]
    for (method, fit) in fits
        base = (experiment="I", p=p, k=k, seed=seed)
        metrics = full_metrics(sim.Omega, fit.Omega, sim.alive)
        g = mean_gamma_by_time(fit.gamma, sim.alive)
        push!(rows, merge(base, (method=method,), metrics, (
            gamma_at_change=g[CHANGEPOINT],
            gamma_elsewhere=nanmean(g[t] for t in 2:T_EXP1 if t != CHANGEPOINT),
            a0=fit.coefficients[1],
            seconds=fit.seconds
        )))
        append!(gamma_rows, gamma_time_rows(base, method, fit.gamma, sim.alive))
    end
    return rows, gamma_rows
end

function experiment1_main()
    options = parse_options(ARGS, Dict("p" => "10", "seeds" => "1:5", "out" => ""))
    p = parse(Int, options["p"])
    out = options["out"]
    for seed in parse_seeds(options["seeds"])
        rows, gamma_rows = experiment1_one(p, seed)
        for r in rows
            @printf("I p=%d seed=%d %-4s RV=%.3f RFE=%.3f offRV=%.3f gamma(t*)=%.3f other=%.3f (%.0fs)\n",
                p, seed, r.method, r.rv, r.rfe, r.offdiag_rv, r.gamma_at_change,
                r.gamma_elsewhere, r.seconds)
        end
        flush(stdout)
        append_rows(out, rows)
        isempty(out) || append_rows(replace(out, ".csv" => "_gamma.csv"), gamma_rows)
    end
end

experiment1_main()

#
# Experiments III and IV (paper, Sections 7.5 and 3.x): endogenous structural
# absence and contemporaneous selection.
#
# Data (p = 10, T = 80, k = 3): Omega_t = Lambda_t Lambda_t' + Diag(Delta) on
# the active set, Delta = 0.2, loadings N(0, 1.5^2/k I) at the start of each
# active spell, and a continuation shock N(0, 1.5^2/k I) with probability
# 0.06. All variables are active at t = 1. A variable absent for at least two
# periods re-enters with probability 0.2, with new loadings. With
# z_tj = Y_tj / sqrt(Sigma_t,jj) the standardized observation and c = -1.25,
# a variable active at t - 1 exits at t under
#
#   history          (Experiment III; IV, mechanism i)  z_{t-1,j} <= c;
#   stochastic       (IV, mechanism ii) with probability
#                    logistic(alpha + beta z_{t-1,j}), beta = -2, alpha set so
#                    that the marginal exit rate equals Phi(c), as in (i);
#   contemporaneous  (IV, mechanism iii) z_tj <= c, where Y_t is drawn from
#                    the Gaussian on all candidate variables and the
#                    variables that exit are dropped from the observed vector.
#
# An exit that would leave fewer than k + 1 active variables is cancelled
# (the variables with the largest z are kept), so that p_t > k throughout.
#
# Fits: k = 3, tau_a = 0.75, as in Experiment II. Precision recovery is
# reported on the full active sets, and on a common node-time block (at least
# five variables over at least ten time points, relaxed when no such block
# exists) against both the principal block and the marginal precision of the
# block variables. For Experiment IV the bias of Delta-hat and the recovery
# among variables far from the selection boundary are also reported: at each
# t, the principal block of the active variables with z_tj > c + 1.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/paper/experiment3.jl mechanism=history \
#       seeds=1:5 out=results/paper/experiment3_history.csv
#

include(joinpath(@__DIR__, "common.jl"))
using Distributions: Normal, cdf

const P_EXP3 = 10
const T_EXP3 = 80
const K_EXP3 = 3
const THRESHOLD = -1.25
const STOCHASTIC_BETA = -2.0
const SAFE_MARGIN = 1.0

logistic(x) = 1 / (1 + exp(-x))

"""Intercept alpha with E[logistic(alpha + beta Z)] = Phi(c), Z ~ N(0, 1)."""
function calibrate_alpha(c, beta)
    target = cdf(Normal(), c)
    grid = range(-8, 8; length=4001)
    weights = exp.(-grid .^ 2 ./ 2)
    weights ./= sum(weights)
    rate(alpha) = sum(weights .* logistic.(alpha .+ beta .* grid))
    low, high = -20.0, 20.0
    for _ in 1:100
        mid = (low + high) / 2
        rate(mid) < target ? (low = mid) : (high = mid)
    end
    return (low + high) / 2
end

const STOCHASTIC_ALPHA = calibrate_alpha(THRESHOLD, STOCHASTIC_BETA)

"""Draw y on the index set S and return it with its standardized values."""
function draw_standardized(Lambda_t, Delta, S, rng)
    Omega = Symmetric(Lambda_t[S, :] * Lambda_t[S, :]' + Diagonal(Delta[S]))
    y = cholesky(Omega).U \ randn(rng, length(S))
    z = y ./ sqrt.(diag(inv(Omega)))
    return y, z
end

function simulate_endogenous(seed; mechanism)
    rng = MersenneTwister(seed)
    p, T_, k = P_EXP3, T_EXP3, K_EXP3
    sd0, sdc = TAU0 / sqrt(k), TAU_CHANGE / sqrt(k)
    Delta = fill(DELTA, p)
    alive = falses(T_, p)
    Lambda = [zeros(p, k) for _ in 1:T_]
    s_true = falses(p, T_)
    Y = zeros(T_, p)
    z = fill(NaN, T_, p)
    absent_for = zeros(Int, p)
    exits_attempted = 0
    exits_cancelled = 0
    at_risk = 0

    alive[1, :] .= true
    for j in 1:p
        Lambda[1][j, :] .= sd0 .* randn(rng, k)
    end
    S = collect(1:p)
    y, zs = draw_standardized(Lambda[1], Delta, S, rng)
    Y[1, S] .= y
    z[1, S] .= zs

    for t in 2:T_
        continuing = [j for j in 1:p if alive[t - 1, j]]
        reentering = [j for j in 1:p if !alive[t - 1, j] && absent_for[j] >= 2 &&
                      rand(rng) < 0.20]
        for j in continuing
            Lambda[t][j, :] .= Lambda[t - 1][j, :]
            if rand(rng) < 0.06
                Lambda[t][j, :] .+= sdc .* randn(rng, k)
                s_true[j, t] = true
            end
        end
        for j in reentering
            Lambda[t][j, :] .= sd0 .* randn(rng, k)
        end
        candidates = sort(vcat(continuing, reentering))
        at_risk += length(continuing)

        if mechanism == :contemporaneous
            y, zs = draw_standardized(Lambda[t], Delta, candidates, rng)
            z_now = Dict(zip(candidates, zs))
            exits = [j for j in continuing if z_now[j] <= THRESHOLD]
            keep_score = j -> z_now[j]
        else
            exits = if mechanism == :history
                [j for j in continuing if z[t - 1, j] <= THRESHOLD]
            else
                [j for j in continuing
                 if rand(rng) < logistic(STOCHASTIC_ALPHA + STOCHASTIC_BETA * z[t - 1, j])]
            end
            keep_score = j -> z[t - 1, j]
        end
        exits_attempted += length(exits)
        # Cancel exits, highest z first, until k + 1 variables stay active.
        sort!(exits; by=keep_score, rev=true)
        while length(candidates) - length(exits) < k + 1 && !isempty(exits)
            popfirst!(exits)
            exits_cancelled += 1
        end
        active = setdiff(candidates, exits)
        alive[t, active] .= true

        if mechanism == :contemporaneous
            for (i, j) in enumerate(candidates)
                j in active || continue
                Y[t, j] = y[i]
                z[t, j] = zs[i]
            end
        else
            y, zs = draw_standardized(Lambda[t], Delta, active, rng)
            Y[t, active] .= y
            z[t, active] .= zs
        end
        for j in 1:p
            absent_for[j] = alive[t, j] ? 0 : absent_for[j] + 1
            alive[t, j] || (s_true[j, t] = false; fill!(view(Lambda[t], j, :), 0.0))
        end
    end

    return (
        Y=Y, alive=BitMatrix(alive), Lambda=Lambda, Delta=Delta, s_true=s_true, z=z,
        exit_rate=(exits_attempted - exits_cancelled) / max(at_risk, 1),
        exits_cancelled=exits_cancelled
    )
end

"""Active variables at t whose standardized value is above c + SAFE_MARGIN."""
safe_at(sim, t) = [
    j for j in 1:size(sim.alive, 2)
    if sim.alive[t, j] && sim.z[t, j] > THRESHOLD + SAFE_MARGIN
]

function experiment3_one(seed; mechanism)
    sim = simulate_endogenous(seed; mechanism=mechanism)
    k = K_EXP3
    Omega_true = [L * L' + Diagonal(sim.Delta) for L in sim.Lambda]
    fits = fit_birthdeath(sim.Y, sim.alive, k, seed; tau_active=0.75)
    block = adaptive_block(sim.alive; min_vars=5, min_length=10)
    ever_active = [j for j in 1:size(sim.alive, 2) if any(view(sim.alive, :, j))]
    prefixed(prefix, nt) = NamedTuple{Tuple(Symbol(prefix, key) for key in keys(nt))}(values(nt))

    rows = NamedTuple[]
    for (method, fit) in fits
        safe_metrics = average_metrics(
            t -> Omega_true[t][safe_at(sim, t), safe_at(sim, t)],
            t -> fit.Omega[t][safe_at(sim, t), safe_at(sim, t)],
            1:size(sim.alive, 1); min_size=2
        )
        push!(rows, merge(
            (experiment=mechanism == :history ? "III" : "IV", mechanism=string(mechanism),
             p=P_EXP3, k=k, seed=seed, method=method,
             exit_rate=sim.exit_rate, exits_cancelled=sim.exits_cancelled,
             mean_active=mean(sum(sim.alive; dims=2))),
            full_metrics(Omega_true, fit.Omega, sim.alive),
            prefixed("block_", block_metrics(Omega_true, fit.Omega, sim.alive, block)),
            prefixed("blockmarg_", block_metrics(Omega_true, fit.Omega, sim.alive, block; marginal=true)),
            block_summary(block),
            prefixed("safe_", safe_metrics),
            (mean_safe=mean(length(safe_at(sim, t)) for t in 1:size(sim.alive, 1)),
             delta_bias=mean(fit.Delta[ever_active] .- sim.Delta[ever_active]),
             delta_relative_bias=mean(fit.Delta[ever_active] ./ sim.Delta[ever_active] .- 1)),
            shock_metrics_paper(fit.gamma, sim.s_true, sim.alive),
            (seconds=fit.seconds,)
        ))
    end
    return rows
end

function experiment3_main()
    options = parse_options(ARGS, Dict("mechanism" => "history", "seeds" => "1:5", "out" => ""))
    mechanism = Symbol(options["mechanism"])
    mechanism in (:history, :stochastic, :contemporaneous) ||
        error("mechanism must be history, stochastic or contemporaneous")
    for seed in parse_seeds(options["seeds"])
        rows = experiment3_one(seed; mechanism=mechanism)
        for r in rows
            @printf("%s %s seed=%d %-4s exit=%.3f RV=%.3f offRV=%.3f blockmargRV=%.3f safeRV=%.3f dBias=%.3f AUROC=%.3f (%.0fs)\n",
                r.experiment, r.mechanism, seed, r.method, r.exit_rate, r.rv, r.offdiag_rv,
                r.blockmarg_rv, r.safe_rv, r.delta_bias, r.auroc, r.seconds)
        end
        flush(stdout)
        append_rows(options["out"], rows)
    end
end

experiment3_main()

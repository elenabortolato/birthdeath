#
# Experiment II: dynamic birth-death factor model (paper, Section 7.4).
#
# Data are generated from the low-rank-plus-diagonal model with node births
# and deaths, and the estimator is fitted with the paper's settings. One CSV
# row is written per simulated trajectory.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/run_experiment2.jl \
#       p=10 seeds=1:20 tau_active=0.75 minimum_active=6 out=results/experiment2/run.csv
#
# Arguments (all optional):
#   p               number of variables, 10 or 50                (default 10)
#   seeds           a:b range or comma-separated list             (default 1:5)
#   tau_active      fitted active innovation scale tau_a          (default 0.75)
#   minimum_active  lower bound on p_t; k+1 respects (S4)         (default K+1)
#   out             output CSV path                               (default stdout only)
#

using LinearAlgebra, Random, Statistics, Printf

include(joinpath(@__DIR__, "..", "..", "src", "BirthDeathGraphicalModels.jl"))
using .BirthDeathGraphicalModels
const BGM = BirthDeathGraphicalModels

# Data-generating values (paper, Section 7.4).
const K = 5
const T = 50
const TAU0 = 1.5
const TAU_CHANGE = 1.5
const DELTA = 0.2
const A0 = -2.5
const A_EVENT = 0.8

function parse_args(args)
    options = Dict(
        "p" => "10", "seeds" => "1:5", "tau_active" => "0.75",
        "minimum_active" => string(K + 1), "out" => ""
    )
    for arg in args
        key, value = split(arg, "="; limit=2)
        haskey(options, key) || error("unknown argument $key")
        options[key] = value
    end
    seeds = occursin(":", options["seeds"]) ?
        collect(range(parse.(Int, split(options["seeds"], ":"))...)) :
        parse.(Int, split(options["seeds"], ","))
    return (
        p=parse(Int, options["p"]),
        seeds=seeds,
        tau_active=parse(Float64, options["tau_active"]),
        minimum_active=parse(Int, options["minimum_active"]),
        out=options["out"]
    )
end

function auroc(scores, labels)
    pos = scores[labels]
    neg = scores[.!labels]
    (isempty(pos) || isempty(neg)) && return NaN
    return mean(x > y ? 1.0 : x == y ? 0.5 : 0.0 for x in pos, y in neg)
end

nanmean(x) = (v = filter(isfinite, collect(Float64, x)); isempty(v) ? NaN : mean(v))
ratio(a, b) = b == 0 ? NaN : a / b

function run_one(p, seed; tau_active, minimum_active)
    data_rng = MersenneTwister(seed)
    fit_rng = MersenneTwister(10_000 + seed)

    # Lower event probability at p = 50 keeps the number of events per time
    # point comparable across dimensions.
    event_probability = p == 10 ? 0.05 : 0.01

    birth_time, death_time = generate_birth_death_times(
        p, T;
        birth_probability=event_probability,
        death_probability=event_probability,
        initial_alive_fraction=0.7,
        minimum_lifetime=15,
        maximum_lifetime=50,
        allow_birth_at_T=false,
        minimum_active=minimum_active,
        rng=data_rng
    )
    sim = generate_birth_death_data(
        p, K, T, birth_time, death_time;
        a0=A0, a_event=A_EVENT,
        initial_sd=TAU0 / sqrt(K), change_sd=TAU_CHANGE / sqrt(K),
        Delta_value=DELTA, rng=data_rng
    )
    alive = sim.alive

    Lambda_init = initialize_loading_paths(alive, K; scale=0.2, rng=fit_rng)
    fit = EM_time_varying(
        sim.Y, alive, K;
        Lambda_init=Lambda_init,
        changepoint_init=true,
        min_segment_length=1,
        bic_threshold=0.0,
        changepoint_max_depth=6,
        changepoint_maxiter=500,
        changepoint_prior_var=1.0,
        changepoint_learning_rate=0.5,
        gamma_active_init=0.8,
        gamma_inactive_init=0.4,
        gamma_floor=1e-10,
        a0_init=-2.0,
        a_birth_init=1.0,
        a_death_init=1.0,
        estimate_shock_coefficients=true,
        coefficient_prior_mean=[0.0, 0.5, 0.5],
        coefficient_prior_sd=[2.0, 0.5, 0.5],
        coefficient_warmup=1,
        coefficient_update_interval=1,
        coefficient_damping=0.05,
        sigma_initial=1.5 / sqrt(K),
        sigma_active=tau_active / sqrt(K),
        sigma_inactive=0.075 / sqrt(K),
        sigma_inactive_start=0.075 / sqrt(K),
        annealed_iterations=1,
        a_delta=1.0,
        b_delta=2.0,
        Delta_init_value=0.10,
        step_delta=1.0,
        step_lambda=1.0,
        niter=500,
        tol=1e-6,
        verbose=false,
        rng=fit_rng
    )

    errors = active_precision_errors(
        sim.Lambda, sim.Delta, fit.Lambda, fit.Delta, alive, 1:T
    )

    # Shock detection is evaluated on continuation pairs only.
    continuations = [
        (j, t) for t in 2:T for j in 1:p if alive[t, j] && alive[t - 1, j]
    ]
    gamma = [fit.gamma[j, t] for (j, t) in continuations]
    truth = BitVector([sim.s_true[j, t] for (j, t) in continuations])
    called = gamma .> 0.5
    tp = count(called .& truth)
    fp = count(called .& .!truth)
    fn = count(.!called .& truth)
    tn = count(.!called .& .!truth)

    births, deaths = BGM.birth_death_counts(alive)
    active_counts = vec(sum(alive, dims=2))

    return (
        p=p,
        seed=seed,
        tau_active=tau_active,
        minimum_active=minimum_active,
        births=sum(births),
        deaths=sum(deaths),
        min_active=minimum(active_counts),
        rv=nanmean(x.rv for x in errors),
        offdiag_rv=nanmean(x.offdiagonal_rv for x in errors),
        partial_rv=nanmean(x.partial_correlation_rv for x in errors),
        rfe=nanmean(x.relative_frobenius_error for x in errors),
        auroc=auroc(gamma, truth),
        precision=ratio(tp, tp + fp),
        recall=ratio(tp, tp + fn),
        specificity=ratio(tn, tn + fp),
        shock_rate=mean(truth),
        n_continuations=length(continuations),
        gamma_true_shocks=nanmean(gamma[truth]),
        gamma_non_shocks=nanmean(gamma[.!truth]),
        a0=fit.coefficient_vector[1],
        a_birth=fit.coefficient_vector[2],
        a_death=fit.coefficient_vector[3]
    )
end

format_value(x::AbstractFloat) = @sprintf("%.6g", x)
format_value(x) = string(x)

function main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        row = run_one(
            options.p, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        @printf(
            "p=%d seed=%d  RV=%.3f  offdiag RV=%.3f  partial RV=%.3f  RFE=%.3f  AUROC=%.3f  (%.0fs)\n",
            row.p, row.seed, row.rv, row.offdiag_rv, row.partial_rv,
            row.rfe, row.auroc, time() - started
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

main()

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
#       p=10 k=5 seeds=1:20 tau_active=0.75 minimum_active=6 \
#       out=results/experiment2/run.csv
#
# Arguments (all optional):
#   p               number of variables, 10 or 50                (default 10)
#   k               latent rank, simulated and fitted             (default 5)
#   seeds           a:b range or comma-separated list             (default 1:5)
#   tau_active      fitted active innovation scale tau_a          (default 0.75)
#   minimum_active  lower bound on p_t; k+1 respects (S4)         (default k+1)
#   out             output CSV path                               (default stdout only)
#
# The functions simulate_experiment2, fit_experiment2 and
# evaluate_experiment2 are also used by oracle_check.jl.
#

using LinearAlgebra, Random, Statistics, Printf

include(joinpath(@__DIR__, "..", "..", "src", "BirthDeathGraphicalModels.jl"))
using .BirthDeathGraphicalModels
const BGM = BirthDeathGraphicalModels

# Data-generating values (paper, Section 7.4).
const T = 50
const TAU0 = 1.5
const TAU_CHANGE = 1.5
const DELTA = 0.2
const A0 = -2.5
const A_EVENT = 0.8

# Fitted prior and initialization settings (paper, Section 7.4).
const TAU_INACTIVE = 0.075
const A_DELTA = 1.0
const B_DELTA = 2.0
const COEFFICIENT_PRIOR_MEAN = [0.0, 0.5, 0.5]
const COEFFICIENT_PRIOR_SD = [2.0, 0.5, 0.5]
const COEFFICIENT_DAMPING = 0.05

function parse_args(args)
    options = Dict(
        "p" => "10", "k" => "5", "seeds" => "1:5", "tau_active" => "0.75",
        "minimum_active" => "", "out" => ""
    )
    for arg in args
        key, value = split(arg, "="; limit=2)
        haskey(options, key) || error("unknown argument $key")
        options[key] = value
    end
    k = parse(Int, options["k"])
    seeds = occursin(":", options["seeds"]) ?
        collect(range(parse.(Int, split(options["seeds"], ":"))...)) :
        parse.(Int, split(options["seeds"], ","))
    minimum_active = isempty(options["minimum_active"]) ?
        k + 1 : parse(Int, options["minimum_active"])
    return (
        p=parse(Int, options["p"]),
        k=k,
        seeds=seeds,
        tau_active=parse(Float64, options["tau_active"]),
        minimum_active=minimum_active,
        out=options["out"]
    )
end

auroc(scores, labels) = begin
    pos = scores[labels]
    neg = scores[.!labels]
    (isempty(pos) || isempty(neg)) ? NaN :
        mean(x > y ? 1.0 : x == y ? 0.5 : 0.0 for x in pos, y in neg)
end

nanmean(x) = (v = filter(isfinite, collect(Float64, x)); isempty(v) ? NaN : mean(v))
ratio(a, b) = b == 0 ? NaN : a / b

"""
Simulate one Experiment II trajectory. `shock_rule = :network` makes each
variable's shock probability depend on its partial correlations with the
variables born and dying at t, with slopes `a_birth` and `a_death`.
"""
function simulate_experiment2(
    p, k, seed;
    minimum_active,
    shock_rule=:counts,
    a_birth=A_EVENT,
    a_death=A_EVENT
)
    data_rng = MersenneTwister(seed)

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
    return generate_birth_death_data(
        p, k, T, birth_time, death_time;
        a0=A0, a_event=A_EVENT,
        initial_sd=TAU0 / sqrt(k), change_sd=TAU_CHANGE / sqrt(k),
        Delta_value=DELTA,
        shock_rule=shock_rule, a_birth=a_birth, a_death=a_death,
        rng=data_rng
    )
end

"""
Fit the estimator to a simulated trajectory with the paper's settings.
`shock_covariates = :network` fits the network-informed shock probabilities.
"""
function fit_experiment2(
    sim, k, seed;
    tau_active,
    shock_covariates=:counts,
    coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
    coefficient_prior_sd=COEFFICIENT_PRIOR_SD
)
    fit_rng = MersenneTwister(10_000 + seed)
    Lambda_init = initialize_loading_paths(sim.alive, k; scale=0.2, rng=fit_rng)
    return EM_time_varying(
        sim.Y, sim.alive, k;
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
        shock_covariates=shock_covariates,
        coefficient_prior_mean=coefficient_prior_mean,
        coefficient_prior_sd=coefficient_prior_sd,
        coefficient_warmup=1,
        coefficient_update_interval=1,
        coefficient_damping=COEFFICIENT_DAMPING,
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_INACTIVE / sqrt(k),
        sigma_inactive_start=TAU_INACTIVE / sqrt(k),
        annealed_iterations=1,
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        Delta_init_value=0.10,
        step_delta=1.0,
        step_lambda=1.0,
        niter=500,
        tol=1e-6,
        verbose=false,
        rng=fit_rng
    )
end

"""
Precision-recovery and shock-detection metrics. `Lambda`, `Delta` and `gamma`
are the fitted loadings, diagonal precisions and shock probabilities.
"""
function evaluate_experiment2(sim, Lambda, Delta, gamma)
    alive = sim.alive
    T_, p = size(alive)
    errors = active_precision_errors(sim.Lambda, sim.Delta, Lambda, Delta, alive, 1:T_)

    # Shock detection is evaluated on continuation pairs only.
    continuations = [
        (j, t) for t in 2:T_ for j in 1:p if alive[t, j] && alive[t - 1, j]
    ]
    scores = [gamma[j, t] for (j, t) in continuations]
    truth = BitVector([sim.s_true[j, t] for (j, t) in continuations])
    called = scores .> 0.5
    tp = count(called .& truth)
    fp = count(called .& .!truth)
    fn = count(.!called .& truth)
    tn = count(.!called .& .!truth)

    return (
        rv=nanmean(x.rv for x in errors),
        offdiag_rv=nanmean(x.offdiagonal_rv for x in errors),
        partial_rv=nanmean(x.partial_correlation_rv for x in errors),
        rfe=nanmean(x.relative_frobenius_error for x in errors),
        auroc=auroc(scores, truth),
        precision=ratio(tp, tp + fp),
        recall=ratio(tp, tp + fn),
        specificity=ratio(tn, tn + fp),
        shock_rate=mean(truth),
        n_continuations=length(continuations),
        gamma_true_shocks=nanmean(scores[truth]),
        gamma_non_shocks=nanmean(scores[.!truth])
    )
end

function run_one(p, k, seed; tau_active, minimum_active)
    sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active)
    fit = fit_experiment2(sim, k, seed; tau_active=tau_active)
    metrics = evaluate_experiment2(sim, fit.Lambda, fit.Delta, fit.gamma)

    births, deaths = BGM.birth_death_counts(sim.alive)
    active_counts = vec(sum(sim.alive, dims=2))

    return merge(
        (
            p=p, k=k, seed=seed,
            tau_active=tau_active,
            minimum_active=minimum_active,
            births=sum(births),
            deaths=sum(deaths),
            min_active=minimum(active_counts)
        ),
        metrics,
        (
            a0=fit.coefficient_vector[1],
            a_birth=fit.coefficient_vector[2],
            a_death=fit.coefficient_vector[3]
        )
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
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        @printf(
            "p=%d k=%d seed=%d  RV=%.3f  offdiag RV=%.3f  partial RV=%.3f  RFE=%.3f  AUROC=%.3f  (%.0fs)\n",
            row.p, row.k, row.seed, row.rv, row.offdiag_rv, row.partial_rv,
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

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

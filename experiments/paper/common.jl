#
# Shared code for the paper's simulation study (Section 7): fitting, metrics,
# common node-time blocks and CSV output. Each experiment script includes it.
#
# Every trajectory is fitted twice:
#
#   em    generalized EM with the paper's settings (point estimate);
#   mcmc  one Metropolis-within-Gibbs chain (MCMC_time_varying) started from
#         the EM fit, MCMC_BURNIN burn-in and MCMC_SAMPLES sampling sweeps
#         (default 500 and 1000), reporting posterior means of Omega_t and of
#         the shock indicators.
#

include(joinpath(@__DIR__, "..", "birthdeath", "run_experiment2.jl"))

const N_BURNIN = parse(Int, get(ENV, "MCMC_BURNIN", "500"))
const N_SAMPLES = parse(Int, get(ENV, "MCMC_SAMPLES", "1000"))
const RUN_MCMC = get(ENV, "PAPER_MCMC", "1") == "1"

"""Parse key=value arguments over a dictionary of defaults (all strings)."""
function parse_options(args, defaults::Dict{String, String})
    options = copy(defaults)
    for arg in args
        key, value = split(arg, "="; limit=2)
        haskey(options, key) || error("unknown argument $key")
        options[key] = value
    end
    return options
end

parse_seeds(text) = occursin(":", text) ?
    collect(range(parse.(Int, split(text, ":"))...)) :
    parse.(Int, split(text, ","))

# ----------------------------------------------------------------------------
# Fitting
# ----------------------------------------------------------------------------

"""
Fit EM with the paper's settings and, if RUN_MCMC, a plain MCMC chain from the
EM fit. Returns the precision estimates and shock probabilities of each.

`fixed_event_effects = true` is the ablation of Section 7.6: a_b = a_d = 0,
so the shock probability is constant over time.
"""
function fit_birthdeath(
    Y, alive, k, seed;
    tau_active=0.75,
    min_segment_length=1,
    Delta_init_value=0.10,
    fixed_event_effects=false,
    coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
    coefficient_prior_sd=COEFFICIENT_PRIOR_SD
)
    fit_rng = MersenneTwister(10_000 + seed)
    Lambda_init = initialize_loading_paths(alive, k; scale=0.2, rng=fit_rng)
    prior_mean = copy(coefficient_prior_mean)
    prior_sd = copy(coefficient_prior_sd)
    event_init = 1.0
    if fixed_event_effects
        # A prior this tight pins a_b and a_d at zero in the EM update.
        prior_mean[2:3] .= 0.0
        prior_sd[2:3] .= 1e-6
        event_init = 0.0
    end
    started = time()
    em = EM_time_varying(
        Y, alive, k;
        Lambda_init=Lambda_init,
        changepoint_init=true,
        min_segment_length=min_segment_length,
        bic_threshold=0.0,
        changepoint_max_depth=6,
        changepoint_maxiter=500,
        changepoint_prior_var=1.0,
        changepoint_learning_rate=0.5,
        gamma_active_init=0.8,
        gamma_inactive_init=0.4,
        gamma_floor=1e-10,
        a0_init=-2.0,
        a_birth_init=event_init,
        a_death_init=event_init,
        estimate_shock_coefficients=true,
        coefficient_prior_mean=prior_mean,
        coefficient_prior_sd=prior_sd,
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
        Delta_init_value=Delta_init_value,
        step_delta=1.0,
        step_lambda=1.0,
        niter=500,
        tol=1e-6,
        verbose=false,
        rng=fit_rng
    )
    em_seconds = time() - started
    fits = Pair{String, Any}[
        "em" => (
            Omega=[L * L' + Diagonal(em.Delta) for L in em.Lambda],
            Delta=copy(em.Delta),
            gamma=em.gamma,
            coefficients=copy(em.coefficient_vector),
            seconds=em_seconds
        )
    ]
    RUN_MCMC || return fits

    started = time()
    chain = MCMC_time_varying(
        Y, alive, k;
        Lambda_init=[copy(L) for L in em.Lambda],
        Delta_init=copy(em.Delta),
        coefficients_init=copy(em.coefficient_vector),
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_INACTIVE / sqrt(k),
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        coefficient_prior_mean=prior_mean,
        coefficient_prior_sd=prior_sd,
        coefficient_mask=fixed_event_effects ? [true, false, false] : [true, true, true],
        n_burnin=N_BURNIN,
        n_samples=N_SAMPLES,
        rng=MersenneTwister(30_000 + seed)
    )
    push!(fits, "mcmc" => (
        Omega=chain.Omega_mean,
        Delta=chain.Delta_mean,
        gamma=chain.gamma,
        coefficients=chain.coefficient_mean,
        seconds=time() - started
    ))
    return fits
end

# ----------------------------------------------------------------------------
# Precision-recovery metrics (Section 7.2)
# ----------------------------------------------------------------------------

"""
Metrics comparing one estimated precision matrix with the truth: RV
(eq. rv_metric), RFE, matrix MSE, norm ratio q_t, scale-adjusted RFE, and the
RV of the off-diagonal entries and of the partial correlations.
"""
function matrix_metrics(truth::AbstractMatrix, estimate::AbstractMatrix)
    difference = estimate - truth
    c = dot(estimate, truth) / max(dot(estimate, estimate), eps())
    n = size(truth, 1)
    return (
        rv=rv_coefficient(truth, estimate),
        rfe=norm(difference) / norm(truth),
        mse=sum(abs2, difference) / n^2,
        norm_ratio=norm(estimate) / norm(truth),
        rfe_adjusted=norm(c * estimate - truth) / norm(truth),
        offdiag_rv=n >= 2 ? offdiagonal_rv(truth, estimate) : NaN,
        partial_rv=n >= 2 ? partial_correlation_rv(truth, estimate) : NaN
    )
end

"""
Average of `matrix_metrics` over time points, comparing `truth(t)` and
`estimate(t)` (matrices on the same index set) at each t in `times`. Time
points with fewer than `min_size` variables are skipped.
"""
function average_metrics(truth, estimate, times; min_size=1)
    per_time = [
        matrix_metrics(truth(t), estimate(t))
        for t in times if size(truth(t), 1) >= min_size
    ]
    isempty(per_time) && return (
        rv=NaN, rfe=NaN, mse=NaN, norm_ratio=NaN, rfe_adjusted=NaN,
        offdiag_rv=NaN, partial_rv=NaN, n_matrices=0
    )
    keys_ = keys(first(per_time))
    means = NamedTuple{keys_}(Tuple(nanmean(getfield(m, key) for m in per_time) for key in keys_))
    return merge(means, (n_matrices=length(per_time),))
end

"""Metrics on the full time-varying active sets."""
full_metrics(Omega_true, Omega_est, alive) = average_metrics(
    t -> Omega_true[t][findall(view(alive, t, :)), findall(view(alive, t, :))],
    t -> Omega_est[t][findall(view(alive, t, :)), findall(view(alive, t, :))],
    [t for t in 1:size(alive, 1) if any(view(alive, t, :))]
)

"""
Metrics on a fixed node-time block (a, b, S). With `marginal = false` the
target is the principal block Omega_t[S, S]; with `marginal = true` it is the
marginal precision [(Omega_t^{-1})_{S,S}]^{-1} of the variables in S, the
other active variables being integrated out (Experiment III).
"""
function block_metrics(Omega_true, Omega_est, alive, block; marginal=false)
    block === nothing && return average_metrics(t -> zeros(0, 0), t -> zeros(0, 0), Int[])
    function restrict(Omega, t)
        marginal || return Omega[block.S, block.S]
        active = findall(view(alive, t, :))
        position = [findfirst(==(j), active) for j in block.S]
        Sigma = inv(Symmetric(Omega[active, active]))
        return Matrix(inv(Symmetric(Sigma[position, position])))
    end
    return average_metrics(
        t -> restrict(Omega_true[t], t),
        t -> restrict(Omega_est[t], t),
        block.a:block.b
    )
end

"""
Longest contiguous interval [a, b] with at least `min_vars` variables active
throughout; ties go to the larger node set. Returns (a, b, S) or nothing.
"""
function longest_block(alive; min_vars=5)
    T_, p = size(alive)
    best = nothing
    for a in 1:T_, b in a:T_
        S = [j for j in 1:p if all(view(alive, a:b, j))]
        length(S) >= min_vars || break
        if best === nothing || b - a > best.b - best.a ||
           (b - a == best.b - best.a && length(S) > length(best.S))
            best = (a=a, b=b, S=S)
        end
    end
    return best
end

"""
Experiment III block: at least `min_vars` common variables over at least
`min_length` time points, relaxing the number of variables down to 2 when no
such block exists, and otherwise the longest block with at least 2 variables.
"""
function adaptive_block(alive; min_vars=5, min_length=10)
    for m in min_vars:-1:2
        block = longest_block(alive; min_vars=m)
        block !== nothing && block.b - block.a + 1 >= min_length && return block
    end
    return longest_block(alive; min_vars=2)
end

block_summary(block) = block === nothing ?
    (block_start=0, block_end=0, block_vars=0) :
    (block_start=block.a, block_end=block.b, block_vars=length(block.S))

# ----------------------------------------------------------------------------
# Shock-recovery metrics (Sections 7.2, 7.6, 7.7)
# ----------------------------------------------------------------------------

continuation_pairs(alive; times=2:size(alive, 1)) = [
    (j, t) for t in times for j in 1:size(alive, 2) if alive[t, j] && alive[t - 1, j]
]

"""
Shock classification at gamma > 0.5 among continuation pairs (optionally
restricted to `times`): precision, recall, specificity, F1, AUROC, Brier
score, mean Bernoulli entropy of gamma, the fraction of gamma below 0.05 or
above 0.95, and the mean gamma among true shocks and non-shocks.
"""
function shock_metrics_paper(gamma, s_true, alive; times=2:size(alive, 1))
    pairs = continuation_pairs(alive; times=times)
    isempty(pairs) && return (
        n_pairs=0, n_shocks=0, precision=NaN, recall=NaN, specificity=NaN,
        f1=NaN, auroc=NaN, brier=NaN, entropy=NaN, fraction_extreme=NaN,
        gamma_true_shocks=NaN, gamma_non_shocks=NaN
    )
    scores = [clamp(gamma[j, t], 0.0, 1.0) for (j, t) in pairs]
    truth = BitVector([s_true[j, t] for (j, t) in pairs])
    called = scores .> 0.5
    tp = count(called .& truth)
    fp = count(called .& .!truth)
    fn = count(.!called .& truth)
    tn = count(.!called .& .!truth)
    precision = ratio(tp, tp + fp)
    recall = ratio(tp, tp + fn)
    entropy(g) = g <= 0 || g >= 1 ? 0.0 : -(g * log(g) + (1 - g) * log(1 - g))
    return (
        n_pairs=length(pairs),
        n_shocks=count(truth),
        precision=precision,
        recall=recall,
        specificity=ratio(tn, tn + fp),
        f1=isfinite(precision) && isfinite(recall) && precision + recall > 0 ?
            2 * precision * recall / (precision + recall) : NaN,
        auroc=auroc(scores, truth),
        brier=mean((scores .- truth) .^ 2),
        entropy=mean(entropy.(scores)),
        fraction_extreme=mean((scores .< 0.05) .| (scores .> 0.95)),
        gamma_true_shocks=nanmean(scores[truth]),
        gamma_non_shocks=nanmean(scores[.!truth])
    )
end

"""Mean posterior shock probability over continuing variables at each t."""
function mean_gamma_by_time(gamma, alive)
    T_ = size(alive, 1)
    return [
        t == 1 ? NaN : nanmean(gamma[j, t] for j in 1:size(alive, 2) if alive[t, j] && alive[t - 1, j])
        for t in 1:T_
    ]
end

# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------

"""Append rows (NamedTuples with the same keys) to a CSV, writing a header once."""
function append_rows(path, rows)
    isempty(path) && return
    new_file = !isfile(path) || filesize(path) == 0
    open(path, "a") do io
        new_file && println(io, join(string.(keys(first(rows))), ","))
        for row in rows
            println(io, join(format_value.(values(row)), ","))
        end
    end
end

"""Rows of mean gamma by time, for shock-probability plots."""
gamma_time_rows(base, method, gamma, alive) = [
    merge(base, (method=method, t=t, mean_gamma=g))
    for (t, g) in enumerate(mean_gamma_by_time(gamma, alive))
]

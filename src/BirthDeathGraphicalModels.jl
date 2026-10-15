module BirthDeathGraphicalModels

using LinearAlgebra
using Distributions
using Random
using Statistics
using Plots

export TensorPrecisionModel,
       generate_birth_death_times,
       generate_birth_death_data,
       initialize_loading_paths,
       EM_time_varying,
       MCMC_time_varying,
       parallel_tempering_time_varying,
       parallel_tempered_EM,
       active_precision,
       active_precision_errors,
       compare_active_precision_matrices,
       rv_coefficient,
       offdiagonal_rv,
       partial_correlation_rv,
       run_example

# ============================================================================
# Numerical helpers
# ============================================================================

"""Numerically stable inverse-logit."""
function inverse_logit(x::Real)
    if x >= 0
        return 1.0 / (1.0 + exp(-x))
    end
    ex = exp(x)
    return ex / (1.0 + ex)
end

"""Numerically stable log(1 + exp(x))."""
log1pexp_stable(x::Real) = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))

"""Numerically stable log(exp(a) + exp(b))."""
logsumexp2(a::Real, b::Real) =
    a > b ? a + log1p(exp(b - a)) : b + log1p(exp(a - b))

# ============================================================================
# Activity helpers
# ============================================================================

"""Return true when variable j starts a new active spell at time t."""
is_birth(alive::AbstractMatrix{Bool}, t::Int, j::Int) =
    alive[t, j] && (t == 1 || !alive[t - 1, j])

"""Return true when variable j is alive at both t-1 and t."""
is_continuation(alive::AbstractMatrix{Bool}, t::Int, j::Int) =
    t >= 2 && alive[t, j] && alive[t - 1, j]

"""Return the global indices of variables alive at time t."""
active_indices(alive::AbstractMatrix{Bool}, t::Int) =
    findall(view(alive, t, :))

"""Return the last time in the active spell containing (s,j)."""
function active_spell_end(alive::AbstractMatrix{Bool}, s::Int, j::Int)
    alive[s, j] || return s - 1
    T = size(alive, 1)
    last_t = s
    while last_t < T && alive[last_t + 1, j]
        last_t += 1
    end
    return last_t
end

# ============================================================================
# Birth/death simulation
# ============================================================================

"""
Generate one active spell per variable.

`birth_time[j]` is the first active time and `death_time[j]` is the final
active time. A birth always receives a new loading row. During a continuation,
the row changes with a probability set by `shock_rule`:

- `:counts` (default): logistic(a0 + a_event * (births[t] + deaths[t])),
  shared by all continuing variables;
- `:network`: logistic(a0 + a_birth * xb[j,t] + a_death * xd[j,t]), where
  xb and xd are the summed absolute partial correlations of j with the
  variables born and dying at t (see `network_event_covariates_at`;
  `network_link` chooses partial correlations or loading inner products). The
  returned `pi_time` is then a p x T matrix.

Inactive entries of `Y` and `Lambda` are stored as zero and must be interpreted
together with `alive`.
"""
function generate_birth_death_data(
    p::Int,
    k::Int,
    T::Int,
    birth_time::AbstractVector{<:Integer},
    death_time::AbstractVector{<:Integer};
    a0::Float64=-1.0,
    a_event::Float64=0.5,
    initial_sd::Float64=1.0 / sqrt(k),
    change_sd::Float64=0.5 / sqrt(k),
    Delta_value::Float64=0.2,
    Δ_value::Union{Nothing, Float64}=nothing,
    probability_min::Float64=1e-6,
    probability_max::Float64=1.0 - 1e-6,
    shock_rule::Symbol=:counts,
    network_link::Symbol=:partial_correlation,
    a_birth::Float64=a_event,
    a_death::Float64=a_event,
    rng::AbstractRNG=Random.default_rng()
)
    shock_rule in (:counts, :network) ||
        throw(ArgumentError("shock_rule must be :counts or :network"))
    p >= 1 || throw(ArgumentError("p must be positive"))
    k >= 1 || throw(ArgumentError("k must be positive"))
    T >= 1 || throw(ArgumentError("T must be positive"))
    length(birth_time) == p ||
        throw(DimensionMismatch("birth_time must have length p"))
    length(death_time) == p ||
        throw(DimensionMismatch("death_time must have length p"))
    initial_sd > 0 || throw(ArgumentError("initial_sd must be positive"))
    change_sd > 0 || throw(ArgumentError("change_sd must be positive"))
    effective_Delta_value = isnothing(Δ_value) ? Delta_value : Δ_value
    effective_Delta_value > 0 || throw(ArgumentError("Delta_value must be positive"))
    0 < probability_min < probability_max < 1 ||
        throw(ArgumentError("probability bounds must satisfy 0 < min < max < 1"))

    alive = falses(T, p)
    for j in 1:p
        first_t = max(1, Int(birth_time[j]))
        last_t = min(T, Int(death_time[j]))
        if first_t <= last_t
            alive[first_t:last_t, j] .= true
        end
    end

    births, deaths = birth_death_counts(alive)
    event_counts = births .+ deaths

    pi_time = Vector{Float64}(undef, T)
    pi_time[1] = clamp(inverse_logit(a0), probability_min, probability_max)
    for t in 2:T
        eta = a0 + a_event * event_counts[t]
        pi_time[t] = clamp(inverse_logit(eta), probability_min, probability_max)
    end

    Lambda = [zeros(Float64, p, k) for _ in 1:T]
    Delta = fill(Float64(effective_Delta_value), p)
    s_true = falses(p, T)

    if shock_rule == :network
        pi_time = fill(
            clamp(inverse_logit(a0), probability_min, probability_max), p, T
        )
        for t in 1:T
            # Births get new rows, continuing rows keep their previous value;
            # this is the pre-shock state the covariates are computed from.
            for j in 1:p
                alive[t, j] || continue
                if is_birth(alive, t, j)
                    view(Lambda[t], j, :) .= initial_sd .* randn(rng, k)
                    s_true[j, t] = true
                else
                    view(Lambda[t], j, :) .= view(Lambda[t - 1], j, :)
                end
            end
            t == 1 && continue

            birth_x, death_x = network_event_covariates_at(
                Lambda[t - 1], Lambda[t], Delta, alive, t; link=network_link
            )
            for j in 1:p
                is_continuation(alive, t, j) || continue
                pi_time[j, t] = clamp(
                    inverse_logit(a0 + a_birth * birth_x[j] + a_death * death_x[j]),
                    probability_min,
                    probability_max
                )
                if rand(rng) < pi_time[j, t]
                    view(Lambda[t], j, :) .+= change_sd .* randn(rng, k)
                    s_true[j, t] = true
                end
            end
        end
    end

    for t in 1:T
        shock_rule == :counts || break
        for j in 1:p
            alive[t, j] || continue

            if is_birth(alive, t, j)
                view(Lambda[t], j, :) .= initial_sd .* randn(rng, k)
                s_true[j, t] = true
            else
                # A continuing row keeps its previous value unless a shock occurs.
                view(Lambda[t], j, :) .= view(Lambda[t - 1], j, :)
                if rand(rng) < pi_time[t]
                    view(Lambda[t], j, :) .+= change_sd .* randn(rng, k)
                    s_true[j, t] = true
                end
            end
        end
    end

    Y = zeros(Float64, T, p)
    for t in 1:T
        idx = active_indices(alive, t)
        isempty(idx) && continue

        Lambda_obs = Lambda[t][idx, :]
        Omega_t = Symmetric(Lambda_obs * Lambda_obs' + Diagonal(Delta[idx]))
        chol_Omega = cholesky(Omega_t)

        # If Omega = U'U and z ~ N(0,I), then U\z ~ N(0,Omega^{-1}).
        Y[t, idx] .= chol_Omega.U \ randn(rng, length(idx))
    end

    return (
        Y=Y,
        alive=BitMatrix(alive),
        Lambda=Lambda,
        Delta=Delta,
        s_true=s_true,
        births=births,
        deaths=deaths,
        event_counts=event_counts,
        pi_time=pi_time,
        Λ=Lambda,
        Δ=Delta,
        π_time=pi_time
    )
end

"""
Generate nonconstant birth and death times from shifted geometric laws.

`minimum_active` is the smallest number of variables allowed to be active at
any time. A death that would take the active count below it is postponed one
period at a time until a birth makes room. Use `k + 1` to respect assumption
(S4), p_t > k. The default of 1 only rules out time points with no active
variables. The adjustment uses no random numbers, so the other draws are
unchanged, and exit then depends only on past activity, so the mechanism is
still predetermined (S-AR).
"""
function generate_birth_death_times(
    p::Int,
    T::Int;
    birth_probability::Float64=0.10,
    death_probability::Float64=0.05,
    initial_alive_fraction::Float64=0.25,
    minimum_lifetime::Int=5,
    maximum_lifetime::Int=T,
    allow_birth_at_T::Bool=true,
    shuffle_variables::Bool=true,
    minimum_active::Int=1,
    rng::AbstractRNG=Random.default_rng()
)
    p >= 1 || throw(ArgumentError("p must be at least 1"))
    T >= 1 || throw(ArgumentError("T must be at least 1"))
    0 < birth_probability <= 1 ||
        throw(ArgumentError("birth_probability must be in (0,1]"))
    0 < death_probability <= 1 ||
        throw(ArgumentError("death_probability must be in (0,1]"))
    0 <= initial_alive_fraction <= 1 ||
        throw(ArgumentError("initial_alive_fraction must be in [0,1]"))
    minimum_lifetime >= 1 ||
        throw(ArgumentError("minimum_lifetime must be at least 1"))
    maximum_lifetime >= minimum_lifetime ||
        throw(ArgumentError("maximum_lifetime must be at least minimum_lifetime"))
    0 <= minimum_active <= p ||
        throw(ArgumentError("minimum_active must be in 0:p"))

    birth_time = Vector{Int}(undef, p)
    death_time = Vector{Int}(undef, p)

    n_initial = round(Int, initial_alive_fraction * p)
    n_initial >= minimum_active || throw(ArgumentError(
        "initial_alive_fraction * p = $n_initial is below minimum_active = $minimum_active"
    ))
    initial_variables = falses(p)
    if n_initial > 0
        initial_variables[randperm(rng, p)[1:n_initial]] .= true
    end

    latest_birth = allow_birth_at_T ? T : max(T - 1, 1)
    birth_distribution = Geometric(birth_probability)
    lifetime_distribution = Geometric(death_probability)

    for j in 1:p
        birth_time[j] = if initial_variables[j]
            1
        else
            clamp(2 + rand(rng, birth_distribution), 1, latest_birth)
        end

        lifetime = clamp(
            minimum_lifetime + rand(rng, lifetime_distribution),
            minimum_lifetime,
            maximum_lifetime
        )

        death_time[j] = min(birth_time[j] + lifetime - 1, T)
    end

    # Postpone deaths that would leave fewer than minimum_active variables.
    # The count at t-1 satisfies the bound, so a shortfall at t can only come
    # from variables whose last active time is t-1.
    for t in 2:T
        n_active = count(j -> birth_time[j] <= t <= death_time[j], 1:p)
        while n_active < minimum_active
            candidates = [
                j for j in 1:p
                if death_time[j] == t - 1 &&
                   death_time[j] - birth_time[j] + 1 < maximum_lifetime
            ]
            isempty(candidates) && throw(ArgumentError(
                "cannot keep $minimum_active variables active at t=$t " *
                "without exceeding maximum_lifetime"
            ))
            death_time[first(candidates)] = t
            n_active += 1
        end
    end

    if shuffle_variables
        permutation = randperm(rng, p)
        birth_time = birth_time[permutation]
        death_time = death_time[permutation]
    end

    return birth_time, death_time
end

# ============================================================================
# State and difference parameterization
# ============================================================================

mutable struct TensorPrecisionModel
    Lambda::Vector{Matrix{Float64}}
    DiffLambda::Vector{Matrix{Float64}}
    Delta::Vector{Float64}

    y::Matrix{Float64}
    alive::BitMatrix
    obs_idx::Vector{Vector{Int}}
    local_pos::Matrix{Int}

    logdets::Vector{Float64}
    yOmegays::Vector{Float64}
    logliks::Vector{Float64}
    Sigma_diag::Vector{Vector{Float64}}
    SigmaLambda::Vector{Matrix{Float64}}

    iter::Int
end

const _STATE_PROPERTY_ALIASES = Dict{Symbol, Symbol}(
    :Λ => :Lambda,
    :ΔΛ => :DiffLambda,
    :Δ => :Delta,
    :yΩys => :yOmegays,
    :Σdiag => :Sigma_diag,
    :ΣΛ => :SigmaLambda
)

function Base.getproperty(state::TensorPrecisionModel, name::Symbol)
    actual_name = get(_STATE_PROPERTY_ALIASES, name, name)
    return getfield(state, actual_name)
end

function Base.setproperty!(state::TensorPrecisionModel, name::Symbol, value)
    actual_name = get(_STATE_PROPERTY_ALIASES, name, name)
    return setfield!(state, actual_name, value)
end

function Base.propertynames(state::TensorPrecisionModel, private::Bool=false)
    base_names = fieldnames(TensorPrecisionModel)
    alias_names = Tuple(keys(_STATE_PROPERTY_ALIASES))
    return (base_names..., alias_names...)
end

"""
Construct loading differences. Births reset the path; inactive rows are zero.
"""
function make_diff_lambda(
    Lambda::Vector{Matrix{Float64}},
    alive::AbstractMatrix{Bool}
)
    T = length(Lambda)
    p, k = size(Lambda[1])
    size(alive) == (T, p) ||
        throw(DimensionMismatch("alive must have size (T,p)"))

    DiffLambda = [zeros(Float64, p, k) for _ in 1:T]
    for t in 1:T, j in 1:p
        alive[t, j] || continue
        if is_birth(alive, t, j)
            view(DiffLambda[t], j, :) .= view(Lambda[t], j, :)
        else
            view(DiffLambda[t], j, :) .=
                view(Lambda[t], j, :) .- view(Lambda[t - 1], j, :)
        end
    end
    return DiffLambda
end

"""Reconstruct all loading paths from difference loadings."""
function reconstruct_lambda_from_diff!(state::TensorPrecisionModel)
    T = length(state.Lambda)
    p, _ = size(state.Lambda[1])

    for t in 1:T
        fill!(state.Lambda[t], 0.0)
        for j in 1:p
            state.alive[t, j] || continue
            if is_birth(state.alive, t, j)
                view(state.Lambda[t], j, :) .= view(state.DiffLambda[t], j, :)
            else
                view(state.Lambda[t], j, :) .=
                    view(state.Lambda[t - 1], j, :) .+
                    view(state.DiffLambda[t], j, :)
            end
        end
    end
    return state
end

function make_tensor_precision_model(
    Lambda_init::Vector{Matrix{Float64}},
    Delta_init::Vector{Float64},
    y::Matrix{Float64},
    alive_input::AbstractMatrix{Bool}
)
    T, p = size(y)
    size(alive_input) == (T, p) ||
        throw(DimensionMismatch("alive must have the same size as y"))
    length(Lambda_init) == T ||
        throw(DimensionMismatch("Lambda_init must contain T matrices"))
    length(Delta_init) == p ||
        throw(DimensionMismatch("Delta_init must have length p"))
    all(Delta_init .> 0) || throw(ArgumentError("Delta_init must be positive"))

    k = size(Lambda_init[1], 2)
    for t in 1:T
        size(Lambda_init[t]) == (p, k) ||
            throw(DimensionMismatch("Lambda_init[$t] must have size ($p,$k)"))
    end

    alive = BitMatrix(alive_input)
    Lambda = [copy(Lambda_init[t]) for t in 1:T]
    for t in 1:T, j in 1:p
        alive[t, j] || fill!(view(Lambda[t], j, :), 0.0)
    end

    DiffLambda = make_diff_lambda(Lambda, alive)
    Delta = copy(Delta_init)
    obs_idx = [active_indices(alive, t) for t in 1:T]
    local_pos = zeros(Int, T, p)

    for t in 1:T
        for (local_j, global_j) in enumerate(obs_idx[t])
            local_pos[t, global_j] = local_j
        end
    end

    state = TensorPrecisionModel(
        Lambda,
        DiffLambda,
        Delta,
        y,
        alive,
        obs_idx,
        local_pos,
        zeros(Float64, T),
        zeros(Float64, T),
        zeros(Float64, T),
        [Float64[] for _ in 1:T],
        [zeros(Float64, 0, k) for _ in 1:T],
        1
    )

    refresh_all_caches!(state)
    return state
end

# Backward-compatible constructor name.
makeTensorPrecisionModel(args...) = make_tensor_precision_model(args...)

# ============================================================================
# Cache updates
# ============================================================================

function refresh_time_cache!(state::TensorPrecisionModel, t::Int)
    idx = state.obs_idx[t]
    m = length(idx)
    k = size(state.Lambda[t], 2)

    if m == 0
        state.logdets[t] = 0.0
        state.yOmegays[t] = 0.0
        state.logliks[t] = 0.0
        state.Sigma_diag[t] = Float64[]
        state.SigmaLambda[t] = zeros(Float64, 0, k)
        return state
    end

    Lambda_obs = state.Lambda[t][idx, :]
    Delta_obs = state.Delta[idx]
    y_obs = view(state.y, t, idx)

    inv_Delta = 1.0 ./ Delta_obs
    DinvLambda = inv_Delta .* Lambda_obs

    P = Symmetric(Matrix{Float64}(I, k, k) + Lambda_obs' * DinvLambda)
    chol_P = cholesky(P)

    # Sigma*Lambda = D^{-1}Lambda(I + Lambda'D^{-1}Lambda)^{-1}
    SigmaLambda = (chol_P \ DinvLambda')'

    logdet_Omega =
        sum(log, Delta_obs) + 2.0 * sum(log, diag(chol_P.L))

    Lambda_ty = Lambda_obs' * y_obs
    y_Omega_y =
        dot(Delta_obs .* y_obs, y_obs) + dot(Lambda_ty, Lambda_ty)

    correction = vec(sum(SigmaLambda .* DinvLambda, dims=2))
    Sigma_diag = inv_Delta .- correction

    state.logdets[t] = logdet_Omega
    state.yOmegays[t] = y_Omega_y
    state.logliks[t] =
        0.5 * (logdet_Omega - y_Omega_y - m * log(2pi))
    state.SigmaLambda[t] = Matrix(SigmaLambda)
    state.Sigma_diag[t] = Vector{Float64}(Sigma_diag)

    return state
end

function refresh_all_caches!(state::TensorPrecisionModel)
    reconstruct_lambda_from_diff!(state)
    for t in eachindex(state.Lambda)
        refresh_time_cache!(state, t)
    end
    return state
end

# ============================================================================
# Event covariates and shock probabilities
# ============================================================================

function birth_death_counts(alive::AbstractMatrix{Bool})
    T, p = size(alive)
    births = zeros(Int, T)
    deaths = zeros(Int, T)

    for t in 2:T, j in 1:p
        if !alive[t - 1, j] && alive[t, j]
            births[t] += 1
        elseif alive[t - 1, j] && !alive[t, j]
            deaths[t] += 1
        end
    end

    return births, deaths
end

function shock_prior_probabilities(
    alive::AbstractMatrix{Bool};
    a0::Float64=-3.0,
    a_birth::Float64=0.5,
    a_death::Float64=0.5,
    standardize::Bool=false,
    probability_min::Float64=1e-12,
    probability_max::Float64=1.0 - 1e-12
)
    T, _ = size(alive)
    births, deaths = birth_death_counts(alive)
    birth_covariate = Float64.(births)
    death_covariate = Float64.(deaths)

    if standardize && T > 1
        for covariate in (birth_covariate, death_covariate)
            mu = mean(view(covariate, 2:T))
            sigma = std(view(covariate, 2:T); corrected=false)
            if sigma > eps(Float64)
                covariate[2:T] .= (covariate[2:T] .- mu) ./ sigma
            else
                covariate[2:T] .= 0.0
            end
            covariate[1] = 0.0
        end
    end

    coefficients = [a0, a_birth, a_death]
    pi_time = shock_probabilities_from_coefficients(
        coefficients,
        birth_covariate,
        death_covariate;
        probability_min=probability_min,
        probability_max=probability_max
    )

    return (
        pi=pi_time,
        π=pi_time,
        births=births,
        deaths=deaths,
        birth_covariate=birth_covariate,
        death_covariate=death_covariate,
        coefficients=(a0=a0, a_birth=a_birth, a_death=a_death)
    )
end

function shock_probabilities_from_coefficients(
    coefficients::AbstractVector,
    birth_covariate::AbstractVector,
    death_covariate::AbstractVector;
    probability_min::Float64=1e-12,
    probability_max::Float64=1.0 - 1e-12
)
    length(coefficients) == 3 ||
        throw(DimensionMismatch("coefficients must have length 3"))
    T = length(birth_covariate)
    length(death_covariate) == T ||
        throw(DimensionMismatch("birth and death covariates must have equal lengths"))

    pi_time = Vector{Float64}(undef, T)
    for t in 1:T
        eta = coefficients[1] +
              coefficients[2] * birth_covariate[t] +
              coefficients[3] * death_covariate[t]
        pi_time[t] = clamp(
            inverse_logit(eta),
            probability_min,
            probability_max
        )
    end
    return pi_time
end

"""
Per-variable shock probabilities from p x T covariate matrices: the
network-informed parameterization, where each continuing variable has its own
covariates.
"""
function shock_probabilities_from_coefficients(
    coefficients::AbstractVector,
    birth_covariate::AbstractMatrix,
    death_covariate::AbstractMatrix;
    probability_min::Float64=1e-12,
    probability_max::Float64=1.0 - 1e-12
)
    length(coefficients) == 3 ||
        throw(DimensionMismatch("coefficients must have length 3"))
    size(birth_covariate) == size(death_covariate) ||
        throw(DimensionMismatch("birth and death covariates must have equal sizes"))

    eta = coefficients[1] .+
          coefficients[2] .* birth_covariate .+
          coefficients[3] .* death_covariate
    return clamp.(inverse_logit.(eta), probability_min, probability_max)
end

"""
Shock probability of variable j at time t. `pi` is either a length-T vector
shared by all variables or a p x T matrix with one probability per variable.
"""
shock_probability(pi::AbstractVector, j::Int, t::Int) = pi[t]
shock_probability(pi::AbstractMatrix, j::Int, t::Int) = pi[j, t]

covariate_value(x::AbstractVector, j::Int, t::Int) = x[t]
covariate_value(x::AbstractMatrix, j::Int, t::Int) = x[j, t]

function check_time_array_size(x::AbstractVector, p::Int, T::Int, name::String)
    length(x) == T || throw(DimensionMismatch("$name must have length $T"))
end
function check_time_array_size(x::AbstractMatrix, p::Int, T::Int, name::String)
    size(x) == (p, T) || throw(DimensionMismatch("$name must have size ($p,$T)"))
end

"""
Absolute off-diagonal partial correlations of the precision matrix built from
the given loading rows and diagonal precisions on `idx`. Returns a p x p
matrix indexed by global variable, zero outside `idx`.
"""
function abs_partial_correlations(
    Lambda_rows::AbstractMatrix,
    Delta::AbstractVector,
    idx::AbstractVector{Int}
)
    p = size(Lambda_rows, 1)
    R = zeros(Float64, p, p)
    isempty(idx) && return R
    L = Lambda_rows[idx, :]
    Omega = L * L' + Diagonal(Delta[idx])
    d = sqrt.(diag(Omega))
    R[idx, idx] .= abs.(Omega ./ (d * d'))
    for j in idx
        R[j, j] = 0.0
    end
    return R
end

"""
Absolute inner products |<lambda_j, lambda_l>| of the given loading rows on
`idx`, as a p x p matrix indexed by global variable with a zero diagonal.
"""
function abs_loading_inner_products(
    Lambda_rows::AbstractMatrix,
    idx::AbstractVector{Int}
)
    p = size(Lambda_rows, 1)
    R = zeros(Float64, p, p)
    isempty(idx) && return R
    L = Lambda_rows[idx, :]
    R[idx, idx] .= abs.(L * L')
    for j in idx
        R[j, j] = 0.0
    end
    return R
end

"""
Strength of the link between active variables, from loading rows: absolute
partial correlations (`link = :partial_correlation`) or absolute loading inner
products (`link = :loadings`).
"""
function event_link_strengths(Lambda_rows, Delta, idx, link::Symbol)
    link == :partial_correlation &&
        return abs_partial_correlations(Lambda_rows, Delta, idx)
    link == :loadings && return abs_loading_inner_products(Lambda_rows, idx)
    throw(ArgumentError("link must be :partial_correlation or :loadings"))
end

"""
Network-informed event covariates, as p x T matrices. For a variable j that
continues from t-1 to t:

    death covariate = sum over l dying at t of link_{t-1}(j, l),
    birth covariate = sum over l born at t of link~_t(j, l).

Both use j's loading row at t-1, so a variable's own shock at t does not enter
its shock probability. For deaths, l's row is also taken at t-1; for births,
l has no earlier row and its row at birth is used.

With `link = :partial_correlation` (default) the link is |rho_jl|, computed
from the precision on the active set (at t-1 for deaths; for births, the
pre-shock precision at t: continuing rows at their t-1 values, newly born
rows at their time-t values). With `link = :loadings` it is the absolute
inner product |sum_h lambda_{t-1,jh} lambda_{lh}|. Entries for non-continuing
pairs are zero.
"""
function network_event_covariates(
    Lambda::Vector{<:AbstractMatrix},
    Delta::AbstractVector,
    alive::AbstractMatrix{Bool};
    link::Symbol=:partial_correlation
)
    T, p = size(alive)
    birth_covariate = zeros(Float64, p, T)
    death_covariate = zeros(Float64, p, T)
    for t in 2:T
        birth_t, death_t = network_event_covariates_at(
            Lambda[t - 1], Lambda[t], Delta, alive, t; link=link
        )
        birth_covariate[:, t] .= birth_t
        death_covariate[:, t] .= death_t
    end
    return birth_covariate, death_covariate
end

"""
Network-informed event covariates at a single time t >= 2, as length-p
vectors. Only the rows of `Lambda_t` belonging to variables born at t are
used, so `Lambda_t` may hold either pre-shock or post-shock values for the
continuing variables.
"""
function network_event_covariates_at(
    Lambda_previous::AbstractMatrix,
    Lambda_t::AbstractMatrix,
    Delta::AbstractVector,
    alive::AbstractMatrix{Bool},
    t::Int;
    link::Symbol=:partial_correlation
)
    p = size(alive, 2)
    birth_covariate = zeros(Float64, p)
    death_covariate = zeros(Float64, p)

    continuing = [j for j in 1:p if alive[t, j] && alive[t - 1, j]]
    isempty(continuing) && return birth_covariate, death_covariate
    born = [l for l in 1:p if alive[t, l] && !alive[t - 1, l]]
    dying = [l for l in 1:p if alive[t - 1, l] && !alive[t, l]]

    if !isempty(dying)
        R = event_link_strengths(
            Lambda_previous, Delta, findall(view(alive, t - 1, :)), link
        )
        for j in continuing
            death_covariate[j] = sum(R[j, l] for l in dying)
        end
    end

    if !isempty(born)
        pre_shock = copy(Lambda_previous)
        for l in born
            pre_shock[l, :] .= view(Lambda_t, l, :)
        end
        R = event_link_strengths(pre_shock, Delta, findall(view(alive, t, :)), link)
        for j in continuing
            birth_covariate[j] = sum(R[j, l] for l in born)
        end
    end

    return birth_covariate, death_covariate
end

# ============================================================================
# Posterior objectives
# ============================================================================

function log_posterior(
    state::TensorPrecisionModel,
    gamma::AbstractMatrix,
    pi_time::AbstractArray;
    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    a_delta::Float64=1.0,
    b_delta::Float64=1.0,
    likelihood_temperature::Float64=1.0
)
    T = length(state.Lambda)
    p, k = size(state.Lambda[1])
    size(gamma) == (p, T) ||
        throw(DimensionMismatch("gamma must have size ($p,$T)"))
    check_time_array_size(pi_time, p, T, "pi_time")

    inv_initial = 1.0 / sigma_initial^2
    inv_active = 1.0 / sigma_active^2
    inv_inactive = 1.0 / sigma_inactive^2

    # The likelihood is raised to the power likelihood_temperature (1 for the
    # posterior; below 1 for the tempered objectives of tempered EM).
    lp = likelihood_temperature * sum(state.logliks)
    logprior = 0.0

    for t in 1:T, j in 1:p
        state.alive[t, j] || continue

        if is_birth(state.alive, t, j)
            for h in 1:k
                delta = state.DiffLambda[t][j, h]
                value = -0.5 * (
                    delta^2 * inv_initial + log(2pi) + 2log(sigma_initial)
                )
                lp += value
                logprior += value
            end
        else
            gamma_jt = clamp(Float64(gamma[j, t]), 1e-12, 1.0 - 1e-12)
            for h in 1:k
                delta = state.DiffLambda[t][j, h]
                active_value = -0.5 * (
                    delta^2 * inv_active + log(2pi) + 2log(sigma_active)
                )
                inactive_value = -0.5 * (
                    delta^2 * inv_inactive + log(2pi) + 2log(sigma_inactive)
                )
                value = gamma_jt * active_value +
                        (1.0 - gamma_jt) * inactive_value
                lp += value
                logprior += value
            end

            pi_t = clamp(Float64(shock_probability(pi_time, j, t)), 1e-12, 1.0 - 1e-12)
            shock_value =
                gamma_jt * log(pi_t) +
                (1.0 - gamma_jt) * log1p(-pi_t)
            lp += shock_value
            logprior += shock_value
        end
    end

    ever_alive = vec(any(state.alive, dims=1))
    if any(ever_alive)
        delta_value = sum(
            (a_delta - 1.0) .* log.(state.Delta[ever_alive]) .-
            b_delta .* state.Delta[ever_alive]
        )
        lp += delta_value
        logprior += delta_value
    end

    return lp, logprior
end

function shock_coefficient_logprior(
    coefficients::AbstractVector;
    prior_mean::AbstractVector=fill(0.0, 3),
    prior_sd::AbstractVector=fill(2.0, 3)
)
    length(coefficients) == 3 ||
        throw(DimensionMismatch("coefficients must have length 3"))
    length(prior_mean) == 3 ||
        throw(DimensionMismatch("prior_mean must have length 3"))
    length(prior_sd) == 3 ||
        throw(DimensionMismatch("prior_sd must have length 3"))
    all(prior_sd .> 0) ||
        throw(ArgumentError("prior standard deviations must be positive"))

    value = 0.0
    for r in 1:3
        z = (coefficients[r] - prior_mean[r]) / prior_sd[r]
        value += -0.5 * (z^2 + log(2pi) + 2log(prior_sd[r]))
    end
    return value
end

function full_log_posterior(
    state::TensorPrecisionModel,
    gamma::AbstractMatrix,
    pi_time::AbstractArray,
    coefficients::AbstractVector;
    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    a_delta::Float64=1.0,
    b_delta::Float64=1.0,
    coefficient_prior_mean::AbstractVector=fill(0.0, 3),
    coefficient_prior_sd::AbstractVector=fill(2.0, 3),
    likelihood_temperature::Float64=1.0
)
    lp, logprior = log_posterior(
        state,
        gamma,
        pi_time;
        sigma_initial=sigma_initial,
        sigma_active=sigma_active,
        sigma_inactive=sigma_inactive,
        a_delta=a_delta,
        b_delta=b_delta,
        likelihood_temperature=likelihood_temperature
    )

    coefficient_prior = shock_coefficient_logprior(
        coefficients;
        prior_mean=coefficient_prior_mean,
        prior_sd=coefficient_prior_sd
    )

    return lp + coefficient_prior, logprior + coefficient_prior
end

# ============================================================================
# E-step
# ============================================================================

function E_step_independent(
    state::TensorPrecisionModel,
    pi_time::AbstractArray;
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01
)
    T = length(state.DiffLambda)
    p, k = size(state.DiffLambda[1])
    check_time_array_size(pi_time, p, T, "pi_time")

    gamma = zeros(Float64, p, T)
    active_constant = -0.5 * k * log(2pi * sigma_active^2)
    inactive_constant = -0.5 * k * log(2pi * sigma_inactive^2)
    inv_active = 1.0 / sigma_active^2
    inv_inactive = 1.0 / sigma_inactive^2

    for t in 1:T
        for j in 1:p
            pi_t = clamp(Float64(shock_probability(pi_time, j, t)), 1e-12, 1.0 - 1e-12)
            log_pi = log(pi_t)
            log_one_minus_pi = log1p(-pi_t)

            if !state.alive[t, j]
                gamma[j, t] = 0.0
            elseif is_birth(state.alive, t, j)
                gamma[j, t] = 1.0
            else
                ss = sum(abs2, view(state.DiffLambda[t], j, :))
                log_active = log_pi + active_constant - 0.5 * ss * inv_active
                log_inactive =
                    log_one_minus_pi + inactive_constant - 0.5 * ss * inv_inactive
                gamma[j, t] = exp(
                    log_active - logsumexp2(log_active, log_inactive)
                )
            end
        end
    end

    return gamma
end

"""
Regularized E-step.

For alpha=beta=1 this is exactly the ordinary mixture posterior. For
alpha,beta >= 1, the continuation probability maximizes the ordinary
Bernoulli variational objective plus a Beta(alpha,beta) regularizer on gamma.
A safeguarded bisection solves the one-dimensional concave problem.
"""
function E_step_independent_advanced(
    state::TensorPrecisionModel,
    pi_time::AbstractArray;
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    gamma_prior_alpha::Float64=1.0,
    gamma_prior_beta::Float64=1.0,
    gamma_min::Float64=1e-10,
    gamma_max::Float64=1.0 - 1e-10,
    max_iter::Int=60,
    tol::Float64=1e-10
)
    gamma_prior_alpha >= 1.0 ||
        throw(ArgumentError("gamma_prior_alpha must be at least 1"))
    gamma_prior_beta >= 1.0 ||
        throw(ArgumentError("gamma_prior_beta must be at least 1"))
    0.0 < gamma_min < gamma_max < 1.0 ||
        throw(ArgumentError("gamma bounds must satisfy 0 < min < max < 1"))

    T = length(state.DiffLambda)
    p, k = size(state.DiffLambda[1])
    check_time_array_size(pi_time, p, T, "pi_time")

    gamma = zeros(Float64, p, T)
    active_constant = -0.5 * k * log(2pi * sigma_active^2)
    inactive_constant = -0.5 * k * log(2pi * sigma_inactive^2)
    inv_active = 1.0 / sigma_active^2
    inv_inactive = 1.0 / sigma_inactive^2

    alpha_minus_one = gamma_prior_alpha - 1.0
    beta_minus_one = gamma_prior_beta - 1.0

    for t in 1:T
        for j in 1:p
            pi_t = clamp(Float64(shock_probability(pi_time, j, t)), 1e-12, 1.0 - 1e-12)
            log_pi = log(pi_t)
            log_one_minus_pi = log1p(-pi_t)

            if !state.alive[t, j]
                gamma[j, t] = 0.0
                continue
            elseif is_birth(state.alive, t, j)
                gamma[j, t] = 1.0
                continue
            end

            ss = sum(abs2, view(state.DiffLambda[t], j, :))
            log_active = log_pi + active_constant - 0.5 * ss * inv_active
            log_inactive =
                log_one_minus_pi + inactive_constant - 0.5 * ss * inv_inactive

            if alpha_minus_one == 0.0 && beta_minus_one == 0.0
                gamma[j, t] = clamp(
                    exp(log_active - logsumexp2(log_active, log_inactive)),
                    gamma_min,
                    gamma_max
                )
                continue
            end

            log_odds = log_active - log_inactive

            derivative(q) =
                log_odds + log1p(-q) - log(q) +
                alpha_minus_one / q - beta_minus_one / (1.0 - q)

            d_low = derivative(gamma_min)
            d_high = derivative(gamma_max)

            if d_low <= 0.0
                gamma[j, t] = gamma_min
            elseif d_high >= 0.0
                gamma[j, t] = gamma_max
            else
                low = gamma_min
                high = gamma_max
                midpoint = 0.5 * (low + high)

                for _ in 1:max_iter
                    midpoint = 0.5 * (low + high)
                    d_mid = derivative(midpoint)

                    if abs(d_mid) <= tol || high - low <= tol
                        break
                    elseif d_mid > 0.0
                        low = midpoint
                    else
                        high = midpoint
                    end
                end

                gamma[j, t] = midpoint
            end
        end
    end

    return gamma
end

# ============================================================================
# M-step: Delta
# ============================================================================

function update_Delta_bayes!(
    state::TensorPrecisionModel,
    gamma::AbstractMatrix,
    pi_time::AbstractArray;
    a::Float64=1.0,
    b::Float64=1.0,
    step::Float64=0.1,
    Delta_min::Float64=1e-6,
    Delta_max::Float64=1e6,
    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    max_backtracks::Int=15,
    min_step::Float64=1e-6,
    likelihood_temperature::Float64=1.0
)
    T, p = size(state.y)
    Delta_old = copy(state.Delta)
    theta_old = log.(Delta_old)
    direction = zeros(Float64, p)

    old_lp, _ = log_posterior(
        state,
        gamma,
        pi_time;
        sigma_initial=sigma_initial,
        sigma_active=sigma_active,
        sigma_inactive=sigma_inactive,
        a_delta=a,
        b_delta=b,
        likelihood_temperature=likelihood_temperature
    )

    for j in 1:p
        Delta_j = Delta_old[j]
        gradient = 0.0
        information = 0.0
        nobs = 0

        for t in 1:T
            local_j = state.local_pos[t, j]
            local_j == 0 && continue

            sigma_jj = state.Sigma_diag[t][local_j]
            y_jt = state.y[t, j]
            score = Delta_j * (sigma_jj - y_jt^2)

            gradient += 0.5 * likelihood_temperature * score
            information +=
                likelihood_temperature *
                (0.5 * Delta_j^2 * sigma_jj^2 + 0.5 * abs(score))
            nobs += 1
        end

        nobs == 0 && continue
        gradient += (a - 1.0) - b * Delta_j
        information += b * Delta_j
        direction[j] = gradient / max(information, 1e-10)
    end

    direction .= clamp.(direction, -5.0, 5.0)
    local_step = step

    for _ in 1:max_backtracks
        theta_proposed = theta_old .+ local_step .* direction
        state.Delta .= clamp.(exp.(theta_proposed), Delta_min, Delta_max)

        proposal_valid = try
            refresh_all_caches!(state)
            true
        catch
            false
        end

        if proposal_valid
            new_lp, _ = log_posterior(
                state,
                gamma,
                pi_time;
                sigma_initial=sigma_initial,
                sigma_active=sigma_active,
                sigma_inactive=sigma_inactive,
                a_delta=a,
                b_delta=b,
                likelihood_temperature=likelihood_temperature
            )

            if isfinite(new_lp) && new_lp >= old_lp
                return true
            end
        end

        state.Delta .= Delta_old
        refresh_all_caches!(state)
        local_step *= 0.5
        local_step < min_step && break
    end

    state.Delta .= Delta_old
    refresh_all_caches!(state)
    return false
end

# ============================================================================
# M-step: difference loadings
# ============================================================================

function lambda_likelihood_grad_hess(
    state::TensorPrecisionModel,
    t::Int,
    j::Int,
    h::Int
)
    local_j = state.local_pos[t, j]
    local_j == 0 && return 0.0, 0.0

    idx = state.obs_idx[t]
    y_obs = view(state.y, t, idx)
    Lambda_obs = state.Lambda[t][idx, :]

    covariance_term = state.SigmaLambda[t][local_j, h]
    data_term = y_obs[local_j] * dot(y_obs, view(Lambda_obs, :, h))

    gradient = covariance_term - data_term
    curvature = -state.Sigma_diag[t][local_j]
    return gradient, curvature
end

function apply_diff_lambda_increment!(
    state::TensorPrecisionModel,
    s::Int,
    j::Int,
    h::Int,
    increment::Float64
)
    last_t = active_spell_end(state.alive, s, j)
    for t in s:last_t
        state.Lambda[t][j, h] += increment
        refresh_time_cache!(state, t)
    end
    return state
end

function diff_lambda_grad_hess(
    state::TensorPrecisionModel,
    gamma::AbstractMatrix,
    s::Int,
    j::Int,
    h::Int;
    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    likelihood_temperature::Float64=1.0
)
    state.alive[s, j] || return 0.0, -Inf

    gradient = 0.0
    hessian = 0.0
    last_t = active_spell_end(state.alive, s, j)

    for t in s:last_t
        g_lik, h_lik = lambda_likelihood_grad_hess(state, t, j, h)
        gradient += likelihood_temperature * g_lik
        hessian += likelihood_temperature * h_lik
    end

    delta = state.DiffLambda[s][j, h]
    if is_birth(state.alive, s, j)
        precision = 1.0 / sigma_initial^2
    else
        gamma_s = clamp(Float64(gamma[j, s]), 0.0, 1.0)
        precision =
            gamma_s / sigma_active^2 +
            (1.0 - gamma_s) / sigma_inactive^2
    end

    gradient -= precision * delta
    hessian -= precision
    return gradient, hessian
end

function local_diff_lambda_objective(
    state::TensorPrecisionModel,
    gamma::AbstractMatrix,
    s::Int,
    j::Int,
    h::Int;
    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    likelihood_temperature::Float64=1.0
)
    state.alive[s, j] || return -Inf

    last_t = active_spell_end(state.alive, s, j)
    value = likelihood_temperature * sum(view(state.logliks, s:last_t))
    delta = state.DiffLambda[s][j, h]

    if is_birth(state.alive, s, j)
        value += -0.5 * (
            delta^2 / sigma_initial^2 + log(2pi * sigma_initial^2)
        )
    else
        gamma_s = clamp(Float64(gamma[j, s]), 1e-12, 1.0 - 1e-12)
        active_value = -0.5 * (
            delta^2 / sigma_active^2 + log(2pi * sigma_active^2)
        )
        inactive_value = -0.5 * (
            delta^2 / sigma_inactive^2 + log(2pi * sigma_inactive^2)
        )
        value += gamma_s * active_value + (1.0 - gamma_s) * inactive_value
    end

    return value
end

function update_DiffLambda_hess!(
    state::TensorPrecisionModel,
    gamma::AbstractMatrix;
    step::Float64=1.0,
    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    max_backtracks::Int=20,
    min_step::Float64=1e-12,
    max_jump::Float64=1.0,
    use_newton::Bool=true,
    shuffle_coords::Bool=true,
    verbose::Bool=false,
    likelihood_temperature::Float64=1.0,
    rng::AbstractRNG=Random.default_rng()
)
    T = length(state.DiffLambda)
    p, k = size(state.DiffLambda[1])
    size(gamma) == (p, T) ||
        throw(DimensionMismatch("gamma must have size ($p,$T)"))

    coords = [
        (s, j, h)
        for s in 1:T
        for j in 1:p
        for h in 1:k
        if state.alive[s, j]
    ]
    shuffle_coords && Random.shuffle!(rng, coords)

    accepted_count = 0
    rejected_count = 0
    total_gain = 0.0

    for (s, j, h) in coords
        old_value = state.DiffLambda[s][j, h]
        old_objective = local_diff_lambda_objective(
            state,
            gamma,
            s,
            j,
            h;
            sigma_initial=sigma_initial,
            sigma_active=sigma_active,
            sigma_inactive=sigma_inactive,
            likelihood_temperature=likelihood_temperature
        )

        gradient, hessian = diff_lambda_grad_hess(
            state,
            gamma,
            s,
            j,
            h;
            sigma_initial=sigma_initial,
            sigma_active=sigma_active,
            sigma_inactive=sigma_inactive,
            likelihood_temperature=likelihood_temperature
        )

        if !isfinite(gradient) || abs(gradient) < 1e-14
            rejected_count += 1
            continue
        end

        direction = if use_newton && isfinite(hessian) && hessian < -1e-10
            -gradient / hessian
        else
            gradient
        end
        direction = clamp(direction, -max_jump, max_jump)

        if !isfinite(direction) || abs(direction) < 1e-14
            rejected_count += 1
            continue
        end

        local_step = step
        accepted = false

        for _ in 1:max_backtracks
            increment = local_step * direction
            state.DiffLambda[s][j, h] = old_value + increment

            proposal_valid = true
            propagation_completed = false
            new_objective = -Inf

            try
                apply_diff_lambda_increment!(state, s, j, h, increment)
                propagation_completed = true
                new_objective = local_diff_lambda_objective(
                    state,
                    gamma,
                    s,
                    j,
                    h;
                    sigma_initial=sigma_initial,
                    sigma_active=sigma_active,
                    sigma_inactive=sigma_inactive,
                    likelihood_temperature=likelihood_temperature
                )
            catch
                proposal_valid = false
            end

            if proposal_valid &&
               isfinite(new_objective) &&
               new_objective >= old_objective
                accepted = true
                accepted_count += 1
                total_gain += new_objective - old_objective
                break
            end

            state.DiffLambda[s][j, h] = old_value
            if propagation_completed
                try
                    apply_diff_lambda_increment!(state, s, j, h, -increment)
                catch
                    refresh_all_caches!(state)
                end
            else
                refresh_all_caches!(state)
            end

            local_step *= 0.5
            local_step < min_step && break
        end

        accepted || (rejected_count += 1)
    end

    refresh_all_caches!(state)
    state.iter += 1

    if verbose
        println(
            "DiffLambda update | accepted=", accepted_count,
            " rejected=", rejected_count,
            " total_gain=", round(total_gain, digits=6)
        )
    end

    return state
end

# ============================================================================
# M-step: shock coefficients
# ============================================================================

function update_shock_coefficients!(
    coefficients::Vector{Float64},
    gamma::AbstractMatrix,
    alive::AbstractMatrix{Bool},
    birth_covariate::AbstractArray,
    death_covariate::AbstractArray;
    prior_mean::AbstractVector=fill(0.0, 3),
    prior_sd::AbstractVector=fill(2.0, 3),
    maxiter::Int=25,
    tolerance::Float64=1e-8,
    max_backtracks::Int=20,
    min_step::Float64=1e-12,
    verbose::Bool=false
)
    T, p = size(alive)
    size(gamma) == (p, T) ||
        throw(DimensionMismatch("gamma must have size ($p,$T)"))
    length(coefficients) == 3 ||
        throw(DimensionMismatch("coefficients must have length 3"))
    check_time_array_size(birth_covariate, p, T, "birth_covariate")
    check_time_array_size(death_covariate, p, T, "death_covariate")
    length(prior_mean) == 3 ||
        throw(DimensionMismatch("prior_mean must have length 3"))
    length(prior_sd) == 3 ||
        throw(DimensionMismatch("prior_sd must have length 3"))
    all(prior_sd .> 0) ||
        throw(ArgumentError("prior standard deviations must be positive"))

    rows = NTuple{3, Float64}[]
    responses = Float64[]

    for t in 2:T
        for j in 1:p
            if is_continuation(alive, t, j)
                push!(rows, (
                    1.0,
                    Float64(covariate_value(birth_covariate, j, t)),
                    Float64(covariate_value(death_covariate, j, t))
                ))
                push!(responses, clamp(Float64(gamma[j, t]), 0.0, 1.0))
            end
        end
    end

    n = length(responses)
    if n == 0
        return (
            coefficients=copy(coefficients),
            pi_time=shock_probabilities_from_coefficients(
                coefficients,
                birth_covariate,
                death_covariate
            ),
            accepted=false,
            objective=NaN
        )
    end

    X = Matrix{Float64}(undef, n, 3)
    for i in 1:n
        X[i, 1], X[i, 2], X[i, 3] = rows[i]
    end

    y = responses
    prior_precision = 1.0 ./ prior_sd.^2

    function coefficient_objective(a::AbstractVector)
        eta = X * a
        value = sum(y[i] * eta[i] - log1pexp_stable(eta[i]) for i in eachindex(y))
        difference = a .- prior_mean
        return value - 0.5 * sum(prior_precision .* difference.^2)
    end

    a = copy(coefficients)
    current_objective = coefficient_objective(a)
    any_accepted = false

    for iteration in 1:maxiter
        eta = X * a
        pi = inverse_logit.(eta)
        gradient = X' * (y .- pi) .- prior_precision .* (a .- prior_mean)
        weights = pi .* (1.0 .- pi)
        information = X' * (weights .* X) + Diagonal(prior_precision)

        direction = try
            information \ gradient
        catch
            pinv(Matrix(information)) * gradient
        end

        all(isfinite, direction) || break
        norm(direction) <= tolerance * (1.0 + norm(a)) && break

        step_size = 1.0
        accepted = false
        for _ in 1:max_backtracks
            proposal = a .+ step_size .* direction
            proposal_objective = coefficient_objective(proposal)
            if isfinite(proposal_objective) && proposal_objective >= current_objective
                a .= proposal
                current_objective = proposal_objective
                accepted = true
                any_accepted = true
                break
            end
            step_size *= 0.5
            step_size < min_step && break
        end

        if verbose
            println(
                "  coefficient iteration=", iteration,
                " objective=", round(current_objective, digits=8),
                " accepted=", accepted,
                " coefficients=", round.(a, digits=5)
            )
        end

        accepted || break
    end

    coefficients .= a
    return (
        coefficients=copy(coefficients),
        pi_time=shock_probabilities_from_coefficients(
            coefficients,
            birth_covariate,
            death_covariate
        ),
        accepted=any_accepted,
        objective=current_objective
    )
end

# ============================================================================
# Annealing and initialization
# ============================================================================

function annealed_sigma_inactive(
    iter::Int,
    n_warmup::Int,
    sigma_start::Float64,
    sigma_target::Float64
)
    n_warmup >= 1 || throw(ArgumentError("n_warmup must be at least one"))
    sigma_start > 0 || throw(ArgumentError("sigma_start must be positive"))
    sigma_target > 0 || throw(ArgumentError("sigma_target must be positive"))
    n_warmup == 1 && return sigma_target

    fraction = clamp((iter - 1) / (n_warmup - 1), 0.0, 1.0)
    return exp((1.0 - fraction) * log(sigma_start) + fraction * log(sigma_target))
end

"""Initialize a constant loading row within each uninterrupted active spell."""
function initialize_loading_paths(
    alive::AbstractMatrix{Bool},
    k::Int;
    scale::Float64=0.5,
    rng::AbstractRNG=Random.default_rng()
)
    T, p = size(alive)
    k >= 1 || throw(ArgumentError("k must be positive"))
    scale > 0 || throw(ArgumentError("scale must be positive"))

    Lambda_init = [zeros(Float64, p, k) for _ in 1:T]
    for t in 1:T, j in 1:p
        alive[t, j] || continue
        if is_birth(alive, t, j)
            view(Lambda_init[t], j, :) .= scale .* randn(rng, k)
        else
            view(Lambda_init[t], j, :) .= view(Lambda_init[t - 1], j, :)
        end
    end
    return Lambda_init
end

# ============================================================================
# Mask-aware binary-splitting initialization
# ============================================================================

function masked_segment_loglikelihood(
    Y::AbstractMatrix,
    alive::AbstractMatrix{Bool},
    Lambda::AbstractMatrix,
    Delta::AbstractVector
)
    Tseg, p = size(Y)
    size(alive) == (Tseg, p) ||
        throw(DimensionMismatch("alive and Y must have the same size"))
    size(Lambda, 1) == p ||
        throw(DimensionMismatch("Lambda has the wrong number of rows"))
    length(Delta) == p ||
        throw(DimensionMismatch("Delta has the wrong length"))

    value = 0.0
    effective_dimension = 0

    for t in 1:Tseg
        idx = active_indices(alive, t)
        isempty(idx) && continue

        Lambda_obs = Lambda[idx, :]
        y_obs = view(Y, t, idx)
        Omega_t = Symmetric(Lambda_obs * Lambda_obs' + Diagonal(Delta[idx]))
        chol_Omega = cholesky(Omega_t)
        logdet_Omega = 2.0 * sum(log, diag(chol_Omega.L))

        value += 0.5 * (
            logdet_Omega - dot(y_obs, Omega_t * y_obs) - length(idx) * log(2pi)
        )
        effective_dimension += length(idx)
    end

    return value, effective_dimension
end

function fit_masked_segment_simple(
    Y::AbstractMatrix,
    alive::AbstractMatrix{Bool},
    k::Int;
    Lambda_init=nothing,
    Delta_init=nothing,
    prior_var::Float64=10.0,
    a_delta::Float64=1.0,
    b_delta::Float64=1.0,
    learning_rate::Float64=0.01,
    maxiter::Int=500,
    tol::Float64=1e-10,
    Delta_min::Float64=1e-10,
    Delta_max::Float64=1e3,
    rng::AbstractRNG=Random.default_rng()
)
    Tseg, p = size(Y)
    size(alive) == (Tseg, p) ||
        throw(DimensionMismatch("alive and Y must have the same size"))

    ever_alive = vec(any(alive, dims=1))
    Lambda = isnothing(Lambda_init) ?
        0.05 .* randn(rng, p, k) : copy(Lambda_init)
    Delta = isnothing(Delta_init) ? fill(0.1, p) : copy(Delta_init)
    Lambda[.!ever_alive, :] .= 0.0
    Delta .= clamp.(Delta, Delta_min, Delta_max)

    function objective(Lambda_x, Delta_x)
        ll, _ = masked_segment_loglikelihood(Y, alive, Lambda_x, Delta_x)
        prior_Lambda = -sum(abs2, Lambda_x[ever_alive, :]) / (2.0 * prior_var)
        prior_Delta = sum(
            (a_delta - 1.0) .* log.(Delta_x[ever_alive]) .-
            b_delta .* Delta_x[ever_alive]
        )
        return ll + prior_Lambda + prior_Delta
    end

    history = Float64[]
    current = objective(Lambda, Delta)
    push!(history, current)
    step = learning_rate

    for _ in 1:maxiter
        grad_Lambda = zeros(Float64, p, k)
        grad_theta = zeros(Float64, p)

        for t in 1:Tseg
            idx = active_indices(alive, t)
            isempty(idx) && continue

            Lambda_obs = Lambda[idx, :]
            y_obs = view(Y, t, idx)
            Omega_t = Symmetric(Lambda_obs * Lambda_obs' + Diagonal(Delta[idx]))
            chol_Omega = cholesky(Omega_t)
            Sigma_t = chol_Omega \ Matrix{Float64}(I, length(idx), length(idx))
            G_t = Sigma_t - y_obs * y_obs'

            grad_Lambda[idx, :] = grad_Lambda[idx, :] + G_t * Lambda_obs
            for (r, j) in enumerate(idx)
                grad_theta[j] += 0.5 * Delta[j] * (Sigma_t[r, r] - y_obs[r]^2)
            end
        end

        grad_Lambda[ever_alive, :] .-= Lambda[ever_alive, :] ./ prior_var
        for j in findall(ever_alive)
            grad_theta[j] += (a_delta - 1.0) - b_delta * Delta[j]
        end

        old_Lambda = copy(Lambda)
        old_Delta = copy(Delta)
        local_step = step
        accepted = false

        for _ in 1:20
            Lambda_proposal = old_Lambda .+ local_step .* grad_Lambda
            Lambda_proposal[.!ever_alive, :] .= 0.0
            theta_proposal = log.(old_Delta) .+ local_step .* grad_theta
            Delta_proposal = clamp.(exp.(theta_proposal), Delta_min, Delta_max)

            proposed = try
                objective(Lambda_proposal, Delta_proposal)
            catch
                -Inf
            end

            if isfinite(proposed) && proposed >= current
                Lambda .= Lambda_proposal
                Delta .= Delta_proposal
                current = proposed
                accepted = true
                step = min(local_step * 1.05, learning_rate)
                break
            end
            local_step *= 0.5
        end

        push!(history, current)
        accepted || break

        if length(history) >= 2
            change = abs(history[end] - history[end - 1])
            change <= tol * (1.0 + abs(history[end - 1])) && break
        end
    end

    return Lambda, Delta, history
end

"""
Pseudo-BIC used by the binary-splitting initialization. This is deliberate:
it is not a standard BIC and should not be "corrected".

It returns -(loglik + k(k+1)/2) / (p*k). The penalty is added rather than
subtracted, has no log(n) factor, and the normalizer p*k is the same for every
segment. In `recursive_binary_split_masked` the split score is the average of
the two halves, so the penalty cancels and

    improvement = (-ll_full + (ll_left + ll_right) / 2) / (p*k).

Because the halves are warm-started from the full fit, this is positive
whenever ll_full < 0, so in practice the recursion is stopped by
`min_segment_length` and `max_depth` rather than by the criterion.
"""
function masked_segment_bic(
    Y::AbstractMatrix,
    alive::AbstractMatrix{Bool},
    Lambda::AbstractMatrix,
    Delta::AbstractVector,
    k::Int
)
    ll, effective_dimension = masked_segment_loglikelihood(Y, alive, Lambda, Delta)
    active = vec(any(alive, dims=1))
   # p_active = count(identity, active)

    loading_parameters = max(
        k * k - div(k * (k - 1), 2),
        0
    )
    n_parameters = max(loading_parameters , 1)
   # sample_size = max(effective_dimension, 2)
    sample_size=length( Lambda)
    return -( ll +n_parameters)/(sample_size)
end

function align_segment_loading(
    current::AbstractMatrix,
    previous::AbstractMatrix,
    common_rows::AbstractVector{Bool}
)
    size(current) == size(previous) ||
        throw(DimensionMismatch("loading matrices must have equal sizes"))
    count(identity, common_rows) == 0 && return Matrix(current)

    C = current[common_rows, :]
    P = previous[common_rows, :]
    decomposition = svd(C' * P)
    rotation = decomposition.U * decomposition.Vt
    return Matrix(current * rotation)
end

function recursive_binary_split_masked(
    Y::AbstractMatrix,
    alive::AbstractMatrix{Bool},
    k::Int,
    start_idx::Int,
    end_idx::Int,
    depth::Int;
    min_segment_length::Int=5,
    bic_threshold::Float64=0.0,
    prior_var::Float64=10.0,
    learning_rate::Float64=0.01,
    maxiter::Int=500,
    max_depth::Int=8,
    a_delta::Float64=1.0,
    b_delta::Float64=1.0,
    Lambda_start=nothing,
    Delta_start=nothing,
    verbose::Bool=true,
    rng::AbstractRNG=Random.default_rng()
)
    Y_segment = Y[start_idx:end_idx, :]
    A_segment = alive[start_idx:end_idx, :]
    segment_length = end_idx - start_idx + 1

    Lambda_full, Delta_full, _ = fit_masked_segment_simple(
        Y_segment,
        A_segment,
        k;
        Lambda_init=Lambda_start,
        Delta_init=Delta_start,
        prior_var=prior_var,
        a_delta=a_delta,
        b_delta=b_delta,
        learning_rate=learning_rate,
        maxiter=maxiter,
        rng=rng
    )
    bic_full = masked_segment_bic(Y_segment, A_segment, Lambda_full, Delta_full, k)

    if verbose
        println(
            repeat("  ", depth),
            "segment [", start_idx, ", ", end_idx, "] BIC=",
            round(bic_full, digits=3)
        )
    end

    if segment_length < 2 * min_segment_length || depth >= max_depth
        return [(
            start=start_idx,
            stop=end_idx,
            Lambda=Lambda_full,
            Delta=Delta_full,
            Λ=Lambda_full,
            Δ=Delta_full,
            bic=bic_full
        )]
    end

    midpoint = start_idx + div(segment_length, 2) - 1
    left_range = start_idx:midpoint
    right_range = (midpoint + 1):end_idx

    Lambda_left, Delta_left, _ = fit_masked_segment_simple(
        Y[left_range, :],
        alive[left_range, :],
        k;
        Lambda_init=Lambda_full,
        Delta_init=Delta_full,
        prior_var=prior_var,
        a_delta=a_delta,
        b_delta=b_delta,
        learning_rate=learning_rate,
        maxiter=maxiter,
        rng=rng
    )

    Lambda_right, Delta_right, _ = fit_masked_segment_simple(
        Y[right_range, :],
        alive[right_range, :],
        k;
        Lambda_init=Lambda_full,
        Delta_init=Delta_full,
        prior_var=prior_var,
        a_delta=a_delta,
        b_delta=b_delta,
        learning_rate=learning_rate,
        maxiter=maxiter,
        rng=rng
    )

    bic_left = masked_segment_bic(
        Y[left_range, :], alive[left_range, :], Lambda_left, Delta_left, k
    )
    bic_right = masked_segment_bic(
        Y[right_range, :], alive[right_range, :], Lambda_right, Delta_right, k
    )
    bic_split = (bic_left + bic_right)/2
    improvement = bic_full - bic_split

    if verbose
        println(
            repeat("  ", depth),
            "candidate split at ", midpoint,
            " | full=", round(bic_full, digits=3),
            " split=", round(bic_split, digits=3),
            " improvement=", round(improvement, digits=3)
        )
    end

    if improvement > bic_threshold
        left_segments = recursive_binary_split_masked(
            Y,
            alive,
            k,
            start_idx,
            midpoint,
            depth + 1;
            min_segment_length=min_segment_length,
            bic_threshold=bic_threshold,
            prior_var=prior_var,
            learning_rate=learning_rate,
            maxiter=maxiter,
            max_depth=max_depth,
            a_delta=a_delta,
            b_delta=b_delta,
            Lambda_start=Lambda_left,
            Delta_start=Delta_left,
            verbose=verbose,
            rng=rng
        )
        right_segments = recursive_binary_split_masked(
            Y,
            alive,
            k,
            midpoint + 1,
            end_idx,
            depth + 1;
            min_segment_length=min_segment_length,
            bic_threshold=bic_threshold,
            prior_var=prior_var,
            learning_rate=learning_rate,
            maxiter=maxiter,
            max_depth=max_depth,
            a_delta=a_delta,
            b_delta=b_delta,
            Lambda_start=Lambda_right,
            Delta_start=Delta_right,
            verbose=verbose,
            rng=rng
        )
        return vcat(left_segments, right_segments)
    end

    return [(
            start=start_idx,
            stop=end_idx,
            Lambda=Lambda_full,
            Delta=Delta_full,
            Λ=Lambda_full,
            Δ=Delta_full,
            bic=bic_full
        )]
end

function detect_changepoints_binary_masked(
    Y::AbstractMatrix,
    alive::AbstractMatrix{Bool},
    k::Int;
    min_segment_length::Int=5,
    bic_threshold::Float64=0.0,
    prior_var::Float64=10.0,
    learning_rate::Float64=0.01,
    maxiter::Int=500,
    max_depth::Int=8,
    a_delta::Float64=1.0,
    b_delta::Float64=1.0,
    Lambda_start=nothing,
    Delta_start=nothing,
    verbose::Bool=true,
    rng::AbstractRNG=Random.default_rng()
)
    T, p = size(Y)
    size(alive) == (T, p) ||
        throw(DimensionMismatch("alive and Y must have the same size"))

    segments = recursive_binary_split_masked(
        Y,
        alive,
        k,
        1,
        T,
        0;
        min_segment_length=min_segment_length,
        bic_threshold=bic_threshold,
        prior_var=prior_var,
        learning_rate=learning_rate,
        maxiter=maxiter,
        max_depth=max_depth,
        a_delta=a_delta,
        b_delta=b_delta,
        Lambda_start=Lambda_start,
        Delta_start=Delta_start,
        verbose=verbose,
        rng=rng
    )
    sort!(segments, by=segment -> segment.start)
    return segments
end

function initialization_from_masked_segments(
    segments,
    alive::AbstractMatrix{Bool},
    k::Int;
    Delta_min::Float64=1e-10,
    Delta_max::Float64=1e3,
    Delta_default::Float64=0.1
)
    T, p = size(alive)
    isempty(segments) && throw(ArgumentError("segments cannot be empty"))

    aligned = Matrix{Float64}[]
    for (ell, segment) in enumerate(segments)
        A = copy(segment.Lambda)
        active = BitVector(vec(any(alive[segment.start:segment.stop, :], dims=1)))
        A[.!active, :] .= 0.0

        if ell > 1
            boundary_t = segment.start
            common = BitVector(
                view(alive, boundary_t - 1, :) .& view(alive, boundary_t, :)
            )
            A = align_segment_loading(A, aligned[end], common)
            A[.!active, :] .= 0.0
        end
        push!(aligned, A)
    end

    Lambda_init = [zeros(Float64, p, k) for _ in 1:T]
    Delta_sum = zeros(Float64, p)
    Delta_weight = zeros(Float64, p)

    for (ell, segment) in enumerate(segments)
        A = aligned[ell]
        for t in segment.start:segment.stop, j in 1:p
            alive[t, j] || continue
            view(Lambda_init[t], j, :) .= view(A, j, :)
        end

        for j in 1:p
            nobs = count(identity, view(alive, segment.start:segment.stop, j))
            if nobs > 0
                Delta_sum[j] += nobs * segment.Delta[j]
                Delta_weight[j] += nobs
            end
        end
    end

    Delta_init = fill(Delta_default, p)
    for j in 1:p
        if Delta_weight[j] > 0
            Delta_init[j] = clamp(
                Delta_sum[j] / Delta_weight[j],
                Delta_min,
                Delta_max
            )
        end
    end

    return Lambda_init, Delta_init
end

# ============================================================================
# Main EM/GEM algorithm
# ============================================================================

function _EM_time_varying_core(
    Y::Matrix{Float64},
    alive_input::AbstractMatrix{Bool},
    k::Int;
    Lambda_init::Vector{Matrix{Float64}},

    changepoint_init::Bool=true,
    min_segment_length::Int=6,
    bic_threshold::Float64=0.0,
    changepoint_prior_var::Float64=10.0,
    changepoint_learning_rate::Float64=0.01,
    changepoint_maxiter::Int=500,
    changepoint_max_depth::Int=8,

    gamma_inactive_init::Float64=0.50,
    gamma_active_init::Float64=0.90,
    gamma_floor::Float64=1e-10,
    gamma_prior_alpha::Float64=1.0,
    gamma_prior_beta::Float64=1.0,
    gamma_min::Float64=1e-10,
    gamma_max::Float64=1.0 - 1e-10,

    a0_init::Float64=0.0,
    a_birth_init::Float64=0.0,
    a_death_init::Float64=0.0,
    estimate_shock_coefficients::Bool=true,
    coefficient_prior_mean::Vector{Float64}=[0.0, 0.0, 0.0],
    coefficient_prior_sd::Vector{Float64}=[2.0, 2.0, 2.0],
    coefficient_maxiter::Int=25,
    coefficient_tolerance::Float64=1e-8,
    coefficient_warmup::Int=1,
    coefficient_damping::Float64=0.10,
    coefficient_update_interval::Int=1,
    standardize_event_counts::Bool=false,
    shock_covariates::Symbol=:counts,
    network_link::Symbol=:partial_correlation,

    niter::Int=100,
    step_delta::Float64=0.01,
    step_lambda::Float64=1.0,

    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    sigma_inactive_start::Float64=0.8 * sigma_active,
    annealed_iterations::Int=30,

    a_delta::Float64=1.0,
    b_delta::Float64=1.0,
    Delta_init_value::Float64=0.1,

    tol::Float64=1e-6,
    verbose::Bool=true,
    rng::AbstractRNG=Random.default_rng()
)
    T, p = size(Y)
    size(alive_input) == (T, p) ||
        throw(DimensionMismatch("alive must have the same size as Y"))
    length(Lambda_init) == T ||
        throw(DimensionMismatch("Lambda_init must contain T matrices"))

    for t in 1:T
        size(Lambda_init[t]) == (p, k) ||
            throw(DimensionMismatch("Lambda_init[$t] must have size ($p,$k)"))
    end

    length(coefficient_prior_mean) == 3 ||
        throw(DimensionMismatch("coefficient_prior_mean must have length 3"))
    length(coefficient_prior_sd) == 3 ||
        throw(DimensionMismatch("coefficient_prior_sd must have length 3"))
    all(coefficient_prior_sd .> 0) ||
        throw(ArgumentError("coefficient prior standard deviations must be positive"))

    sigma_initial > 0 || throw(ArgumentError("sigma_initial must be positive"))
    sigma_active > 0 || throw(ArgumentError("sigma_active must be positive"))
    sigma_inactive > 0 || throw(ArgumentError("sigma_inactive must be positive"))
    sigma_inactive_start > 0 ||
        throw(ArgumentError("sigma_inactive_start must be positive"))
    sigma_active > sigma_inactive ||
        @warn "sigma_active should normally be larger than sigma_inactive"

    0 <= gamma_inactive_init <= 1 ||
        throw(ArgumentError("gamma_inactive_init must be in [0,1]"))
    0 <= gamma_active_init <= 1 ||
        throw(ArgumentError("gamma_active_init must be in [0,1]"))
    0 <= gamma_floor < 0.5 ||
        throw(ArgumentError("gamma_floor must be in [0,0.5)"))
    0 < coefficient_damping <= 1 ||
        throw(ArgumentError("coefficient_damping must be in (0,1]"))
    coefficient_warmup >= 0 ||
        throw(ArgumentError("coefficient_warmup cannot be negative"))
    coefficient_update_interval >= 1 ||
        throw(ArgumentError("coefficient_update_interval must be at least one"))
    annealed_iterations >= 1 ||
        throw(ArgumentError("annealed_iterations must be at least one"))
    niter >= 1 || throw(ArgumentError("niter must be at least one"))

    shock_covariates in (:counts, :network) ||
        throw(ArgumentError("shock_covariates must be :counts or :network"))

    alive = BitMatrix(alive_input)

    active_counts = vec(sum(alive, dims=2))
    n_empty = count(==(0), active_counts)
    n_low_rank = count(m -> 0 < m <= k, active_counts)
    n_empty > 0 &&
        @warn "Some time points have no active variables" n_empty
    n_low_rank > 0 &&
        @warn "Some time points have p_t <= k, so the factor term is not low rank there (assumption S4)" n_low_rank k

    function clamp_continuation_gamma!(gamma_matrix::AbstractMatrix)
        for t in 2:T, j in 1:p
            if is_continuation(alive, t, j)
                gamma_matrix[j, t] = clamp(
                    gamma_matrix[j, t],
                    max(gamma_floor, gamma_min),
                    min(1.0 - gamma_floor, gamma_max)
                )
            end
        end
        return gamma_matrix
    end

    function continuation_gamma_summary(gamma_matrix::AbstractMatrix)
        values = Float64[
            gamma_matrix[j, t]
            for t in 2:T
            for j in 1:p
            if is_continuation(alive, t, j)
        ]
        isempty(values) && return (
            mean=NaN,
            maximum=NaN,
            count=0,
            count_above_half=0
        )
        return (
            mean=mean(values),
            maximum=maximum(values),
            count=length(values),
            count_above_half=count(>(0.5), values)
        )
    end

    Lambda = [copy(Lambda_init[t]) for t in 1:T]
    Delta = fill(Float64(Delta_init_value), p)
    segments = nothing

    if changepoint_init
        verbose && println("Running mask-aware binary-splitting initialization")

        # Start the splitting from the supplied initial values: each variable's
        # Lambda_init row at its first active time, and Delta_init_value. The
        # segment fits then refine them under the model's Gamma(a_delta, b_delta)
        # prior.
        Lambda_start = zeros(Float64, p, k)
        for j in 1:p
            first_t = findfirst(view(alive, :, j))
            isnothing(first_t) && continue
            view(Lambda_start, j, :) .= view(Lambda_init[first_t], j, :)
        end
        # Lambda = 0 is a stationary point of the segment objective, so an
        # all-zero start would never move; fall back to the random start.
        all(iszero, Lambda_start) && (Lambda_start = nothing)
        Delta_start = fill(Float64(Delta_init_value), p)

        segments = detect_changepoints_binary_masked(
            Y,
            alive,
            k;
            min_segment_length=min_segment_length,
            bic_threshold=bic_threshold,
            prior_var=changepoint_prior_var,
            learning_rate=changepoint_learning_rate,
            maxiter=changepoint_maxiter,
            max_depth=changepoint_max_depth,
            a_delta=a_delta,
            b_delta=b_delta,
            Lambda_start=Lambda_start,
            Delta_start=Delta_start,
            verbose=verbose,
            rng=rng
        )
        Lambda, Delta = initialization_from_masked_segments(
            segments,
            alive,
            k;
            Delta_default=Float64(Delta_init_value)
        )
    end

    if verbose
        initial_differences = make_diff_lambda(Lambda, alive)
        continuation_norms = Float64[
            norm(view(initial_differences[t], j, :))
            for t in 2:T
            for j in 1:p
            if is_continuation(alive, t, j)
        ]
        println("Number of continuations: ", length(continuation_norms))
        if isempty(continuation_norms)
            println("No continuation differences are available.")
        else
            println("Nonzero differences: ", count(>(1e-10), continuation_norms))
            println("Mean difference norm: ", mean(continuation_norms))
            println("Maximum difference norm: ", maximum(continuation_norms))
        end
    end

    shock_prior = shock_prior_probabilities(
        alive;
        a0=a0_init,
        a_birth=a_birth_init,
        a_death=a_death_init,
        standardize=standardize_event_counts
    )
    birth_covariate = shock_prior.birth_covariate
    death_covariate = shock_prior.death_covariate

    # :counts uses the birth and death counts at t, shared by all variables.
    # :network gives each continuing variable its own covariates, the summed
    # link strengths with the variables born and dying at t, using partial
    # correlations or loading inner products (network_link; see
    # network_event_covariates). They depend on the parameters, so they
    # are recomputed from the current estimate at every iteration and held
    # fixed within it.
    if shock_covariates == :network
        birth_covariate, death_covariate =
            network_event_covariates(Lambda, Delta, alive; link=network_link)
    end

    shock_coefficients = [a0_init, a_birth_init, a_death_init]
    pi_time = shock_probabilities_from_coefficients(
        shock_coefficients,
        birth_covariate,
        death_covariate
    )

    gamma = zeros(Float64, p, T)
    for t in 1:T, j in 1:p
        if is_birth(alive, t, j)
            gamma[j, t] = 1.0
        elseif is_continuation(alive, t, j)
            gamma[j, t] = gamma_inactive_init
        end
    end

    if changepoint_init && segments !== nothing && length(segments) > 1
        for ell in 2:length(segments)
            changepoint = segments[ell].start
            for j in 1:p
                if is_continuation(alive, changepoint, j)
                    gamma[j, changepoint] = gamma_active_init
                end
            end
        end
    end
    clamp_continuation_gamma!(gamma)

    state = make_tensor_precision_model(Lambda, Delta, Y, alive)

    logpost_history = Float64[]
    loglik_history = Float64[]
    logprior_history = Float64[]
    Delta_mean_history = Float64[]
    mstep_gain_history = Float64[]
    delta_accepted_history = Bool[]
    sigma_inactive_history = Float64[]
    coefficient_history = Vector{Vector{Float64}}()
    pi_time_history = Vector{Array{Float64}}()
    coefficient_objective_history = Float64[]
    coefficient_accepted_history = Bool[]
    gamma_mean_history = Float64[]
    gamma_max_history = Float64[]

    push!(coefficient_history, copy(shock_coefficients))
    push!(pi_time_history, copy(pi_time))

    sigma_current = annealed_sigma_inactive(
        1,
        annealed_iterations,
        sigma_inactive_start,
        sigma_inactive
    )

    initial_objective, initial_prior = full_log_posterior(
        state,
        gamma,
        pi_time,
        shock_coefficients;
        sigma_initial=sigma_initial,
        sigma_active=sigma_active,
        sigma_inactive=sigma_current,
        a_delta=a_delta,
        b_delta=b_delta,
        coefficient_prior_mean=coefficient_prior_mean,
        coefficient_prior_sd=coefficient_prior_sd
    )

    push!(logpost_history, initial_objective)
    push!(logprior_history, initial_prior)
    push!(loglik_history, initial_objective - initial_prior)
    push!(Delta_mean_history, mean(state.Delta))
    push!(sigma_inactive_history, sigma_current)

    initial_gamma_summary = continuation_gamma_summary(gamma)
    push!(gamma_mean_history, initial_gamma_summary.mean)
    push!(gamma_max_history, initial_gamma_summary.maximum)

    coefficient_start = max(
        annealed_iterations,
        estimate_shock_coefficients ? coefficient_warmup : 0
    )

    for iter in 1:niter
        sigma_current = annealed_sigma_inactive(
            iter,
            annealed_iterations,
            sigma_inactive_start,
            sigma_inactive
        )

        if shock_covariates == :network && iter > 1
            birth_covariate, death_covariate =
                network_event_covariates(state.Lambda, state.Delta, alive; link=network_link)
            pi_time = shock_probabilities_from_coefficients(
                shock_coefficients,
                birth_covariate,
                death_covariate
            )
        end

        if iter > 1
            gamma = E_step_independent_advanced(
                state,
                pi_time;
                sigma_active=sigma_active,
                sigma_inactive=sigma_current,
                gamma_prior_alpha=gamma_prior_alpha,
                gamma_prior_beta=gamma_prior_beta,
                gamma_min=gamma_min,
                gamma_max=gamma_max
            )
            clamp_continuation_gamma!(gamma)
        end

        objective_before, _ = full_log_posterior(
            state,
            gamma,
            pi_time,
            shock_coefficients;
            sigma_initial=sigma_initial,
            sigma_active=sigma_active,
            sigma_inactive=sigma_current,
            a_delta=a_delta,
            b_delta=b_delta,
            coefficient_prior_mean=coefficient_prior_mean,
            coefficient_prior_sd=coefficient_prior_sd
        )

        accepted_delta = update_Delta_bayes!(
            state,
            gamma,
            pi_time;
            a=a_delta,
            b=b_delta,
            step=step_delta,
            Delta_min=1e-10,
            Delta_max=1e6,
            sigma_initial=sigma_initial,
            sigma_active=sigma_active,
            sigma_inactive=sigma_current
        )

        update_DiffLambda_hess!(
            state,
            gamma;
            step=step_lambda,
            sigma_initial=sigma_initial,
            sigma_active=sigma_active,
            sigma_inactive=sigma_current,
            max_backtracks=20,
            min_step=1e-12,
            max_jump=1.0,
            use_newton=true,
            shuffle_coords=true,
            verbose=verbose,
            rng=rng
        )

        coefficient_accepted = false
        coefficient_objective = NaN

        eligible_iteration = iter - (coefficient_start + 1)
        should_update_coefficients =
            estimate_shock_coefficients &&
            iter > coefficient_start &&
            mod(eligible_iteration, coefficient_update_interval) == 0

        if should_update_coefficients
            old_coefficients = copy(shock_coefficients)
            proposed_coefficients = copy(shock_coefficients)

            coefficient_update = update_shock_coefficients!(
                proposed_coefficients,
                gamma,
                alive,
                birth_covariate,
                death_covariate;
                prior_mean=coefficient_prior_mean,
                prior_sd=coefficient_prior_sd,
                maxiter=coefficient_maxiter,
                tolerance=coefficient_tolerance,
                verbose=false
            )
            coefficient_objective = coefficient_update.objective

            if coefficient_update.accepted
                shock_coefficients .=
                    (1.0 - coefficient_damping) .* old_coefficients .+
                    coefficient_damping .* coefficient_update.coefficients
                pi_time .= shock_probabilities_from_coefficients(
                    shock_coefficients,
                    birth_covariate,
                    death_covariate
                )
                coefficient_accepted = true
            end
        end

        objective_after, prior_after = full_log_posterior(
            state,
            gamma,
            pi_time,
            shock_coefficients;
            sigma_initial=sigma_initial,
            sigma_active=sigma_active,
            sigma_inactive=sigma_current,
            a_delta=a_delta,
            b_delta=b_delta,
            coefficient_prior_mean=coefficient_prior_mean,
            coefficient_prior_sd=coefficient_prior_sd
        )

        likelihood_after = objective_after - prior_after
        mstep_gain = objective_after - objective_before
        gamma_summary = continuation_gamma_summary(gamma)

        push!(logpost_history, objective_after)
        push!(logprior_history, prior_after)
        push!(loglik_history, likelihood_after)
        push!(Delta_mean_history, mean(state.Delta))
        push!(mstep_gain_history, mstep_gain)
        push!(delta_accepted_history, accepted_delta)
        push!(sigma_inactive_history, sigma_current)
        push!(coefficient_history, copy(shock_coefficients))
        push!(pi_time_history, copy(pi_time))
        push!(coefficient_objective_history, coefficient_objective)
        push!(coefficient_accepted_history, coefficient_accepted)
        push!(gamma_mean_history, gamma_summary.mean)
        push!(gamma_max_history, gamma_summary.maximum)

        if verbose
            mean_pi = pi_time isa AbstractMatrix ?
                finite_mean(pi_time[j, t] for t in 2:T for j in 1:p if is_continuation(alive, t, j)) :
                T > 1 ? mean(view(pi_time, 2:T)) : pi_time[1]
            println(
                "iter=", iter,
                " Q=", round(objective_after, digits=6),
                " gain=", round(mstep_gain, digits=8),
                " mean_gamma=", round(gamma_summary.mean, sigdigits=5),
                " max_gamma=", round(gamma_summary.maximum, sigdigits=5),
                " gamma>0.5=", gamma_summary.count_above_half,
                "/", gamma_summary.count,
                " mean_pi=", round(mean_pi, sigdigits=5),
                " a0=", round(shock_coefficients[1], digits=5),
                " sigma_inactive=", round(sigma_current, digits=5),
                " coefficient_update=", coefficient_accepted,
                " Delta_update=", accepted_delta
            )
        end

        if mstep_gain < -tol
            @warn "Conditional M-step objective decreased" iteration=iter gain=mstep_gain
            step_lambda *= 0.5
            step_delta *= 0.5
        end

        convergence_threshold = tol * (1.0 + abs(objective_before))
        warmup_finished = iter > coefficient_start
        if warmup_finished && abs(mstep_gain) <= convergence_threshold
            verbose && println("Converged at iteration ", iter)
            break
        end
    end

    if shock_covariates == :network
        birth_covariate, death_covariate =
            network_event_covariates(state.Lambda, state.Delta, alive; link=network_link)
    end
    pi_time = shock_probabilities_from_coefficients(
        shock_coefficients,
        birth_covariate,
        death_covariate
    )

    gamma_final = E_step_independent_advanced(
        state,
        pi_time;
        sigma_active=sigma_active,
        sigma_inactive=sigma_inactive,
        gamma_prior_alpha=gamma_prior_alpha,
        gamma_prior_beta=gamma_prior_beta,
        gamma_min=gamma_min,
        gamma_max=gamma_max
    )
    clamp_continuation_gamma!(gamma_final)

    final_objective, final_prior = full_log_posterior(
        state,
        gamma_final,
        pi_time,
        shock_coefficients;
        sigma_initial=sigma_initial,
        sigma_active=sigma_active,
        sigma_inactive=sigma_inactive,
        a_delta=a_delta,
        b_delta=b_delta,
        coefficient_prior_mean=coefficient_prior_mean,
        coefficient_prior_sd=coefficient_prior_sd
    )

    return (
        Lambda=state.Lambda,
        DiffLambda=state.DiffLambda,
        Delta=state.Delta,
        gamma=gamma_final,
        pi_time=pi_time,
        births=shock_prior.births,
        deaths=shock_prior.deaths,
        birth_covariate=birth_covariate,
        death_covariate=death_covariate,
        shock_coefficients=(
            a0=shock_coefficients[1],
            a_birth=shock_coefficients[2],
            a_death=shock_coefficients[3]
        ),
        coefficient_vector=copy(shock_coefficients),
        segments=segments,
        alive=state.alive,
        state=state,
        logpost_history=logpost_history,
        ll_history=loglik_history,
        lprior_history=logprior_history,
        Delta_mean_history=Delta_mean_history,
        mstep_gain_history=mstep_gain_history,
        delta_accepted_history=delta_accepted_history,
        sigma_inactive_history=sigma_inactive_history,
        coefficient_history=coefficient_history,
        pi_time_history=pi_time_history,
        coefficient_objective_history=coefficient_objective_history,
        coefficient_accepted_history=coefficient_accepted_history,
        gamma_mean_history=gamma_mean_history,
        gamma_max_history=gamma_max_history,
        final_logposterior=final_objective,
        final_loglikelihood=final_objective - final_prior,
        final_logprior=final_prior
    )
end

# User-facing wrapper. It accepts both the original Unicode keywords and the
# ASCII keywords used internally.
function EM_time_varying(
    Y::Matrix{Float64},
    alive_input::AbstractMatrix{Bool},
    k::Int;
    kwargs...
)
    normalized = Dict{Symbol, Any}(pairs(kwargs))

    aliases = Dict(
        :Λ_init => :Lambda_init,
        :γ_inactive_init => :gamma_inactive_init,
        :γ_active_init => :gamma_active_init,
        :γ_floor => :gamma_floor,
        :γ_prior_alpha => :gamma_prior_alpha,
        :γ_prior_beta => :gamma_prior_beta,
        :γ_min => :gamma_min,
        :γ_max => :gamma_max,
        :σ_initial => :sigma_initial,
        :σ_active => :sigma_active,
        :σ_inactive => :sigma_inactive,
        :σ_inactive_start => :sigma_inactive_start,
        :Δ_init_value => :Delta_init_value
    )

    for (old_name, new_name) in aliases
        if haskey(normalized, old_name)
            haskey(normalized, new_name) &&
                throw(ArgumentError(
                    "Specify only one of $old_name and $new_name"
                ))
            normalized[new_name] = pop!(normalized, old_name)
        end
    end

    result = _EM_time_varying_core(
        Y,
        alive_input,
        k;
        normalized...
    )

    return merge(
        result,
        (
            Λ=result.Lambda,
            ΔΛ=result.DiffLambda,
            Δ=result.Delta,
            γ=result.gamma,
            π_time=result.pi_time,
            Δmean_history=result.Delta_mean_history,
            σ_inactive_history=result.sigma_inactive_history,
            π_time_history=result.pi_time_history
        )
    )
end

# ============================================================================
# Parallel-tempered EM
# ============================================================================

"""One EM replica at likelihood temperature `beta`."""
mutable struct EMReplica
    beta::Float64
    state::TensorPrecisionModel
    coefficients::Vector{Float64}
    pi_time::Vector{Float64}
end

"""
One generalized EM iteration on the tempered objective
beta * loglik + log prior: E-step, Delta step, loading-difference step and a
damped update of the shock coefficients, as in `EM_time_varying`.
"""
function tempered_em_step!(replica::EMReplica, settings, birth_x, death_x, rng)
    gamma = E_step_independent_advanced(
        replica.state, replica.pi_time;
        sigma_active=settings.sigma_active,
        sigma_inactive=settings.sigma_inactive
    )
    update_Delta_bayes!(
        replica.state, gamma, replica.pi_time;
        a=settings.a_delta, b=settings.b_delta, step=1.0,
        Delta_min=1e-10, Delta_max=1e6,
        sigma_initial=settings.sigma_initial,
        sigma_active=settings.sigma_active,
        sigma_inactive=settings.sigma_inactive,
        likelihood_temperature=replica.beta
    )
    update_DiffLambda_hess!(
        replica.state, gamma;
        step=1.0, max_jump=1.0,
        sigma_initial=settings.sigma_initial,
        sigma_active=settings.sigma_active,
        sigma_inactive=settings.sigma_inactive,
        likelihood_temperature=replica.beta,
        rng=rng
    )
    update = update_shock_coefficients!(
        copy(replica.coefficients), gamma, replica.state.alive, birth_x, death_x;
        prior_mean=settings.coefficient_prior_mean,
        prior_sd=settings.coefficient_prior_sd
    )
    if update.accepted
        replica.coefficients .=
            (1 - settings.coefficient_damping) .* replica.coefficients .+
            settings.coefficient_damping .* update.coefficients
        replica.pi_time = shock_probabilities_from_coefficients(
            replica.coefficients, birth_x, death_x
        )
    end
    return replica
end

"""
Shock probabilities from an E-step, and the untempered log posterior used by
`EM_time_varying` for its `final_logposterior`.
"""
function replica_posterior(replica::EMReplica, settings)
    gamma = E_step_independent_advanced(
        replica.state, replica.pi_time;
        sigma_active=settings.sigma_active,
        sigma_inactive=settings.sigma_inactive
    )
    value, _ = full_log_posterior(
        replica.state, gamma, replica.pi_time, replica.coefficients;
        sigma_initial=settings.sigma_initial,
        sigma_active=settings.sigma_active,
        sigma_inactive=settings.sigma_inactive,
        a_delta=settings.a_delta,
        b_delta=settings.b_delta,
        coefficient_prior_mean=settings.coefficient_prior_mean,
        coefficient_prior_sd=settings.coefficient_prior_sd
    )
    return gamma, value
end

"""
Parallel-tempered (replica-exchange) generalized EM.

One replica per likelihood temperature in `betas` runs generalized EM on the
tempered objective beta * loglik + log prior, all starting from the same
values (typically an `EM_time_varying` fit). Every `swap_interval`
iterations, neighbouring replicas propose to exchange states. With
`exchange = :metropolis` (default) the exchange is accepted with probability
min(1, exp((beta_i - beta_j) (loglik_j - loglik_i))), as in parallel
tempering. With `:greedy` it is accepted exactly when that ratio exceeds 1,
that is when it increases the sum of the replicas' tempered log likelihoods;
with `:always` every proposed exchange is made. Low-temperature replicas move between shock
configurations more easily, and exchanges carry their states towards
beta = 1. After `n_tempered` iterations, the replica state with the highest
untempered log posterior is kept and refined by untempered generalized EM
until the gain falls below `tol` or `n_polish` iterations.

This is a search for a higher mode of the posterior; it is not a sampler.
Shock probabilities use the birth and death counts at each time.
"""
function parallel_tempered_EM(
    Y::Matrix{Float64},
    alive_input::AbstractMatrix{Bool},
    k::Int;
    Lambda_init::Vector{Matrix{Float64}},
    Delta_init::Vector{Float64},
    coefficients_init::Vector{Float64},
    betas::Vector{Float64}=collect(0.2 .^ ((0:8) ./ 8)),
    swap_interval::Int=5,
    exchange::Symbol=:metropolis,
    n_tempered::Int=200,
    n_polish::Int=500,
    tol::Float64=1e-6,
    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    a_delta::Float64=1.0,
    b_delta::Float64=1.0,
    coefficient_prior_mean::Vector{Float64}=[0.0, 0.0, 0.0],
    coefficient_prior_sd::Vector{Float64}=[2.0, 2.0, 2.0],
    coefficient_damping::Float64=0.05,
    rng::AbstractRNG=Random.default_rng()
)
    exchange in (:metropolis, :greedy, :always) ||
        throw(ArgumentError("exchange must be :metropolis, :greedy or :always"))
    alive = BitMatrix(alive_input)
    births, deaths = birth_death_counts(alive)
    birth_x = Float64.(births)
    death_x = Float64.(deaths)
    settings = (
        sigma_initial=sigma_initial, sigma_active=sigma_active,
        sigma_inactive=sigma_inactive, a_delta=a_delta, b_delta=b_delta,
        coefficient_prior_mean=coefficient_prior_mean,
        coefficient_prior_sd=coefficient_prior_sd,
        coefficient_damping=coefficient_damping
    )
    replicas = [
        EMReplica(
            beta,
            make_tensor_precision_model(
                [copy(L) for L in Lambda_init], copy(Delta_init), Y, alive
            ),
            copy(coefficients_init),
            shock_probabilities_from_coefficients(coefficients_init, birth_x, death_x)
        )
        for beta in betas
    ]
    n = length(replicas)
    swap_accept = zeros(Int, max(n - 1, 0))
    swap_tries = zeros(Int, max(n - 1, 0))

    for iter in 1:n_tempered
        for replica in replicas
            tempered_em_step!(replica, settings, birth_x, death_x, rng)
        end
        if iter % swap_interval == 0
            first_pair = isodd(div(iter, swap_interval)) ? 1 : 2
            for i in first_pair:2:(n - 1)
                a, b = replicas[i], replicas[i + 1]
                log_ratio = (a.beta - b.beta) *
                    (sum(b.state.logliks) - sum(a.state.logliks))
                swap_tries[i] += 1
                accept = exchange == :metropolis ? log(rand(rng)) < log_ratio :
                         exchange == :greedy ? log_ratio > 0 : true
                if accept
                    a.state, b.state = b.state, a.state
                    a.coefficients, b.coefficients = b.coefficients, a.coefficients
                    a.pi_time, b.pi_time = b.pi_time, a.pi_time
                    swap_accept[i] += 1
                end
            end
        end
    end

    # Keep the replica state with the highest untempered log posterior.
    scores = [replica_posterior(replica, settings)[2] for replica in replicas]
    chosen = argmax(scores)
    best = replicas[chosen]
    best.beta = 1.0

    previous = scores[chosen]
    polish_iterations = 0
    for _ in 1:n_polish
        tempered_em_step!(best, settings, birth_x, death_x, rng)
        polish_iterations += 1
        _, current = replica_posterior(best, settings)
        abs(current - previous) <= tol * (1 + abs(previous)) && break
        previous = current
    end

    gamma, final_logposterior = replica_posterior(best, settings)
    refresh_all_caches!(best.state)
    return (
        Lambda=best.state.Lambda,
        Delta=best.state.Delta,
        gamma=gamma,
        coefficient_vector=copy(best.coefficients),
        final_logposterior=final_logposterior,
        chosen_beta=betas[chosen],
        replica_logposteriors=scores,
        polish_iterations=polish_iterations,
        swap_acceptance=swap_accept ./ max.(swap_tries, 1)
    )
end

# ============================================================================
# MCMC
# ============================================================================

"""Gaussian log likelihood of one time point, by the determinant lemma."""
function time_loglikelihood(
    Lambda_obs::AbstractMatrix,
    Delta_obs::AbstractVector,
    y_obs::AbstractVector
)
    k = size(Lambda_obs, 2)
    DinvLambda = Lambda_obs ./ Delta_obs
    P = Symmetric(Matrix{Float64}(I, k, k) + Lambda_obs' * DinvLambda)
    C = cholesky(P; check=false)
    issuccess(C) || return -Inf
    logdet_Omega = sum(log, Delta_obs) + 2.0 * sum(log, diag(C.U))
    Lty = Lambda_obs' * y_obs
    quad = dot(Delta_obs .* y_obs, y_obs) + dot(Lty, Lty)
    return 0.5 * (logdet_Omega - quad - length(y_obs) * log(2pi))
end

"""
Change in the log likelihood when row j of the loadings moves by `increment`
over the active spell starting at s. Returns the change and the new per-time
log likelihoods.
"""
function row_increment_loglik_change(state, s, j, increment)
    last_t = active_spell_end(state.alive, s, j)
    new_values = Vector{Float64}(undef, last_t - s + 1)
    change = 0.0
    for (i, t) in enumerate(s:last_t)
        idx = state.obs_idx[t]
        Lambda_obs = state.Lambda[t][idx, :]
        view(Lambda_obs, state.local_pos[t, j], :) .+= increment
        new_values[i] = time_loglikelihood(
            Lambda_obs, state.Delta[idx], view(state.y, t, idx)
        )
        change += new_values[i] - state.logliks[t]
    end
    return change, new_values
end

function apply_row_increment!(state, s, j, increment, new_values)
    last_t = active_spell_end(state.alive, s, j)
    view(state.DiffLambda[s], j, :) .+= increment
    for (i, t) in enumerate(s:last_t)
        view(state.Lambda[t], j, :) .+= increment
        state.logliks[t] = new_values[i]
    end
    return state
end

gaussian_row_logdensity(d, sigma) =
    -length(d) * log(sigma) - sum(abs2, d) / (2 * sigma^2)

"""Settings shared by all chains of a sampler run."""
struct MCMCSettings
    sigma_initial::Float64
    sigma_active::Float64
    sigma_inactive::Float64
    a_delta::Float64
    b_delta::Float64
    coefficient_prior_mean::Vector{Float64}
    coefficient_prior_sd::Vector{Float64}
    birth_x::Vector{Float64}
    death_x::Vector{Float64}
    rows::Vector{Tuple{Int, Int}}
    continuations::Vector{Tuple{Int, Int}}
end

"""
One Markov chain. `beta` is its likelihood temperature: the chain targets
p(y | theta)^beta p(theta), and beta = 1 is the posterior.
"""
mutable struct TemperedChain
    beta::Float64
    state::TensorPrecisionModel
    z::BitMatrix
    coefficients::Vector{Float64}
    pi_time::Vector{Float64}
    rw_scale::Dict{Symbol, Float64}
    rw_accept::Dict{Symbol, Int}
    rw_tries::Dict{Symbol, Int}
    delta_scale::Float64
    delta_accept::Int
    delta_tries::Int
    coefficient_scale::Float64
    coefficient_accept::Int
    coefficient_tries::Int
    flip_accept::Int
    flip_tries::Int
end

function make_tempered_chain(beta, Lambda_init, Delta_init, coefficients_init, Y, alive, settings)
    state = make_tensor_precision_model(
        [copy(L) for L in Lambda_init], copy(Delta_init), Y, alive
    )
    coefficients = copy(coefficients_init)
    pi_time = shock_probabilities_from_coefficients(
        coefficients, settings.birth_x, settings.death_x
    )

    # Start each indicator at its more probable state given the initial D.
    T, p = size(alive)
    z = falses(p, T)
    for (t, j) in settings.continuations
        d = view(state.DiffLambda[t], j, :)
        z[j, t] = log(pi_time[t]) + gaussian_row_logdensity(d, settings.sigma_active) >
                  log1p(-pi_time[t]) + gaussian_row_logdensity(d, settings.sigma_inactive)
    end

    kinds = (:birth, :active, :inactive)
    return TemperedChain(
        beta, state, z, coefficients, pi_time,
        Dict(key => 0.5 for key in kinds),
        Dict(key => 0 for key in kinds),
        Dict(key => 0 for key in kinds),
        0.3, 0, 0, 0.2, 0, 0, 0, 0
    )
end

function coefficient_logtarget(chain, settings, a)
    pis = shock_probabilities_from_coefficients(a, settings.birth_x, settings.death_x)
    value = 0.0
    for (t, j) in settings.continuations
        value += chain.z[j, t] ? log(pis[t]) : log1p(-pis[t])
    end
    return value -
        sum(abs2, (a .- settings.coefficient_prior_mean) ./ settings.coefficient_prior_sd) / 2
end

"""One Metropolis-within-Gibbs sweep at the chain's temperature."""
function mcmc_sweep!(chain::TemperedChain, settings::MCMCSettings, rng::AbstractRNG)
    state = chain.state
    alive = state.alive
    T, p = size(alive)
    k = size(state.Lambda[1], 2)
    beta = chain.beta
    sigma_of(z) = z ? settings.sigma_active : settings.sigma_inactive

    for (t, j) in settings.rows
        d = Vector(view(state.DiffLambda[t], j, :))

        if is_birth(alive, t, j)
            kind = :birth
            sigma = settings.sigma_initial
        else
            pi_t = chain.pi_time[t]
            # Exact conditional of the indicator given D.
            log_active = log(pi_t) + gaussian_row_logdensity(d, settings.sigma_active)
            log_inactive = log1p(-pi_t) + gaussian_row_logdensity(d, settings.sigma_inactive)
            probability = exp(log_active - logsumexp2(log_active, log_inactive))
            chain.z[j, t] = rand(rng) < probability

            # Flip move: switch component, draw D from its prior.
            proposed_z = !chain.z[j, t]
            proposed_d = sigma_of(proposed_z) .* randn(rng, k)
            change, new_values = row_increment_loglik_change(state, t, j, proposed_d .- d)
            log_ratio = beta * change + (proposed_z ?
                log(pi_t) - log1p(-pi_t) :
                log1p(-pi_t) - log(pi_t))
            chain.flip_tries += 1
            if log(rand(rng)) < log_ratio
                apply_row_increment!(state, t, j, proposed_d .- d, new_values)
                chain.z[j, t] = proposed_z
                d = proposed_d
                chain.flip_accept += 1
            end

            kind = chain.z[j, t] ? :active : :inactive
            sigma = sigma_of(chain.z[j, t])
        end

        # Random walk within the current component.
        step = chain.rw_scale[kind] * sigma .* randn(rng, k)
        change, new_values = row_increment_loglik_change(state, t, j, step)
        log_ratio = beta * change +
            gaussian_row_logdensity(d .+ step, sigma) -
            gaussian_row_logdensity(d, sigma)
        chain.rw_tries[kind] += 1
        if log(rand(rng)) < log_ratio
            apply_row_increment!(state, t, j, step, new_values)
            chain.rw_accept[kind] += 1
        end
    end

    # Diagonal precisions: random walk on log Delta_j.
    for j in 1:p
        times = [t for t in 1:T if alive[t, j]]
        isempty(times) && continue
        old = state.Delta[j]
        proposal = old * exp(chain.delta_scale * randn(rng))
        state.Delta[j] = proposal
        new_values = [
            time_loglikelihood(
                state.Lambda[t][state.obs_idx[t], :],
                state.Delta[state.obs_idx[t]],
                view(state.y, t, state.obs_idx[t])
            ) for t in times
        ]
        state.Delta[j] = old
        change = sum(new_values) - sum(state.logliks[times])
        # Gamma prior on Delta plus the Jacobian of the log scale.
        log_ratio = beta * change +
            settings.a_delta * (log(proposal) - log(old)) -
            settings.b_delta * (proposal - old)
        chain.delta_tries += 1
        if log(rand(rng)) < log_ratio
            state.Delta[j] = proposal
            state.logliks[times] .= new_values
            chain.delta_accept += 1
        end
    end

    # Shock-probability coefficients; their conditional has no likelihood term.
    proposal = chain.coefficients .+ chain.coefficient_scale .* randn(rng, 3)
    chain.coefficient_tries += 1
    if log(rand(rng)) <
       coefficient_logtarget(chain, settings, proposal) -
       coefficient_logtarget(chain, settings, chain.coefficients)
        chain.coefficients .= proposal
        chain.pi_time = shock_probabilities_from_coefficients(
            chain.coefficients, settings.birth_x, settings.death_x
        )
        chain.coefficient_accept += 1
    end
    return chain
end

"""Adapt random-walk scales towards acceptance rates near 0.3 (burn-in only)."""
function adapt_chain!(chain::TemperedChain)
    for key in keys(chain.rw_scale)
        chain.rw_tries[key] > 0 || continue
        rate = chain.rw_accept[key] / chain.rw_tries[key]
        chain.rw_scale[key] *= exp(rate - 0.3)
        chain.rw_accept[key] = 0
        chain.rw_tries[key] = 0
    end
    chain.delta_scale *= exp(chain.delta_accept / max(chain.delta_tries, 1) - 0.3)
    chain.coefficient_scale *=
        exp(chain.coefficient_accept / max(chain.coefficient_tries, 1) - 0.3)
    chain.delta_accept = chain.delta_tries = 0
    chain.coefficient_accept = chain.coefficient_tries = 0
    chain.flip_accept = chain.flip_tries = 0
    return chain
end

"""Exchange the model states of two chains; temperatures and scales stay."""
function swap_chain_states!(a::TemperedChain, b::TemperedChain)
    a.state, b.state = b.state, a.state
    a.z, b.z = b.z, a.z
    a.coefficients, b.coefficients = b.coefficients, a.coefficients
    a.pi_time, b.pi_time = b.pi_time, a.pi_time
    return nothing
end

"""
Parallel-tempering sampler for the birth-death factor model.

Runs one chain per likelihood temperature in `betas` (default: nine
temperatures, geometric from 1 to 0.2; neighbouring temperatures must be close
enough for exchanges to be accepted, and a ladder of 1, 0.6, 0.35, 0.2 gave
exchange rates near zero in Experiment II); the chain at
`betas[1]` (normally 1.0) samples the posterior targeted by `EM_time_varying`
with counts-based shock probabilities. Each sweep updates every chain with
`mcmc_sweep!`, then proposes exchanges between neighbouring temperatures,
alternating even and odd pairs, accepted with probability
min(1, exp((beta_i - beta_j) (loglik_j - loglik_i))). Hotter chains move
between shock configurations more easily, and the exchanges carry those moves
to the posterior chain.

The sweep updates, for every active (t, j): z_jt from its exact conditional
given D_jt (continuations); (z_jt, D_jt) by a flip move that switches the
component and draws D_jt from that component's prior, so the acceptance ratio
is the tempered likelihood ratio times the prior odds of z; and D_jt by a
Gaussian random walk scaled to the current component. Then each Delta_j by a
random walk on log Delta_j, and the coefficients by a Gaussian random walk.
Random-walk scales adapt during burn-in only.

Returns posterior means from the first chain (`gamma`, with births set to 1
as in `EM_time_varying`; `Omega_mean`, full p x p per time and meaningful on
the active set; `Delta_mean`; `coefficient_mean`), the mean shock indicators
at every temperature (`gamma_by_beta`), acceptance rates and exchange rates.
"""
function parallel_tempering_time_varying(
    Y::Matrix{Float64},
    alive_input::AbstractMatrix{Bool},
    k::Int;
    Lambda_init::Vector{Matrix{Float64}},
    Delta_init::Vector{Float64},
    coefficients_init::Vector{Float64}=[-2.0, 0.0, 0.0],
    betas::Vector{Float64}=collect(0.2 .^ ((0:8) ./ 8)),
    sigma_initial::Float64=1.0,
    sigma_active::Float64=1.0,
    sigma_inactive::Float64=0.01,
    a_delta::Float64=1.0,
    b_delta::Float64=1.0,
    coefficient_prior_mean::Vector{Float64}=[0.0, 0.0, 0.0],
    coefficient_prior_sd::Vector{Float64}=[2.0, 2.0, 2.0],
    n_burnin::Int=1000,
    n_samples::Int=2000,
    rng::AbstractRNG=Random.default_rng()
)
    alive = BitMatrix(alive_input)
    T, p = size(Y)
    births, deaths = birth_death_counts(alive)
    settings = MCMCSettings(
        sigma_initial, sigma_active, sigma_inactive, a_delta, b_delta,
        coefficient_prior_mean, coefficient_prior_sd,
        Float64.(births), Float64.(deaths),
        [(t, j) for t in 1:T for j in 1:p if alive[t, j]],
        [(t, j) for t in 2:T for j in 1:p if is_continuation(alive, t, j)]
    )
    chains = [
        make_tempered_chain(beta, Lambda_init, Delta_init, coefficients_init, Y, alive, settings)
        for beta in betas
    ]
    n_chains = length(chains)
    swap_accept = zeros(Int, max(n_chains - 1, 0))
    swap_tries = zeros(Int, max(n_chains - 1, 0))

    z_sums = [zeros(Float64, p, T) for _ in 1:n_chains]
    Omega_sum = [zeros(Float64, p, p) for _ in 1:T]
    Delta_sum = zeros(Float64, p)
    coefficient_sum = zeros(Float64, 3)

    for sweep in 1:(n_burnin + n_samples)
        burning = sweep <= n_burnin

        for chain in chains
            mcmc_sweep!(chain, settings, rng)
        end

        for i in (isodd(sweep) ? 1 : 2):2:(n_chains - 1)
            a, b = chains[i], chains[i + 1]
            log_ratio = (a.beta - b.beta) *
                (sum(b.state.logliks) - sum(a.state.logliks))
            swap_tries[i] += 1
            if log(rand(rng)) < log_ratio
                swap_chain_states!(a, b)
                swap_accept[i] += 1
            end
        end

        if burning && sweep % 50 == 0
            foreach(adapt_chain!, chains)
        end

        if !burning
            for (i, chain) in enumerate(chains)
                z_sums[i] .+= chain.z
            end
            cold = chains[1]
            for t in 1:T
                idx = cold.state.obs_idx[t]
                isempty(idx) && continue
                L = cold.state.Lambda[t][idx, :]
                Omega_sum[t][idx, idx] .+= L * L' + Diagonal(cold.state.Delta[idx])
            end
            Delta_sum .+= cold.state.Delta
            coefficient_sum .+= cold.coefficients
        end
    end

    gamma_by_beta = [zs ./ n_samples for zs in z_sums]
    for gamma in gamma_by_beta, t in 1:T, j in 1:p
        is_birth(alive, t, j) && (gamma[j, t] = 1.0)
    end
    cold = chains[1]
    refresh_all_caches!(cold.state)

    return (
        gamma=gamma_by_beta[1],
        gamma_by_beta=gamma_by_beta,
        betas=copy(betas),
        Omega_mean=[O ./ n_samples for O in Omega_sum],
        Delta_mean=Delta_sum ./ n_samples,
        coefficient_mean=coefficient_sum ./ n_samples,
        state=cold.state,
        acceptance=(
            flip=cold.flip_accept / max(cold.flip_tries, 1),
            random_walk=Dict(key => cold.rw_accept[key] / max(cold.rw_tries[key], 1)
                             for key in keys(cold.rw_accept)),
            delta=cold.delta_accept / max(cold.delta_tries, 1),
            coefficients=cold.coefficient_accept / max(cold.coefficient_tries, 1)
        ),
        swap_acceptance=swap_accept ./ max.(swap_tries, 1)
    )
end

"""
Metropolis-within-Gibbs sampler for the birth-death factor model: the
single-chain case of `parallel_tempering_time_varying`, at temperature 1.
"""
MCMC_time_varying(Y::Matrix{Float64}, alive::AbstractMatrix{Bool}, k::Int; kwargs...) =
    parallel_tempering_time_varying(Y, alive, k; kwargs..., betas=[1.0])

# ============================================================================
# Precision diagnostics
# ============================================================================

function active_precision(
    Lambda_t::AbstractMatrix,
    Delta::AbstractVector,
    alive_t::AbstractVector{Bool}
)
    idx = findall(alive_t)
    if isempty(idx)
        return (indices=idx, Omega=zeros(Float64, 0, 0))
    end
    Omega = Lambda_t[idx, :] * Lambda_t[idx, :]' + Diagonal(Delta[idx])
    return (indices=idx, Omega=Matrix(Omega))
end

"""
RV/Frobenius cosine coefficient for two same-sized matrices.
For symmetric cross-product matrices this is the usual RV coefficient.
"""
function rv_coefficient(A::AbstractMatrix, B::AbstractMatrix)
    size(A) == size(B) ||
        throw(DimensionMismatch("A and B must have the same dimensions"))
    denominator = norm(A) * norm(B)
    denominator > eps(Float64) || return NaN
    return dot(A, B) / denominator
end

"""Copy of A with the diagonal set to zero."""
function offdiagonal_part(A::AbstractMatrix)
    B = Matrix{Float64}(A)
    B[diagind(B)] .= 0.0
    return B
end

"""
Off-diagonal partial correlations -Omega_ij / sqrt(Omega_ii Omega_jj), with
the diagonal set to zero.
"""
function offdiagonal_partial_correlation(Omega::AbstractMatrix)
    d = sqrt.(diag(Omega))
    return offdiagonal_part(-Matrix{Float64}(Omega) ./ (d * d'))
end

"""
RV computed on off-diagonal entries only. The full-matrix RV is dominated by
the diagonal, so it can be high when the off-diagonal structure is poorly
recovered. NaN when either matrix has no nonzero off-diagonal entries.
"""
offdiagonal_rv(A::AbstractMatrix, B::AbstractMatrix) =
    rv_coefficient(offdiagonal_part(A), offdiagonal_part(B))

"""Off-diagonal RV on the partial-correlation scale."""
partial_correlation_rv(A::AbstractMatrix, B::AbstractMatrix) =
    rv_coefficient(
        offdiagonal_partial_correlation(A),
        offdiagonal_partial_correlation(B)
    )

function active_precision_errors(
    Lambda_true::Vector{Matrix{Float64}},
    Delta_true::AbstractVector,
    Lambda_est::Vector{Matrix{Float64}},
    Delta_est::AbstractVector,
    alive::AbstractMatrix{Bool},
    times
)
    results = NamedTuple[]

    for t in times
        true_t = active_precision(Lambda_true[t], Delta_true, view(alive, t, :))
        est_t = active_precision(Lambda_est[t], Delta_est, view(alive, t, :))

        if isempty(true_t.indices)
            push!(results, (
                time=t,
                dimension=0,
                frobenius_error=0.0,
                relative_frobenius_error=NaN,
                mean_absolute_error=0.0,
                mean_relative_error=NaN,
                root_mean_squared_error=0.0,
                rv=NaN,
                offdiagonal_rv=NaN,
                partial_correlation_rv=NaN
            ))
            continue
        end

        difference = est_t.Omega - true_t.Omega
        frobenius_error = norm(difference)
        relative_frobenius_error =
            frobenius_error / max(norm(true_t.Omega), eps(Float64))
        mean_absolute_error = mean(abs.(difference))
        mean_relative_error =
            mean_absolute_error / max(mean(abs.(true_t.Omega)), eps(Float64))
        root_mean_squared_error = sqrt(mean(abs2.(difference)))
        rv = rv_coefficient(true_t.Omega, est_t.Omega)
        offdiag_rv = offdiagonal_rv(true_t.Omega, est_t.Omega)
        partial_rv = partial_correlation_rv(true_t.Omega, est_t.Omega)

        push!(results, (
            time=t,
            dimension=length(true_t.indices),
            frobenius_error=frobenius_error,
            relative_frobenius_error=relative_frobenius_error,
            mean_absolute_error=mean_absolute_error,
            mean_relative_error=mean_relative_error,
            root_mean_squared_error=root_mean_squared_error,
            rv=rv,
            offdiagonal_rv=offdiag_rv,
            partial_correlation_rv=partial_rv
        ))
    end

    return results
end

function compare_active_precision_matrices(
    Lambda_true::Vector{Matrix{Float64}},
    Delta_true::AbstractVector,
    Lambda_est::Vector{Matrix{Float64}},
    Delta_est::AbstractVector,
    alive::AbstractMatrix{Bool},
    times;
    color=:viridis
)
    panels = Any[]

    for t in times
        true_t = active_precision(Lambda_true[t], Delta_true, view(alive, t, :))
        est_t = active_precision(Lambda_est[t], Delta_est, view(alive, t, :))

        if isempty(true_t.indices)
            push!(panels, plot(title="No active variables at t=$t", framestyle=:none))
            push!(panels, plot(title="No active variables at t=$t", framestyle=:none))
            continue
        end

        pair_clim = (
            min(minimum(true_t.Omega), minimum(est_t.Omega)),
            max(maximum(true_t.Omega), maximum(est_t.Omega))
        )
        labels = string.(true_t.indices)

        push!(panels, heatmap(
            true_t.Omega;
            title="True Omega (t=$t, p_t=$(length(labels)))",
            clim=pair_clim,
            color=color,
            colorbar=true,
            aspect_ratio=:equal,
            yflip=true,
            xticks=(1:length(labels), labels),
            yticks=(1:length(labels), labels),
            titlefontsize=9,
            tickfontsize=6
        ))

        push!(panels, heatmap(
            est_t.Omega;
            title="Estimated Omega (t=$t)",
            clim=pair_clim,
            color=color,
            colorbar=true,
            aspect_ratio=:equal,
            yflip=true,
            xticks=(1:length(labels), labels),
            yticks=(1:length(labels), labels),
            titlefontsize=9,
            tickfontsize=6
        ))
    end

    return plot(
        panels...;
        layout=(length(times), 2),
        size=(800, 280 * length(times))
    )
end

# ============================================================================
# Reproducible example
# ============================================================================

function finite_mean(values)
    kept = [Float64(x) for x in values if isfinite(x)]
    return isempty(kept) ? NaN : mean(kept)
end

function run_example(; seed::Int=20260722)
    rng = MersenneTwister(seed)
    p, k, T = 20, 2, 50

    birth_time, death_time = generate_birth_death_times(
        p,
        T;
        birth_probability=0.01,
        death_probability=0.01,
        initial_alive_fraction=0.95,
        minimum_lifetime=8,
        maximum_lifetime=50,
        allow_birth_at_T=false,
        shuffle_variables=true,
        rng=rng
    )

    sim = generate_birth_death_data(
        p,
        k,
        T,
        birth_time,
        death_time;
        a0=-2.5,
        a_event=.30,
        initial_sd=0.50 / sqrt(k),
        change_sd=0.5 / sqrt(k),
        Delta_value=0.2,
        rng=rng
    )

    Y = sim.Y
    alive = sim.alive
    Lambda_true = sim.Lambda
    Delta_true = sim.Delta

    Lambda_init = initialize_loading_paths(alive, k; scale=0.5, rng=rng)

    result = EM_time_varying(
        Y,
        alive,
        k;
        Lambda_init=Lambda_init,
        changepoint_init=true,
        min_segment_length=2,
        bic_threshold=0.0,
        changepoint_prior_var=1.0,
        changepoint_learning_rate=0.05,
        changepoint_maxiter=500,
        changepoint_max_depth=8,
        gamma_inactive_init=0.65,
        gamma_active_init=0.9,
        niter=p * 10,
        step_delta=1.0,
        step_lambda=1.0,
        a0_init=4.0,
        a_birth_init=4.0,
        a_death_init=4.0,
        gamma_prior_alpha=1.0,
        gamma_prior_beta=1.0,
        estimate_shock_coefficients=true,
        coefficient_warmup=1,
        coefficient_prior_mean=[0.0, 0.0, 0.0],
        coefficient_prior_sd=fill(0.5 , 3),
        sigma_initial=1.0 /  (k),
        sigma_active=1.0 /  (k),
        sigma_inactive=0.1/(k) ,
        sigma_inactive_start=0.1 /(k),
        annealed_iterations=1,
        Delta_init_value=0.1,
        verbose=true,
        rng=rng
    )

    diagnostics = active_precision_errors(
        Lambda_true,
        Delta_true,
        result.Lambda,
        result.Delta,
        alive,
        1:T
    )

    println("Mean absolute error: ",
        finite_mean(x.mean_absolute_error for x in diagnostics))
    println("Mean relative error: ",
        finite_mean(x.mean_relative_error for x in diagnostics))
    println("Mean relative Frobenius error: ",
        finite_mean(x.relative_frobenius_error for x in diagnostics))
    println("Root mean squared error: ",
        finite_mean(x.root_mean_squared_error for x in diagnostics))
    println("Mean RV coefficient: ",
        finite_mean(x.rv for x in diagnostics))
    println("Mean off-diagonal RV: ",
        finite_mean(x.offdiagonal_rv for x in diagnostics))
    println("Mean partial-correlation RV: ",
        finite_mean(x.partial_correlation_rv for x in diagnostics))

    return (
        simulation=sim,
        result=result,
        diagnostics=diagnostics,
        activity_plot=heatmap(
            alive';
            xlabel="Time",
            ylabel="Variable",
            title="Variable activity",
            colorbar=false
        ),
        true_changepoint_plot=heatmap(
            sim.s_true;
            xlabel="Time",
            ylabel="Variable",
            title="True changepoints"
        ),
        posterior_changepoint_plot=heatmap(
            result.gamma;
            xlabel="Time",
            ylabel="Variable",
            title="Posterior changepoint probabilities"
        ),
        precision_plot=compare_active_precision_matrices(
            Lambda_true,
            Delta_true,
            result.Lambda,
            result.Delta,
            alive,
            [1, 10, 20, 40, 50]
        )
    )
end

end # module

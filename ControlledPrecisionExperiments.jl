
module ControlledPrecisionExperiments

using LinearAlgebra
using Random
using Statistics
using Distributions

export ExperimentConfig,
       PrecisionScenario,
       constant_precision_scenario,
       piecewise_precision_scenario,
       sparse_precision,
       sparse_graph_scenario,
       factor_shock_scenario,
       simulate_from_precision,
       fit_one_experiment,
       run_replicates,
       summarize_replicates,
       offdiag_rv,
       precision_to_partial_correlation

# `rv_coefficient` is deliberately not exported. BirthDeathGraphicalModels
# exports a function of the same name, and scripts load both modules, so
# exporting it here makes every unqualified call ambiguous:
#   UndefVarError: two or more modules export different bindings with this name
# Use ControlledPrecisionExperiments.rv_coefficient when the version defined
# in this module is the one wanted.

"""
Configuration shared by all controlled experiments.

`nrep_time` is the number of independent Gaussian observations generated at
each time. The current BirthDeathGraphicalModels estimator accepts one row per
time, so `nrep_time=1` is directly compatible. Larger values are useful for
oracle/information experiments and require aggregation or an estimator extension.
"""
Base.@kwdef struct ExperimentConfig
    p::Int = 20
    T::Int = 50
    k_fit::Int = 2
    nrep_time::Int = 1
    seed::Int = 1234
    all_alive::Bool = true
end

"""
A data-generating scenario defined directly through precision matrices.

Fields
------
Omega:
    Vector of full p x p positive-definite precision matrices.
alive:
    T x p activity mask.
changepoints:
    Times at which the precision regime changes.
metadata:
    Scenario-specific information.
"""
struct PrecisionScenario
    Omega::Vector{Matrix{Float64}}
    alive::BitMatrix
    changepoints::Vector{Int}
    metadata::NamedTuple
end

function check_spd(A::AbstractMatrix; atol::Float64=1e-10)
    size(A, 1) == size(A, 2) ||
        throw(DimensionMismatch("precision matrix must be square"))
    isapprox(A, A'; atol=atol, rtol=1e-8) ||
        throw(ArgumentError("precision matrix must be symmetric"))
    cholesky(Symmetric(Matrix(A)))
    return true
end

function repeat_precision(Omega::AbstractMatrix, T::Int)
    check_spd(Omega)
    return [Matrix{Float64}(Omega) for _ in 1:T]
end

"""
Constant precision matrix supplied by the user.

This is the first experiment to run. It tests whether the estimator invents
temporal changes when the truth is constant.
"""
function constant_precision_scenario(
    Omega::AbstractMatrix,
    T::Int;
    alive::Union{Nothing,AbstractMatrix{Bool}}=nothing
)
    p = size(Omega, 1)
    A = isnothing(alive) ? trues(T, p) : BitMatrix(alive)
    size(A) == (T, p) ||
        throw(DimensionMismatch("alive must have size (T,p)"))

    return PrecisionScenario(
        repeat_precision(Omega, T),
        A,
        Int[],
        (kind=:constant,)
    )
end

"""
Piecewise-constant precision matrices supplied directly.

Example:
    piecewise_precision_scenario(
        [Omega1, Omega2, Omega3],
        [1:20, 21:35, 36:50]
    )
"""
function piecewise_precision_scenario(
    regime_matrices::Vector{<:AbstractMatrix},
    regimes::Vector{<:UnitRange};
    alive::Union{Nothing,AbstractMatrix{Bool}}=nothing
)
    length(regime_matrices) == length(regimes) ||
        throw(DimensionMismatch("one precision matrix is required per regime"))
    isempty(regimes) && throw(ArgumentError("regimes cannot be empty"))

    p = size(regime_matrices[1], 1)
    T = maximum(last, regimes)
    Omega = Vector{Matrix{Float64}}(undef, T)

    covered = falses(T)
    for (M, r) in zip(regime_matrices, regimes)
        size(M) == (p, p) ||
            throw(DimensionMismatch("all precision matrices must be p x p"))
        check_spd(M)
        for t in r
            covered[t] && throw(ArgumentError("regimes overlap at time $t"))
            Omega[t] = Matrix{Float64}(M)
            covered[t] = true
        end
    end
    all(covered) || throw(ArgumentError("regimes must cover every time 1:T"))

    A = isnothing(alive) ? trues(T, p) : BitMatrix(alive)
    size(A) == (T, p) ||
        throw(DimensionMismatch("alive must have size (T,p)"))

    cps = [first(regimes[i]) for i in 2:length(regimes)]
    return PrecisionScenario(
        Omega,
        A,
        cps,
        (kind=:piecewise, regimes=regimes)
    )
end

"""
Construct a sparse, diagonally dominant precision matrix.

Positive definiteness is guaranteed by strict diagonal dominance.
"""
function sparse_precision(
    p::Int;
    edge_probability::Float64=0.10,
    edge_scale::Float64=0.20,
    diagonal_buffer::Float64=0.30,
    rng::AbstractRNG=Random.default_rng()
)
    0.0 <= edge_probability <= 1.0 ||
        throw(ArgumentError("edge_probability must be in [0,1]"))
    diagonal_buffer > 0 ||
        throw(ArgumentError("diagonal_buffer must be positive"))

    Omega = zeros(Float64, p, p)
    for i in 1:(p - 1), j in (i + 1):p
        if rand(rng) < edge_probability
            value = edge_scale * randn(rng)
            Omega[i, j] = value
            Omega[j, i] = value
        end
    end

    for i in 1:p
        Omega[i, i] = sum(abs, view(Omega, i, :)) + diagonal_buffer
    end

    check_spd(Omega)
    return Omega
end

"""
Sparse-graph scenario generated outside the fitted low-rank model.

`change_fraction` controls the fraction of off-diagonal pairs reconsidered at
each changepoint. The matrix is re-diagonalized after every graph change.
"""
function sparse_graph_scenario(
    p::Int,
    T::Int;
    changepoints::Vector{Int}=[div(T, 2) + 1],
    edge_probability::Float64=0.10,
    edge_scale::Float64=0.20,
    change_fraction::Float64=0.10,
    diagonal_buffer::Float64=0.30,
    rng::AbstractRNG=Random.default_rng()
)
    all(t -> 2 <= t <= T, changepoints) ||
        throw(ArgumentError("changepoints must lie in 2:T"))
    cps = sort(unique(changepoints))

    current = sparse_precision(
        p;
        edge_probability=edge_probability,
        edge_scale=edge_scale,
        diagonal_buffer=diagonal_buffer,
        rng=rng
    )

    Omega = Vector{Matrix{Float64}}(undef, T)
    cp_set = Set(cps)

    for t in 1:T
        if t in cp_set
            for i in 1:(p - 1), j in (i + 1):p
                if rand(rng) < change_fraction
                    value = rand(rng) < edge_probability ?
                        edge_scale * randn(rng) : 0.0
                    current[i, j] = value
                    current[j, i] = value
                end
            end
            # Recompute the diagonal from the off-diagonal entries only.
            # Summing the whole row would add the diagonal left over from the
            # previous regime, so the diagonal would grow at every changepoint
            # and the partial correlations would shrink as an artifact.
            for i in 1:p
                current[i, i] = 0.0
                current[i, i] =
                    sum(abs, view(current, i, :)) + diagonal_buffer
            end
            check_spd(current)
        end
        Omega[t] = copy(current)
    end

    return PrecisionScenario(
        Omega,
        trues(T, p),
        cps,
        (
            kind=:sparse_graph,
            edge_probability=edge_probability,
            change_fraction=change_fraction
        )
    )
end

"""
Correctly specified low-rank-plus-diagonal scenario.

Rows of Lambda change independently at the specified global changepoints.
This scenario returns the actual factor parameters in metadata.
"""
function factor_shock_scenario(
    p::Int,
    k::Int,
    T::Int;
    changepoints::Vector{Int}=[div(T, 2) + 1],
    row_change_probability::Float64=0.30,
    initial_sd::Float64=0.5 / sqrt(k),
    change_sd::Float64=0.5 / sqrt(k),
    delta_value::Float64=0.2,
    rng::AbstractRNG=Random.default_rng()
)
    cps = sort(unique(changepoints))
    cp_set = Set(cps)

    Lambda = [zeros(Float64, p, k) for _ in 1:T]
    Lambda[1] .= initial_sd .* randn(rng, p, k)
    s_true = falses(p, T)
    s_true[:, 1] .= true

    for t in 2:T
        Lambda[t] .= Lambda[t - 1]
        if t in cp_set
            for j in 1:p
                if rand(rng) < row_change_probability
                    Lambda[t][j, :] .+= change_sd .* randn(rng, k)
                    s_true[j, t] = true
                end
            end
        end
    end

    Delta = fill(delta_value, p)
    Omega = [
        Matrix(Lambda[t] * Lambda[t]' + Diagonal(Delta))
        for t in 1:T
    ]

    return PrecisionScenario(
        Omega,
        trues(T, p),
        cps,
        (
            kind=:factor_shock,
            Lambda=Lambda,
            Delta=Delta,
            s_true=s_true
        )
    )
end

"""
Generate Gaussian observations directly from a sequence of precision matrices.

Returns:
- `Y`: T x p when nrep_time == 1;
- `Yrep`: T x nrep_time x p for all settings.

Inactive variables are stored as zero, matching the estimator's convention.
"""
function simulate_from_precision(
    scenario::PrecisionScenario;
    nrep_time::Int=1,
    rng::AbstractRNG=Random.default_rng()
)
    nrep_time >= 1 || throw(ArgumentError("nrep_time must be positive"))

    T = length(scenario.Omega)
    p = size(scenario.Omega[1], 1)
    size(scenario.alive) == (T, p) ||
        throw(DimensionMismatch("scenario activity mask has the wrong size"))

    Yrep = zeros(Float64, T, nrep_time, p)

    for t in 1:T
        idx = findall(view(scenario.alive, t, :))
        isempty(idx) && continue

        Omega_active = scenario.Omega[t][idx, idx]
        check_spd(Omega_active)
        U = cholesky(Symmetric(Omega_active)).U

        for r in 1:nrep_time
            Yrep[t, r, idx] .= U \ randn(rng, length(idx))
        end
    end

    Y = nrep_time == 1 ? Matrix(Yrep[:, 1, :]) : nothing
    return (Y=Y, Yrep=Yrep)
end

function rv_coefficient(A::AbstractMatrix, B::AbstractMatrix)
    size(A) == size(B) ||
        throw(DimensionMismatch("matrices must have equal size"))
    denominator = norm(A) * norm(B)
    denominator > eps(Float64) || return NaN
    return dot(A, B) / denominator
end

function offdiagonal(A::AbstractMatrix)
    B = Matrix(A)
    B[diagind(B)] .= 0.0
    return B
end

offdiag_rv(A::AbstractMatrix, B::AbstractMatrix) =
    rv_coefficient(offdiagonal(A), offdiagonal(B))

function precision_to_partial_correlation(Omega::AbstractMatrix)
    d = sqrt.(diag(Omega))
    P = -Matrix(Omega) ./ (d * d')
    P[diagind(P)] .= 1.0
    return P
end

function compare_precision_sequences(
    Omega_true::Vector{<:AbstractMatrix},
    Omega_est::Vector{<:AbstractMatrix},
    alive::AbstractMatrix{Bool}
)
    T = length(Omega_true)
    length(Omega_est) == T ||
        throw(DimensionMismatch("true and estimated sequences differ in length"))

    rows = NamedTuple[]
    for t in 1:T
        idx = findall(view(alive, t, :))
        isempty(idx) && continue

        A = Matrix(Omega_true[t][idx, idx])
        B = Matrix(Omega_est[t][idx, idx])
        D = B - A

        partial_A = precision_to_partial_correlation(A)
        partial_B = precision_to_partial_correlation(B)

        push!(rows, (
            time=t,
            dimension=length(idx),
            rv=rv_coefficient(A, B),
            offdiag_rv=offdiag_rv(A, B),
            partial_rv=rv_coefficient(partial_A, partial_B),
            relative_frobenius=norm(D) / max(norm(A), eps(Float64)),
            rmse=sqrt(mean(abs2, D))
        ))
    end
    return rows
end

"""
Mean of the finite values of an iterator, NaN when there are none.

Every aggregate below has to tolerate empty or non-finite input. A true
precision matrix with no edges at all — common at small p, for instance
sparse_precision(5; edge_probability=0.10) — makes every off-diagonal RV
NaN, and `mean` over the resulting empty generator throws
"reducing over an empty collection of unknown element type is not allowed".
"""
function safe_mean(values)
    total = 0.0
    count_finite = 0
    for v in values
        x = Float64(v)
        if isfinite(x)
            total += x
            count_finite += 1
        end
    end
    return count_finite == 0 ? NaN : total / count_finite
end

"""
Fit one experiment with BirthDeathGraphicalModels.EM_time_varying.

The BirthDeathGraphicalModels module must already be loaded, for example:

    include("BirthDeathGraphicalModels.jl")
    using .BirthDeathGraphicalModels
    include("ControlledPrecisionExperiments.jl")
    using .ControlledPrecisionExperiments
"""
function fit_one_experiment(
    scenario::PrecisionScenario,
    k_fit::Int;
    seed::Int=1,
    estimator_kwargs...
)
    isdefined(Main, :BirthDeathGraphicalModels) ||
        throw(ArgumentError(
            "Load BirthDeathGraphicalModels before calling fit_one_experiment"
        ))

    BGM = getfield(Main, :BirthDeathGraphicalModels)
    rng = MersenneTwister(seed)
    sim = simulate_from_precision(scenario; nrep_time=1, rng=rng)
    Y = sim.Y
    alive = scenario.alive

    Lambda_init = BGM.initialize_loading_paths(
        alive,
        k_fit;
        scale=0.2,
        rng=rng
    )

    result = BGM.EM_time_varying(
        Y,
        alive,
        k_fit;
        Lambda_init=Lambda_init,
        estimator_kwargs...,
        rng=rng
    )

    Omega_est = [
        Matrix(result.Lambda[t] * result.Lambda[t]' +
               Diagonal(result.Delta))
        for t in eachindex(result.Lambda)
    ]

    diagnostics = compare_precision_sequences(
        scenario.Omega,
        Omega_est,
        alive
    )

    continuation_gamma = [
        result.gamma[j, t]
        for t in 2:size(alive, 1), j in 1:size(alive, 2)
        if alive[t, j] && alive[t - 1, j]
    ]

    return (
        seed=seed,
        scenario=scenario,
        Y=Y,
        fit=result,
        Omega_est=Omega_est,
        diagnostics=diagnostics,
        mean_rv=safe_mean(x.rv for x in diagnostics),
        early_rv=safe_mean(
            x.rv for x in diagnostics if x.time <= min(10, length(diagnostics))
        ),
        mean_offdiag_rv=safe_mean(x.offdiag_rv for x in diagnostics),
        mean_partial_rv=safe_mean(x.partial_rv for x in diagnostics),
        mean_relative_frobenius=safe_mean(x.relative_frobenius for x in diagnostics),
        gamma_mean=isempty(continuation_gamma) ? NaN : mean(continuation_gamma),
        gamma_above_half=isempty(continuation_gamma) ?
            0 : count(>(0.5), continuation_gamma)
    )
end

function run_replicates(
    scenario_builder::Function,
    k_fit::Int;
    seeds=1:20,
    estimator_kwargs...
)
    results = NamedTuple[]
    for seed in seeds
        rng = MersenneTwister(seed)
        scenario = scenario_builder(rng)
        push!(
            results,
            fit_one_experiment(
                scenario,
                k_fit;
                seed=seed,
                estimator_kwargs...
            )
        )
    end
    return results
end

function finite_values(values)
    return [Float64(x) for x in values if isfinite(x)]
end

function mean_sd(values)
    x = finite_values(values)
    isempty(x) && return (mean=NaN, sd=NaN, n=0)
    return (
        mean=mean(x),
        sd=length(x) > 1 ? std(x) : 0.0,
        n=length(x)
    )
end

function summarize_replicates(results)
    return (
        mean_rv=mean_sd(r.mean_rv for r in results),
        early_rv=mean_sd(r.early_rv for r in results),
        mean_offdiag_rv=mean_sd(r.mean_offdiag_rv for r in results),
        mean_partial_rv=mean_sd(r.mean_partial_rv for r in results),
        relative_frobenius=mean_sd(
            r.mean_relative_frobenius for r in results
        ),
        gamma_mean=mean_sd(r.gamma_mean for r in results),
        gamma_above_half=mean_sd(r.gamma_above_half for r in results)
    )
end

end # module


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------
#
# include("BirthDeathGraphicalModels.jl")
# using .BirthDeathGraphicalModels
#
# include("ControlledPrecisionExperiments.jl")
# using .ControlledPrecisionExperiments
#
# p, T, k_fit = 20, 50, 2
#
# # 1. Constant precision outside the true factor model.
# rng = MersenneTwister(10)
# Omega0 = sparse_precision(
#     p;
#     edge_probability=0.10,
#     edge_scale=0.15,
#     rng=rng
# )
# constant_scenario = constant_precision_scenario(Omega0, T)
#
# constant_fit = fit_one_experiment(
#     constant_scenario,
#     k_fit;
#     seed=1,
#     changepoint_init=false,
#     estimate_shock_coefficients=false,
#     a0_init=-3.0,
#     gamma_prior_alpha=1.0,
#     gamma_prior_beta=1.0,
#     sigma_initial=0.5 / sqrt(k_fit),
#     sigma_active=0.5 / sqrt(k_fit),
#     sigma_inactive_start=0.3 / sqrt(k_fit),
#     sigma_inactive=0.1 / sqrt(k_fit),
#     annealed_iterations=50,
#     Delta_init_value=0.2,
#     niter=200,
#     verbose=false
# )
#
# # 2. Piecewise sparse precision outside the true factor model.
# scenario_builder = rng -> sparse_graph_scenario(
#     p,
#     T;
#     changepoints=[21, 36],
#     edge_probability=0.10,
#     change_fraction=0.15,
#     rng=rng
# )
#
# replicates = run_replicates(
#     scenario_builder,
#     k_fit;
#     seeds=1:20,
#     changepoint_init=false,
#     estimate_shock_coefficients=false,
#     a0_init=-2.5,
#     gamma_prior_alpha=1.0,
#     gamma_prior_beta=1.0,
#     sigma_initial=0.5 / sqrt(k_fit),
#     sigma_active=0.5 / sqrt(k_fit),
#     sigma_inactive_start=0.3 / sqrt(k_fit),
#     sigma_inactive=0.1 / sqrt(k_fit),
#     annealed_iterations=50,
#     Delta_init_value=0.2,
#     niter=200,
#     verbose=false
# )
#
# println(summarize_replicates(replicates))

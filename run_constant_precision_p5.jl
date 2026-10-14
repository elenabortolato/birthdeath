#
# Non-time-varying precision matrix, p = 5.
#
# The truth is constant over all T time points. This is the diagnostic
# experiment described in ControlledPrecisionExperiments.constant_precision_scenario:
# it tests whether the estimator invents temporal changes when there are none.
#
# Three constant truths are used:
#   1. constant factor precision, correctly specified (Lambda Lambda' + Delta);
#   2. constant sparse precision, outside the fitted model, sparse graph;
#   3. constant sparse precision, outside the fitted model, denser graph.
#
# Because p = 5 there are only 10 off-diagonal pairs, so edge_probability = 0.10
# leaves most replicates with a single edge. Scenario 3 raises it to 0.30 so the
# graph is non-trivial.
#
# Outputs, written to constant_precision_p5_output/:
#   constant_precision_p5_summary.csv        one row per scenario, over seeds
#   constant_precision_p5_by_time.csv        metrics per time point, per seed
#   constant_precision_p5_stability.csv      how much the estimate moves in time
#   julia_export_constant_p5/                observations + true precisions,
#                                            for the R competitor scripts
#

using LinearAlgebra
using Random
using Statistics
using Printf
using DelimitedFiles

include("BirthDeathGraphicalModels.jl")
using .BirthDeathGraphicalModels

include("ControlledPrecisionExperiments.jl")
using .ControlledPrecisionExperiments

const P = 5
const T = 50
const K_FIT = 2
const SEEDS = 1:10
const OUTPUT_DIRECTORY = joinpath(@__DIR__, "constant_precision_p5_output")
const EXPORT_DIRECTORY = joinpath(OUTPUT_DIRECTORY, "julia_export_constant_p5")

# Set P5_EXPORT_ONLY=1 to regenerate the CSV bundles without refitting. The
# scenarios and the observations are functions of the seed alone, so the data
# written this way is identical to the data the fits saw.
const EXPORT_ONLY = get(ENV, "P5_EXPORT_ONLY", "0") == "1"

# Same estimator settings as run_experiments.jl, rescaled to K_FIT.
const COMMON_FIT = (
    changepoint_init=false,
    estimate_shock_coefficients=false,
    a0_init=-2.5,
    a_birth_init=0.0,
    a_death_init=0.0,
    gamma_inactive_init=0.50,
    gamma_active_init=0.90,
    gamma_prior_alpha=1.0,
    gamma_prior_beta=1.0,
    sigma_initial=0.50 / sqrt(K_FIT),
    sigma_active=0.50 / sqrt(K_FIT),
    sigma_inactive_start=0.30 / sqrt(K_FIT),
    sigma_inactive=0.10 / sqrt(K_FIT),
    annealed_iterations=50,
    Delta_init_value=0.20,
    step_lambda=0.50,
    step_delta=0.10,
    niter=200,
    verbose=false
)

# ---------------------------------------------------------------------------
# Scenario builders. Every one of them is constant in time.
# ---------------------------------------------------------------------------

const SCENARIOS = [
    (
        name="constant factor",
        note="correctly specified, k_true = $(K_FIT)",
        builder=rng -> factor_shock_scenario(
            P, K_FIT, T;
            changepoints=Int[],
            row_change_probability=0.0,
            rng=rng
        )
    ),
    (
        name="constant sparse p=0.10",
        note="outside the fitted model, edge_probability = 0.10",
        builder=rng -> constant_precision_scenario(
            sparse_precision(
                P;
                edge_probability=0.10,
                edge_scale=0.20,
                diagonal_buffer=0.30,
                rng=rng
            ),
            T
        )
    ),
    (
        name="constant sparse p=0.30",
        note="outside the fitted model, edge_probability = 0.30",
        builder=rng -> constant_precision_scenario(
            sparse_precision(
                P;
                edge_probability=0.30,
                edge_scale=0.20,
                diagonal_buffer=0.30,
                rng=rng
            ),
            T
        )
    )
]

# ---------------------------------------------------------------------------
# How much does the estimate move over time when the truth does not?
# ---------------------------------------------------------------------------

function temporal_stability(result)
    Omega = result.Omega_est
    steps = [norm(Omega[t] - Omega[t - 1]) for t in 2:length(Omega)]
    scale = mean(norm(O) for O in Omega)
    gamma = result.fit.gamma
    return (
        mean_step=mean(steps),
        max_step=maximum(steps),
        relative_mean_step=mean(steps) / scale,
        total_drift=norm(Omega[end] - Omega[1]) / scale,
        mean_gamma=mean(gamma[:, 2:end]),
        gamma_above_half=count(>(0.5), gamma[:, 2:end]) / length(gamma[:, 2:end])
    )
end

function scenario_slug(name)
    slug = replace(lowercase(name), r"[^a-z0-9]+" => "_")
    return strip(slug, '_')
end

# Reproduce exactly what fit_one_experiment() feeds the estimator: the scenario
# is built from MersenneTwister(seed), then the observations are drawn from a
# fresh MersenneTwister(seed).
function regenerate_data(builder, seed)
    scenario = builder(MersenneTwister(seed))
    sim = simulate_from_precision(scenario; nrep_time=1, rng=MersenneTwister(seed))
    return (scenario=scenario, Y=sim.Y)
end

function export_scenario(directory, name, builder, seeds)
    mkpath(directory)
    metadata_rows = Vector{Any}()
    observation_rows = Vector{Any}()
    precision_rows = Vector{Any}()

    for (experiment, seed) in enumerate(seeds)
        data = regenerate_data(builder, seed)
        push!(metadata_rows, (experiment, P, T, K_FIT, K_FIT, seed))
        for t in 1:T, j in 1:P
            push!(observation_rows, (experiment, t, j, data.Y[t, j]))
        end
        for t in 1:T, i in 1:P, j in 1:P
            push!(precision_rows, (experiment, t, i, j, data.scenario.Omega[t][i, j]))
        end
    end

    write_table(joinpath(directory, "metadata.csv"), metadata_rows,
                ["experiment", "p", "T", "K_TRUE", "K_FIT", "data_seed"])
    write_table(joinpath(directory, "observations.csv"), observation_rows,
                ["experiment", "time", "variable", "value"])
    write_table(joinpath(directory, "precisions.csv"), precision_rows,
                ["experiment", "time", "row", "col", "value"])
    println("exported \"", name, "\" to ", directory)
end

function mean_sd_pair(values)
    x = [Float64(v) for v in values if isfinite(v)]
    isempty(x) && return (NaN, NaN)
    return (mean(x), length(x) > 1 ? std(x) : 0.0)
end

# ---------------------------------------------------------------------------

function main()
    mkpath(EXPORT_DIRECTORY)

    # Every scenario gets its own bundle. The p=0.10 truth is often edgeless at
    # p = 5, which leaves off-diagonal RV undefined, so the R comparison should
    # use the p=0.30 bundle for off-diagonal metrics.
    for scenario in SCENARIOS
        export_scenario(
            joinpath(EXPORT_DIRECTORY, scenario_slug(scenario.name)),
            scenario.name,
            scenario.builder,
            SEEDS
        )
    end

    if EXPORT_ONLY
        println("\nP5_EXPORT_ONLY=1: bundles regenerated, no fitting performed.")
        return nothing
    end

    summary_rows = Vector{Any}()
    time_rows = Vector{Any}()
    stability_rows = Vector{Any}()

    for (scenario_index, scenario) in enumerate(SCENARIOS)
        @printf("\n[%d/%d] %s (%s)\n",
                scenario_index, length(SCENARIOS), scenario.name, scenario.note)

        results = run_replicates(
            scenario.builder,
            K_FIT;
            seeds=SEEDS,
            COMMON_FIT...
        )

        s = summarize_replicates(results)

        stability = [temporal_stability(r) for r in results]
        mean_step, sd_step = mean_sd_pair(x.relative_mean_step for x in stability)
        drift, sd_drift = mean_sd_pair(x.total_drift for x in stability)

        push!(summary_rows, (
            scenario=scenario.name,
            note=scenario.note,
            p=P, T=T, k_fit=K_FIT, n_seeds=length(SEEDS),
            mean_rv=s.mean_rv.mean, sd_rv=s.mean_rv.sd,
            early_rv=s.early_rv.mean,
            offdiag_rv=s.mean_offdiag_rv.mean, sd_offdiag_rv=s.mean_offdiag_rv.sd,
            partial_rv=s.mean_partial_rv.mean,
            relative_frobenius=s.relative_frobenius.mean,
            sd_relative_frobenius=s.relative_frobenius.sd,
            gamma_mean=s.gamma_mean.mean,
            relative_step=mean_step, sd_relative_step=sd_step,
            total_drift=drift, sd_total_drift=sd_drift
        ))

        for (r, st) in zip(results, stability)
            push!(stability_rows, (
                scenario=scenario.name, seed=r.seed,
                mean_step=st.mean_step, max_step=st.max_step,
                relative_mean_step=st.relative_mean_step,
                total_drift=st.total_drift,
                mean_gamma=st.mean_gamma,
                fraction_gamma_above_half=st.gamma_above_half
            ))

            for d in r.diagnostics
                push!(time_rows, (
                    scenario=scenario.name, seed=r.seed, time=d.time,
                    rv=d.rv, offdiag_rv=d.offdiag_rv, partial_rv=d.partial_rv,
                    relative_frobenius=d.relative_frobenius, rmse=d.rmse
                ))
            end
        end

        @printf("     RV %.4f (sd %.4f)   offdiag RV %.4f   partial RV %.4f   relFrob %.4f\n",
                s.mean_rv.mean, s.mean_rv.sd, s.mean_offdiag_rv.mean,
                s.mean_partial_rv.mean, s.relative_frobenius.mean)
        @printf("     temporal step %.5f of ||Omega||   total drift %.5f   mean gamma %.4f\n",
                mean_step, drift, s.gamma_mean.mean)

    end

    mkpath(OUTPUT_DIRECTORY)
    write_table(joinpath(OUTPUT_DIRECTORY, "constant_precision_p5_summary.csv"), summary_rows)
    write_table(joinpath(OUTPUT_DIRECTORY, "constant_precision_p5_by_time.csv"), time_rows)
    write_table(joinpath(OUTPUT_DIRECTORY, "constant_precision_p5_stability.csv"), stability_rows)

    print_summary(summary_rows)
    return summary_rows
end

function write_table(path, rows, header=nothing)
    isempty(rows) && return
    open(path, "w") do io
        if header === nothing
            println(io, join(string.(keys(rows[1])), ","))
        else
            println(io, join(header, ","))
        end
        for r in rows
            println(io, join([value_to_field(v) for v in values(r)], ","))
        end
    end
    println("wrote ", path, " (", length(rows), " rows)")
end

value_to_field(v::AbstractString) = occursin(",", v) ? "\"$v\"" : v
value_to_field(v) = string(v)

function print_summary(rows)
    println("\n", repeat("=", 108))
    println("Non-time-varying precision matrix, p = $P, T = $T, k_fit = $K_FIT, ",
            length(SEEDS), " seeds")
    println(repeat("=", 108))
    @printf("%-24s %8s %8s %10s %10s %10s %12s %10s\n",
            "scenario", "RV", "sd", "offdiagRV", "partialRV", "relFrob", "step/||Om||", "gamma")
    println(repeat("-", 108))
    for r in rows
        @printf("%-24s %8.4f %8.4f %10.4f %10.4f %10.4f %12.5f %10.4f\n",
                r.scenario, r.mean_rv, r.sd_rv, r.offdiag_rv, r.partial_rv,
                r.relative_frobenius, r.relative_step, r.gamma_mean)
    end
    println(repeat("-", 108))
    println("step/||Om|| is the mean relative change of the estimated precision")
    println("matrix between consecutive times. The truth is constant, so this")
    println("measures spurious temporal movement; 0 is ideal.")
end

main()

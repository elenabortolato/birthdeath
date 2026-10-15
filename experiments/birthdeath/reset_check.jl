#
# Reset-prior check for Experiment II shock detection.
#
# Compares two definitions of a shock, on data simulated under each:
#
#   jump   (shock_center = :previous, the paper's model)
#          lambda_jt = lambda_{j,t-1} + D_jt,  D_jt ~ N(0, tau_a^2/k I);
#   reset  (shock_center = :zero)
#          lambda_jt ~ N(0, tau_a^2/k I), independent of lambda_{j,t-1}.
#
# Without a shock both use the spike D_jt ~ N(0, tau_i^2/k I). Each trajectory
# is simulated with shock_type = :jump and with shock_type = :reset (same
# seed, so the same births, deaths and shock times), and both models are
# fitted to both, by generalized EM and by plain MCMC started from the EM fit.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/reset_check.jl \
#       p=10 k=2 seeds=1:5 tau_active=1.5 minimum_active=6 \
#       out=results/experiment2/reset_p10_k2_tau150_min6.csv
#
# MCMC_BURNIN and MCMC_SAMPLES set the number of sweeps (default 1000, 2000).
#

include(joinpath(@__DIR__, "mcmc_check.jl"))

const DATA_TYPES = (:jump, :reset)
const MODELS = (jump=:previous, reset=:zero)

function reset_one(p, k, seed; tau_active, minimum_active)
    rows = NamedTuple[]
    for data in DATA_TYPES
        sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active, shock_type=data)
        for (model, shock_center) in pairs(MODELS)
            fit = fit_experiment2(sim, k, seed; tau_active=tau_active, shock_center=shock_center)
            chain = MCMC_time_varying(
                sim.Y, sim.alive, k;
                Lambda_init=[copy(L) for L in fit.Lambda],
                Delta_init=copy(fit.Delta),
                coefficients_init=copy(fit.coefficient_vector),
                sigma_initial=TAU0 / sqrt(k),
                sigma_active=tau_active / sqrt(k),
                sigma_inactive=TAU_INACTIVE / sqrt(k),
                a_delta=A_DELTA,
                b_delta=B_DELTA,
                coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
                coefficient_prior_sd=COEFFICIENT_PRIOR_SD,
                n_burnin=N_BURNIN,
                n_samples=N_SAMPLES,
                shock_center=shock_center,
                rng=MersenneTwister(30_000 + seed)
            )
            base = (p=p, k=k, seed=seed, tau_active=tau_active,
                    minimum_active=minimum_active, data=string(data), model=string(model))
            push!(rows, merge(
                base, (method="em",),
                shock_metrics(sim, fit.gamma),
                precision_metrics(sim, em_precisions(fit)),
                (flip_acceptance=NaN, a0=fit.coefficient_vector[1],
                 logposterior=fit.final_logposterior)
            ))
            push!(rows, merge(
                base, (method="mcmc",),
                shock_metrics(sim, chain.gamma),
                precision_metrics(sim, chain.Omega_mean),
                (flip_acceptance=chain.acceptance.flip, a0=chain.coefficient_mean[1],
                 logposterior=NaN)
            ))
        end
    end
    return rows
end

function reset_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        rows = reset_one(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        @printf("p=%d k=%d seed=%d  (%.0fs)\n", options.p, options.k, seed, time() - started)
        for r in rows
            @printf("   data=%-5s model=%-5s %-4s AUROC=%.3f recall=%.2f precision=%.2f offRV=%.3f\n",
                r.data, r.model, r.method, r.auroc, r.recall, r.precision, r.offdiag_rv)
        end
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

reset_main()

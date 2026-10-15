#
# MCMC check for Experiment II shock detection.
#
# The fitting-strategy comparison showed that the posterior mode found by
# generalized EM favours solutions without shocks. This script asks whether
# the posterior mass, rather than its mode, identifies the true shocks. For
# each trajectory it compares:
#
#   em          the default generalized EM fit;
#   mcmc_em     MCMC_time_varying started from the EM fit;
#   mcmc_truth  MCMC_time_varying started from the true parameters.
#
# The two chains target the same posterior. If they agree, the chain has
# mixed well enough for the comparison with EM to be meaningful.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/mcmc_check.jl \
#       p=10 k=2 seeds=1:20 tau_active=1.5 minimum_active=6 \
#       out=results/experiment2/mcmc_p10_k2.csv
#

include(joinpath(@__DIR__, "run_experiment2.jl"))

const N_BURNIN = parse(Int, get(ENV, "MCMC_BURNIN", "1000"))
const N_SAMPLES = parse(Int, get(ENV, "MCMC_SAMPLES", "2000"))

"""Shock-detection metrics from shock probabilities, on continuation pairs."""
function shock_metrics(sim, gamma)
    alive = sim.alive
    T_, p = size(alive)
    continuations = [
        (j, t) for t in 2:T_ for j in 1:p if alive[t, j] && alive[t - 1, j]
    ]
    scores = [gamma[j, t] for (j, t) in continuations]
    truth = BitVector([sim.s_true[j, t] for (j, t) in continuations])
    called = scores .> 0.5
    tp = count(called .& truth)
    fp = count(called .& .!truth)
    fn = count(.!called .& truth)
    return (
        auroc=auroc(scores, truth),
        recall=ratio(tp, tp + fn),
        precision=ratio(tp, tp + fp),
        gamma_true_shocks=nanmean(scores[truth]),
        gamma_non_shocks=nanmean(scores[.!truth])
    )
end

"""RV and off-diagonal RV of estimated precision matrices on the active sets."""
function precision_metrics(sim, Omega_est)
    rv = Float64[]
    offdiag = Float64[]
    for t in 1:size(sim.alive, 1)
        idx = findall(view(sim.alive, t, :))
        isempty(idx) && continue
        L = sim.Lambda[t][idx, :]
        truth = L * L' + Diagonal(sim.Delta[idx])
        estimate = Omega_est[t][idx, idx]
        push!(rv, rv_coefficient(truth, estimate))
        push!(offdiag, offdiagonal_rv(truth, estimate))
    end
    return (rv=nanmean(rv), offdiag_rv=nanmean(offdiag))
end

em_precisions(fit) = [L * L' + Diagonal(fit.Delta) for L in fit.Lambda]

function mcmc_one(p, k, seed; tau_active, minimum_active)
    sim = simulate_experiment2(p, k, seed; minimum_active=minimum_active)
    fit = fit_experiment2(sim, k, seed; tau_active=tau_active)

    mcmc_settings = (
        sigma_initial=TAU0 / sqrt(k),
        sigma_active=tau_active / sqrt(k),
        sigma_inactive=TAU_INACTIVE / sqrt(k),
        a_delta=A_DELTA,
        b_delta=B_DELTA,
        coefficient_prior_mean=COEFFICIENT_PRIOR_MEAN,
        coefficient_prior_sd=COEFFICIENT_PRIOR_SD,
        n_burnin=N_BURNIN,
        n_samples=N_SAMPLES
    )
    chains = [
        "mcmc_em" => MCMC_time_varying(
            sim.Y, sim.alive, k;
            Lambda_init=[copy(L) for L in fit.Lambda],
            Delta_init=copy(fit.Delta),
            coefficients_init=copy(fit.coefficient_vector),
            mcmc_settings...,
            rng=MersenneTwister(30_000 + seed)
        ),
        "mcmc_truth" => MCMC_time_varying(
            sim.Y, sim.alive, k;
            Lambda_init=[copy(L) for L in sim.Lambda],
            Delta_init=copy(sim.Delta),
            coefficients_init=[A0, A_EVENT, A_EVENT],
            mcmc_settings...,
            rng=MersenneTwister(40_000 + seed)
        )
    ]

    base = (p=p, k=k, seed=seed, tau_active=tau_active, minimum_active=minimum_active)
    rows = NamedTuple[]
    push!(rows, merge(
        base, (method="em",),
        shock_metrics(sim, fit.gamma),
        precision_metrics(sim, em_precisions(fit)),
        (flip_acceptance=NaN, a0=fit.coefficient_vector[1],
         a_birth=fit.coefficient_vector[2], a_death=fit.coefficient_vector[3])
    ))
    for (name, chain) in chains
        push!(rows, merge(
            base, (method=name,),
            shock_metrics(sim, chain.gamma),
            precision_metrics(sim, chain.Omega_mean),
            (flip_acceptance=chain.acceptance.flip, a0=chain.coefficient_mean[1],
             a_birth=chain.coefficient_mean[2], a_death=chain.coefficient_mean[3])
        ))
    end
    return rows
end

function mcmc_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        rows = mcmc_one(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active
        )
        summary = join(
            [@sprintf("%s AUROC=%.3f recall=%.2f offRV=%.3f", r.method, r.auroc, r.recall, r.offdiag_rv)
             for r in rows], " | "
        )
        @printf("p=%d k=%d seed=%d  %s  (%.0fs)\n",
            options.p, options.k, seed, summary, time() - started)
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

if abspath(PROGRAM_FILE) == @__FILE__
    mcmc_main()
end

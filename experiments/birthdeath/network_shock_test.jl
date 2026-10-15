#
# Test of the network-informed shock probability.
#
# Shocks are simulated with
#
#   logit pi_jt = a0 + a_birth * sum_{l born at t} |rho~_{t,jl}|
#                    + a_death * sum_{l dying at t} |rho_{t-1,jl}|,
#
# so a continuing variable is more likely to be shocked when it is strongly
# linked to the variables entering or leaving. With link=loadings the link is
# the absolute loading inner product |<lambda_{t-1,j}, lambda_l>| instead of
# the partial correlation. Each trajectory is fitted with
#
#   network  the network-informed shock probability (the true model);
#   counts   the current model, one probability per time from birth and
#            death counts, with the paper's priors;
#   oracle   the network model started from the true parameters.
#
# Usage (from the repository root):
#
#   julia --project=. experiments/birthdeath/network_shock_test.jl \
#       p=10 k=5 seeds=1:20 tau_active=1.5 minimum_active=6 link=loadings \
#       out=results/experiment2/network_p10_k5.csv
#

include(joinpath(@__DIR__, "oracle_check.jl"))

# True slopes on the summed link strengths, chosen so that a typical link adds
# about 1.2 to 1.4 to the log odds of a shock. Nonzero summed |rho| average
# about 0.3 to 0.4; nonzero summed loading inner products average about 2.7.
const NETWORK_SLOPE = Dict(:partial_correlation => 4.0, :loadings => 0.5)

# Weakly informative prior for the network model's coefficients.
const NETWORK_PRIOR_MEAN = [0.0, 0.0, 0.0]
const NETWORK_PRIOR_SD = [2.0, 2.0, 2.0]

function network_one(p, k, seed; tau_active, minimum_active, link)
    slope = NETWORK_SLOPE[link]
    sim = simulate_experiment2(
        p, k, seed;
        minimum_active=minimum_active,
        shock_rule=:network,
        network_link=link,
        a_birth=slope,
        a_death=slope
    )
    network_priors = (
        coefficient_prior_mean=NETWORK_PRIOR_MEAN,
        coefficient_prior_sd=NETWORK_PRIOR_SD
    )

    # How much the true shock probabilities vary across variables.
    continuations = [
        (j, t) for t in 2:size(sim.alive, 1) for j in 1:p
        if sim.alive[t, j] && sim.alive[t - 1, j]
    ]
    true_pi = [sim.pi_time[j, t] for (j, t) in continuations]

    network_fit = fit_experiment2(
        sim, k, seed;
        tau_active=tau_active, shock_covariates=:network, network_link=link,
        network_priors...
    )
    counts_fit = fit_experiment2(sim, k, seed; tau_active=tau_active)
    oracle = oracle_em(
        sim, k, seed;
        tau_active=tau_active,
        true_coefficients=[A0, slope, slope],
        shock_covariates=:network,
        network_link=link,
        network_priors...
    )

    network = evaluate_experiment2(sim, network_fit.Lambda, network_fit.Delta, network_fit.gamma)
    counts = evaluate_experiment2(sim, counts_fit.Lambda, counts_fit.Delta, counts_fit.gamma)
    oracle_metrics = evaluate_experiment2(
        sim, oracle.state.Lambda, oracle.state.Delta, oracle.gamma
    )

    return (
        p=p, k=k, seed=seed, tau_active=tau_active, minimum_active=minimum_active,
        link=link, true_slope=slope,
        shock_rate=network.shock_rate,
        true_pi_mean=mean(true_pi),
        true_pi_max=maximum(true_pi),
        network_auroc=network.auroc,
        network_recall=network.recall,
        network_precision=network.precision,
        network_offdiag_rv=network.offdiag_rv,
        network_logposterior=network_fit.final_logposterior,
        network_a0=network_fit.coefficient_vector[1],
        network_a_birth=network_fit.coefficient_vector[2],
        network_a_death=network_fit.coefficient_vector[3],
        counts_auroc=counts.auroc,
        counts_recall=counts.recall,
        counts_precision=counts.precision,
        counts_offdiag_rv=counts.offdiag_rv,
        oracle_auroc=oracle_metrics.auroc,
        oracle_recall=oracle_metrics.recall,
        oracle_precision=oracle_metrics.precision,
        oracle_offdiag_rv=oracle_metrics.offdiag_rv,
        oracle_logposterior=oracle.logposterior,
        oracle_a_birth=oracle.coefficients[2],
        oracle_a_death=oracle.coefficients[3]
    )
end

function network_main()
    options = parse_args(ARGS)
    io = isempty(options.out) ? nothing : open(options.out, "w")
    header_written = false

    for seed in options.seeds
        started = time()
        row = network_one(
            options.p, options.k, seed;
            tau_active=options.tau_active,
            minimum_active=options.minimum_active,
            link=options.link
        )
        @printf(
            "p=%d k=%d seed=%d  AUROC network=%.3f counts=%.3f oracle=%.3f  a_b=%.2f a_d=%.2f  (%.0fs)\n",
            row.p, row.k, row.seed, row.network_auroc, row.counts_auroc,
            row.oracle_auroc, row.network_a_birth, row.network_a_death,
            time() - started
        )
        flush(stdout)
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

network_main()

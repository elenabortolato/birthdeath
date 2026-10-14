# TIME VARYING GRAPHICAL MODEL 

using LinearAlgebra, Distributions, Random, Statistics

Random.seed!(12345)
# ==================================================
# Utilities
# ==================================================
logsumexp(v) = maximum(v) + log(sum(exp.(v .- maximum(v))))

# Safe cholesky
function safe_cholesky(Λ, Δ; jitter=1e-8)
    Ω = Λ*Λ' + Diagonal(Δ)
    return cholesky(Symmetric(Ω + jitter*I); check=false)
end

# ==================================================
# Forward-Backward for HMM
# ==================================================
function forward_backward(loglik, logA, init_logp)
    T = size(loglik, 2)
    α = fill(-Inf, 2, T)
    β = fill(0.0, 2, T)
    
    α[:, 1] .= init_logp .+ loglik[:, 1]
    
    for t in 2:T
        for s in 1:2
            α[s, t] = logsumexp([α[sp, t-1] + logA[sp, s] + loglik[s, t] for sp in 1:2])
        end
    end
    for t in T-1:-1:1
        for s in 1:2
            β[s, t] = logsumexp([β[sn, t+1] + logA[s, sn] + loglik[sn, t+1] for sn in 1:2])
        end
    end

    γ = zeros(2, T)
    for t in 1:T
        logZ = logsumexp(α[:, t] + β[:, t])
        γ[:, t] .= exp.(α[:, t] + β[:, t] .- logZ)
    end
    
    # Pairwise ξ_t(s_prev, s_curr) for π update
    ξ = zeros(2,2,T-1)
    for t in 1:T-1
        denom = logsumexp([α[sp,t] + logA[sp,sn] + loglik[sn,t+1] + β[sn,t+1] for sp in 1:2, sn in 1:2])
        for sp in 1:2, sn in 1:2
            ξ[sp,sn,t] = exp(α[sp,t] + logA[sp,sn] + loglik[sn,t+1] + β[sn,t+1] - denom)
        end
    end
    
    return γ, ξ
end 
 
# ==================================================
# E-step
# ==================================================
function E_step_HMM_OLD(Λlist, π_; σ_active=2.0, σ_inactive=0.01)
    T = length(Λlist)
    p, k = size(Λlist[1])
    γ_active = zeros(p, k, T)
    ξ_active = zeros(p, k) # sum of ξ for inactive->active
    
    # Transition matrix
    logA = [
        log(1-π_) log(π_);
        log(1-π_)   log(π_)
    ]
    init_logp = [log(0.99), log(0.01)]
    
    for j in 1:p, h in 1:k
        loglik = zeros(2, T)
        for t in 1:T
            loglik[1,t] = logpdf(Normal(0, σ_active), Λlist[t][j,h])
        end
        loglik[2,1] = loglik[1,1]
        for t in 2:T
            loglik[2,t] = logpdf(Normal(Λlist[t-1][j,h], σ_inactive), Λlist[t][j,h])
        end
        
        γ, ξ = forward_backward(loglik, logA, init_logp)
        γ_active[j,h,:] .= γ[1,:]
        ξ_active[j,h] = sum(ξ[2,1,:]) # inactive -> active
    end
    
    return γ_active, ξ_active
end


# Vector log-sum-exp (do NOT shadow SpecialFunctions.logsumexp)
logsumexp_vec(v) = maximum(v) + log(sum(exp.(v .- maximum(v))))

# Elementwise pairwise log-sum-exp
function logsumexp2(a, b)
    m = max.(a, b)
    return m .+ log.(exp.(a .- m) .+ exp.(b .- m))
end

# ==================================================
# Build log-likelihood tensor
# ==================================================
# loglik :: 2 × T × p × k
function build_loglik(Λlist, σ_active, σ_inactive)
    T = length(Λlist)
    p, k = size(Λlist[1])

    loglik = zeros(2, T, p, k)

    for t in 1:T
        Λt = Λlist[t]
        loglik[1,t,:,:] .= logpdf.(Normal(0, σ_active), Λt)

        if t == 1
            loglik[2,t,:,:] .= loglik[1,t,:,:]
        else
            loglik[2,t,:,:] .= logpdf.(Normal.(Λlist[t-1], σ_inactive), Λt)
        end
    end

    return loglik
end

 
# ==================================================
# E-step
# ==================================================
function E_step_HMM(Λlist, π_; σ_active=2.0, σ_inactive=0.01)
    T = length(Λlist)
    p, k = size(Λlist[1])
    γ_active = zeros(p, k, T)
    ξ_active = zeros(p, k) # sum of ξ for inactive->active
    
    
    init_logp = [log(0.99), log(0.01)]
    logA = [
        log(1-π_)  log(π_);
        log(1-π_)  log(π_)
    ]

    init_logp = [log(0.99), log(0.01)]

    loglik = build_loglik(Λlist, σ_active, σ_inactive)
    γ, ξ = forward_backward_3d(loglik, logA, init_logp)

    return γ, ξ
end

function E_step_HMM(Λlist, π_; σ_active=2.0, σ_inactive=0.01)

    T = length(Λlist)
    p, k = size(Λlist[1])

    # outputs
    γ = zeros(p, k, T)
    ξ = zeros(p, k, 2, 2, T-1)

    # HMM parameters
    logA = [
        log(1 - π_)   log(π_);
        log(1 - π_)   log(π_)
    ]

    init_logp = [log(0.99), log(0.01)]

    for j in 1:p, h in 1:k

        # build log-likelihood for THIS (j,h)
        loglik = zeros(2, T)

        for t in 1:T
            loglik[1,t] = logpdf(Normal(0, σ_active), Λlist[t][j,h])
        end

        loglik[2,1] = loglik[1,1]
        for t in 2:T
            loglik[2,t] =
                logpdf(
                    Normal(Λlist[t-1][j,h], σ_inactive),
                    Λlist[t][j,h]
                )
        end

        # run HMM
        γ_jh, ξ_jh = forward_backward(loglik, logA, init_logp)

        # store
        γ[j,h,:] .= γ_jh[1,:]               # active state
        ξ[j,h,:,:,:] .= ξ_jh
    end

    return γ, ξ
end


# ==================================================
# Δ update (Bayesian)
# ==================================================
function update_Delta_bayes_old(Y, Λlist, Δ; a=1.0, b=1.0, step=0.1)
    T = length(Y)
    p = length(Y[1])
    θ = log.(Δ)
    θ_new = copy(θ)
    
    UnitVector(j,p) = let v=zeros(p); v[j]=1.0; v end
    
    for j in 1:p
        g = 0.0
        I = 0.0
        for t in 1:T
            Λ = Λlist[t]
            y = Y[t]
            ch = safe_cholesky(Λ, Δ)
            z = ch \ y
            Σjj = sum((ch \ UnitVector(j,p)).^2)
            g += 0.5*Δ[j]*(Σjj - z[j]^2)
            I += 0.5*Δ[j]^2*Σjj^2
        end
        g += (a-1) - b*Δ[j]
        I += b*Δ[j]
        θ_new[j] += step*g/I
    end
    
    return exp.(θ_new)
end

function update_Delta_bayes(
    Y, Λlist, Δ;
    a=1.0,
    b=1.0,
    step=0.1,
    Δ_min=1e-6,
    Δ_max=1e6
)
    T = length(Y)
    p = length(Y[1])

    θ = log.(Δ)          # unconstrained parameters
    θ_new = copy(θ)

    UnitVector(j,p) = (v=zeros(p); v[j]=1.0; v)

    for j in 1:p
        g = 0.0
        I = 0.0

        for t in 1:T
            Λ = Λlist[t]
            y = Y[t]

            ch = safe_cholesky(Λ, Δ)
            z = ch \ y

            ej = UnitVector(j,p)
            v = ch \ ej
            Σjj = dot(v, v)

            g += 0.5 * Δ[j] * (Σjj - z[j]^2)
            I += 0.5 * Δ[j]^2 * Σjj^2
        end

        # Gamma prior
        g += (a - 1) - b * Δ[j]
        I += b * Δ[j]

        θ_new[j] += step * g / I
    end

    # Map back and clamp
    Δ_new = exp.(θ_new)
    Δ_new .= clamp.(Δ_new, Δ_min, Δ_max)

    return Δ_new
end

# ==================================================
# Soft EM update for Λ
# ==================================================
function update_Lambda_softEM(Y, Λlist, Δ, γ; step=1e-3, σ_active=2.0, σ_inactive=0.01)
    T = length(Λlist)
    p, k = size(Λlist[1])
    Λnew = deepcopy(Λlist)
    
    for t in 1:T
        Λt = Λnew[t]
        y = Y[t]
        Ω = Λt*Λt' + Diagonal(Δ)
        ch = cholesky(Symmetric(Ω))
        ΣΛ = ch \ (ch' \ Λt)
        α = ch \ y
        G_data = ΣΛ - α*(α'*Λt)
        
        for j in 1:p, h in 1:k
            γt = γ[j,h,t]
            grad_prior = -γt*Λt[j,h]/σ_active^2
            if t>1
                grad_prior -= (1-γt)*(Λt[j,h]-Λnew[t-1][j,h])/σ_inactive^2
            end
            if t<T
                grad_prior -= (1-γt)*(Λt[j,h]-Λnew[t+1][j,h])/σ_inactive^2
            end
            Λnew[t][j,h] += step*(G_data[j,h]+grad_prior) 
        end
    end
    
    return Λnew
end

function update_Lambda_softEM_Laplace(
    Y, Λlist, Δ, γ;
    step=1e-3,
    σ_active=2.0,
    σ_inactive=0.01,
    eps=1e-6)
    T = length(Λlist)
    p, k = size(Λlist[1])
    Λnew = deepcopy(Λlist)
    
    b_active=σ_active
    b_inactive= σ_active
    for t in 1:T
        Λt = Λnew[t]
        y = Y[t]

        Ω = Λt * Λt' + Diagonal(Δ)
        ch = cholesky(Symmetric(Ω))

        ΣΛ = ch \ (ch' \ Λt)
        α = ch \ y
        G_data = ΣΛ - α * (α' * Λt)

        for j in 1:p, h in 1:k
            γt = γ[j,h,t]

            # -----------------------------
            # Laplace prior gradients
            # -----------------------------
            # Active: |Λ|
            grad_prior =
                -γt * (Λt[j,h] / sqrt(Λt[j,h]^2 + eps)) / b_active

            # Inactive: |Λ_t - Λ_{t-1}|
            if t > 1
                diff = Λt[j,h] - Λnew[t-1][j,h]
                grad_prior -=
                    (1 - γt) * (diff / sqrt(diff^2 + eps)) / b_inactive
            end

            # Inactive: |Λ_{t+1} - Λ_t|
            if t < T
                diff = Λt[j,h] - Λnew[t+1][j,h]
                grad_prior -=
                    (1 - γt) * (diff / sqrt(diff^2 + eps)) / b_inactive
            end

            # Gradient ascent on posterior
            Λnew[t][j,h] += step * (G_data[j,h] + grad_prior)
        end
    end

    return Λnew
end

function log_posterior_old(
    Y, Λ, Δ, γ, ξ, π_;
    σ_active=2.0,
    σ_inactive=0.01,
    a=1.0,
    b=1.0
)
    T = length(Y)
    p, k = size(Λ[1])

    lp = 0.0

    # --------------------------------------------------
    # (1) Data term (Gaussian precision parameterization)
    # y ~ N(0, Ω^{-1}), Ω = ΛΛ' + diag(Δ)
    # --------------------------------------------------
    for t in 1:T
        Ω = Λ[t] * Λ[t]' + Diagonal(Δ)
        ch = cholesky(Symmetric(Ω))
        lp += logdet(ch)              # +½ log|Ω|
        lp -= 0.5 * dot(Y[t], Ω * Y[t])
    end

    # --------------------------------------------------
    # (2) Λ prior (HMM-weighted Gaussian)
    # --------------------------------------------------
    for j in 1:p, h in 1:k
        for t in 1:T
            γt = γ[j,h,t]

            # Active state
            lp -= 0.5 * γt * (Λ[t][j,h]^2) / σ_active^2
            lp -= 0.5 * γt * log(2π * σ_active^2)

            # Inactive (random walk)
            if t > 1
                diff = Λ[t][j,h] - Λ[t-1][j,h]
                lp -= 0.5 * (1 - γt) * (diff^2) / σ_inactive^2
                lp -= 0.5 * (1 - γt) * log(2π * σ_inactive^2)
        
            end
        end
    end

    # --------------------------------------------------
    # (3) HMM transition prior
    # --------------------------------------------------
   
    lp += sum(ξ[:, :, 2, 1, :]) * log(π_)       # inactive → active
    lp += sum(ξ[:, :, 2, 2, :]) * log(1 - π_)   # inactive → inactive
    
    # --------------------------------------------------
    # (4) Δ prior (Gamma)
    # --------------------------------------------------
    lp += sum((a - 1) * log.(Δ) - b .* Δ)

    return lp
end

function log_posterior(
    Y, Λ, Δ, γ, π_;
    σ_active=2.0,
    σ_inactive=0.01,
    a=1.0,
    b=1.0
)
    T = length(Y)
    p, k = size(Λ[1])

    lp = 0.0

    # --------------------------------------------------
    # (1) Likelihood: y_t ~ N(0, Ω_t^{-1})
    # --------------------------------------------------
    for t in 1:T
        Ω = Λ[t] * Λ[t]' + Diagonal(Δ)
        ch = cholesky(Symmetric(Ω))
        lp += logdet(ch)                    # +½ log|Ω|
        lp -= 0.5 * dot(Y[t], Ω * Y[t])     # −½ y'Ωy
    end

    # --------------------------------------------------
    # (2) Λ prior: spike-and-slab with RW slab
    # --------------------------------------------------
    for j in 1:p, h in 1:k
        for t in 1:T
            γt = γ[j,h,t]

            # ---- active ----
            lp -= 0.5 * γt * (Λ[t][j,h]^2) / σ_active^2
            lp -= 0.5 * γt * log(2π * σ_active^2)

            # ---- inactive ----
            if t == 1
                lp -= 0.5 * (1 - γt) * (Λ[t][j,h]^2) / σ_inactive^2
                lp -= 0.5 * (1 - γt) * log(2π * σ_inactive^2)
            else
                diff = Λ[t][j,h] - Λ[t-1][j,h]
                lp -= 0.5 * (1 - γt) * (diff^2) / σ_inactive^2
                lp -= 0.5 * (1 - γt) * log(2π * σ_inactive^2)
            end
        end
    end

    # --------------------------------------------------
    # (3) Bernoulli prior on γ (NOT transitions!)
    # --------------------------------------------------
    lp += sum(γ) * log(π_) + sum(1 .- γ) * log(1 - π_)
     # Bernoulli prior on γ
 
    # --------------------------------------------------
    # (4) Δ prior: Gamma(a,b)
    # --------------------------------------------------
    lp += sum((a - 1) * log.(Δ) - b .* Δ)

    return lp
end

# ==================================================
# EM loop
# ==================================================
function stochastic_EM_time_varying_old(Y, k; niter=500, step=0.01, π_init=0.1, a_delta=1.0, b_delta=1.0, verbose=true)
    T = length(Y)
    p = length(Y[1])
    
    Λ = [randn(p,k)*0.1 for _ in 1:T]
    for t in 2:T
        Λ[t] = copy(Λ[t-1])
    end
    Δ = ones(p)*0.1
    π_ = π_init
    ll_history = Float64[]
    
    for iter in 1:niter
        # E-step
        γ, ξ = E_step_HMM(Λ, π_; σ_active=2.0, σ_inactive=0.01)
        size(γ)
        size(s_true)
        # ==================================================
 
        # M-step: π update
        π_ = mean(γ)
 
     #   π_ = sum(ξ) / sum(1 .- γ[:, :, 1:end-1]) # soft inactive->active
       # π = sum(ξ[:, :, 2, 1, :]) / sum(1 .- γ[:, :, 1:end-1])
 
        # M-step: Λ update
        Λ = update_Lambda_softEM(Y, Λ, Δ, γ; step=step, σ_active=2.0, σ_inactive=0.01)
        
        # M-step: Δ update
        Δ = update_Delta_bayes(Y, Λ, Δ; a=a_delta, b=b_delta, step=step)
        
        # Log-likelihood
        lp=log_posterior(Y, Λ, Δ, γ, ξ, π_; σ_active=2.0, σ_inactive=0.01, a=a_delta, b=b_delta)
       # ll = sum([-0.5*logdet(Λ[t]*Λ[t]'+Diagonal(Δ)) - 0.5*(y'*inv(Λ[t]*Λ[t]'+Diagonal(Δ))*y) for (t,y) in enumerate(Y)])
        push!(ll_history,lp)
         
        if verbose && iter % 20 == 0
            println("Iter $iter: LL=$(round(lp,digits=3)), π_=$(round(π_,digits=3))")
        end
    end
    
    return (Λ=Λ, Δ=Δ, π_=π_, γ=γ, ll_history=ll_history)
end

function stochastic_EM_time_varying_check(
    Y, k;
    niter=200,
    step=0.01,
    π_init=0.1,
    a_delta=1.0,
    b_delta=1.0,
    tol=1e-6,
    verbose=true)

    T = length(Y)
    p = length(Y[1])

    # -------------------------------
    # Initialization
    # -------------------------------
    Λ = [randn(p,k)*0.1 for _ in 1:T]
    for t in 2:T
        Λ[t] = copy(Λ[t-1])
    end

    Δ = ones(p) * 0.1
    π_ = π_init

    ll_history = Float64[]

    # initial E-step
    γ, ξ = E_step_HMM(Λ, π_; σ_active=2.0, σ_inactive=0.01)

    lp_old = log_posterior(
        Y, Λ, Δ, γ, ξ, π_;
        σ_active=2.0,
        σ_inactive=0.01,
        a=a_delta,
        b=b_delta
    )

    push!(ll_history, lp_old)

    # -------------------------------
    # EM loop
    # -------------------------------
    for iter in 1:niter

        # =============================
        # E-step
        # =============================
        γ, ξ = E_step_HMM(Λ, π_; σ_active=2.0, σ_inactive=0.01)

        lp_E = log_posterior(
            Y, Λ, Δ, γ, ξ, π_;
            σ_active=2.0,
            σ_inactive=0.01,
            a=a_delta,
            b=b_delta
        )

        if lp_E < lp_old - tol
            @warn "Log-posterior decreased after E-step at iter $iter: Δ=$(lp_E - lp_old)"
        end

        lp_old = lp_E

        # =============================
        # M-step: π
        # =============================
      #  π_ = sum(ξ[:,:,2,1,:]) / sum(1 .- γ[:,:,1:end-1])
        π_ = mean( γ)
        π_ = clamp(π_, 1e-6, 1 - 1e-6)

        lp_π = log_posterior(
            Y, Λ, Δ, γ, ξ, π_;
            σ_active=2.0,
            σ_inactive=0.01,
            a=a_delta,
            b=b_delta
        )

        if lp_π < lp_old - tol
            @warn "Log-posterior decreased after π update at iter $iter: Δ=$(lp_π - lp_old)"
        end

        lp_old = lp_π

        # =============================
        # M-step: Λ
        # =============================
        Λ = update_Lambda_softEM(
            Y, Λ, Δ, γ;
            step=step,
            σ_active=2.0,
            σ_inactive=0.01
        )

        lp_Λ = log_posterior(
            Y, Λ, Δ, γ, ξ, π_;
            σ_active=2.0,
            σ_inactive=0.01,
            a=a_delta,
            b=b_delta
        )

        if lp_Λ < lp_old - tol
            @warn "Log-posterior decreased after Λ update at iter $iter: Δ=$(lp_Λ - lp_old)"
        end

        lp_old = lp_Λ

        # =============================
        # M-step: Δ
        # =============================
        Δ = update_Delta_bayes(
            Y, Λ, Δ;
            a=a_delta,
            b=b_delta,
            step=step
        )

        lp_new = log_posterior(
            Y, Λ, Δ, γ, ξ, π_;
            σ_active=2.0,
            σ_inactive=0.01,
            a=a_delta,
            b=b_delta
        )

        if lp_new < lp_old - tol
            @warn "Log-posterior decreased after Δ update at iter $iter: Δ=$(lp_new - lp_old)"
        end

        lp_old = lp_new
        push!(ll_history, lp_new)

        if verbose && iter % 20 == 0
            println(
                "Iter $iter | LP = $(round(lp_new,digits=3)) | π = $(round(π_,digits=3))"
            )
        end
    end

    return (
        Λ = Λ,
        Δ = Δ,
        π_ = π_,
        γ = γ,
        ll_history = ll_history
    )
end

function stochastic_EM_time_varying(
    Y, k;
    niter=200,
    step=0.01,
    π_init=0.1,
    a_delta=1.0,
    b_delta=1.0,
    tol=1e-8,
    verbose=true
)

    T = length(Y)
    p = length(Y[1])

    # -------------------------------
    # Initialization
    # -------------------------------
    Λ = [randn(p,k)*0.1 for _ in 1:T]
    for t in 2:T
        Λ[t] = copy(Λ[t-1])
    end

    Δ = ones(p) * 0.1
    π_ = π_init

    γ, ξ = E_step_HMM(Λ, π_; σ_active=2.0, σ_inactive=0.01)

    lp = log_posterior(
        Y, Λ, Δ, γ, ξ, π_;
        σ_active=2.0,
        σ_inactive=0.01,
        a=a_delta,
        b=b_delta
    )

    ll_history = [lp]

    # -------------------------------
    # EM loop
    # -------------------------------
    for iter in 1:niter

        # =============================
        # E-step (always accept)
        # =============================
        γ, ξ = E_step_HMM(Λ, π_; σ_active=2.0, σ_inactive=0.01)

        # =============================
        # π update (accept/reject)
        # =============================
        π_candidate = mean(γ)#sum(ξ[:,:,2,1,:]) / sum(1 .- γ[:,:,1:end-1])
        π_candidate = clamp(π_candidate, 1e-6, 1 - 1e-6)

        lp_candidate = log_posterior(
            Y, Λ, Δ, γ, ξ, π_candidate;
            σ_active=2.0,
            σ_inactive=0.01,
            a=a_delta,
            b=b_delta
        )

        if lp_candidate ≥ lp - tol 
            π_ = π_candidate
            lp = lp_candidate
        end

        # =============================
        # Λ update (accept/reject)
        # =============================
        Λ_candidate = update_Lambda_softEM(
            Y, Λ, Δ, γ;
            step=step,
            σ_active=2.0,
            σ_inactive=0.01
        )

        lp_candidate = log_posterior(
            Y, Λ_candidate, Δ, γ, ξ, π_;
            σ_active=2.0,
            σ_inactive=0.01,
            a=a_delta,
            b=b_delta
        )

        if iter <20
        if lp_candidate ≥ lp - tol  
            Λ = Λ_candidate
            lp = lp_candidate
        end
         end
        # =============================
        # Δ update (accept/reject)
        # =============================
        Δ_candidate = update_Delta_bayes(
            Y, Λ, Δ;
            a=a_delta,
            b=b_delta,
            step=step
        )

        lp_candidate = log_posterior(
            Y, Λ, Δ_candidate, γ, ξ, π_;
            σ_active=2.0,
            σ_inactive=0.01,
            a=a_delta,
            b=b_delta
        )

        if lp_candidate ≥ lp - tol  
            Δ = Δ_candidate
            lp = lp_candidate
        end

        push!(ll_history, lp)

        if verbose && iter % 20 == 0
            println(
                "Iter $iter | LP=$(round(lp,digits=3)) | π=$(round(π_,digits=3))"
            )
        end
    end

    return (
        Λ = Λ,
        Δ = Δ,
        π_ = π_,
        γ = γ,
        ll_history = ll_history
    )
end

# ==================================================
# Generate time-varying data
# ==================================================
function generate_time_varying_data(p,k,T; π_change=0.05)
    Λ = Vector{Matrix{Float64}}(undef,T)
    Δ = ones(p)*0.1
    Λ[1] = randn(p,k).*rand(Bernoulli(0.5),p,k)
    s_true = falses(p,k,T)
    s_true[:,:,1] .= true
    
    for t in 2:T
        Λ[t] = copy(Λ[t-1])+randn(p,k)*0.01
        for j in 1:p, h in 1:k
            if rand() < π_change
                Λ[t][j,h] = randn()*rand(Bernoulli(0.5))
                s_true[j,h,t] = true
            end
        end
    end
    
    Y = Vector{Vector{Float64}}(undef,T)
    for t in 1:T
        Ω = (Λ[t]*Λ[t]' + Diagonal(Δ))
        Ω = Symmetric( Ω )
        SIGMA=inv(Ω)
        Y[t] = rand(MvNormal(zeros(p),SIGMA))
    end
    
    return Y, Λ, Δ, s_true
end

# ==================================================
# Example run
# ==================================================
p,k,T = 10,2,40
changerate = 1/T/k
Y, Λ_true, Δ_true, s_true = generate_time_varying_data(p,k,T; π_change=changerate)
s_true
sum(s_true)-p*k
maxiter=2000
result = stochastic_EM_time_varying_check(Y,k; niter=maxiter, step=0.0001, π_init=0.3, verbose=true)


result.π_
size(result.γ)
println("Done. Final π = $(result.π_)")

times=[1,15,40]
function compare_precision_matrices(Λ_true_list, Δ_true, Λ_est_list, Δ_est, times)
    n = length(times)
    plots = []
    
    for (idx, t) in enumerate(times)
        Ω_true = Λ_true_list[t] * Λ_true_list[t]' + Diagonal(Δ_true)
        Ω_est = Λ_est_list[t] * Λ_est_list[t]' + Diagonal(Δ_est)
        
        p1 = heatmap(Ω_true, 
                     title="True Ω (t=$t)",
                   #  clim=(-8,8),
                     color=:viridis)
        
        p2 = heatmap(Ω_est,
                     title="Estimated Ω  (t=$t)",
                  #   clim=(-8,8),
                     color=:viridis)
        
        push!(plots, p1, p2)
    end
    
    return plot(plots..., layout=(n, 2), size=(800, 300*n))
end

function compare_covariance_matrices(Λ_true_list, Δ_true, Λ_est_list, Δ_est, times)
    n = length(times)
    plots = []
    
    for (idx, t) in enumerate(times)
        Ω_true = Λ_true_list[t] * Λ_true_list[t]' + Diagonal(Δ_true)
        Ω_est = Λ_est_list[t] * Λ_est_list[t]' + Diagonal(Δ_est)
        
        p1 = heatmap(inv(Ω_true), 
                     title="True Ω^-1 (t=$t)",
                  #   clim=(-8,8),
                     color=:viridis)
        
        p2 = heatmap(inv(Ω_est),
                     title="Estimated Ω^-1  (t=$t)",
                   #     clim=(-8,8),
                     color=:viridis)
        
        push!(plots, p1, p2)
    end
    
    return plot(plots..., layout=(n, 2), size=(800, 300*n))
end
compare_precision_matrices(Λ_true, Δ_true, result.Λ, result.Δ, times)
compare_covariance_matrices(Λ_true, Δ_true, result.Λ, result.Δ, times)

plot(result.ll_history[1:maxiter], 
     xlabel="Iteration", 
     ylabel="Log-likelihood",
     title="EM Convergence",
     legend=false)


 # RV Coefficient
 rv_coeff = zeros(T)  # Initialize RV coefficient array
 for time_ in 1:T
     true_prec=Λ_true[ time_]*Λ_true[ time_]'+Diagonal(Δ_true[1:p])
     hat_prec=result.Λ[time_]*result.Λ[time_]'+Diagonal(result.Δ )
     rv_coeff[time_] = (tr((true_prec' * hat_prec))) / sqrt(tr(true_prec' * true_prec) * tr(hat_prec' * hat_prec))
  end
  rv_coeff
  plot(rv_coeff, title = "RV Coefficient", ylabel = "RV Coefficient", xlabel = "Time"
  ,
  ylim=(0,1)
  )
 
s_true
result.γ
p 
k
T

size(s_true)
size(result.γ)


using Statistics

"""
    roc_auc(scores, labels)

Compute ROC curve and AUC.
- scores :: AbstractVector{<:Real}  (predicted probabilities)
- labels :: AbstractVector{Bool} or {0,1}

Returns: fpr, tpr, auc
"""
function roc_auc(scores, labels)
    @assert length(scores) == length(labels)

    labels = Bool.(labels)
    P = sum(labels)
    N = length(labels) - P
    @assert P > 0 && N > 0 "Need both positive and negative samples"

    # Sort by descending score
    order = sortperm(scores, rev=true)
    scores = scores[order]
    labels = labels[order]

    tpr = Float64[]
    fpr = Float64[]

    tp = 0
    fp = 0
    last_score = Inf

    for i in eachindex(scores)
        if scores[i] != last_score
            push!(tpr, tp / P)
            push!(fpr, fp / N)
            last_score = scores[i]
        end

        if labels[i]
            tp += 1
        else
            fp += 1
        end
    end

    push!(tpr, tp / P)
    push!(fpr, fp / N)

    # Trapezoidal AUC
    auc = sum((fpr[i] - fpr[i-1]) * (tpr[i] + tpr[i-1]) / 2
              for i in 2:length(fpr))

    return fpr, tpr, auc
end


fpr, tpr, auc = roc_auc(vec(result.γ), vec(s_true))
histogram(vec(result.γ))

plot(fpr,tpr)
heatmap(result.γ[9,:,:])
heatmap(s_true[9,:,:])
γ_sum = [mean(result.γ[j, :, :], dims=1)  for j in 1:p]  # sums over h
s_sum = [mean(s_true[j, :, :], dims=1) for j in 1:p]

heatmap(s_sum)
heatmap(γ_sum)

# Sum (or mean) over h dimension for all j
γ_sum = mean(result.γ, dims=2)  # average over h
γ_sum = dropdims(γ_sum; dims=2) # remove the singleton dimension -> p × T

s_sum = mean(s_true, dims=2)
s_sum = dropdims(s_sum; dims=2) # p × T

# Plot
heatmap(s_sum, title="s_true summed/averaged over h", xlabel="Time", ylabel="j")
heatmap(γ_sum, title="γ summed/averaged over h", xlabel="Time", ylabel="j")

auc

#Pkg.add("GraphicalLasso")
using LinearAlgebra, GraphicalLasso, Random, Distributions
 

# Select the optimal tuning parameter from a range
λvalues = [0.1]
size(Y)
Ymat = reduce(hcat, Y)   # size: 40 × p
SIGMA = Ymat * Ymat' / T
optimalλ = tuningselect(SIGMA, obs, λvalues, verbose=false)
println("Optimal λ: ", optimalλ)

# Apply the graphical lasso algorithm
resultGLASSO = glasso(SIGMA, T, optimalλ)
resultGLASSO.W
# Apply the graphical lasso algorithm.

heatmap(resultGLASSO.W,  color=:viridis)

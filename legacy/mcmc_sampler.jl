# MCMC Sampler for birth-death factor graphical models
using Distributions
using MultivariateStats
using LinearAlgebra , TensorCore
using Random 
using Statistics
using StatsFuns 
using Plots
using StatsPlots  # Import StatsPlots for boxplot functionality
using Base.Threads
using SparseArrays
 
p, d, T = 10, 3, 100

NSIM, ISIM=50, 29
RV= zeros(T, NSIM) 

Random.seed!(ISIM)
# p=variables , d=latent dimension, T=time steps
Λ = randn(p, d, T) .* rand(Bernoulli(0.5), p, d, T) # array of latent factor loadings
Λ_sparse = [sparse(Λ[:, :, t]) for t in 1:T]  # each Λ[:, :, t] is a 2D matrix
 
 # Sparse matrix multiplication
Δ = rand(p) # Diagonal precision component, constant over time (Ω_t=Λ_tΛ_t'+Δ)
y= randn(T, p) # data
ϵ= randn(T, p) # errors
η= randn(T, d)# latent factors
 
Ω = Array{Float64, 3}(undef, p, p, T)  # Initialize precision tensors

# Latent state indicator (shocks)  
probs = [1,0,0]  # Example: 30% for -1, 40% for 0, 30% for 1
# Sample with probabilities
s = reshape(rand(Categorical(probs), T * d * p) .- 2, (p, d, T))
#note -2 maps 1:3 to -1:1

s[:, :, 1] .= 1  # Set all elements in the first slice to 1
s[:, :, T] .= -1  # Set all elements in the last slice to 0
 
# Big shock that acts on all factors
#s[:, 1:3, 20] .= 1
for t in 2:(T)  # Start from t=2 to handle t-1 case
    for j in 1:p
        for h in 1:d
            if s[j, h, t] == 0
                Λ_sparse[t][j, h]  = 0.0
            elseif s[j, h, t] == -1
                Λ_sparse[t][j, h] = Λ_sparse[t-1][j, h]
            elseif s[j, h, t] == 1
                Λ_sparse[t][j, h] = rand(Normal(0, 1))  # Draw from a normal distribution
            end
        end
    end
end

Δ_true = copy(Δ)
Ω_true = [Λ_sparse[t]*Λ_sparse[t]'+Diagonal(Δ[1:p]) for t in 1:T]   # Initialize precision tensors
heatmap(Ω_true[1])
  
 

using BenchmarkTools
#@benchmark compute_inv_ch(Ω_true[1]) better with cholesky
#@benchmark compute_inv(Ω_true[1])
t=1
for t in 1:T
    A=Ω_true[t] # Compute the precision matrix
    C = cholesky(Symmetric(A))
    cov_y = inv(Matrix(A))  # A^{-1}
   # cov_y = compute_inv_ch(A)
    y[t, :] = rand(MvNormal(Hermitian(cov_y)))
    Pt = I(d) + Λ_sparse[t]' * Diagonal(Δ[1:p]) * Λ_sparse[t]
    Pt = (Pt+Pt')/2
    Pt_dense = Matrix(Pt) 
    η[t, :] = rand(MvNormal(Pt_dense)) 
    K =  (Diagonal(Δ[1:p]) * Λ_sparse[t] * inv(Pt_dense))'
    ϵ[t, :] = η[t, :]' * K + y[t, :]'
end

  # Store the true Δ for comparison later
# latent states
heatmap(s[:, 1,:], title = "s[1]", xlabel = "time")
heatmap(s[:, 2,:], title = "s[2]", xlabel = "time")
heatmap(s[:, 3,:], title = "s[3]", xlabel = "time")
 
# Covariate matrix: three columns (intercepts, .... birth, death) for each time 1:T
F=1
XBD = zeros(T, F)
XBD[:, 1] .= 1.0
 
# Initialize BETA as a 4D tensor (d × 3 × 3 × F)
# factors (d), previous states (3), current states (3), features (F) (columns in XBD)
BETA = randn(d, 3,3,F) * 0.01
BETA[:,:,2:3,:].= -2 # Set the birth and death coefficients

# prior variance for Lambda
σλ =1 # Standard deviation for the prior distribution of λ_{hj,t}
Plam = σλ^2 * I(d)
invPlam = inv(Plam)
# Prior for error variance
a_delta, b_delta=2,2
 
 # Define the TensorPrecisionModel structure
mutable struct TensorPrecisionModel
    Λ::Vector{SparseMatrixCSC{Float64, Int}}# Dimensions: [p × d × T]
    Δ::Vector{Float64}     # Diagonal precision, length p
    x::Matrix{Float64}     # Observations: [T × p]
    Ω::Vector{SparseMatrixCSC{Float64, Int}}   # Precision tensors: [p × p x T]
    C::Vector{Cholesky{Float64, Matrix{Float64}}} 
    logliks::Vector{Float64}
    logdets::Vector{Float64}
    yΩys::Vector{Float64}
    # Workspace buffers
    #u_old::Vector{Float64}
    #u_new::Vector{Float64}
    #Ω_inv_u::Vector{Float64}
end



function TensorPrecisionModel(Λ_sparse::Vector{SparseMatrixCSC{Float64, Int}},
   Δ::Vector{Float64}, y::Matrix{Float64})
   logliks = zeros(T)
   logdets=zeros(T)
   yΩys = zeros(T)
   Ω  =   [sparse(zeros(p,p)) for t in 1:T]    
   C= [cholesky(zeros(p,p)+I(p))  for t in 1:T] # Cholesky factors
   for t in 1:T
       Ω[t] = sparse( Λ_sparse[t]) * sparse(Λ_sparse[t]') +  Diagonal(Δ[1:p])
      # Ω[t] = (Ω[t] + Ω[t]') / 2  # Symmetrize
       C[t] = cholesky(Matrix(Ω[t])) # Symmetrize Ωt
       logdets[t]=2*sum(log.(diag(C[t].L))) # Determinant of Ωt
     #  logdet( Ω[t] )
       yΩys[t]=y[t,:]' * Ω[t] * y[t,:]
       logliks[t] = 0.5*( logdets[t] -  yΩys[t] )
   end
   TensorPrecisionModel(Λ_sparse, Δ, y, Ω ,  C,logliks, logdets, yΩys)
end

model = TensorPrecisionModel(Λ_sparse, Δ, y)
  

# jitter
jitter=zeros(p, d)
jitter .=  rand(Normal(1,0.01),p, d) # Small random perturbation
model.Λ[1] =  sparse((Λ_sparse[1].*jitter ))  
 
for tt in 2:T
    jitter=zeros(p, d)
jitter .=  rand(Normal(1,0.01),p, d)  # Small random perturbation
model.Λ[tt] =  sparse((Λ_sparse[tt].*jitter ))  
end

s[:, 1:3, 2:(T-1)] .= -1 
 
heatmap(s[:, 1,:], title = "s[1]", xlabel = "time")
heatmap(s[:, 2,:], title = "s[2]", xlabel = "time")
heatmap(s[:, 3,:], title = "s[3]", xlabel = "time")
 # Set all elements in the first slice to 1
# Reinitialize the model with the new Λ
model = TensorPrecisionModel(model.Λ, Δ, y)   
 
  
#@benchmark likelihood_s_0(y[5,:], model,     5,2,1)
#@benchmark likelihood_s_0_opt(y[5,:], model, 5,2,1)
 
 
function compute_transition_probs(X::Matrix{Float64}, BETA::Array{Float64,4})
    T, _ = size(X)
    d, _, _, _ = size(BETA)
    #  factors (d), previous states (3), current states (3), features (F) (columns in XBD)
    P = zeros(T, d, 3, 3)  # States are now properly mapped to 1,2,3
    
    for t in 1:T
        x = X[t, :]
        for h in 1:d
            for w in 1:3  # previous state (1=-1, 2=0, 3=1)
                # Get linear predictors for all possible next states
                predlin = [dot(BETA[h, w, v, :], x) for v in 1:3]
                
                # Apply constraints:
                # If current state is 0 (w=2), cannot stay at 0
                if w == 2
                    predlin[2] = -Inf  # Make probability of staying at 0 zero
                end

            # QUESTION : AVOID THE -1 -> -1 Transition???
            # If current state is 1 (w=3), cannot transition to 1
               # if w == 3 # If the previous state is 1, set the probability of transitioning to state 1 to 0
                 #   predlin[3] = -Inf 
                     # Renormalize to ensure probs sum to 1
                #end
                probs = softmax(predlin)
                P[t, h, w, :] .= probs
                
            end
        end
    end
    return P
end
# Compute transition probabilities
P = compute_transition_probs(XBD, BETA)

 
function update_λ_opt!(model::TensorPrecisionModel, t::Int, j::Int, h::Int, newval::Float64)        # [T × d]
    xt= model.x[t, :]  # Get the observation at time t
    uold=copy(model.Λ[t][:,h])
    unew=copy(uold)
    unew[j] = newval
    C=copy(model.C[t])  
    FLAG=0
try
    C = lowrankupdate(C,   Vector(unew))     # add   u_new u_new'
    C = lowrankdowndate(C, Vector(uold))     # remove u_old u_old'
catch e
    @warn "Cholesky low-rank step failed: $e"
    println("flag")
    FLAG=1
end
    if FLAG== 0
      # Efficient determinant update using matrix determinant lemma
      CMAT=Matrix(C)  # Convert to dense matrix for further calculations
      # 2×2 matrix
      ytT_u = dot(xt, unew)
      ytT_v = dot(xt, uold)
      yupdatey=ytT_u^2 - ytT_v^2 
      model.yΩys[t] += yupdatey # Update the quadratic form
      model.C[t] = C  # Update the Cholesky factor in the model
      update=unew*unew' - uold*uold'  # Rank-2 update
      model.Ω[t] += update
      model.logliks[t]=   0.5 * (2*sum(log.(diag(CMAT)))-model.yΩys[t])
      model.Λ[t][j,h]=newval
      else
        # If downdate fails, fall back to the old method    
        Ωt = copy(model.Ω[t])
        Ωt += unew * unew' - uold * uold'  # Rank-2 update    
        model.Ω[t] = Ωt
        yt_Ωt_yt = xt' * Ωt * xt
        model.yΩys[t] = yt_Ωt_yt
        model.Λ[t][j,h]=newval
        model.logliks[t]= 0.5 * (logdet(cholesky(Ωt)) - yt_Ωt_yt)
     end
end
  
  
 
function log_posterior_beta(BETA::Array{Float64,4}, X::Matrix{Float64}, 
    s::Array{Int64,3}, prior_std::Float64=1.0)
    P = compute_transition_probs(X, BETA)
    loglik = 0.0
    # Likelihood term using actual transitions in s
    @inbounds for t in 2:size(X,1), h in 1:size(s,2),  j in 1:size(s,1)
    prev_state =   s[j,h,t-1]+2 # Previous state (1-3)
    curr_state =   s[j,h,t]+2  # Current state
    loglik += log(P[t,h,prev_state,curr_state] )
    end
    # Prior term (Gaussian)
    logprior = -sum(BETA.^2) / (2*prior_std^2)
return loglik + logprior
end

# Function to update BETA using a random walk Metropolis-Hastings step
# This function perturbs one element of BETA at a time
function update_beta(BETA_current::Array{Float64,4}, X::Matrix{Float64},
   s::Array{Int64,3}, proposal_std::Float64=0.1)
 
current_logpost = log_posterior_beta(BETA_current, X, s)
acc=zeros(d*3*3*F) # Accumulator for accepted proposals 
accum=1
# Random walk proposal (only perturb one element per update)
for h in 1:d
    for w in 1:3
        for v in 1:3
            for f in 1:F
                BETA_proposed = copy(BETA_current)
                BETA_proposed[h, w, v, f] += randn() * proposal_std
                proposed_logpost = log_posterior_beta(BETA_proposed, X, s)
                log_accept_ratio = proposed_logpost - current_logpost

                if log(rand()) < log_accept_ratio
                    BETA_current[h,w,v,f] = BETA_proposed[h,w,v,f]
                     current_logpost = proposed_logpost
                    acc[accum]=1
                end
                accum=accum+1
            end

        end
    end
   
end
return BETA_current, mean(acc)
end


   
function update_eta_eps(eta::Matrix{Float64}, ϵ::Matrix{Float64},
    y::Matrix{Float64}, model::TensorPrecisionModel)
    
    invDelta = Diagonal(1.0 ./ model.Δ[1:p])  # not 1/model.Δ
    d=size(model.Λ[1])[2]# Inverse of the diagonal matrix
    I_d = I(d)  # Identity matrix of size d
    @inbounds for t in 1:T
        Pt =  (Matrix(I_d + model.Λ[t]' * invDelta * model.Λ[t])) # Compute Pt
        Pt = (Pt + Pt') / 2  # Ensure symmetry
         # Compute the inverse of Pt
        CPt = cholesky(Pt)  # Cholesky decomposition of Pt
        invCPt = inv(CPt) 
        # Sample from standard Normal
        η_sample = rand(MvNormal(zeros(d), I_d))  # Sample from standard Normal
        # Scale by Cholesky factor
        ηt = CPt.L * η_sample  # Scale by the Cholesky factor of Pt
        eta[t, :] = ηt  # Update eta for time step t
        # Compute K = invPt * Λt' * invDelta (transposed)
        invPt = inv(Pt) 
        K = (invDelta * model.Λ[t] * invPt)'  # Shape: d × p
        # set the errors deterministially 
        ϵ[t, :] = y[t, :] - (ηt' * K)'  # Update eps using y and η
    end
return eta, ϵ
end
  

iter=1
 

for h in 1:d, v in 1:3, w in 1:3 BETA[h,v,w,:] .= [-1 ]/2 end
for h in 1:d, w in 1:1 BETA[h,1,w,:] .= [1 ]/2 end
for h in 1:d, w in 1:1 BETA[h,:,1,:] .= [1,0,0]/2 end
P=compute_transition_probs(XBD, BETA)
 
 
σλ=1.0
s
maxiter=2000
# Initialize Λ with the true values
Ω_hat=zeros(p,p,T)  # Initialize Ω_hat with the same structure as Ω
BETA_mean=copy(BETA)

n_store = Int(maxiter)  # Number of samples to store
Λ_samples = Array{Float64, 4}(undef, p, d, T, n_store)
s_samples = Array{Int64, 4}(undef, p, d, T, n_store)
loglik_store = Vector{Float64}(undef, n_store)


acc_beta= zeros(maxiter)  # Store acceptance rates for BETA updates
  # Current BETA for updates
  
 
using Profile
t=1
for iter in 1:maxiter
    # Update η and ϵ
    η, ϵ = update_eta_eps(η, ϵ, y, model)
    #  time step t from 1 to T-1
    t = mod(t, T-1) + 1# t=rand(1:(T-1))#  select a time step t
    println(t)
        if t == 1  # at time t=1  s[j,h,1] == 1 fixed
            @inbounds for h in 1:d 
                @inbounds for j in 1:p
                    s_prop =   1 # Always 1 at t=1
                    old_val = model.Λ[1][j, h]
                    prop_val=  old_val + rand(Normal(0, 0.15)) 
                    
                    modelprop= deepcopy(model)
                    loglik_proposed=0;loglik_old=0
                    # Compute log-prior and log-likelihoods for the proposed and old values
                    logpriorlambda_old= logpdf(Normal(0, σλ/d), old_val)
                    logpriorlambda_proposed= logpdf(Normal(0, σλ/d), prop_val)
                 #  update_λ_opt!(model, 1, j, h, old_val)  # n)  # this is already stored
                    update_λ_opt!(modelprop, 1, j, h, prop_val) 
                    loglik_proposed+= modelprop.logliks[1]
                    loglik_old+=  model.logliks[1]
                    τ= 2
                    while τ <=T  && s[j,h,τ] == -1
                        #update_λ_opt!(model, τ, j, h, old_val ) # no need to update the old 
                        update_λ_opt!(modelprop, τ, j, h, prop_val)
                        loglik_proposed+=  modelprop.logliks[τ]
                        loglik_old+= model.logliks[τ]
                        τ += 1
                    end
                    if log(rand())<  (loglik_proposed +0+logpriorlambda_proposed - 
                            loglik_old -logpriorlambda_old)
                            # Accept the proposed value
                            model=modelprop
                            println("Accepted new value for λ_{", j, ",", h, ",", t, "}: ", prop_val)
                    else    # reject
                            println("Rejected new value for λ_{", j, ",", h, ",", t, "}: ", prop_val)
                    end    #  end of accept-reject step for lambda proposed
                end # end loop in j
            end # end loop in h
        else # if t is not 1: update s
            for h in 1:d
                for j in 1:p
                    vals = filter(!=(-1), s[j,h,1:(t-1)])
                    future_state = s[j,h,t+1] + 2
                    prev_state   = s[j,h,t-1] + 2
            
                    # Identify times τ where this entry remains undefined
                    τs = [t]

# Forward in time
τ = t + 1
while τ <= T && s[j,h,τ] == -1
    push!(τs, τ)
    τ += 1
end

# Backward in time
τ = t-1
push!(τs, τ)
while τ >1 && s[j,h,τ] == -1
    push!(τs, τ-1)
    τ -= 1
end

                    
                    τs = unique(sort(τs))
                    # ---- Step 1: proposal distribution for s[j,h,t] ----
                    pi = P[t+1,h,:,future_state] .* P[t,h,prev_state,:]
            
                    # Apply constraints (e.g. forbid 0→0)
                    if model.Λ[t-1][j,h] == 0.0
                        pi[2] = 0.0  # Forbid transition to 0 if previous λ was 0
                    end
                    if s[j,h, t+1] == 0.0
                        pi[2] = 0.0  # Forbid transition to 0 if previous λ was 0
                    end
            
                    # Normalize to valid probabilities
                    probs = pi ./ sum(pi)
            
                    # ---- Step 2: sample proposal once ----
                    s_prop = rand(Categorical(probs)) - 2  # map 1:3 → -1:1
            
                    # ---- Step 3: build proposed model ----
                    modelprop = deepcopy(model)
                    loglik_proposed = 0.0
                    loglik_old      = 0.0
                    logpriorlambda_proposed = 0.0
                    logpriorlambda_old      = 0.0
            
                    if s_prop == -1
                        previous_val = model.Λ[t-1][j,h]
                        for update_t in τs
                            update_λ_opt!(modelprop, update_t, j, h, previous_val)
                            loglik_proposed += modelprop.logliks[update_t]
                            loglik_old      += model.logliks[update_t]
                        end
            
                    elseif s_prop == 0
                        for update_t in τs
                            update_λ_opt!(modelprop, update_t, j, h, 0.0)
                            loglik_proposed += modelprop.logliks[update_t]
                            loglik_old      += model.logliks[update_t]
                        end
            
                    elseif s_prop == 1
                        old_val  = model.Λ[t][j,h]
                        prop_val = old_val + rand(Normal(0,.150))
                         
                        for update_t in τs
                            update_λ_opt!(modelprop, update_t, j, h, prop_val)
                            loglik_proposed += modelprop.logliks[update_t]
                            loglik_old      += model.logliks[update_t]
                        end
                    end

                    if s[j,h,t] > 0
                        # If the next state is not defined, use the current value
                         logpriorlambda_old=logpdf(Normal(0, σλ/d), model.Λ[t][j, h])
                           end
                      if s_prop > 0
                        # If the next state is not defined, use the current value
                        logpriorlambda_proposed=logpdf(Normal(0, σλ/d), modelprop.Λ[t][j, h])
                      end
                    # ---- Step 4: MH acceptance ----
                    q_old = probs[s[j,h,t] + 2]
                    q_new = probs[s_prop + 2]
            
                    if log(rand())<  (loglik_proposed + logpriorlambda_proposed - 
                        loglik_old - logpriorlambda_old)
                    
                        # Accept
                        model = modelprop
                        s[j,h,t] = s_prop
                    end
                end
            end
        end
        
  
 

    # Update Δ DONT' DO THAT FOR NOW AS 
    shape = a_delta + 0.5 * T
    ϵ² = ϵ .^ 2
    rates = b_delta .+ 0.5 .* vec(sum(ϵ², dims=1))
    Delta_old=copy(model.Δ)
   # model.Δ .= rand.(Gamma.(shape, 1.0 ./ rates))
      
   # model = TensorPrecisionModel(model.Λ, model.Δ , y)
    proposal_std = 0.5
    BETA, acc_beta[iter] = update_beta(BETA, XBD, s, proposal_std)
    BETA_mean= BETA_mean + BETA
    println("Iteration: ", iter, " - Log-likelihood: ", sum(model.logliks)) 
    
   
    s_samples[:, :, :, iter] = (s)
    loglik_store[iter] = sum(model.logliks)
    tt=1
    @inbounds for tt in 1:T
        Λ_samples[:, :, tt, iter] = (model.Λ[tt])
        Ω_hat[:, :, tt] =  Ω_hat[:, :, tt] +  model.Ω[tt]# Λ_samples[:, :, t, iter] *Λ_samples[:, :, t, iter]'+ Diagonal(Δ[1:p])
    end
    
end


 

plot(loglik_store, title = "Log-likelihood over iterations", xlabel = "Iteration", ylabel = "Log-likelihood")
 
using BenchmarkTools
#@benchmark update_λ_opt!(model, 4, 1, 2,0.3)    
#@benchmark update_λ!(model, update_t, j, h, old_val)    
mean(acc_beta    )                       
 
BETA_mean = BETA_mean ./ maxiter  # Average over iterations
 
heatmap(s[:, 1,:], title = "s[1]", xlabel = "time")
heatmap(s[:, 2,:], title = "s[2]", xlabel = "time")
heatmap(s[:, 3,:], title = "s[3]", xlabel = "time")
 
 
shat=s_samples[:, :, :, 1]
for t in 1:maxiter
    shat=shat+s_samples[:, :, :, t]
end

heatmap(shat[:, 1, :]/maxiter)
 
 # Create custom color scheme
state_colors = cgrad([:red, :white, :blue], [0.0, 0.5, 1.0], categorical=false)
#maxiter=20000
time_=1
 
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
Λ_true=Λ_sparse
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")
heatmap(Λ_true[time_]*Λ_true[ time_]'+Diagonal(Δ),   title =  "True Ω at t=$time_", clim=(-a_Ω_hat, a_Ω_hat),    c=state_colors )
 
time_=2
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")
heatmap(Λ_true[ time_]*Λ_true[ time_]'+Diagonal(Δ),  c=state_colors,
clim=(-a_Ω_hat, a_Ω_hat), title =  "True Ω at t=$time_",  )



time_=5
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")
heatmap(Λ_true[ time_]*Λ_true[ time_]'+Diagonal(Δ),  c=state_colors,
clim=(-a_Ω_hat, a_Ω_hat), title =  "True Ω at t=$time_",  )
 




iter


time_=10
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")


        heatmap(Λ_true[ time_]*Λ_true[ time_]'+Diagonal(Δ_true),  c=state_colors,
clim=(-a_Ω_hat, a_Ω_hat), title =  "True Ω at t=$time_" )


time_=15
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
 
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")
heatmap(Λ_true[ time_]*Λ_true[time_]'+Diagonal(Δ),  c=state_colors,
clim=(-a_Ω_hat, a_Ω_hat), title =  "True Ω at t=$time_",  )

time_=25
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")
heatmap(Λ_true[ time_]*Λ_true[time_]'+Diagonal(Δ_true),  c=state_colors,
clim=(-a_Ω_hat, a_Ω_hat), title =  "True Ω at t=$time_" )
 


time_=30
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")
heatmap(Λ_true[time_]*Λ_true[ time_]'+Diagonal(Δ_true),  c=state_colors,
clim=(-a_Ω_hat, a_Ω_hat), title =  "True Ω at t=$time_" )

time_=35
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")
heatmap(Λ_true[ time_]*Λ_true[time_]'+Diagonal(Δ_true),  c=state_colors,
clim=(-a_Ω_hat, a_Ω_hat), title =  "True Ω at t=$time_" )


time_=40
a_Ω_hat = maximum(abs, Ω_hat[:,:,time_]/maxiter)
# Plot with centered color gradient
heatmap(Ω_hat[:,:,time_]/maxiter,  c=state_colors,  clim=(-a_Ω_hat, a_Ω_hat), 
        title="Estimated Ω at t=$time_")
heatmap(Λ_true[time_]*Λ_true[time_]'+Diagonal(Δ_true),  c=state_colors,
clim=(-a_Ω_hat, a_Ω_hat), title =  "True Ω at t=$time_" )
 
 

 # RV Coefficient
rv_coeff = zeros(T)  # Initialize RV coefficient array
for time_ in 1:T
    true_prec=Matrix(Λ_true[ time_]*Λ_true[ time_]'+Diagonal(Δ_true))
    hat_prec=Matrix(Ω_hat[:,:,time_]/maxiter)
    rv_coeff[time_] = (tr((true_prec' * hat_prec))) / sqrt(tr(true_prec' * true_prec) * tr(hat_prec' * hat_prec))
 end
 rv_coeff


boxplot(rv_coeff, title = "RV Coefficient", ylabel = "RV Coefficient", xlabel = "Time")
 


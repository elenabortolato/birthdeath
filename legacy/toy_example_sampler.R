# Gibbs sampler for birth-death graphical model
#  


# library import 
if(1==1){
library(pgdraw)
library(MASS)
library(Rcpp)
library(statmod)
library(mnormt)
library(mvtnorm)
}
#polya gamma data augmentation
cppFunction('

List pgg_m_sigma_(const Eigen::Map<Eigen::MatrixXd>  & omega,
                       const Eigen::Map<Eigen::MatrixXd>  & X,
                       const Eigen::Map<Eigen::MatrixXd>  & invB,
                       const Eigen::Map<Eigen::VectorXd>  & KTkappaplusinvBtimesb){
  int n = X.rows();
  int p = X.cols();
  // The matrix A stores XT Omega X + B^{-1}, that is, Sigma^{-1}
  Eigen::MatrixXd A(p,p);
  for (int j1 = 0; j1 < p; j1 ++){
    for (int j2 = j1; j2 < p; j2 ++){
      A(j1,j2) = invB(j1, j2);
      for (int i = 0; i < n; i++){
        A(j1,j2) = A(j1,j2) + X(i,j1) * X(i,j2) * omega(i);
      }
      A(j2,j1) = A(j1,j2);
    }
  }
  Eigen::LLT<Eigen::MatrixXd> lltofA(A);
  Eigen::MatrixXd lower = lltofA.matrixL();
  Eigen::VectorXd x = lltofA.solve(KTkappaplusinvBtimesb);
  return List::create(Named("m")=x,
                      Named("Sigma_inverse") = A,
                      Named("Cholesky_inverse") = lower,
                      Named("Cholesky") = lower.inverse());
}', depends = "RcppEigen")


#initialize everything     
set.seed(1411)


if(1==1){  
  # dimension of the data
  p=50#ncol(y)
  d=3# fixed
  T=25#nrow(y)
  
  
  a_delta=1;b_delta=1# variance of error covariance 
  a_lam=1;b_lam=1# governing loadings variance
  
  #furthermore
  # covariate for the probability of Lambda=0|s=1
  q=1 # for now intercept
  
  # initialize precision matrix for lambda
  Plam=2*diag(d)
  invPlam=solve(Plam)
  
  # initailize precision matrix of the "errors"
  deltas = 0.1*rep(1,p)
  Delta=diag(deltas)
  invDelta=diag(1/diag((Delta)))
  eps=mvnfast::rmvn(T,sigma =diag(p), mu=rep(0,p))
  # this is a copy of errors useful for the gibbs sampler
  
  # we can omit this parameter- it is a constant
  p_constant =1 #2*exp(1)*log(p)/p      
  
  # initialize LAMBDA CONSTANT over time
  Lambda=array(rnorm(p*d,0,1)*rbinom(p*d,size = 1,prob = 10/p),
               dim = c(p,d,T))
  
  
  #latent states ("shocks") - just initialized
  s=array(rbinom(d*T,1,0.0 ),dim = c(d,T))
  
  #probability of Lambda!=0 | s=1  
  gamma=array(runif(d*T*p,0.999,1),dim = c(p,d,T)) # probabilities (x%*% alpha)
  #pre-store array used for backward-forward sampling
  fhnorm=array(NA,dim = c(2,d,T))
  Psi=array(rbinom(p*d,1,gamma),dim = c(p,d,T))
  #Lambda=Lambda_*Psi
  
  #"big shock that acts on all factors"
  Lambda[,1:3,1:12]=  0.0*Lambda[,1:3,1:12]+  
    array(runif(3*p, -0.5,+0.5)*rbinom(p*3,size = 1,prob = 15/p),dim = c(p,3,12))
  
  
  
  par(mfrow=c(2,3))
  image(Lambda[,1,T:1],main=expression(Lambda[1]), ylab = "time")
  image(Lambda[,2,T:1],main=expression(Lambda[2]))
  image(Lambda[,3,T:1],main=expression(Lambda[3]))
   
  Lambda_=Lambda
  Lambda_true=Lambda  #save the true matrix of loadings
  Delta_true=Delta # save true Delta
  
  # initialize factors and y 
  eta=matrix(0, nrow=T, ncol = d)
  y=eps # initialize
 
  # and y coherent to the loadings defined
  for(t in 1:T){
    Pt=diag(d)+ t(Lambda[,,t])%*%invDelta%*%(Lambda[,,t])
     
    eps[t,]=mvnfast::rmvn(1,sigma =invDelta,mu=rep(0, p))  
    eta[t,] =  mvnfast::rmvn(1,sigma =diag(d),mu=rep(0, d ))   
    
    vart=solve(tcrossprod(Lambda[,,t])+Delta)
    y[t,]= mvnfast::rmvn(1,sigma = vart,mu=rep(0, p ))   
  }
  
  precision= ginv(cov(y))
  image(precision) # this is the empirical precision matrix,  
  #not accounting for time variation
  # beta birth and death 
  BETA0 =matrix(rep(c(-4.61,-2.1,-2.2),d), ncol=d)
  BETA1 =matrix(rep(c(-1.93,-2.1,-2.2),d), ncol=d)
  rownames(BETA0)=rownames(BETA1)=c("beta0", "beta_birth", "beta_death")
  BETA0
  XBD=matrix(0, nrow=T, ncol=3)
  XBD[,1]=1
  #intercepts, n.biths n. deaths
  XBD[12,c(2:3)]=2 # IMPRECISE INFORMATION OF SHOCKS HAPPENING (13-14-15-16-20)
  XBD[2,c(2)]=8
  XBD[13,1]=1
  # additional covariates 
  # (for prediction of some lambda == 0 / individual information on p variables ! )
  # "dropping probability"
  X=array(1, dim=c(p,1,T))
  # coefficients for modeling the dropping probability (PR(\LAMDBA=0|S=1))--> psi
  BETA_DROP=matrix(rep(3,d), ncol=d)
  rownames(BETA_DROP)[1]="Intercept"
  

}

Pt=(diag(d)+t(Lambda[,,1])%*%invDelta%*%Lambda[,,1])
invPt=solve(Pt)

# note that we can write y as:
-invDelta%*%Lambda[,,1]%*%invPt%*%eta[1,]+eps[1,]
#or 
-solve(Delta+Lambda[,,1]%*%t(Lambda[,,1]))%*%(Lambda[,,1])%*%eta[1,]+eps[1,]

# range of the COV matrix
range(solve((tcrossprod(Lambda_true[,,t])+Delta_true)))

par(mfrow=c(1,2))
par(mar=c(3,3,3,3))
t=1
image(tcrossprod(Lambda_true[,,t]), main="t=1")
t=25
image(tcrossprod(Lambda_true[,,t]), main="t=20")



# the same
Lambda=jitter(Lambda_true,amount = 0.1)
eta=jitter(eta,amount = .1)
Lambda-Lambda_true
 
s[,1]=1 
s[,13]=1


eps# check that the Lambdas are coherent with the shocks
for(t in 1:T){
  W0=which(s[,t]==0)
  if(t>1) Lambda[,W0,t]=Lambda[,W0,t-1]
}

t=1
image(tcrossprod(Lambda[,,t]), main="t=1")
t=27
image(tcrossprod(Lambda[,,t]), main="t=27")


 
p#s[,4]=1
# I will store just the Monte Carlo mean of the 
# "time varying - precision matrix" Lambda_t Lambda_t^T 
# in this tensor
LL=array(0, c(p,p,T))
Deltas=Delta
 


T
llik_hist=rep(0,10000)

llik_lambda_hist=rep(0,10000)
llik_hist_ext=rep(0,10000)



iter=1


# ETA INDIPEND STAND 1.44
# U 1.86
# MY UPDATE 0.87


### gibbs 
system.time(
for( iter in iter:10000){
  
  # eta and epsilon
   
   
  for(t in 1:T){
    I=diag(d)
    Pt <- I + t(Lambda[,,t]) %*% invDelta %*% Lambda[,,t]
    invPt=solve(Pt)
    
    CPt=chol(Pt)
    invCPt=chol(invPt)
    A <- t(-invDelta %*% Lambda[,,t] %*% t(invCPt))
    var_update = solve(I+A%*%invDelta%*% t(A) )
    mean_update = var_update%*% ( A%*%invDelta%*%(y[t,])) 
    var_update=0.5*(var_update+t( var_update))
  #eta[t,]= rmvnorm(1,rep(0,d), var_update )+c(mean_update)
    eta[t,]=rmvnorm(1,rep(0,d),sigma = diag(d))  
    eta[t,]=eta[t,]-mean(eta[t,])
    eta[t,]=eta[t,]/sd(eta[t,]) 
    cvar=chol(var_update)
    eta[t,]=cvar%*%eta[t,]+c(mean_update)
    eta[t,]=(CPt)%*%eta[t,]
    #eta[t,]=rmvnorm(1,rep(0,d),sigma = Pt) 
    K=t(invDelta%*%(Lambda[,,t])%*%invPt) 
    eps[t,]=y[t,]+c(eta[t,])%*%K
  }
   
  
  tt=sample(2:(T-1),1)
  
  for(t in tt) {
    for (h in 1:d) {
      # Compute probabilities
      if (s[h,t-1]==1) p1 = plogis(sum(BETA1[,h] * XBD[t,]))  # P(s_t=1 | s_{t-1}=1)
      if (s[h,t-1]==0) p1 = plogis(sum(BETA0[,h] * XBD[t,]))  # P(s_t=1 | s_{t-1}=0)
      
      pfuture_1 = plogis(sum(BETA1[,h] * XBD[t+1,]))  # P(s_{t+1} = 1 | s_t = 1)
      pfuture_0 = plogis(sum(BETA0[,h] * XBD[t+1,]))  # P(s_{t+1} = 1 | s_t = 0)
      
      future_1=dbinom(s[h,t+1],1,pfuture_1)
      future_0=dbinom(s[h,t+1],1,pfuture_0)
      # Compute likelihood contributions
      likelihood_1 = exp(sum(dnorm(Lambda[,h,t], 0, sqrt(Plam[h,h]), log=T)))
      likelihood_0 = exp(sum(dnorm(Lambda[,h,t-1], 0, sqrt(Plam[h,h]), log=T)))
      
      # Since likelihood_0 involves equality, we simplify by considering it as part of the model structure
      is_equal <- all(Lambda[, h, t] == Lambda[, h, t - 1])
      likelihood_0 <- if (is_equal) 1 else 0
      # Compute posterior probability
      num = p1 * future_1 * likelihood_1
      den = num + (1 - p1) * (future_0) * likelihood_0
      
      # Sample s_t
      
      s[h,t] = rbinom(1, 1, num / den)
      if(den==0) s[h,t] =0
    }
  }
  if(2==2){ 
    deltas=rgamma(p, a_delta + 0.5*T, b_delta+0.5*colSums(eps^2,na.rm = T))
    Delta=diag(deltas)
    invDelta=diag(1/deltas)
  }

    
  t=1
  
    for(t in 1:T){
      w1=which(s[,t]==1)
      w0=which(s[,t]==0)
      
      if(t==1)  w1=c(1:d) #all have shock
      if(t==1)  w0=which(1:d>Inf) #none has previous value
     
      # define time slots between shocks for each latent factor (h=1...d)
      # and save them in a list of length d
      # i.e. if the series of latent (s_h,t ... s_h,t+4) is (1,0,0,0,1) 
      # time points t, t+1, t+2, t+3, before the next shock for s_h (t+4)
      
      if(t==T)  {seq1=list();seq1[1:length(w1)]=T}
      else seq1=lapply(w1, function (w) c(t:(t+min(T,which( s[w,(t+1):T]==1))[1]-1)))
      
      # if at time t at least one sh among d undergoes a shock then update
      if(length(w1)>0){
        Lambda[,w0 ,t]=Lambda[,w0 ,t-1]
        #trim to T
        seq1=lapply(1:length(seq1),function (w)  seq1[[w]][seq1[[w]]<=T])
        longest_len_seq=max(sapply(1:length(seq1), function(x) length(seq1[[x]])))
        longest_seq=which.max(sapply(1:length(seq1), function(x) length(seq1[[x]])))
        # which s_h (state) has the longest sequence of non-shocks
        len_seqs=(sapply(1:length(seq1), function(x) length(seq1[[x]])))
        
      
        
        Pt=diag(d)+ t(Lambda[,,t])%*%(invDelta )%*%(Lambda[,,t])
        invPt=solve(Pt)
         
        
        
    for(j in 1:p){
    #Compute residuals u_i^{(j)} = u_i - sum_{r ≠ j} λ_r * v_{r,i}
    lambda_j_excluded <- Lambda[-j, , t]
    eps_j_excluded <- eps[ ,-j , drop = FALSE]
    eps_j <- eps[t:(t+longest_len_seq-1) , j]
   # maxlen=length(t:(t+longest_len_seq-1))
    eta_j <- matrix(0, nrow = longest_len_seq, ncol = d)
    for (i in t:(t+longest_len_seq-1)) {
      lambda_j_excluded <- Lambda[-j, , t]
      eta_j[i-t+1 , ] <- eta[i, ] - (t(lambda_j_excluded) %*% eps_j_excluded[i , ])
    }
    
     
    
    w_j <- t(eta_j) %*% eps_j
    D_j_inv <- diag(1 / diag(Plam))
    epsj_sq <- sum(eps_j^2)
    
    post_prec <- D_j_inv + epsj_sq * diag(d)
    post_cov <- solve(post_prec)
    post_mean <- post_cov %*% w_j
    
    Lambda[j, ,t]  <- MASS::mvrnorm(1, mu = post_mean, Sigma = post_cov)[ ]
     
  }

      }
      else{Lambda[,,t]=Lambda[,,t-1]}
    }
   
  
  if(8==8){
    pgg_m_and_sigma <- function(omega, precomputed){
      return(pgg_m_sigma_(omega, precomputed$X, precomputed$invB, precomputed$KTkappaplusinvBtimesb))
    }
    pgg_precomputation <- function(Y, X, b, B){
      invB <- MASS::ginv(B)
      invBtimesb <- invB %*% (b)
      Ykappa <- matrix(Y - rep(0.5, length(Y)), ncol=1)
      XTkappa <- t(X) %*% Ykappa
      KTkappaplusinvBtimesb <- XTkappa + (invBtimesb)
      return(list(n=nrow(X), p=ncol(X), X=X, Y=Y, b=b, B=B,
                  invB=invB, invBtimesb=c(invBtimesb), KTkappaplusinvBtimesb=
                    c(KTkappaplusinvBtimesb)))
    } 
    # state 1
    
    #variance prior for betas
    scale_beta=1/T
    B=diag(3)*scale_beta# -0.0095/state$S*(state$n^-1) 
    
    for(h in 1:d){
      betas=(BETA0[,h])
      ones_indices <- which(s[h,] == 0)
      ones_indices=setdiff(ones_indices,T)
      which1=ones_indices+1 # which prediction can be made with this set of beta
      
      pred =  XBD[ones_indices,1:3]%*%betas
      logit_phi = plogis(pred)
      s_= matrix(1, nrow = length(which1), ncol = 1)
      logit_phi0 = logit_phi[which(s[h,which1]==0)]
      p_constant=1 
      which_zero = which(runif(length(logit_phi0))<
                           ((1-logit_phi0)/(1- logit_phi0*p_constant)))
      s_[  which(s[h,which1]==0)[which_zero] ] = 0
      #s_=as.vector(s_)
      XBDH=matrix(XBD[ones_indices,1:3], ncol=3)
      pgg_precomputed <- pgg_precomputation(s_, X = XBDH, b = betas, B = B)
      pgg_kernel <- function(beta){
        zs <- (pgg_precomputed$X %*% beta)
        w <- (pgdraw::pgdraw(1, zs))
        res <- pgg_m_and_sigma(w, precomputed = pgg_precomputed)
        beta_ <- unbiasedmcmc:::fast_rmvnorm_chol(1, res$m, res$Cholesky)[1,]
        return(list(beta = beta_))
      }
      beta_ <- pgg_kernel(betas)
      BETA0[,h]=beta_$beta
    }
    
    
    for(h in 1:d){
      betas=(BETA1[,h])
      ones_indices <- which(s[h,] == 1)
      ones_indices= setdiff(ones_indices,T)
      which1=ones_indices+1 # which prediction can be made with this set of beta
      
      pred =  XBD[ones_indices,1:3]%*%betas
      logit_phi = plogis(pred)
      s_= matrix(1, nrow = length(which1), ncol = 1)
      logit_phi0 = logit_phi[which(s[h,which1]==0)]
      p_constant=1 
      which_zero = which(runif(length(logit_phi0))<
                           ((1-logit_phi0)/(1- logit_phi0*p_constant)))
      s_[  which(s[h,which1]==0)[which_zero] ] = 0
      #s_=as.vector(s_)
      XBDH=matrix(XBD[ones_indices,1:3], ncol=3)
      
      pgg_precomputed <- pgg_precomputation(s_, X = XBDH, b = betas, B = B)
      pgg_kernel <- function(beta){
        zs <- (pgg_precomputed$X %*% beta)
        w <- (pgdraw::pgdraw(1, zs))
        res <- pgg_m_and_sigma(w, precomputed = pgg_precomputed)
        beta_ <- unbiasedmcmc:::fast_rmvnorm_chol(1, res$m, res$Cholesky)[1,]
        return(list(beta = beta_))
      }
      beta_ <- pgg_kernel(betas)
      BETA1[,h]=beta_$beta
    }
  }
  print(iter)
  
  llik=0
  llik_ext=0
  llik_lambda=0
  if(iter>1) Deltas=Deltas+Delta
  for(t in 1:T){
    LL_iter=Lambda[,,t]%*%t(Lambda[,,t]) 
    if(iter>0) LL[,,t]=LL[,,t]+LL_iter
    
    
    Peta=diag(d)+ t(Lambda[,,t])%*%invDelta%*%(Lambda[,,t])
    invPeta<-ginv(Peta)
    K=t(invDelta%*%(Lambda[,,t])%*%invPeta)
    
    
    llik_lambda= llik_lambda+sum(dnorm(Lambda[,,t], 0, sd =   Plam[1,1], log=T))
    llik_ext= llik_ext+dmvnorm(eta[t,], rep(0,d), sigma = Peta, log=T)
    
    
    PREC=(LL_iter+Delta)
    SIG=solve(PREC)
    SIG=0.5*(SIG+t(SIG))
    isSymmetric(SIG)
    llik=llik+dmvnorm(y[t,],mean = rep(0,p), sigma = (SIG),
                      log=T)
    
  }
  llik_lambda_hist[iter]= llik_lambda
  llik_hist[iter]=llik
  llik_hist_ext[iter]=llik_ext
  plot_fig=1
  #plot results (new Lambda)
  if(plot_fig==1){
   # par(mfrow=c(3,3))
  #  image(Lambda[,1,T:1],main=expression(Lambda[1]), ylab = "time")
  #  image(Lambda[,2,T:1],main=expression(Lambda[2]))
  #  image(Lambda[,3,T:1],main=expression(Lambda[3]))
    #show probability of shock
    #image(fhnorm[c(2,2),1,T:1], main=expression(pr[1]))
    #image(fhnorm[c(2,2),2,T:1], main=expression(pr[2]))
    #image(fhnorm[c(2,2),3,T:1], main=expression(pr[3]))
    image(s[,T:1], main="s")
  }
})
iter
par(mfrow=c(1, 1))
par(mar=c(3,3,3,3))
par(mfrow=c(3,1))

iter=iter-3

Lambda[,,14]

plot(llik_hist[10:iter], main="p=100 - eta ", type="l") 
#abline(v=5)
#points(llik_hist_ext[1:iter], col=2, pch=2) 
#points(llik_lambda_hist[1:iter], col=3, pch=3) 
#legend("center", col=c(1,2,3), pch=1:3,legend = c("loglik", "logprior eta | lambda", "logprior lambda"))

 
#abline(v=5)
iter=iter-10

plot(llik_hist_ext[1:iter], col=2, pch=2, type="l") 
plot(llik_lambda_hist[1:iter], col=3, pch=3, type="l") 
legend("topright", col=c(1,2,3), pch=1:3,
       cex=0.8,legend = c("loglik", "logprior eta | lambda", "logprior lambda"))

  



LL=LL/iter
Deltas=Deltas/iter
par(mar=c(3,3,3,3)-c(2.53,2.53,1,2.53))



par(mfrow=c(2,7))

pp=p
for(i in c(1,5,9, 13,17,21,25)){ 
# LL[,,i][which.min(LL[,,i])]=- max(LL[,,i])
  image(((LL[,,i]    )[1:pp,1:pp] ) ,axes=F,
        col = cm.colors(7),
        main = bquote(hat(Omega)[.(i)] )) }


for(i in c(1,5,9,13,17,21,25)){ 
#  LL[which.min(LL)]=- max(LL)
#  Lambda_true[,,i][which.min( Lambda_true[,,i])]=-max(( Lambda_true[,,i]))
  image(((tcrossprod(Lambda_true[,,i]) )[1:pp,1:pp]  )   ,
        axes=F,  col = cm.colors(20), main = bquote( Omega[.(i)]^{-1})) }

# Gibbs sampler for birth-death graphical model
#  


# library import 

library(pgdraw)
library(MASS)
library(Rcpp)
library(statmod)
library(mnormt)
library(mvtnorm)

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
set.seed(12525)


if(1==1){  
  # dimension of the data
  p=30#ncol(y)
  d=3# fixed
  T=25#nrow(y)
  
  
  a_delta=2;b_delta=1# variance of error covariance 
  a_lam=2;b_lam=1# governing loadings variance
  
  #furthermore
  # covariate for the probability of Lambda=0|s=1
  q=1 # for now intercept
  
  # initialize precision matrix for lambda
  Plam=2*diag(d)
  invPlam=solve(Plam)
  
  # initailize precision matrix of the "errors"
  deltas = 0.2*rep(1,p)
  Delta=diag(deltas)
  invDelta=diag(1/diag((Delta)))
  eps=mvnfast::rmvn(T,sigma = invDelta, mu=rep(0,p))
  # this is a copy of errors useful for the gibbs sampler
  
  # we can omit this parameter- it is a constant
  p_constant =1 #2*exp(1)*log(p)/p      
  
  # initialize LAMBDA CONSTANT over time
  Lambda=array(rnorm(p*d,0,1)*rbinom(p*d,size = 1,prob = 0.9999),
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
  Lambda[,1:3,1:12]=  0.2*Lambda[,1:3,1:12]+  
    array(runif(3*p, -0.5,+0.5)*rbinom(p*3,size = 1,prob = 0.99),dim = c(p,3,12))
  
  
  
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
    invPeta<-ginv(Pt)
    Peta
    eps[t,]=mvnfast::rmvn(1,sigma =invDelta,mu=rep(0, p))  
    eta[t,] =  mvnfast::rmvn(1,sigma =Pt,mu=rep(0, d ))   
    
    vart=solve(tcrossprod(Lambda[,,t])+Delta)
    y[t,]= mvnfast::rmvn(1,sigma = vart,mu=rep(0, p ))   
  }
  
  precision= ginv(cov(y))
  image(precision) # this is the empirical precision matrix,  
  #not accounting for time variation
  # beta birth and death 
  BETA0 =matrix(rep(c(-4.61,02.1,02.2),d), ncol=d)
  BETA1 =matrix(rep(c(-1.93,-02.1,02.2),d), ncol=d)
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
t=20
image(tcrossprod(Lambda_true[,,t]), main="t=20")


Lambda=jitter(Lambda_true,amount = 0.5)
 
Lambda-Lambda_true
 

# check that the Lambdas are coherent with the shocks
for(t in 1:T){
  W0=which(s[,t]==0)
  if(t>1) Lambda[,W0,t]=Lambda[,W0,t-1]
}

 

s[,13]=1
#s[,4]=1
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
### gibbs 
for( iter in iter:1000){
  
  # eta and epsilon
   
   
  for(t in 1:T){
    I=diag(d)
    Pt <- I + t(Lambda[,,t]) %*% invDelta %*% Lambda[,,t]
    invPt=solve(Pt)
    dim(Pt)
    CP=chol(invPt)
    A <- t(-invDelta %*% Lambda[,,t] %*% t(CP))
    var_update = solve(I+A%*%Delta%*% t(A) )
    mean_update = var_update%*% ( A%*%Delta%*%(y[t,])) 
    var_update=0.5*(var_update+t( var_update))
   # eta[t,]= rmvnorm(1,rep(0,d), var_update )+c(mean_update)
    eta[t,]=rmvnorm(1,rep(0,d),sigma = diag(d))  
    eta[t,]=eta[t,]-mean(eta[t,])
    eta[t,]=eta[t,]/sd(eta[t,])
    cvar=chol(var_update)
    eta[t,]=cvar%*%eta[t,]+c(mean_update)
    
    K=t(invDelta%*%(Lambda[,,t])%*%invPt) 
    eps[t,]=y[t,]+(eta[t,])%*%K
  }
   
    
    t=1
    for(j in 1:p){
    # Step 1: Compute residuals u_i^{(j)} = u_i - sum_{r ≠ j} λ_r * v_{r,i}
    lambda_j_excluded <- Lambda[-j, , t]
    eps_j_excluded <- eps[,-j , drop = FALSE]
    eta_j <- matrix(0, nrow = T, ncol = d)
    for (i in 1:T) {
      eta_j[i, ] <- eta[i, ] - t(lambda_j_excluded) %*% eps_j_excluded[i, ]
    }
    
    # Step 2: Compute v_j and w_j
    eps_j <- eps[, j]
    w_j <- t(eta_j) %*% eps_j
    
    # Step 3: Compute D_j
   # D_j_diag <- tau^2 * psi[j, ] * phi[j, ]^2
    D_j_inv <- diag(1 / diag(Plam))
    
    # Step 4: Compute posterior precision and covariance
    epsj_sq <- sum(eps_j^2)
    post_prec <- D_j_inv + epsj_sq * diag(d)
    post_cov <- solve(post_prec)
    
    # Step 5: Compute posterior mean
    post_mean <- post_cov %*% w_j
    
    # Step 6: Sample from multivariate normal
    Lambda[j,,t]  <- MASS::mvrnorm(1, mu = post_mean, Sigma = post_cov)
     
  }

  t
   
  
  for(t in 2:25) Lambda[,  ,t]=Lambda[,  ,t-1]
  
  
  
  
  
  print(iter)
  
  llik=0
  llik_ext=0
  llik_lambda=0
  if(iter>1) Deltas=Deltas+Delta
  for(t in 1:T){
    LL_iter=Lambda[,,t]%*%t(Lambda[,,t]) 
    if(iter>1) LL[,,t]=LL[,,t]+LL_iter
    
    
    Peta=diag(d)+ t(Lambda[,,t])%*%invDelta%*%(Lambda[,,t])
    invPeta<-ginv(Peta)
    K=t(invDelta%*%(Lambda[,,t])%*%invPeta)
    
    
    llik_lambda= llik_lambda+sum(dnorm(Lambda[,,t], 0, sd =   Plam[1,1], log=T))
    llik_ext= llik_ext+dmvnorm(eta[t,], rep(0,d), sigma = Peta, log=T)
    
    llik=llik+dmvnorm(y[t,],mean = rep(0,p), sigma = solve(LL_iter+Delta),
                      log=T)
    
  }
  llik_lambda_hist[iter]= llik_lambda
  llik_hist[iter]=llik
  llik_hist_ext[iter]=llik_ext
  plot_fig=1
  #plot results (new Lambda)
  if(plot_fig==1){
    par(mfrow=c(3,3))
    image(Lambda[,1,T:1],main=expression(Lambda[1]), ylab = "time")
    image(Lambda[,2,T:1],main=expression(Lambda[2]))
    image(Lambda[,3,T:1],main=expression(Lambda[3]))
    #show probability of shock
    #image(fhnorm[c(2,2),1,T:1], main=expression(pr[1]))
    #image(fhnorm[c(2,2),2,T:1], main=expression(pr[2]))
    #image(fhnorm[c(2,2),3,T:1], main=expression(pr[3]))
    image(s[,T:1], main="s")
  }
}

par(mfrow=c(1, 1))
par(mar=c(3,3,3,3))
par(mfrow=c(3,1))

iter=iter


plot(llik_hist[1:iter], main="p=100 ") 
#abline(v=5)
points(llik_hist_ext[1:iter], col=2, pch=2) 
points(llik_lambda_hist[1:iter], col=3, pch=3) 
#legend("center", col=c(1,2,3), pch=1:3,legend = c("loglik", "logprior eta | lambda", "logprior lambda"))

plot(llik_hist[1:iter], main="p=100") 
#abline(v=5)
plot(llik_hist_ext[1:iter], col=2, pch=2) 
plot(llik_lambda_hist[1:iter], col=3, pch=3) 
legend("topright", col=c(1,2,3), pch=1:3,legend = c("loglik", "logprior eta | lambda", "logprior lambda"))

par(mfrow=c(1,1))
plot(llik_hist[1:iter]+llik_hist_ext[1:iter]+llik_lambda_hist[1:iter], main="logposterior") 
abline(v=5)
points(llik_hist_ext[1:iter], col=2, pch=2) 
points(llik_lambda_hist[1:iter], col=3, pch=3) 
legend("bottomleft", col=c(1,2,3), pch=1:3,legend = c("loglik", "logprior eta | lambda", "logprior lambda"))

title("non problematic, p=10")

etaiter 
(eta)
Lambda

iter

LL=LL/iter
Deltas=Deltas/iter
par(mar=c(3,3,3,3)-c(2.53,2.53,1,2.53))

par(mfrow=c(1,2))


for(i in c(3)){ 
 # LL[which.min(LL)]=- max(LL)
  image(((LL[,,i]       ) ) ,axes=F,
        col = cm.colors(8),
        main = bquote(hat(Omega)[.(i)] )) }
for(i in c(3)){
#  LL[which.min(LL)]=- max(LL)
#  Lambda_true[,,i][which.min( Lambda_true[,,i])]=-max(( Lambda_true[,,i]))
  image(((tcrossprod(Lambda_true[,,i]) )   )   ,
        axes=F,  col = cm.colors(8), main = bquote( Omega[.(i)]^{-1})) }

range( tcrossprod(Lambda[,,1])- tcrossprod(Lambda_true[,,1]))
Delta
Delta_true
Lambda
Lambda_true[,,13]
'invPlam=  diag(1/diag(Plam))
EPS_LAMBDA=(sapply(1:25, function (t) eps[t,]%*%((Lambda_[, , t]))))

for ( t in 1:25){
  for (j in 1:p){
    # j fixed, i have Lambda_t (time fixed  and d components)
    epsj_act=eps[ ,j]
    
    epsj2I=c(crossprod(epsj_act))*diag(d)
    
    #EPS_LAMBDA=(sapply(1:25, function (t) eps[t,]%*%((Lambda[, , t]))))
    etaj= eta - t(EPS_LAMBDA) +t(sapply(1:25, function (t) eps[t,j]%*%((Lambda[j, , t]))))
   ## nb different than   t((eps[,j])*(Lambda[j, , ]))
    wj=(eps[ ,j])%*% etaj
    var_L=solve(invPlam+epsj2I)
    mu_L=var_L%*%t(wj)
    lj=mvnfast::rmvn(1, mu = mu_L ,sigma  =  var_L)
    #L_[j,,t] = lj
    #L[j,]=t(t(L_[j,])) * Theta[j,]  #sparse
    Lambda_[j, ,t]=lj
    #update the t-th row of the bog matrix
    EPS_LAMBDA[,t]=eps[t,]%*%((Lambda_[, , t]))
  }
}

Lambda=Lambda_*Psi
'

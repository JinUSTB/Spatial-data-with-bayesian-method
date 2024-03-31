########################################################
####----------------------------------####
####--Full conditional distributions--####
####----------------------------------####

#### LIKELIHOOD ####
likelihood <- function(beta0, beta, W){
  W.repmat = rep(W, times=ni)
  temp = beta0+Z%*%t(beta)+W.repmat
  ttemp = temp*temp1+temp2*exp(temp)
}

#### BETA0 ####
pro.beta0 <- function(beta0, beta, W){
  sum(likelihood(beta0, beta, W))-beta0^2/(2*sigma2)}

#### BETA ####
pro.beta <- function(beta0, beta, W){
  sum(likelihood(beta0, beta, W))-sum(beta^2)/(2*sigma2)}

#### W ####
# 1.W.ICAR----------------------------------
pro.Wi.ICAR <- function(beta0, beta, W, tau2, i){
  temp = likelihood(beta0, beta, W)
  k=0; for (j in 1:i){k = k+ni[j]}
  likelihood.i = sum(temp[(k-ni[i]+1):k])
  likelihood.i-(tau2*D[i,i]*(W[i]-sum(E[i,]*W)/D[i,i])^2/2)
}
# 2. W.GRF-----------------------------------
pro.Wi.GRF <- function(beta0, beta, W, tau2, i, Rmminv){
  temp = likelihood(beta0, beta, W)
  k=0; for (j in 1:i){k = k+ni[j]}
  likelihood.i = sum(temp[(k-ni[i]+1):k])
  likelihood.i-(tau2*Rmminv[i,i]*(W[i]+(sum(Rmminv[i,]*W)-Rmminv[i,i]*W[i])/Rmminv[i,i])^2/2)
} 


#### PHI ####
pro.phi <- function(W, phi, tau2){
  Rmm = exp(-(phi*Dmm)^nu)
  Rmminv = solve(diag(rep(1e-10,m))+Rmm)
  phi.temp = phi^(aphi-1)*exp(-bphi*phi)
  exp(-tau2*(W%*%Rmminv%*%t(W))/2)*phi.temp/sqrt(abs(det(Rmm)))
} 


########################################################
####----------------####
####--update steps--####
####----------------####

#### Update beta0 ####
beta0.update <- function(beta0,beta,W,sd.beta0){
    beta0.can <- rnorm(1, mean = beta0, sd = 2.4*sqrt(sd.beta0+10^(-10)))
    accept.prob <- exp(pro.beta0(beta0.can,beta,W)-pro.beta0(beta0,beta,W))
    if (runif(1) >= accept.prob) {
        beta0
    } else {
        beta0.can
    }
}


#### Update beta ####
beta.update <- function(beta0,beta,W,cov.beta){
    beta.can = rmvnorm(1, mean = beta, sigma = cov.beta)
    accept.prob <- exp(pro.beta(beta0, beta.can, W)-pro.beta(beta0, beta, W))
    if (runif(1) >= accept.prob) {
        beta 
    } else {
        beta.can
    }
}

#### 1. Update W.ICAR ####
Wi.update.ICAR <- function(beta0, beta, W, tau2, i){
  Wi.can <- rnorm(1, mean = W[i], sd = sqrt(1/(tau2*D[i,i])))
  W.can = replace(W, i, Wi.can)
  accept.prob <- exp(pro.Wi.ICAR(beta0, beta, W.can, tau2, i)-pro.Wi.ICAR(beta0, beta, W, tau2, i))
  if (runif(1) >= accept.prob) {
    W[i]
  } else {
    Wi.can
  }
}

#### 2. Update W.GRF ####
Wi.update.GRF <- function(beta0, beta, W, tau2, Rmminv, i){
  Wi.can <- rnorm(1, mean = W[i], sd = sqrt(1/(tau2*Rmminv[i,i])))
  W.can = replace(W, i, Wi.can)
  accept.prob <- exp(pro.Wi.GRF(beta0, beta, W.can, tau2, i, Rmminv)-pro.Wi.GRF(beta0, beta, W, tau2, i, Rmminv))
  if (runif(1) >= accept.prob){
    W[i]
  } else {
    Wi.can
  }
}



#### Update tau2 ####
tau2.update <- function(W, Sigma){
  atau.star = atau+qr(Sigma)$rank/2
  btau.star = btau+W%*%Sigma%*%t(W)/2
  tau2  = rgamma(1, shape = atau.star, scale = 1/btau.star)
}


#### Update phi ####
phi.update <- function(W, phi, tau2){
  phi.can <- rnorm(1, mean = phi, sd = 0.4)
  accept.prob <- pro.phi(W, phi.can, tau2)/pro.phi(W, phi, tau2)
  if (is.finite(accept.prob) == FALSE | runif(1) >= accept.prob) {
    phi
  } else {
    phi.can
  }
}



###################################################################################
# 1.ICAR
#### ------------------------- ####
#### --- mh.ICAR.algorithm --- ####
#### ------------------------- ####
mh.recurrent.spatial.ICAR <- function(n.sims, beta0, beta, W, tau2, L.para0, burnin){
  # start value ------------------------------------------------- 
  beta0.l=beta0; beta.l=beta; W.l=matrix(0,1,m); tau2.l=tau2
  # mcmc -------------------------------------------------------- 
  draws <- matrix(NA, nrow=n.sims, ncol=L.para0+2)
  cpoi <- matrix(NA, nrow = n, ncol = n.sims)
  sd.beta0 = 0.4
  cov.beta = 0.16*diag(length(beta))
  for (i in 1:n.sims){
    beta0.l <- beta0.update(beta0.l, beta.l, W.l,sd.beta0)
    beta.l <- beta.update(beta0.l, beta.l, W.l, cov.beta)
    for (j in 1:m){
      W[j] <- Wi.update.ICAR(beta0.l, beta.l, W.l, tau2.l, j)
    }
    W.l = W-mean(W)
    tau2.l <- tau2.update(W.l,Sigma)
    cpoi[,i] = likelihood(beta0.l, beta.l, W.l)
    draws[i,] <- cbind(beta0.l, beta.l, W.l, tau2.l)
    if (i>=5000){
        sd.beta0 = (sd(draws[1:i,1]))^2
        sd_beta = matrix(NA,nrow=length(beta),ncol=length(beta))
        for (k in (2:length(beta)+1)){
            for (j in 2:k){
                sd_beta[k-1,j-1]=sd_beta[j-1,k-1]=cov(draws[1:i,k],draws[1:i,j])}
        }
        cov.beta=2.4^2*(sd_beta+10^(-10)*diag(length(beta)))/2
    }
 }
  result = list(draws[(burnin + 1):n.sims, ],cpoi[,(burnin + 1):n.sims])
  return(result) 
}


# 2.GRF
#### ------------------------ ####
#### --- mh.GRF.algorithm --- ####
#### ------------------------ ####
mh.recurrent.spatial.GRF <- function(n.sims, beta0, beta, W, tau2, phi0, L.para0, burnin){
  ##start value-----------------------------------------------------------------
  beta0.l=beta0; beta.l=beta; W.l=matrix(0,1,m); tau2.l=tau2; phi.l=phi0;
  ##mcmc----------------------------------------------------------------------
  draws <- matrix(NA, nrow = n.sims, ncol = L.para0+3)
  cpoi <- matrix(NA, nrow = n, ncol = n.sims)
  sd.beta0 = 0.4
  cov.beta = 0.16*diag(length(beta))
  for (i in 1:n.sims){
    sill=0.9999
    Rmm.l = sill*exp(-(phi.l*Dmm)^nu)+diag(1-sill, m, m)
    Rmminv.l = solve(diag(rep(1e-10,m))+Rmm.l)
    
    beta0.l <- beta0.update(beta0.l, beta.l, W.l,sd.beta0)
    beta.l <- beta.update(beta0.l, beta.l, W.l,cov.beta)
    for (j in 1:m){
      W[j] <- Wi.update.GRF(beta0.l, beta.l, W.l, tau2.l, Rmminv.l, j)
    }
    W.l = W
    tau2.l <- tau2.update(W.l,Rmminv.l)
    phi.l <- phi.update(W.l, phi.l, tau2.l)
    cpoi[,i] = likelihood(beta0.l, beta.l, W.l)
    draws[i,] <- cbind(beta0.l, beta.l, W.l, tau2.l, phi.l)
    if (i>=500){
        sd.beta0 = (sd(draws[1:i,1]))^2
        sd_beta = matrix(NA,nrow=length(beta),ncol=length(beta))
        for (k in 2:(length(beta)+1)){
            for (j in 2:k){
                sd_beta[k-1,j-1]=sd_beta[j-1,k-1]=cov(draws[1:i,k],draws[1:i,j])}
        }
        cov.beta=2.4^2*(sd_beta+10^(-10)*diag(length(beta)))/2
    }
  }
  result = list(draws[(burnin + 1):n.sims, ],cpoi[,(burnin + 1):n.sims])
  return(result) 
}





# 3.nofrailty
#### ------------------------------ ####
#### --- mh.nofrailty.algorithm --- ####
#### ------------------------------ ####
mh.recurrent.spatial.nofrailty <- function(n.sims, beta0, beta, W, tau2, L.para0, burnin){
  ##start value-----------------------------------------------------------------
  beta0.l=beta0; beta.l=beta; W.l=matrix(0,1,m); 
  ##mcmc----------------------------------------------------------------------
  draws <- matrix(NA, nrow = n.sims, ncol = L.para0+1)
  cpoi <- matrix(NA, nrow = n, ncol = n.sims)
  sd.beta0 = 0.4
  cov.beta = 0.16*diag(length(beta))
  for (i in 1:n.sims){
    beta0.l <- beta0.update(beta0.l, beta.l, W.l,sd.beta0)
    beta.l <- beta.update(beta0.l, beta.l, W.l,cov.beta)
    cpoi[,i] <- likelihood(beta0.l, beta.l, W.l)
    draws[i,] <- cbind(beta0.l, beta.l, W.l)
    if (i>=100000){
        sd.beta0 = (sd(draws[1:i,1]))^2
        sd_beta = matrix(NA,nrow=length(beta),ncol=length(beta))
        for (k in 2:(length(beta)+1)){
            for (j in 2:k){
                sd_beta[k-1,j-1]=sd_beta[j-1,k-1]=cov(draws[1:i,k],draws[1:i,j])}
        }
        cov.beta=2.4^2*(sd_beta+10^(-10)*diag(length(beta)))/2
    }
  }
  result = list(draws[(burnin + 1):n.sims, ],cpoi[,(burnin + 1):n.sims])
  return(result) 
}



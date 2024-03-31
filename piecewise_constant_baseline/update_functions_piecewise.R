########################################################
####----------------------------------####
####--Full conditional distributions--####
####----------------------------------####

#### LIKELIHOOD ####
likelihood <- function(r0, beta, W){
  W.repmat = rep(W, times=ni)
  temp = Z%*%t(beta)+W.repmat
  temp3 = -Y%*%t(T.ind)%*%r0
  ttemp = temp*temp1+temp3*exp(temp)
}
  

#### BETA ####
pro.beta <- function(r0, beta, W){
    sum(likelihood(r0, beta, W))-sum(beta^2)/(2*sigma2)}

#### W ####
# 1.W.ICAR----------------------------------
pro.Wi.ICAR <- function(r0, beta, W, tau2, i){
    temp = likelihood(r0, beta, W)
    likelihood.i = sum(temp[(sum(ni[1:i])-ni[i]+1):(sum(ni[1:i]))])
    likelihood.i-(tau2*D[i,i]*(W[i]-sum(E[i,]*W)/D[i,i])^2/2)
}
# 2. W.GRF-----------------------------------
pro.Wi.GRF <- function(r0, beta, W, tau2, i, Rmminv){
    temp = likelihood(r0, beta, W)
    likelihood.i = sum(temp[(sum(ni[1:i])-ni[i]+1):(sum(ni[1:i]))])
    likelihood.i-(tau2*Rmminv[i,i]*(W[i]+(sum(Rmminv[i,]*W)-Rmminv[i,i]*W[i])/Rmminv[i,i])^2/2)
} 

#### PHI ####
pro.phi <- function(W, phi, tau2){
    Rmm = exp(-(phi*Dmm)^nu)
    Rmminv = ginv(diag(rep(1e-10,m))+Rmm)
    phi.temp = dgamma(phi, shape = aphi, scale = 1/bphi)
    exp(-tau2*(W%*%Rmminv%*%t(W))/2)*phi.temp/sqrt(abs(det(Rmm)))
} 


########################################################
####----------------####
####--update steps--####
####----------------####

#### Update r0 ####
r0.update <- function(beta,W){
    W.repmat = rep(W, times=ni)
    temp = Z%*%t(beta)+W.repmat
    r0scale = colSums(dN%*%t(T.ind))
    r0rate = T.ind%*%t(Y)%*%exp(temp)
    r0[1] = rgamma(1, shape = alambda+r0scale[1], rate = 25+r0rate[1])
    for (l in 2:(length(interval)-1)) {
        r0[l]  = rgamma(1, shape = alambda+r0scale[l], rate = alambda/r0[l-1]+r0rate[l])  
    }
    r0
}

#### Update beta ####
beta.update <- function(r0,beta,W,cov.beta){
    beta.can = rmvnorm(1, mean = beta, sigma = cov.beta)
    accept.prob <- exp(pro.beta(r0, beta.can, W)-pro.beta(r0, beta, W))
    if (runif(1) >= accept.prob) {
        beta 
    } else {
        beta.can
    }
}


#### Update W.ICAR ####
Wi.update.ICAR <- function(r0, beta, W, tau2, i){
    Wi.can <- rnorm(1, mean = W[i], sd = sqrt(1/(tau2*D[i,i])))
    W.can = replace(W, i, Wi.can)
    accept.prob <- exp(pro.Wi.ICAR(r0, beta, W.can, tau2, i)-pro.Wi.ICAR(r0, beta, W, tau2,i))
    if (runif(1) >= accept.prob) {
        W[i]
    } else {
        Wi.can
    }
}


#### Update W.GRF ####
Wi.update.GRF <- function(r0, beta, W, tau2, Rmminv, i){
    Wi.can <- rnorm(1, mean = W[i], sd = sqrt(1/(tau2*Rmminv[i,i])))
    W.can = replace(W, i, Wi.can)
    accept.prob <- exp(pro.Wi.GRF(r0, beta, W.can, tau2, i, Rmminv)-pro.Wi.GRF(r0, beta, W, tau2, i, Rmminv))
    if (is.finite(accept.prob) == FALSE | runif(1) >= accept.prob){
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




####################################################################################
#### ------------------------- ####
#### --- mh.ICAR.algorithm --- ####
#### ------------------------- ####
mh.recurrent.spatial.ICAR <- function(n.sims, beta, W, tau2, L.para0, burnin){
    # start value ------------------------------------------------- 
    beta.l = beta; W.l = W; tau2.l = tau2
    # mcmc -------------------------------------------------------- 
    draws <- matrix(NA, nrow=n.sims, ncol=L.para0+length(r0)+length(W))
    cpoi <- matrix(NA, nrow = n, ncol = n.sims)
    cov.beta = 0.16*diag(length(beta))
    for (i in 1:n.sims){
        r0.l <- r0.update(beta.l,W.l)
        beta.l <- beta.update(r0.l, beta.l, W.l, cov.beta)
        for (j in 1:m){
            W[j] <- Wi.update.ICAR(r0.l, beta.l, W.l, tau2.l, j)
        }
        W.l = W-mean(W)
        tau2.l <- tau2.update(W.l,Sigma)
        if (i>=100000){
          sd_beta = matrix(NA,nrow=length(beta),ncol=length(beta))
          for (k in 1:length(beta)){
            for (j in 1:k){
              sd_beta[k,j]=sd_beta[j,k]=cov(draws[1:i,k],draws[1:i,j])}
          }
          cov.beta=2.4^2*(sd_beta+10^(-10)*diag(length(beta)))/2
        }
    cpoi[,i] = likelihood(r0.l, beta.l, W.l)
    draws[i,] <- cbind(beta.l, tau2.l, t(r0.l), W.l)
}
    result = list(draws[(burnin + 1):n.sims, ],cpoi[,(burnin + 1):n.sims])
    return(result) 
}



####################################################################################
#### ------------------------ ####
#### --- mh.GRF.algorithm --- ####
#### ------------------------ ####
mh.recurrent.spatial.GRF <- function(n.sims, beta, W, tau2, phi0, L.para0, burnin){
    ##start value-----------------------------------------------------------------
    beta.l=beta; W.l=W; tau2.l=tau2; phi.l=phi0;
    ##mcmc----------------------------------------------------------------------
    draws <- matrix(NA, nrow = n.sims, ncol = L.para0+length(r0)+length(W))
    cpoi <- matrix(NA, nrow = n, ncol = n.sims)
    cov.beta = 0.16*diag(length(beta))
    for (i in 2:n.sims){
        Rmm.l = exp(-(phi.l*Dmm)^nu)
        Rmminv.l = solve(diag(rep(1e-10,m))+Rmm.l)
        r0.l <- r0.update(beta.l,W.l)
        beta.l <- beta.update(r0.l, beta.l, W.l, cov.beta)
        for (j in 1:m){
            W[j] <- Wi.update.GRF(r0.l, beta.l, W.l, tau2.l, Rmminv.l, j)
        }
        W.l = W
        tau2.l <- tau2.update(W.l,Rmminv.l)
        phi.l <- phi.update(W.l, phi.l, tau2.l)
        if (i>=100000){
          sd_beta = matrix(NA,nrow=length(beta),ncol=length(beta))
          for (k in 1:length(beta)){
            for (j in 1:k){
              sd_beta[k,j]=sd_beta[j,k]=cov(draws[1:i,k],draws[1:i,j])}
          }
          cov.beta=2.4^2*(sd_beta+10^(-10)*diag(length(beta)))/2
        }
        cpoi[,i] = likelihood(r0.l, beta.l, W.l)
        draws[i,] <- cbind(beta.l, tau2.l, phi.l, t(r0.l), W.l)
    }
    result = list(draws[(burnin + 1):n.sims, ],cpoi[,(burnin + 1):n.sims])
    return(result) 
}



# 4.nofrailty
#### ------------------------------ ####
#### --- mh.nofrailty.algorithm --- ####
#### ------------------------------ ####
mh.recurrent.spatial.nofrailty <- function(n.sims, beta, W, tau2, L.para0, burnin){
    ##start value-----------------------------------------------------------------
    beta.l=beta; W.l=matrix(0,1,m); 
    ##mcmc----------------------------------------------------------------------
    draws <- matrix(0, nrow = n.sims, ncol = L.para0)
    var1 = var2 = matrix(0.16, nrow=n.sims, ncol=1)
    COV = m1 = m2 = matrix(0, nrow=n.sims, ncol=1)
    for (i in 2:n.sims){
        r0.l <- r0.update(beta.l,W.l)
        beta.l <- beta.update(r0.l, beta.l, W.l, var1[i], var2[i], COV[i])
        draws[i,] <- cbind(r0.l, beta.l)
        if (i>=100){
            var1[i+1] = var(draws[1:i,1])
            var2[i+1] = var(draws[1:i,2])
            COV[i+1] = cov(draws[1:i,1],draws[1:i,2])
        }
    }
    return(draws[(burnin + 1):n.sims, ]) 
}




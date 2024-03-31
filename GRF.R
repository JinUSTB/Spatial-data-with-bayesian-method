# data example
#### load data--------------------------------
#GRF_data = as.matrix(read.table(file ="GRFdata.csv", header=TRUE, sep =","))

#### data we need-----------------------------
#locations0 = GRF_data[, c("Long","Lat")]; is locations for all subjects, with longitude and latitude
#and dim(locations0)=(n,2)
#n = length(all subjects) is the number of subjects 
#C = GRF_data[,"C"] is the censoring time for all subjects and dim(C)=(n,1)
#Z = GRF_data[,c("Z1","Z2","Z3")] is p-dimensional covariates, and dim(Z)=(n,3)
#T_hap = GRF_data[,2:41] is recurrent time for all subjects, and here dim(T_hap)=(n,40)

#*Tips: Please let the subjects from the same longitude and latitude next to each other like the example data



####################################################################################
####----------------####
#### ---Packages--- ####
####----------------####
library(survival)
library(coda)
library(fields)
library(maps)
library(ggplot2)
library(MASS)
library(mvtnorm)

####################################################################################
#### -------------------- ####
#### ----input data------ ####
#### -------------------- ####

## read the example data
GRF_data = as.matrix(read.table(file ="GRFdata.csv", header=TRUE, sep =","))
## data we need
locations0 = GRF_data[, c("Long","Lat")];
C = GRF_data[,"C"]
Z = GRF_data[,c("Z1","Z2","Z3")]
T_hap = GRF_data[,2:41]

####################################################################################
#### --------------------- ####
#### ----generate data---- ####
#### --------------------- ####
## Generate spatial data
locations = GRF_data[!duplicated(locations0), c("Long","Lat")];
cor.dist = function(x1, x2) rdist.earth(x1, x2, miles=FALSE, R=6378.137)
Dmm = rdist(locations, locations)
m = nrow(Dmm)
ni = matrix(0,1,m)
for (i in 1:m){
  longsame = as.matrix(locations0[,1]%in%locations[i,1])
  latsame = as.matrix(locations0[,2]%in%locations[i,2])
  ni[i] = sum(longsame*latsame)  
} 

## Generate recurrent event data
n=length(C)
ind = (T_hap <= matrix(C,nrow=n,ncol=ncol(T_hap)))
T = unique(sort(T_hap*ind))[-1]              #all observed time in a row 
dN = matrix(NA, nrow=n, ncol=length(T))
for (i in 1:n) {dN[i,]=is.element(t(T), T_hap[i,])}
Y = matrix(C,n,length(T))>matrix(T,n,length(T))
dt = matrix(diff(append(0,T)),1,length(T))
temp1 = matrix(rowSums(dN),n,1)
temp2 = matrix(-rowSums(Y*matrix(dt,nrow=n,ncol=length(T), byrow=TRUE)),n,1)


####################################################################################
#### -------------------------- ####
#### --- Iterative solution --- ####
#### -------------------------- ####
## load functions 
source('update_functions.R')

## setting hyper-parameters 
L.para0 = ncol(Z)+m
sigma2=1000; atau=0.001; btau=0.001;
tau2=1; nu=1; phi0=-log(0.001)/max(Dmm); aphi=2; bphi=(aphi-1)/phi0;
burnin=5000; nsave=5000; nsims = burnin+nsave
beta0=0; beta=matrix(0,1,ncol(Z)); W = matrix(0,1,m)

#### grf------------------------------
parameter_GRF <- mh.recurrent.spatial.GRF(nsims, beta0, beta, W, tau2, phi0, L.para0, burnin)
mh.draws_GRF = parameter_GRF[[1]]
like_il_GRF = parameter_GRF[[2]]
### nofrailty -------------------------
parameter_NO <- mh.recurrent.spatial.nofrailty(nsims, beta0, beta, W, tau2, L.para0, burnin)
mh.draws_NO = parameter_NO[[1]]
like_il_NO = parameter_NO[[2]]


####################################################################################
####################################################################################
####################################################################################
#### ------------------ ####
#### --- GRF result --- ####
#### ------------------ ####

## Parameters-----------------
beta0_GRF = mean(mh.draws_GRF[,1])
beta_GRF = t(as.matrix(colMeans(mh.draws_GRF[,2:(ncol(Z)+1)])))
W_GRF = colMeans(mh.draws_GRF[,(ncol(Z)+2):(L.para0+1)])
tau2_GRF = median(mh.draws_GRF[,(L.para0+2)])
phi_GRF = mean(mh.draws_GRF[,(L.para0+3)])
Para_GRF = c(beta_GRF,1/tau2_GRF,phi_GRF)
## sd-----------------
sd_beta_GRF = apply(mh.draws_GRF[,2:(ncol(Z)+1)],2,sd)
sd_tau2_GRF = sd(mh.draws_GRF[,(L.para0+2)])
sd_phi_GRF = sd(mh.draws_GRF[,(L.para0+3)])
sd_GRF = c(sd_beta_GRF,sd_tau2_GRF,sd_phi_GRF)
credibleup_GRF = as.matrix(Para_GRF)+1.96*as.matrix(sd_GRF)
crediblelow_GRF = as.matrix(Para_GRF)-1.96*as.matrix(sd_GRF)
result_GRF = cbind(as.matrix(Para_GRF), as.matrix(sd_GRF), crediblelow_GRF,credibleup_GRF)

## DIC_GRF-----------------
like_i_GRF = mean(colSums(like_il_GRF))
like_GRF = sum(likelihood(beta0_GRF, beta_GRF, W_GRF))
DIC_GRF = 2*like_GRF-4*like_i_GRF; 
pD_GRF = 2*like_GRF-2*like_i_GRF;
## LPML_GRF----------------
Like_GRF = exp(like_il_GRF)
omega_bar_GRF = rowMeans(1/Like_GRF)
omega_GRF = 1/Like_GRF
aa_GRF = omega_GRF<=t(matrix(sqrt(nsave)*omega_bar_GRF,nsave,n,byrow = TRUE))
bb_GRF = omega_GRF>=t(matrix(sqrt(nsave)*omega_bar_GRF,nsave,n,byrow = TRUE))
omega_tilde_GRF = aa_GRF*omega_GRF+bb_GRF*t(matrix(sqrt(nsave)*omega_bar_GRF,nsave,n,byrow = TRUE))
CPO_hat_GRF = rowSums(Like_GRF*omega_tilde_GRF)/rowSums(omega_tilde_GRF)
LPML_GRF = sum(log(CPO_hat_GRF))



####################################################################################
#### ------------------------ ####
#### --- nofrailty result --- ####
#### ------------------------ ####

## Parameters-----------------
beta0_NO = mean(mh.draws_NO[,1])
beta_NO = t(as.matrix(colMeans(mh.draws_NO[,2:(ncol(Z)+1)])))
## sd-----------------
sd_beta_NO = apply(mh.draws_NO[,2:(ncol(Z)+1)],2,sd)
credibleup_NO = as.matrix(t(beta_NO))+1.96*as.matrix(sd_beta_NO)
crediblelow_NO = as.matrix(t(beta_NO))-1.96*as.matrix(sd_beta_NO)
result_NO = cbind(t(beta_NO), as.matrix(sd_beta_NO),crediblelow_NO,credibleup_NO)

## DIC_NO-----------------
W_NO = rep(0,1,m)
like_i_NO = mean(colSums(like_il_NO))
like_NO = sum(likelihood(beta0_NO, beta_NO, W_NO))
DIC_NO = 2*like_NO-4*like_i_NO;
pD_NO = 2*like_NO-2*like_i_NO;
## LPML_NO----------------
Like_NO = exp(like_il_NO)
omega_bar_NO = rowMeans(1/Like_NO)
omega_NO = 1/Like_NO
aa_NO = omega_NO<=t(matrix(sqrt(nsave)*omega_bar_NO,nsave,n,byrow = TRUE))
bb_NO = omega_NO>=t(matrix(sqrt(nsave)*omega_bar_NO,nsave,n,byrow = TRUE))
omega_tilde_NO = aa_NO*omega_NO+bb_NO*t(matrix(sqrt(nsave)*omega_bar_NO,nsave,n,byrow = TRUE))
CPO_hat_NO = rowSums(Like_NO*omega_tilde_NO)/rowSums(omega_tilde_NO)
LPML_NO = sum(log(CPO_hat_NO))


####################################################################################
#### -------------------------- ####
#### ------ Print result ------ ####
#### -------------------------- ####
DIC_NO    #is the DIC value for nofrailty
pD_NO     #is the effective number of parameters
LPML_NO   #is the LPML value for nofrailty
result_NO #is the estimates for parameters (covariates effects) and their sd and 95%credible interval
DIC_GRF   #is the DIC value for nofrailty
pD_GRF    #is the effective number of parameters
LPML_GRF  #is the LPML value for nofrailty
result_GRF #is the estimates for parameters (covariates effects+tau2+phi) and their sd and 95%credible interval

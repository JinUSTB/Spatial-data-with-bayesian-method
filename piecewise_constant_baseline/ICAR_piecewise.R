####################################################################################
#### ------------------ ####
#### ----input data---- ####
#### ------------------ ####
#### Load data--------------------------------
#ICAR_data = as.matrix(read.table(file ="ICARdata.csv", header=TRUE, sep =","))

#### Data we need-----------------------------
#n = length(all subjects) is the number of subjects 
#C = ICAR_data[,"C"] is the censoring time for all subjects and dim(C)=(n,1)
#Z = ICAR_data[,c("Z1","Z2","Z3")] is covariates, and dim(Z)=(n,3)
#T_hap = ICAR_data[,2:41] is recurrent event time for all subjects, and here dim(T_hap)=(n,40)
#regions0 = ICAR_data[, "units"] is the area code, and dim(regions0)=(n,1)
#E = as.matrix(read.table(file ="adjacency.csv", header=TRUE, sep =",")) is the adjacency of any two area. 
#dim(E)=(m,m), where m is the number of different area, here m=10
#E[i,j] be 1 if area code i and j are neighbors and 0 otherwise. E[i,i]=0

#*tips: Please let the subjects that from the same area next to each other like the example data.


####################################################################################
#### -------------- ####
#### ---Packages--- ####
#### -------------- ####
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
ICAR_data = as.matrix(read.table(file ="ICARdata.csv", header=TRUE, sep =","))
## data we need
C = ICAR_data[,"C"]
Z = ICAR_data[,c("Z1","Z2","Z3")]
T_hap = ICAR_data[,2:41] 
regions0 = ICAR_data[, "units"];
E = as.matrix(read.table(file ="adjacency.csv", header=TRUE, sep =","))
#------------------------------------------------------------------

####################################################################################
#### ----------------------- ####
#### ----generate data------ ####
#### ----------------------- ####
## Generate spatial data
regions = ICAR_data[!duplicated(regions0), "units"];
D = diag(colSums(E))
Sigma = D-E
m = nrow(E)
ni = matrix(0,1,m)
for (i in 1:m){ni[i] = sum(regions0==regions[i])}


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


##interval for piecewise############
length.inter = 1
inter.temp = seq(0,max(T),by=length.inter)
interval = seq(0,length(inter.temp)*length.inter,by=length.inter)

T.ind = matrix(NA,nrow=length(interval)-1,ncol=length(T))
for (i in 2:length(interval)){
    T.ind[i-1,] = (T >= interval[i-1] & T <= interval[i])
}

r0 <- matrix(NA, nrow=length(interval)-1, ncol=1)


####################################################################################
#### -------------------------- ####
#### --- Iterative solution --- ####
#### -------------------------- ####
## load functions 
source('update_functions_piecewise.R')
set.seed(20240313)

#### default choice for hyper parameters ------------------------------
sigma2=1000; atau=0.001; btau=0.001; tau2=1; alambda=5
burnin=5000; nsave=5000; nsims = burnin+nsave
beta=matrix(0,1,ncol(Z)); 
W = matrix(1,1,m)
L_para0=ncol(Z)+1

#### ICAR ------------------------------
parameter_ICAR <- mh.recurrent.spatial.ICAR(nsims, beta, W, tau2, L_para0, burnin)
mh.draws_ICAR = parameter_ICAR[[1]]
like_il_ICAR = parameter_ICAR[[2]]

####################################################################################
####################################################################################
####################################################################################
#### ------------------- ####
#### --- ICAR result --- ####
#### ------------------- ####

## Parameters-----------------
beta_ICAR = matrix(colMeans(mh.draws_ICAR[,1:ncol(Z)]),1,ncol(Z))
tau2_ICAR = median(mh.draws_ICAR[,(ncol(Z)+1)])
r0_ICAR = matrix(colMeans(mh.draws_ICAR[,(ncol(Z)+2):(ncol(Z)+1+length(r0))]),1,length(r0))
W_ICAR = matrix(colMeans(mh.draws_ICAR[,(ncol(Z)+2+length(r0)):(L_para0+length(r0)+length(W))]),1,length(W))
Para_ICAR = c(beta_ICAR,1/tau2_ICAR)
## PSD-----------------
sd_beta_ICAR = apply(mh.draws_ICAR[,1:ncol(Z)],2,sd)
sd_tau2_ICAR = sd(1/mh.draws_ICAR[,(ncol(Z)+1)])
sd_ICAR = c(sd_beta_ICAR,sd_tau2_ICAR)
credibleup_ICAR = as.matrix(Para_ICAR)+1.96*as.matrix(sd_ICAR)
crediblelow_ICAR = as.matrix(Para_ICAR)-1.96*as.matrix(sd_ICAR)
result_ICAR = cbind(as.matrix(Para_ICAR), as.matrix(sd_ICAR), crediblelow_ICAR,credibleup_ICAR)

## DIC_ICAR-----------------
like_i_ICAR = mean(colSums(like_il_ICAR))
like_ICAR = sum(likelihood(t(r0_ICAR), beta_ICAR, W_ICAR))
DIC_ICAR = 2*like_ICAR-4*like_i_ICAR
pD_ICAR = 2*like_ICAR-2*like_i_ICAR
## LPML_ICAR----------------
Like_ICAR = exp(like_il_ICAR)
omega_bar_ICAR = rowMeans(1/Like_ICAR)
omega_ICAR = 1/Like_ICAR
aa_ICAR = omega_ICAR<=t(matrix(sqrt(nsave)*omega_bar_ICAR,nsave,n,byrow = TRUE))
bb_ICAR = omega_ICAR>=t(matrix(sqrt(nsave)*omega_bar_ICAR,nsave,n,byrow = TRUE))
omega_tilde_ICAR = aa_ICAR*omega_ICAR+bb_ICAR*t(matrix(sqrt(nsave)*omega_bar_ICAR,nsave,n,byrow = TRUE))
CPO_hat_ICAR = rowSums(Like_ICAR*omega_tilde_ICAR)/rowSums(omega_tilde_ICAR)
LPML_ICAR = sum(log(CPO_hat_ICAR))


####################################################################################
#### -------------------------- ####
#### ------ Print result ------ ####
#### -------------------------- ####
DIC_ICAR   #is the DIC value for nofrailty
pD_ICAR    #is the effective number of parameters
LPML_ICAR  #is the LPML value for nofrailty
result_ICAR #is the estimates for parameters (covariates effects+tau2) and their sd and 95%credible interval


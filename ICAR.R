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


####################################################################################
#### -------------------------- ####
#### --- Iterative solution --- ####
#### -------------------------- ####
## load functions 
source('update_functions.R')

#### default choice for hyper parameters ------------------------------
sigma2=1000; atau=0.001; btau=0.001; tau2=1 
burnin=5000; nsave=5000; nsims = burnin+nsave
beta0=0; beta=matrix(0,1,ncol(Z)); 
W = matrix(0,1,m)
L_para0=ncol(Z)+m

#### ICAR ------------------------------
parameter_ICAR <- mh.recurrent.spatial.ICAR(nsims, beta0, beta, W, tau2, L_para0, burnin)
mh.draws_ICAR = parameter_ICAR[[1]]
like_il_ICAR = parameter_ICAR[[2]]

### nofrailty -------------------------
parameter_NO <- mh.recurrent.spatial.nofrailty(nsims, beta0, beta, W, tau2, L_para0, burnin)
mh.draws_NO = parameter_NO[[1]]
like_il_NO = parameter_NO[[2]]


####################################################################################
####################################################################################
####################################################################################
#### ------------------- ####
#### --- ICAR result --- ####
#### ------------------- ####

## Parameters-----------------
beta0_ICAR = mean(mh.draws_ICAR[,1])
beta_ICAR = matrix(colMeans(mh.draws_ICAR[,2:(ncol(Z)+1)]),1,ncol(Z))
W_ICAR = colMeans(mh.draws_ICAR[,(ncol(Z)+2):(L_para0+1)])
tau2_ICAR = median(mh.draws_ICAR[,(L_para0+2)])
Para_ICAR = c(beta_ICAR,1/tau2_ICAR)
## PSD-----------------
sd_beta_ICAR = apply(mh.draws_ICAR[,2:(ncol(Z)+1)],2,sd)
sd_tau2_ICAR = sd(mh.draws_ICAR[,L_para0+2])
sd_ICAR = c(sd_beta_ICAR,sd_tau2_ICAR)
credibleup_ICAR = as.matrix(Para_ICAR)+1.96*as.matrix(sd_ICAR)
crediblelow_ICAR = as.matrix(Para_ICAR)-1.96*as.matrix(sd_ICAR)
result_ICAR = cbind(as.matrix(Para_ICAR), as.matrix(sd_ICAR), crediblelow_ICAR,credibleup_ICAR)

## DIC_ICAR-----------------
like_i_ICAR = mean(colSums(like_il_ICAR))
like_ICAR = sum(likelihood(beta0_ICAR, beta_ICAR, W_ICAR))
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
DIC_NO = 2*like_NO-4*like_i_NO
pD_NO = 2*like_NO-2*like_i_NO
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
DIC_NO     #is the DIC value for nofrailty
pD_NO      #is the effective number of parameters
LPML_NO    #is the LPML value for nofrailty
result_NO  #is the estimates for parameters (covariates effects) and their sd and 95%credible interval
DIC_ICAR   #is the DIC value for nofrailty
pD_ICAR    #is the effective number of parameters
LPML_ICAR  #is the LPML value for nofrailty
result_ICAR #is the estimates for parameters (covariates effects+tau2) and their sd and 95%credible interval

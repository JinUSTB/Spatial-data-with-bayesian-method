
##ICAR#################################


1. This folder contains the ICARdata.csv and adjacency.csv as the example data and the code for analyzing the data with ICAR frailty.
2. Open the file “ICAR_piecewise.R” with R and run.



####################################################################################
#### ------------------- ####
#### ----input data---- ####
#### ------------------- ####

#### load data--------------------------------
#ICAR_data = as.matrix(read.table(file ="ICARdata.csv", header=TRUE, sep =","))

#### data we need-----------------------------
N = length(all subjects) is the number of subjects 
C = ICAR_data[,"C"] is the censoring time for all subjects and dim(C)=(n,1)
Z = ICAR_data[,c("Z1","Z2","Z3")] is covariates, and here dim(Z)=(n,3)
T_hap = ICAR_data[,2:41] is recurrent event time for each subjects, and here dim(T_hap)=(n,40)
regions0 = ICAR_data[, "units"] is the area code, and dim(regions0)=(n,1)
E = as.matrix(read.table(file ="adjacency.csv", header=TRUE, sep =",")) is the adjacency of any two area. dim(E)=(m,m), where m is the number of different area, here m=10
E[i,j] be 1 if area code i and j are neighbors and 0 otherwise. E[i,i]=0


*Tips*: Please set the subjects that from the same area next to each other like the example data(ICARdata)

####################################################################################
#### ------------------- ####
#### ----output data---- ####
#### ------------------- ####

DIC_ICAR   #is the DIC value for nofrailty 
pD_ICAR    #is the effective number of parameters 
LPML_ICAR  #is the LPML value for nofrailty 
result_ICAR #is the estimates for parameters (covariates effects+tau2) and their sd and 95%credible interval





##GRF###################################################
1. This folder contains the GRFdata.csv as the example data and the code for analyzing the data with GRF frailty.
2. Open the file "GRF_piecewise.R" with R and run.


####################################################################################
#### ------------------- ####
#### ----input data---- ####
#### ------------------- ####
#### load data--------------------------------
#GRF_data = as.matrix(read.table(file ="GRFdata.csv", header=TRUE, sep =","))

#### data we need-----------------------------
#locations0 = GRF_data[, c("Long","Lat")]; is locations for all subjects, with longitude and latitude and dim(locations0)=(n,2)
#n = length(all subjects) is the number of subjects 
#C = GRF_data[,"C"] is the censoring time for all subjects and dim(C)=(n,1)
#Z = GRF_data[,c("Z1","Z2","Z3")] is p-dimensional covariates, and dim(Z)=(n,3)
#T_hap = GRF_data[,2:41] is recurrent time for all subjects, and here dim(T_hap)=(n,40)

#*Tips: Please let the subjects from the same longitude and latitude next to each other like the example data



####################################################################################
#### ------------------- ####
#### ----output data---- ####
#### ------------------- ####
DIC_GRF   #is the DIC value for nofrailty
pD_GRF    #is the effective number of parameters
LPML_GRF  #is the LPML value for nofrailty
result_GRF #is the estimates for parameters (covariates effects+tau2+phi) and their sd and 95%credible interval








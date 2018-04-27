rm(list=ls())
library(RcppArmadillo)
library(Rcpp)
source('K:\\Projects\\CEC_Eval\\Code\\CUSUM_remodel\\SRR_cpp\\20180427\\SRR_fixed_effect_lib.R')
source('K:\\Projects\\CEC_Eval\\Code\\CUSUM_remodel\\SRR_cpp\\20180427\\SRR_fixed_effect_lib_cpp.R')
sourceCpp("K:\\Projects\\CEC_Eval\\Code\\CUSUM_remodel\\SRR_cpp\\20180427\\SRR_profile.cpp")

Y_name <- "Y"
z_names <- c("z1","z2")

N<-100000
nF<-1000
size<-N/nF
Fac_name <- "Facility" 
Facility<-rep(1:nF,each=size)
z1<-rnorm(N)
z2<-rnorm(N)
beta<-c(1,-1)
gamma<-rnorm(nF)
gamma_long<-rep(gamma, each=size)
exp_eta<-exp(as.matrix(cbind(z1,z2))%*%beta+gamma_long)
p<-exp_eta/(1+exp_eta)
Y<-rbinom(N,size=1,prob=p)
data=data.frame(Y,z1,z2,Facility)
pt<-proc.time()
SRR_fixed_fit = SRR_fixed(data0 = data, Y_name, z_names, Fac_name, test = "Score", refer_fac = NULL,
                          criterion = 1e-8, max.iter = 1000, bound = 8, size_cut = 10)
proc.time()-pt

pt<-proc.time()
SRR_fixed_fit_cpp = SRR_fixed_cpp(data0 = data, Y_name, z_names, Fac_name, test = "Score", refer_fac = NULL,
                          criterion = 1e-8, max.iter = 1000, bound = 8, size_cut = 10)
proc.time()-pt
max(abs(SRR_fixed_fit_cpp$srr-SRR_fixed_fit$srr))
max(abs(SRR_fixed_fit_cpp$beta-SRR_fixed_fit$beta))
max(abs(SRR_fixed_fit_cpp$gamma-SRR_fixed_fit$gamma))
max(abs(SRR_fixed_fit_cpp$iterN-SRR_fixed_fit$iterN))#sometimes can be not identical due to slightly different convergence
max(abs(SRR_fixed_fit_cpp$Observed-SRR_fixed_fit$Observed))
max(abs(SRR_fixed_fit_cpp$Expected-SRR_fixed_fit$Expected))
identical(SRR_fixed_fit_cpp$Fnames,SRR_fixed_fit$Fnames)
save.image("test.Rdata")

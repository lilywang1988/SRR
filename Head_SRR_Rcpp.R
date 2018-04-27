library(RcppArmadillo)#@
library(Rcpp) #@
source('K:\\Projects\\CEC_Eval\\Code\\CUSUM_remodel\\SRR_cpp\\20180427\\SRR_fixed_effect_lib_cpp.R') #@
sourceCpp("K:\\Projects\\CEC_Eval\\Code\\CUSUM_remodel\\SRR_cpp\\20180427\\SRR_profile.cpp") #@
SRR_fixed_fit_cpp = SRR_fixed_cpp(data0 = read09, Y_name, z_names, Fac_name, test = "Score", refer_fac = NULL,
                                  criterion = 1e-8, max.iter = 1000, bound = 8, size_cut = 10) #@

#SEX-SPECIFIC FEMALE#

#In this script we fit univariate models to all traits. Using AIC to infer
#the best fit model, we use the resulting parameters in the KL divergence and
#MI calculation. Output from the script includes a tile plot of data availability
#exported as a pdf, MI as a function of age per trait exported as a pdf, .RDS
#files of the problem, x solution, y solution, aic, and raw MI values with ordinal
# and continuous separated. 

#Packages 
library(yada)
library(tidyverse)
library(doParallel)
library(foreach)
registerDoParallel(detectCores())

#clear workspace
rm(list=ls())

# Re-direct print statements to a text file for a permanent record of
# processing
sink("results/output_col_F.txt")

# Check that a results folder exists in the working directory
if(! ("results" %in% dir()) ) {
  stop("There is no 'results' folder in the working directory")
}

# The data directory is /results
data_dir <- "results"

#The "ID" that uniquely identifies this analysis (pooled-sex, collapsed stages):
analysis_name <- 'MI_female_col'

###############################################################################

########################
#     Initial Prep     #
########################

set.seed(695432)

# Import data based on var_info file and yada functions

var_info <-  yada::load_var_info('data/US_allvar_var_info.csv')
data_file <- 'data/SVAD_US_F.csv'
dat <- load_cp_data(data_file, var_info)
df <- dat$cp_df

###############################################################################

############################
##    Univariate Models   ##
############################

# Extract the main problem
#
# Subsequent steps rely on the analysis_name, "US", which should be a unique
# analysis "ID" for files in the results folder (data_dir).
main_problem <- dat$problem
save_problem(data_dir, analysis_name, main_problem)

# Load the main problem from file
problem0 <- readRDS(build_file_path(data_dir, analysis_name, "main_problem"))

# Build ordinal problems (main problem and cross-validation folds)
ord_prob_list <- build_univariate_ord_problems(data_dir,
                                               analysis_name,
                                               add_folds=F)

base_seed <- 264528
set.seed(base_seed)
seed_vect <- sample.int(1000000, length(ord_prob_list), replace=F)

# Solve the ordinal problems in parallel and save to user-defined directory
ord_success <-
  foreach::foreach(i=1:length(ord_prob_list), .combine=cbind) %dopar% {
    yada::solve_ord_problem(data_dir,
                            analysis_name,
                            ord_prob_list[[i]],
                            anneal_seed=seed_vect[i])
  }

# Build continuous problems (main problem and cross-validation folds)
cont_prob_list <- build_univariate_cont_problems(data_dir,
                                                 analysis_name,
                                                 add_folds=F)

# Solve the continuous problems in parallel and save to user-defined directory
cont_success <-
  foreach::foreach(i=1:length(cont_prob_list), .combine=cbind) %dopar% {
    yada::solve_cont_problem(data_dir, analysis_name, cont_prob_list[[i]])
  }

###############################################################################

#######################################
# AIC Model Selection and Mutual Info #
#######################################

#####################
#Necessary Functions#
#####################

##  Ord MI Function
calc_univ_ord_mi <- function(th_v, mod_spec, th_x, x0, xcalc) {
  # Calculate the mutual information for a univariate ordinal model given the
  # baseline age, x0, and a vector at which to calculate the prior and
  # posterior densitites, xcalc.
  fprior <- calc_x_density(xcalc,th_x)
  
  M <- mod_spec$M
  pv     <- rep(NA, M+1)
  kl_div <- rep(NA, M+1)
  for (m in 0:M) {
    x_post_obj  <- calc_x_posterior(m,
                                    th_x,
                                    th_v,
                                    mod_spec,
                                    xcalc)
    kl_div[m+1] <- calc_kl_div(x_post_obj, th_x)
    pv[m+1] <- calc_q(x0, th_v, m, mod_spec)
  }
  return(sum(pv*kl_div))
}

## Cont Kl function
calc_cont_kl_div_vect <- function(th_w, mod_spec, xcalc, wcalc) {
  # Calculat a vector of KL divergences for each entry in wcalc
  N1 <- length(wcalc)
  N2 <- length(xcalc)
  kl_div <- rep(NA, N1)
  
  # Calculate the prior
  fprior0 <- calc_x_density(xcalc,th_x)
  
  # Use a single for loop to calcualte the likelihood matrix (rather than
  # unwrapping matrices to make a single call to calc_neg_log_lik_cont).
  dx <- xcalc[2] - xcalc[1]
  for (n1 in 1:N1) {
    eta_vect <- calc_neg_log_lik_vect_cont(th_w, xcalc, rep(wcalc[n1], N2), 
                                           mod_spec)
    lik_vect <- exp(-eta_vect)
    fprior <- fprior0
    fpost <- fprior * lik_vect
    fpost <- fpost / dx / sum(fpost)
    ind_bad <- is.na(fprior) | is.na(fpost)
    fprior <- fprior[!ind_bad]
    fpost <- fpost[!ind_bad]
    ind_bad <- !is.finite(fprior) | !is.finite(fpost)
    fprior <- fprior[!ind_bad]
    fpost <- fpost[!ind_bad]
    ind_bad <- (fprior == 0) | (fpost == 0)
    fprior <- fprior[!ind_bad]
    fpost <- fpost[!ind_bad]
    kl_div[n1] <- sum(fpost * log2(fpost/fprior)) * dx
  }
  
  return(kl_div)
}

## Cont MI function
calc_univ_cont_mi <- function(th_w, mod_spec, x0, wcalc, kl_div) {
  # Calculate the mutual information for a univariate cont. model given the
  # baseline age, x0, and a vector of KL divergences with length(wcalc) that
  # was output by calc_cont_kl_div_vect.
  
  N <- length(wcalc)
  dw <- wcalc[2] - wcalc[1]
  pw <- calc_neg_log_lik_vect_cont(th_w, rep(x0, N), wcalc, mod_spec)
  pw <- exp(-pw)
  pw <- pw / sum(pw) / dw
  return(dw*sum(pw*kl_div))
}

#################
# Solve x for MI#
#################

#offset for the weibull fit
weib_offset <- 0.002

#fit mixed weibull to age data
weib_fit <- mixtools::weibullRMM_SEM(problem0$x + weib_offset,
                                     k=3,
                                     maxit=2000)
#save the results of the fit
theta_x <- list(fit_type='offset_weib_mix',
                fit=weib_fit,
                weib_offset=weib_offset)

#save rds to used later
saveRDS(theta_x,build_file_path(data_dir, analysis_name, "solutionx"))

#Combine AIC model selection with MI calculation
#Outputs the best model and stores the MI vector in a list "MI_ord" or "MI_cont"

# Read the problem in to be used
problem <- readRDS(build_file_path(data_dir, analysis_name, 'main_problem'))

################
# Ordinal MI   #
################

# Pull out ordinal first
j_ord <- problem$mod_spec$J

#vector to calculate prior and posterior densities
xcalc = seq(0,23,by=0.01)

#baseline age
x0 = seq(0,23,by=.1)

#create and empty list to store ordinal results (j = 44)
MI_ord <- setNames(vector("list", j_ord), problem$var_names[1:44])

# Loop over each ord var, select the best model,print best model,do MI based on
#best model and store results in MI_ord
for (i in 1:j_ord) {
  
  var_name <- problem$var_names[i]
  
  print(paste0("AIC model selection for ", var_name))
  
  aic_output <- yada::build_aic_output(data_dir, analysis_name, var_name,
                                       format_df=T, save_file=TRUE)
  
  c <- filter(aic_output, rank == 1)$model
  
  md <- paste0(data_dir,"/solutiony_",analysis_name,"_","ord","_","j_",i,
               "_",var_name,"_",c,".rds")
  
  file <- readRDS(md)
  
  th_v = file[[1]]
  
  mod_spec = file[[2]]
  
  print(c)
  
  MI_ord[[i]] <- foreach(n=1:length(x0), .multicombine = T,
                         .packages = c("yada","doParallel")) %dopar%{
                           calc_univ_ord_mi(th_v, mod_spec, theta_x, x0[n],xcalc)
                         }
  
}

################
# Continuous MI#
################

# Pull out # cont variables for loop below
k_cont <- problem$mod_spec$K

#vector to calculate density over w
wcalc <- seq(0, 500, by=1)

#weibull fit needed in calc_x_density
th_x <- readRDS(build_file_path(data_dir, analysis_name, "solutionx"))

#create an empty list to store continuous results (k = 18)
MI_cont <- setNames(vector("list", k_cont), problem$var_names[45:62])

# Loop over each cont var, select the best model with aic, print best model,  
#do KL div,do MI based on best model and kl results and store results in MI_cont
for (i in 1:k_cont) {
  
  var_name <- problem$var_names[45:62][i]
  
  print(paste0("AIC model selection for ", var_name))
  
  aic_output <- yada::build_aic_output(data_dir, analysis_name, var_name,
                                       format_df=T, save_file=TRUE)
  
  c <- filter(aic_output, rank == 1)$model
  
  md<-paste0(data_dir,"/solutiony_",analysis_name,"_","cont","_","k_",i,
             "_",var_name,"_",c,".rds")
  
  file <- readRDS(md)
  
  th_w = file[[1]]
  
  mod_spec = file[[2]]
  
  print(c)
  
  kl <- calc_cont_kl_div_vect(th_w = th_w, mod_spec = mod_spec, xcalc = xcalc, 
                              wcalc = wcalc)
  
  MI_cont[[i]] <- foreach(n=1:length(x0), .multicombine = T,
                          .packages = c("yada","doParallel")) %dopar%{
                            calc_univ_cont_mi(th_w = th_w, mod_spec = mod_spec, x0 = x0[n], 
                                              wcalc = wcalc, kl_div = kl)
                          }
  
}

# Plot the MI over baseline age
pdf(file.path(data_dir,"MI_col_F.pdf"))
for(i in 1:j_ord){
  var_name <- problem$var_names[i]
  plot(x0, unlist(MI_ord[[i]]), ylim=c(0,max(unlist(MI_ord[[i]]))),
       xlab="Age [years]",ylab="Mutual Information",
       type="l", lwd=2, col="grey", main = var_name) 
}
for(i in 1:k_cont){
  var_name <- problem$var_names[45:62][i]
  plot(x0, unlist(MI_cont[[i]]),ylim = c(0,max(unlist(MI_cont[[i]]))),
       xlab="Age [years]",ylab="Mutual Information",
       type="l", lwd=2, col="grey", main = var_name)
}
dev.off()

#Write MI results to an RDS to be used for later use
write_rds(MI_ord, "results/MI_ord_col_F.rds")
write_rds(MI_cont,"results/MI_cont_col_F.rds")

# Close all clusters used for paralle processing
stopImplicitCluster()

# End the re-directing of print statements to file
sink()
#####################################END########################################

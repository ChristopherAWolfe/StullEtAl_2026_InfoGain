#SEX-SPECIFIC FEMALE#

#In this script we re-run the analyses on ordinal variables only on the uncollapsed 
#stages (pooled-sex). Using AIC to infer the best fit model, we use the resulting  
#parameters in the KL divergence and MI calculation. Output from the script includes 
#MI as a function of age per trait exported as a pdf (ordinal only), .RDS files of the
#problem, x solution, y solution, aic, and raw MI values with new ordinal values.

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
sink("results/output_all_F.txt")

# Check that a results folder exists in the working directory
if(! ("results" %in% dir()) ) {
  stop("There is no 'results' folder in the working directory")
}

# The data directory is /results
data_dir <- "results"

#The "ID" that uniquely identifies this analysis (pooled-sex, uncollapsed stages):
analysis_name <- 'MI_female_all'

###############################################################################

########################
#     Initial Prep     #
########################

set.seed(695432)

# Import data based on var_info file and yada functions

var_info <-  yada::load_var_info('data/US_allvar_var_info_ungroup.csv')
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

# Plot the MI over baseline age
pdf(file.path(data_dir,"MI_ord_uncollapsed_F.pdf"))
for(i in 1:j_ord){
  var_name <- problem$var_names[i]
  plot(x0, unlist(MI_ord[[i]]), ylim=c(0,max(unlist(MI_ord[[i]]))),
       xlab="Age [years]",ylab="Mutual Information",
       type="l", lwd=2, col="grey", main = var_name) 
}
dev.off()

#Write MI results to an RDS to be used for later use
write_rds(MI_ord, "results/MI_ord_uncollapsed_F.rds")


# Close all clusters used for paralle processing
stopImplicitCluster()

# End the re-directing of print statements to file
sink()
#####################################END########################################

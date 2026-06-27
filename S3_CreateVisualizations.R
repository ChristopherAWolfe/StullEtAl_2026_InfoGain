################################################################################
#                                                                              #
#                         Script 3: Create Visualizations                      #
#                                                                              #
################################################################################
#
# The following script imports RESULTS from Script 2 and exports IG visualizations.
# A user can modify characterization of individual plots for desired use.
################################################################################

## Step 1: Package Dependencies
### See Script 1 to ensure all requisite packages are installed.

library(yada)

## Step 2: Load in Main Problem file
### Useful for plot labels

main_problem <- readRDS(
  "results\\problem_MI.rds"
)

## Step 3: Prep Variables

j_ord <- main_problem$mod_spec$J

k_cont <- main_problem$mod_spec$K

ordinal_IG <- readRDS("results\\IG_ord_allVars.rds")
cont_IG <- readRDS("results\\IG_cont_allVars.rds")

data_dir <- "results"

x0 = seq(0,23,by=.1)

## Step 4: Create Visualizations
### Here a user can save directly to a PDF.The data_dir is the same directory
### utilized in Script 2.

pdf(file.path(data_dir,"IG_plot.pdf"))
for(i in 1:j_ord){
  var_name <- main_problem$var_names[i]
  plot(x0, unlist(ordinal_IG[[i]]), ylim=c(0,max(unlist(IG_ord[[i]]))),
       xlab="Age [years]",ylab="Information Gain",
       type="l", lwd=2, col="grey", main = var_name)
}
for(i in 1:k_cont){
  var_name <- main_problem$var_names[45:62][i]
  plot(x0, unlist(cont_IG[[i]]),ylim = c(0,max(unlist(cont_IG[[i]]))),
       xlab="Age [years]",ylab="Information Gain",
       type="l", lwd=2, col="grey", main = var_name)
}
dev.off()

####################################END#########################################

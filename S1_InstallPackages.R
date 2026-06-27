################################################################################
#                                                                              #
#                         Script 1: Install Packages                           #
#                                                                              #
################################################################################
#
# The following script installs the requisite packages necessary to quantify 
# mutual information using the `yada` syntax and MCP model.This step only needs
# to be completed once. Once installed, a user can skip to Script 2 and Script 3 
# for subsequent analyses.
################################################################################

## Step 1: Install Packages

install.packages("devtools")
install.packages("tidyverse")
install.packages("doparallel")
install.packages("foreach")
install.packages("mixtools")
devtools::install_github("MichaelHoltonPrice/yada", ref="dev")

####################################END#########################################

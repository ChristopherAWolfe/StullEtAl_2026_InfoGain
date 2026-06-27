#In this script, we calculate the frequency of 'missing' data by age and plot
#for publication. 

#Necessary packages
library(yada)
library(tidyverse)
library(dplyr)
library(ggplot2)

#clear workspace
rm(list=ls())

# Check that a results folder exists in the working directory
if(! ("results" %in% dir()) ) {
  stop("There is no 'results' folder in the working directory")
}

# The data directory is /results
data_dir <- "results"

###############################################################################

########################
#         Prep         #
########################

set.seed(695432)

# Import data based on var_info file and yada functions

var_info <-  yada::load_var_info('data/US_allvar_var_info.csv')
data_file <- 'data/SVAD_US.csv'
dat <- load_cp_data(data_file, var_info)
df <- dat$cp_df

###############################################################################

########################
##    Missing Data    ##
########################

# Frequency of missing data per response
data_na <- data.frame(matrix(NA, nrow=nrow(df), ncol=ncol(df[4:65])))
for (i in 1:ncol(data_na)) {
  data_na[i] <- is.na(df[i+3])
}

names(data_na) <- names(df[4:65])

missing_df <- data.frame(FDL=numeric(),
                         FMSB=numeric(),
                         FDB=numeric(),
                         TDL=numeric(),
                         TPB=numeric(),
                         TMSB=numeric(),
                         TDB=numeric(),
                         FBDL=numeric(),
                         HDL=numeric(),
                         HPB=numeric(),
                         HMSB=numeric(),
                         HDB=numeric(),
                         RDL=numeric(),
                         RPB=numeric(),
                         RMSB=numeric(),
                         RDB=numeric(),
                         UDL=numeric(),
                         UMSB=numeric(),
                         max_M1=numeric(),
                         max_M2=numeric(),
                         max_M3=numeric(),
                         max_PM1=numeric(),
                         max_PM2=numeric(),
                         max_C=numeric(),
                         max_I1=numeric(),
                         max_I2=numeric(),
                         man_M1=numeric(),
                         man_M2=numeric(),
                         man_M3=numeric(),
                         man_PM1=numeric(),
                         man_PM2=numeric(),
                         man_C=numeric(),
                         man_I1=numeric(),
                         man_I2=numeric(),
                         FH_EF=numeric(),
                         FGT_EF=numeric(),
                         FLT_EF=numeric(),
                         FDE_EF=numeric(),
                         TPE_EF=numeric(),
                         TDE_EF=numeric(),
                         FBPE_EF=numeric(),
                         FBDE_EF=numeric(),
                         HH_Oss=numeric(),
                         HGT_Oss=numeric(),
                         HLT_Oss=numeric(),
                         HPE_EF=numeric(),
                         HC_Oss=numeric(),
                         HT_Oss=numeric(),
                         HLE_Oss=numeric(),
                         HDE_EF=numeric(),
                         HME_EF=numeric(),
                         RPE_EF=numeric(),
                         RDE_EF=numeric(),
                         UPE_EF=numeric(),
                         UDE_EF=numeric(),
                         CT_EF=numeric(),
                         CC_Oss=numeric(),
                         TC_Oss=numeric(),
                         ISPR_EF=numeric(),
                         ILIS_EF=numeric(),
                         PC_Oss = numeric(),
                         IC_EF = numeric())
for(j in 1:ncol(data_na)) {
  pct <- length(which(data_na[j]==TRUE))/nrow(data_na)
  missing_df[1,j] <- pct*100
}

colnames(missing_df) <- colnames(df[4:65])

# Write missing_df to file
print("Frequency of missing data per response variable:")
print(missing_df)

## VISUALIZATION##

# Calculate data availability frequencies
df$agey <- floor(df$agey)  # convert age as integer
age_vec <- sort(unique(df$agey))  # vector of unique ages
wide_df <- NULL
N_vec <- NULL

for(i in 1:length(age_vec)) {
  temp_df <- df[which(df$agey==age_vec[i]),]
  N_vec <- c(N_vec, nrow(temp_df))
  NA_vec <- as.vector(colSums(is.na(temp_df[4:65])))/N_vec[i]#count NA
  wide_df <- cbind(wide_df, (1-NA_vec))
  colnames(wide_df) <- age_vec[1:(i)]
}
wide_df <- as.data.frame(cbind(colnames(df[4:65]),wide_df))# add var as column
colnames(wide_df) <- c('var',age_vec)

# Format data for geom_raster
long_df <- tidyr::gather(wide_df, key="age_int", value="freq", -var)
long_df$var<-factor(long_df$var, levels=c("FDL","FMSB","FDB","TDL","TPB","TMSB",
      "TDB","FBDL","HDL","HPB","HMSB","HDB","RDL","RPB","RMSB","RDB",
      "UDL","UMSB","max_M1","max_M2","max_M3","max_PM1","max_PM2","max_C",
      "max_I1","max_I2","man_M1","man_M2","man_M3","man_PM1","man_PM2",
      "man_C","man_I1","man_I2","FH_EF","FGT_EF","FLT_EF","FDE_EF","TPE_EF",
      "TDE_EF","FBPE_EF","FBDE_EF","HH_Oss","HGT_Oss","HLT_Oss","HPE_EF",
      "HC_Oss","HT_Oss","HLE_Oss","HDE_EF","HME_EF","RPE_EF","RDE_EF",
      "UPE_EF","UDE_EF","CT_EF","CC_Oss","TC_Oss","ISPR_EF",
      "ILIS_EF", "PC_Oss", "IC_EF"))
long_df$age_int <- factor(long_df$age_int, levels=age_vec)
long_df$freq <- as.numeric(long_df$freq)

# Generate plot
pdf(file.path(data_dir,"Missing.pdf"))
print(
  ggplot(long_df, aes(x=age_int, y=var)) + 
    geom_tile(aes(fill=(freq*100)), color="black") + coord_equal() +
    theme_minimal() + scale_fill_gradient(low="white",high="dodgerblue4") + 
    labs(x="Age [years]", y="Variable", fill="% Available") + 
    theme(axis.text.x = element_text(size = 6), 
          axis.text.y = element_text(size = 6))
)
dev.off()
###############################END##############################################

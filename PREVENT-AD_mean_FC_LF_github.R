##########################################################################################################################

#Longitudinal functional connectivity during rest and task is differentially 
#related to Alzheimer’s pathology and episodic memory in older adults

#Analysis of PREVENT-AD data - mean FC
#Larissa Fischer - Multimodal Neuroimaging Lab, DZNE Magdeburg, 2024

#########################################################################################################################


#Import packages:
library(readxl)
library(tidyverse)
library(corrplot)
library(dplyr)
library(rstatix)
library(ggpubr)
library(psych)
library(datasets)
library(ggplot2)
library(writexl)
library(ggplot2)
library(plotly)
library(kableExtra)
library(knitr)
library(DT)
library(matrixStats)

################################################################################
#set working directory:
setwd("/Users/your/path")
################################################################################
#run this for condition 001 to 005, those are the different timepoints
#BL, FU12, FU24, FU36, FU48

data_FC_6 <- read_csv("condition_005_results.csv", na = "NA") 
#these csv files contain the values for functional connectivity strength for all participant 
colnames(data_FC_6) <- gsub("Brainnetome\\.", "", colnames(data_FC_6))
colnames(data_FC_6) <- gsub("CG_R_7_7", "Subg_R_7_7", colnames(data_FC_6))
colnames(data_FC_6) <- gsub("CG_L_7_7", "Subg_L_7_7", colnames(data_FC_6))
colnames(data_FC_6)

# define vectors for meta-ROIs
MTL <- c("PhG_L_6_1", "PhG_R_6_1", "PhG_L_6_2", "PhG_R_6_2", "PhG_L_6_3", "PhG_R_6_3", "PhG_L_6_4", "PhG_R_6_4", 
         "PhG_L_6_6", "PhG_R_6_6", "Hipp_L_2_1", "Hipp_R_2_1", "Hipp_L_2_2", "Hipp_R_2_2")
PMC <- c("PCun_L_4_1", "PCun_R_4_1", "PCun_L_4_2", "PCun_R_4_2", "PCun_L_4_3", "PCun_R_4_3", "PCun_L_4_4", "PCun_R_4_4",
         "CG_L_7_1", "CG_R_7_1", "CG_L_7_4", "CG_R_7_4")
mPFC <- c("SFG_L_7_7", "SFG_R_7_7", "OrG_L_6_1", "OrG_R_6_1", "Subg_L_7_7", "Subg_R_7_7")

# function to classify functional connections
classify_connection <- function(connection) {
  # extract the 2 regions of a connection
  regions <- strsplit(connection, "-")[[1]]
  region1 <- regions[1]
  region2 <- regions[2]
  
  # within?
  within_MTL <- region1 %in% MTL && region2 %in% MTL
  within_PMC <- region1 %in% PMC && region2 %in% PMC
  within_mPFC <- region1 %in% mPFC && region2 %in% mPFC
  
  # between?
  between_MTL_PMC <- (region1 %in% MTL && region2 %in% PMC) || (region1 %in% PMC && region2 %in% MTL)
  between_mPFC_PMC <- (region1 %in% mPFC && region2 %in% PMC) || (region1 %in% PMC && region2 %in% mPFC)
  between_MTL_mPFC <- (region1 %in% MTL && region2 %in% mPFC) || (region1 %in% mPFC && region2 %in% MTL)
  
  # give out classification
  if (within_MTL) {
    return("within_MTL")
  } else if (within_PMC) {
    return("within_PMC")
  } else if (within_mPFC) {
    return("within_mPFC")
  } else if (between_MTL_PMC) {
    return("between_MTL-PMC")
  } else if (between_mPFC_PMC) {
    return("between_mPFC-PMC")
  } else if (between_MTL_mPFC) {
    return("between_MTL-mPFC")
  } else {
    return("Unknown")
  }
}

connections <- names(data_FC_6)
classified_connections <- sapply(connections, classify_connection)
meta_ROIs <- data.frame(Connection = connections, Category = classified_connections)


# Meta-ROI-classification as new row in data_FC_6 
data_FC_6 <- rbind("Meta-ROI" = classified_connections, data_FC_6)
# Make it heading
colnames(data_FC_6) <- unlist(data_FC_6[1, ])
data_FC_6 <- data_FC_6[-1, ]
# Convert non-numeric columns to numeric
data_FC_6[, -c(1, 2)] <- sapply(data_FC_6[, -c(1, 2)], as.numeric)
# Add a suffix with a consecutive number to each column name
colnames(data_FC_6) <- paste0(colnames(data_FC_6), "_", seq_along(colnames(data_FC_6)))


#Means for each Meta-ROI
######
# list to save mean FC for each Meta-ROI (label of classified_connection)
mean_data_list <- list()
#  iteration
for (connection_label in classified_connections) {
  # select columns that start with specific classified_connection
  selected_cols <- colnames(data_FC_6)[startsWith(colnames(data_FC_6), connection_label)]
  # extract and make numeric
  numeric_data <- data_FC_6 %>%
    select(all_of(selected_cols)) %>%
    mutate(across(everything(), as.numeric))
  # mean for each row of the selected columns
  mean_values <- rowMeans(numeric_data, na.rm = TRUE)

  # add mean to list
  mean_data_list[[connection_label]] <- mean_values
}
# dataframe
mean_data <- as.data.frame(mean_data_list)
#write_xlsx(mean_data, "meta_ROI_means_analysis_6_FU48.xlsx")    #BL, FU12, FU24, FU36, FU48
#save(mean_data, file = "meta_ROI_means_analysis_6_FU48.RData")

# list to save standard deviation of FC for each Meta-ROI (label of classified_connection)
sd_data_list <- list()
for (connection_label in classified_connections) {
  selected_cols <- colnames(data_FC_6)[startsWith(colnames(data_FC_6), connection_label)]
    numeric_data <- data_FC_6 %>%
    select(all_of(selected_cols)) %>%
    mutate(across(everything(), as.numeric))
    sd_values <- apply(numeric_data, 1, sd, na.rm = TRUE)
    sd_data_list[[connection_label]] <- sd_values
}
sd_data <- as.data.frame(sd_data_list)

#write_xlsx(sd_data, "meta_ROI_sd_analysis_6_FU48.xlsx")    #BL, FU12, FU24, FU36, FU48
#save(sd_data, file = "meta_ROI_sd_analysis_6_FU48.RData")



################################################################################
#manual sanity check
################################################################################
data_FC_6_orig <- read_csv("condition_001_results.csv", na = "NA")
colnames(data_FC_6_orig) <- gsub("Brainnetome\\.", "", colnames(data_FC_6_orig))
colnames(data_FC_6_orig) <- gsub("CG_R_7_7", "Subg_R_7_7", colnames(data_FC_6_orig))
colnames(data_FC_6_orig) <- gsub("CG_L_7_7", "Subg_L_7_7", colnames(data_FC_6_orig))
colnames(data_FC_6_orig)

# sum
data_FC_6_orig <- data_FC_6_orig %>%
  mutate(
    FC_within_mPFC = `SFG_R_7_7-SFG_L_7_7` + `OrG_L_6_1-SFG_L_7_7` +
      `OrG_L_6_1-SFG_R_7_7` + `OrG_R_6_1-SFG_L_7_7` +
      `OrG_R_6_1-SFG_R_7_7` + `OrG_R_6_1-OrG_L_6_1` +
      `Subg_L_7_7-SFG_L_7_7` + `Subg_L_7_7-SFG_R_7_7` +
      `Subg_L_7_7-OrG_L_6_1` + `Subg_L_7_7-OrG_R_6_1` +
      `Subg_R_7_7-SFG_L_7_7` + `Subg_R_7_7-SFG_R_7_7` +
      `Subg_R_7_7-OrG_L_6_1` + `Subg_R_7_7-OrG_R_6_1` +
      `Subg_R_7_7-Subg_L_7_7`)

# mean
data_FC_6_orig <- data_FC_6_orig %>%
  mutate(mean_FC_within_mPFC = FC_within_mPFC / 15)


    
 
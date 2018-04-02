## TO DO: Write code to identify characteristics of act at which transmission occurred (e.g., ai vs. vi, adh = 0/1)

# options(warn = 2) # Turn warnings into errors
# options(error = browser) # Enter browser on error

if(grepl(pattern = "Linux", Sys.info()['sysname'])) {
  setwd("/gscratch/csde/kpeebles")
} else {
  setwd("~/Documents/code/aspire")
}

## Attach packages
library(readxl)
library(data.table)
library(dplyr)
library(ggplot2)
library(mnormt)

## Source functions
source_files <- list.files(paste0(getwd(), "/fx/"), pattern = "*.R", full.names = T)
sapply(source_files, function(x) source(x))

## Read and load data
suppressWarnings(load_data())

## Load parameters
params_dt <- as.data.table(read_excel(paste0(getwd(), "/parameters/parameters.xlsx"), range = "A1:B20", col_names = T))
params <- lapply(params_dt[, name], function(x) { x = params_dt[name == x, value] })
names(params) <- params_dt[, name]

## Convert joint probabilities in age mixing matrix into probability of male partner age conditional on female partner age
age_mix_mat_cond <- sapply(1:ncol(age_mix_mat), function(x) {
  sapply(1:nrow(age_mix_mat), function(y) round(age_mix_mat[y, x]/sum(age_mix_mat[, x]), 2))
})
rownames(age_mix_mat_cond) <- colnames(age_mix_mat_cond) <- colnames(age_mix_mat)

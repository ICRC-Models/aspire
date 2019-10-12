## TO DO: Write code to identify characteristics of act at which transmission occurred (e.g., ai vs. vi, adh = 0/1)

# options(warn = 2) # Turn warnings into errors
# options(error = browser) # Enter browser on error

if(grepl(pattern = "Linux", Sys.info()['sysname'])) {
  wd <- "/gscratch/csde/kpeebles"
} else {
  wd <- "~/Documents/code/aspire"
}

## Attach packages
library(readxl)
library(data.table)
library(dplyr)
library(ggplot2)
library(mnormt)
library(parallel)
library(zoo)
library(survival)

## Source functions
source_files <- list.files(paste0(wd, "/fx/"), pattern = "*.R", full.names = T)
sapply(source_files, function(x) source(x))

## Read and load data
suppressWarnings(load_data())

## Load parameters
params_dt <- as.data.table(read_excel(paste0(wd, "/parameters/parameters.xlsx"), range = "A1:B97", col_names = T))
params <- lapply(params_dt[, name], function(x) { x = params_dt[name == x, value] })
names(params) <- params_dt[, name]


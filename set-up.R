## Purpose: Source functions and load data needed to run ASPIRE model simulations
## Author:  Kathryn Peebles

if(grepl(pattern = "Linux", Sys.info()['sysname'])) {
  wd <- "/gscratch/csde/kpeebles"
} else {
  wd <- "~/Documents/code/aspire"
}

# If on Hyak, append local R library to default R library
if(wd == "/gscratch/csde/kpeebles") {
  .libPaths(c("/gscratch/csde/kpeebles/R", .libPaths()))
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
params_dt <- as.data.table(read_excel(paste0(wd, "/parameters/parameters.xlsx"), range = "A1:B26", col_names = T))
params <- lapply(params_dt[, name], function(x) { x = params_dt[name == x, value] })
names(params) <- params_dt[, name]


## Purpose: Create predictive models for variables with missing data that are used in simulations.
## Date:    9 November 2018
## Author:  Kathryn Peebles

library(data.table)

setwd("~/Documents/code/aspire/data-private")

# Load datasets
load(paste0(getwd(), "/hiv.RData")) # baseline variables
load(paste0(getwd(), "/ba.RData"))  # follow-up CRF
load(paste0(getwd(), "/bba.RData")) # baseline CRF
acasi <- as.data.table(read.csv(file = paste0(getwd(), "/ACASI.csv"))) # ACASI

# Predictive model for unprotected sex in the last week at baseline
coef <- glm(unplast7base ~ basehivsp + baseSTD + noalc + age + married, data = hiv, family = binomial(link = "logit"))$coef
coef

# Predictive model for engaging in any AI at baseline. Transactional sex was a significant predictor of AI in adjusted analyses in our brief report, but was only collected at baseline. Create predictive model using baseline data only, and use model to impute missing values at both baseline and follow-up in simulations.
acasi <- as.data.table(acasi)
acasi[, ptid := gsub(pattern = "-", replacement = "", x = ptid)]
acasi[, ptid := gsub(pattern = " ", replacement = "", x = ptid)]
acasi[, anal_any := ifelse(QANAL == "Skipped by participant  ", NA, as.numeric(as.character(QANAL)) > 0)]

# Merge data tables. Do not keep participants who have responses in ACASI dataset, but not in follow-up data.
hiv <- as.data.table(hiv)
dt <- merge(x  = hiv[, .(ptid, site, noalc, age, edu, married, bdepo, bneten)], 
            y  = acasi[survey == 1, .(ptid, survey, QEXCH, anal_any)], 
            by = "ptid", 
            all.x = T)

# Create clean variable for transactional sex
dt[, trans_sex := ifelse(QEXCH == "                        " | QEXCH == "Skipped by participant  ", NA_real_, ifelse(QEXCH == "No                      ", 0, 1))]

coef <- glm(anal_any ~ site + edu + bdepo + bneten + trans_sex, data = dt, family = binomial(link = "logit"))$coef
coef

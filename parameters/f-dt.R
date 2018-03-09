## Data from Elizabeth Brown, email dated 2017-11-17. Note that there are 2,614 observations here, while Baeten, 2016 reports 2,629 women in the trial. 15 women did not report a partner in the previous 3 months, and so are excluded from model analyses.

library(data.table)
library(splitstackshape)

load("~/Dropbox/UW Epi Program/RB/sources-eff-dilution-ring/ring-eff/data/forKathrn.RData")
dt <- as.data.table(d1)

f_dt <- dt[, .(Country, Age, arm, N)]

f_dt <- expandRows(f_dt, "N")

f_dt[Country == "Malawi", country := "mal"]
f_dt[Country == "South Africa", country := "sa"]
f_dt[Country == "Uganda", country := "uga"]
f_dt[Country == "Zimbabwe", country := "zim"]

f_dt[, Country := NULL]

f_dt[, arm := ifelse(test = arm == "Placebo",
                     yes  = 0,
                     no   = 1)]

f_dt[, age_cat := c("18-21", "22-26", "27-45")[findInterval(x = f_dt$Age, vec = c(18, 22, 27))]] # age categories (18-21, 22-26, 27-45)

setkey(x = f_dt, country, Age)

f_dt[, id := 1:nrow(f_dt)]

setnames(x = f_dt, old = "Age", new = "age")

setcolorder(x = f_dt, neworder = c("id", "age", "age_cat", "country", "arm"))

save(f_dt, file = "~/Dropbox/UW Epi Program/RB/Effect of AI on ring efficacy/ai-ring/data/f_dt.RData")

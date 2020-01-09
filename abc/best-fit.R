## Purpose: Identify parameter set with highest weighted density in ABC model selection process
## Author:  Kathryn Peebles
## Date:    5 April 2018

library(ggplot2)
library(ggtern)
library(data.table)
library(scales)

setwd("~/Documents/code/aspire/dev/test-sims/abc/hyak")

load(paste0(getwd(), "/t5/particles_acc_5.RDATA"))

particles_acc_5[, mean := a/(a + b)]

hm <- ggtern::kde2d.weighted(x = particles_acc_5$prop_ai, y = particles_acc_5$mean, n = 100, w = particles_acc_5$w)

hm_max <- hm$z == max(hm$z)
hm_row <- which.max(rowSums(hm_max))
hm_col <- which.max(colSums(hm_max))
bestfit <- c(hm$x[hm_row], hm$y[hm_col])
bestfit

setkey(particles_acc_5, prop_ai)

View(particles_acc_5[prop_ai > 0.11 & prop_ai < 0.111])

jpeg(file = "~/Dropbox/UW Epi Program/RB/ASPIRE/Presentation-2018-04-05/prior_particles.jpeg", width = 610, height = 480)
image(x    = hm, 
      col  = heat.colors(1000), 
      xlim = c(0, 1), 
      xlab = "Proportion of women engaging in AI", 
      ylab = "Mean of beta distribution",
      main = "Particles sampled from priors")
dev.off()

for(i in 1:5) {
  load(paste0(getwd(), "/t", i, "/particles_acc_", i, ".RDATA"))
  particles <- eval(parse(text = paste0("particles_acc_", i)))
  
  particles[, mean := a/(a + b)]
  
  print(sum(particles$t == i)/nrow(particles))
  
  hm <- ggtern::kde2d.weighted(x = particles$prop_ai, y = particles$mean, n = 100, w = particles$w)

  jpeg(file = paste0("~/Dropbox/UW Epi Program/RB/ASPIRE/Presentation-2018-04-05/intermed-distribution-", i, ".jpeg"), width = 610, height = 480)
  image(x    = hm,
        col  = heat.colors(1000),
        xlim = c(0, 1),
        xlab = "Proportion of women engaging in AI",
        ylab = "Mean of beta distribution",
        main = paste0("Intermediate approximate posterior at iteration ", i))
  dev.off()
  
  if(i == 5) {
    jpeg(file = "~/Dropbox/UW Epi Program/RB/ASPIRE/Presentation-2018-04-05/bestfit.jpeg", width = 610, height = 480)
    image(x    = hm,
        col  = heat.colors(1000),
        xlim = c(0, 1),
        xlab = "Proportion of women engaging in AI",
        ylab = "Mean of beta distribution",
        main = paste0("Intermediate approximate posterior at iteration ", i))
    points(x = bestfit[1], y = bestfit[2], col = "blue", pch = 19)
    dev.off()
  }
}



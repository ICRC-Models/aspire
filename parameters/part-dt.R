## Email dated 2017-11-16 from Elizabeth Brown has distribution of number of extra partners by country. The same email shows the distribution of the number of coital acts in the previous 7 days. While there are clear differences in the distribution of number of partners by country, coital frequency is similar across countries. It therefore does not appear that women have the same number of sex acts with a third partner as with a primary partner in ASPIRE. We will collapse partner number categories to be one partner and 2+ partners.

part_dt <- data.table(country = c("mal", "sa", "uga", "zim"),
                      one_p   = c(250, 1210, 135, 582),
                      two_p   = c(sum(12, 6, 3),
                                  sum(180, 21, 5),
                                  sum(72, 17, 29),
                                  sum(40, 22, 30)))

save(part_dt, file = "~/Dropbox/UW Epi Program/RB/Effect of AI on ring efficacy/ai-ring/data/part_dt.RData")

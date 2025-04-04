library(pheatmap)
library(ggplot2)

dat <- read.csv("innerjoin.csv", header=T, row.names = 1) #use the fold change compared to 28celcius, only from top 200 degs from each comparison (28 vs 20, 28 vs 37)

a<-pheatmap(dat, fontsize_row = 15, cellwidth = 59)

ggsave("aheatmap.png", plot = a, width = 20, height = 39, dpi = 300, units = "cm", limitsize = FALSE)

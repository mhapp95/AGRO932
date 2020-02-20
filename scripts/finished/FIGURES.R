rm(list = ls())
setwd('/work/soybean/mhapp95/AGRO932/AGRO932/data/')

library(data.table)
library(readr)
library(ggplot2)

thetas <- fread("raw/theta.txt")
png("figures/DistThetas.png", width=400, height=400)
hist(thetas$Pairwise)
dev.off()

fst <- read.table("raw/fst_win.txt", skip=1, header=FALSE)
names(fst)[c(3,5)] <- c("midp", "fst")
png("figures/ScatterPlotFst.png", width=400, height=400)
plot(fst$midp, fst$fst, xlab="Physical position", ylab="Fst", col="#5f9ea0", pch=16)
dev.off()

rm(list = ls())
setwd('/work/soybean/mhapp95/AGRO932/AGRO932/data/')

library(data.table)
library(readr)
library(ggplot2)


s<-scan('raw/out.sfs')
s<-s[-c(1,length(s))]
s<-s/sum(s)
png("figures/SFS.png", width=400, height=400)
barplot(s,names=1:length(s),main='SFS', col="slategray", xlab="Frequency", ylab="Density")
dev.off()

thetas <- fread("raw/theta.txt")
png("figures/DistThetas.png", width=400, height=400)
hist(thetas$Pairwise, col="seagreen, xlab="Pairwise Thetas", main="Distribution of Pairwise Thetas")
dev.off()

fst <- read.table("raw/fst_win.txt", skip=1, header=FALSE)
names(fst)[c(3,5)] <- c("midp", "fst")
png("figures/ScatterPlotFst.png", width=400, height=400)
plot(fst$midp, fst$fst, xlab="Physical position (bp)", ylab="Fst", col="#5f9ea0", pch=16, main="Per Site Fst")
dev.off()

thetastwo <- fread("raw/thetasWindow.pestPG")
png("figures/DistTajima.png", width=400, height=400)
hist(thetastwo$Tajima, breaks=30, col="plum", xlab="Tajima's D", main="Distribution of Tajima's D")
dev.off()

gff <- fread(cmd='grep -v "#" reference/ZeaMays_Chloroplast.gff3', header=FALSE, data.table=FALSE)
names(gff) <- c("seq", "source", "feature", "start", "end", "score", "strand", "phase", "att")
g <- subset(gff, feature %in% "gene")
g$geneid <- gsub(".*gene:|;biotype.*", "", g$att)
### + strand
gp <- subset(g, strand %in% "+") 
# nrow(gp) 75
### get the 5k upstream of the + strand gene model
gp_up <- gp
gp_up$end <- gp_up$start - 1
gp_up$start <- gp_up$end - 5000 
### get the 5k downstream of the + strand gene model
gp_down <- gp
gp_down$start <- gp_down$end + 1
gp_down$end <- gp_down$start + 5000
gm <- subset(g, strand %in% "-")

library(GenomicRanges)
library(plyr)

get_mean_theta <- function(feature){
  # gf_file: gene feature file [chr, ="cache/mt_gene_up5k.txt"]
  theta <- fread("raw/theta.txt", data.table=FALSE)
  names(theta)[1] <- "seq"
  ### define the subject file for theta values
  grc <- with(theta, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))
  ### define the query file for genomic feature
  grf <- with(feature, GRanges(seqnames=seq, IRanges(start=start, end=end), geneid=geneid))
  ### find overlaps between the two
  tb <- findOverlaps(query=grf, subject=grc)
  tb <- as.matrix(tb)
  out1 <- as.data.frame(grf[tb[,1]])
  out2 <- as.data.frame(grc[tb[,2]])
  ### for each genomic feature, find the sites with non-missing data
  out <- cbind(out1, out2[, "start"]) 
  names(out)[ncol(out)] <- "pos"
  #define unique identifier and merge with the thetas
  out$uid <- paste(out$seqnames, out$pos, sep="_")
  theta$uid <- paste(theta$seq, theta$Pos, sep="_")
  df <- merge(out, theta[, c(-1, -2)], by="uid")
# for each upstream 5k region, how many theta values
  mx <- ddply(df, .(geneid), summarise,
            Pairwise = mean(Pairwise, na.rm=TRUE),
            thetaH = mean(thetaH, na.rm=TRUE),
            nsites = length(uid))
}

up5k <- get_mean_theta(gp_up)
down5k <- get_mean_theta(gp_down)

library(ggplot2)
up5k$feature <- "up 5k"
down5k$feature <- "down 5k"
res <- rbind(up5k, down5k)
png("figures/USvsDSTheta.png", width=400, height=400)
ggplot(res, aes(x=feature, y=Pairwise, fill=feature)) + 
  geom_violin(trim=FALSE)+
  labs(title="Theta value", x="", y = "Log10 (theta)")+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_brewer(palette="Blues") + 
  theme_classic()
dev.off()

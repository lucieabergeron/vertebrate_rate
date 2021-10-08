library(ape)
library(RColorBrewer)
library(dplyr)
library(phytools)

tree<-read.tree("tree_rooted.nwk")
plot(tree)
nodelabels(c(75,80,120,96,98,100,135,89,78,105,76,110,140,147),c(75,80,120,96,98,100,135,89,78,105,76,110,140,147),cex=0.5)

# Calibrating table
node<- c(75,80,120,96,98,100,135,89,78,105,76,110,140,147)
age.min <- c(416.00,66.00,48.50,222.80,54.00,160.70,37.10,51.81,260.00,61.70,313.40,23.50,378.20,50.00)
age.max <- c(425.40,66.00,48.50,222.80,54.00,160.70,37.10,51.81,260.00,61.70,313.40,23.50,378.20,50.00)
soft.bounds <- c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)

# Calibration
cal <- data.frame(node, age.min, age.max, soft.bounds)
chr <- chronos(tree,lambda=0, calibration = cal, model="relaxed")

# Extract the substitution rate
rate_sp=cbind(chr$edge[,2],attr(chr, "rate"))
rate_sp_table=cbind(chr$tip.label, rate_sp[which(rate_sp[,1]<75),])
write.table(rate_sp_table, file="sub_rate_uce_morecal.txt")
attr(chr, "rate")
write.tree(chr, file="tree_calibrated_UCE_morecal.nwk")


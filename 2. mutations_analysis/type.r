ata = read.csv("sp_all_type.txt", sep="\t")

all = c(sum(data[,3]), sum(data[,4]),
   sum(data[,5]), sum(data[,6]),
   sum(data[,7]), sum(data[,8]),
   sum(data[,9]), sum(data[,10]), sum(data[,11]))


## Per group
subdata_mammals = data[which(data[,2]=="Mammal"),]
all_mammals = c(sum(subdata_mammals[,3]), sum(subdata_mammals[,4]),
   sum(subdata_mammals[,5]), sum(subdata_mammals[,6]),
   sum(subdata_mammals[,7]), sum(subdata_mammals[,8]),
   sum(subdata_mammals[,9]), sum(subdata_mammals[,10]), sum(subdata_mammals[,11]))

subdata_birds = data[which(data[,2]=="Bird"),]
all_birds = c(sum(subdata_birds[,3]), sum(subdata_birds[,4]),
   sum(subdata_birds[,5]), sum(subdata_birds[,6]),
   sum(subdata_birds[,7]), sum(subdata_birds[,8]),
   sum(subdata_birds[,9]), sum(subdata_birds[,10]), sum(subdata_birds[,11]))

subdata_fish = data[which(data[,2]=="Fish"),]
all_fish = c(sum(subdata_fish[,3]), sum(subdata_fish[,4]),
   sum(subdata_fish[,5]), sum(subdata_fish[,6]),
   sum(subdata_fish[,7]), sum(subdata_fish[,8]),
   sum(subdata_fish[,9]), sum(subdata_fish[,10]), sum(subdata_fish[,11]))

subdata_rep = data[which(data[,2]=="Reptile"),]
all_rep = c(sum(subdata_rep[,3]), sum(subdata_rep[,4]),
   sum(subdata_rep[,5]), sum(subdata_rep[,6]),
   sum(subdata_rep[,7]), sum(subdata_rep[,8]),
   sum(subdata_rep[,9]), sum(subdata_rep[,10]), sum(subdata_rep[,11]))

## Proportion for the various group
prop_mam = (all_mammals/sum(all_mammals))*100
prop_bird = (all_birds/sum(all_birds))*100
prop_fish = (all_fish/sum(all_fish))*100
prop_rep = (all_rep/sum(all_rep))*100

##
mut_type_mam=matrix(c(prop_mam[7],0,prop_mam[8],0,prop_mam[9],0, prop_mam[1], prop_mam[2], prop_mam[3], prop_mam[4], prop_mam[5], prop_mam[6]), nrow=2)
mut_type_bird=matrix(c(prop_bird[7],0,prop_bird[8],0,prop_bird[9],0, prop_bird[1], prop_bird[2], prop_bird[3], prop_bird[4], prop_bird[5], prop_bird[6]), nrow=2)
mut_type_fish=matrix(c(prop_fish[7],0,prop_fish[8],0,prop_fish[9],0, prop_fish[1], prop_fish[2], prop_fish[3], prop_fish[4], prop_fish[5], prop_fish[6]), nrow=2)
mut_type_rep=matrix(c(prop_rep[7],0,prop_rep[8],0,prop_rep[9],0, prop_rep[1], prop_rep[2], prop_rep[3], prop_rep[4], prop_rep[5], prop_rep[6]), nrow=2)



mut_type = cbind(mut_type_mam[,1], mut_type_bird[,1], mut_type_fish[,1], mut_type_rep[,1],
  mut_type_mam[,2], mut_type_bird[,2], mut_type_fish[,2], mut_type_rep[,2],
  mut_type_mam[,3], mut_type_bird[,3], mut_type_fish[,3], mut_type_rep[,3],
  mut_type_mam[,4], mut_type_bird[,4], mut_type_fish[,4], mut_type_rep[,4],
  mut_type_mam[,5], mut_type_bird[,5], mut_type_fish[,5], mut_type_rep[,5],
  mut_type_mam[,6], mut_type_bird[,6], mut_type_fish[,6], mut_type_rep[,6])

png("type_of_mutation_diff.png", width = 1300, height = 850)
par(mar=c(10,11,4,2), mgp=c(3, 1.5, 0))
barplot(mut_type[1,]+mut_type[2,], space = c(0,0,0,0,.5,0,0,0,.5,0,0,0,.5,0,0,0,.5,0,0,0,.5,0,0,0), ylim=c(0,max(mut_type[1,]+mut_type[2,])+3),
  col= c("darksalmon","khaki1","lightsteelblue1","darkolivegreen1"),
  cex.axis = 3, cex.names=3, yaxt="n")
barplot(mut_type[1,], space = c(0,0,0,0,.5,0,0,0,.5,0,0,0,.5,0,0,0,.5,0,0,0,.5,0,0,0), ylim=c(0,max(mut_type[1,]+mut_type[2,])+3), add=T,
col=c("firebrick3","goldenrod2","royalblue4", "limegreen"),
  cex.axis = 3, cex.names=3, yaxt="n")
legend("topleft", bty = "n", fill=c("firebrick3","goldenrod2","royalblue4", "limegreen"),c("Mammals", "Birds", "Fish", "Reptiles"), cex=2)
legend(4.2, 52.5, bty = "n", fill=c("grey72", "gray33"),c("CpG", "non-CpG"), cex=2)
axis(1, at=c(1.99,6.49,10.99,15.49,20.02,24.49), labels =c("A>C","A>G","A>T","C>A","C>G","C>T"), cex.axis=2.75)
mtext(c("(T>G)", "(T>C)", "(T>A)", "(G>T)", "(G>C)", "(G>A)"), 1, at=c(1.99,6.49,10.99,15.49,20.02,24.49), line=3, cex=2)
axis(2, cex.axis=3, las=2)
# significance
lines(c(-0.0045,4.0219),c(13,13),lwd=3)
lines(c(-0.0045,-0.0045),c(13,12),lwd=3)
lines(c(4.0219,4.0219),c(13,12),lwd=3)
text(2,16.5,label="p = 0.001", cex=2, font=2)
text(2,14,label="**", cex=3.5)
#
lines(c(13.491,17.517),c(17,17),lwd=3)
lines(c(13.491,13.491),c(17,16),lwd=3)
lines(c(17.517,17.517),c(17,16),lwd=3)
text(15.5,20.5,label="p = 0.032", cex=2, font=2)
text(15.5,18,label="*", cex=3.5)
#
mtext("Type of mutation", 1, line=6, cex=3)
mtext("Pourcentage of total mutations per group", 2, line=5, cex=3)
box()
dev.off()

## Chi square
mut_type=matrix(c(all_mammals[7],0,all_mammals[8],0,all_mammals[9],0, all_mammals[1], all_mammals[2], all_mammals[3], all_mammals[4], all_mammals[5], all_mammals[6]), nrow=2)
mam_type = mut_type[1,]+mut_type[2,]
mam_cpg = c(sum(mut_type[1,]),sum(mut_type[2,]))
mut_type=matrix(c(all_birds[7],0,all_birds[8],0,all_birds[9],0, all_birds[1], all_birds[2], all_birds[3], all_birds[4], all_birds[5], all_birds[6]), nrow=2)
bir_type = mut_type[1,]+mut_type[2,]
bir_cpg = c(sum(mut_type[1,]),sum(mut_type[2,]))
mut_type=matrix(c(all_fish[7],0,all_fish[8],0,all_fish[9],0, all_fish[1], all_fish[2], all_fish[3], all_fish[4], all_fish[5], all_fish[6]), nrow=2)
fish_type = mut_type[1,]+mut_type[2,]
fish_cpg = c(sum(mut_type[1,]),sum(mut_type[2,]))
mut_type=matrix(c(all_rep[7],0,all_rep[8],0,all_rep[9],0, all_rep[1], all_rep[2], all_rep[3], all_rep[4], all_rep[5], all_rep[6]), nrow=2)
rep_type = mut_type[1,]+mut_type[2,]
rep_cpg = c(sum(mut_type[1,]),sum(mut_type[2,]))

# Chisquare
table=rbind(mam_type,bir_type,fish_type,rep_type)
chisq.test(table)
##X-squared = 32.239, df = 20, p-value = 0.04082
# without turtle
##X-squared = 20.916, df = 15, p-value = 0.1395

## Per type
table1 = rbind(c(table[1,1],sum(table[1,2:6])),
    c(table[2,1],sum(table[2,2:6])),
    c(table[3,1],sum(table[3,2:6])),
    c(table[4,1],sum(table[4,2:6])))
table2 = rbind(c(table[1,2],sum(table[1,c(1,3,4,5,6)])),
    c(table[2,2],sum(table[2,c(1,3,4,5,6)])),
    c(table[3,2],sum(table[3,c(1,3,4,5,6)])),
    c(table[4,2],sum(table[4,c(1,3,4,5,6)])))
table3 = rbind(c(table[1,3],sum(table[1,c(1,2,4,5,6)])),
    c(table[2,3],sum(table[2,c(1,2,4,5,6)])),
    c(table[3,3],sum(table[3,c(1,2,4,5,6)])),
    c(table[4,3],sum(table[4,c(1,2,4,5,6)])))
table4 = rbind(c(table[1,4],sum(table[1,c(1,2,3,5,6)])),
    c(table[2,4],sum(table[2,c(1,2,3,5,6)])),
    c(table[3,4],sum(table[3,c(1,2,3,5,6)])),
    c(table[4,4],sum(table[4,c(1,2,3,5,6)])))
table5 = rbind(c(table[1,5],sum(table[1,c(1,2,3,4,6)])),
    c(table[2,5],sum(table[2,c(1,2,3,4,6)])),
    c(table[3,5],sum(table[3,c(1,2,3,4,6)])),
    c(table[4,5],sum(table[4,c(1,2,3,4,6)])))
table6 = rbind(c(table[1,6],sum(table[1,c(1,2,3,4,5)])),
    c(table[2,6],sum(table[2,c(1,2,3,4,5)])),
    c(table[3,6],sum(table[3,c(1,2,3,4,5)])),
    c(table[4,6],sum(table[4,c(1,2,3,4,5)])))

chisq.test(table1)
chisq.test(table2)
chisq.test(table3)
chisq.test(table4)
chisq.test(table5)
chisq.test(table6)


# CPG
table=rbind(mam_cpg,bir_cpg,fish_cpg,rep_cpg)
chisq.test(table)
##X-squared = 4.6411, df = 3, p-value = 0.2
# %
# Mam 23.4852
# birds 21.88552
# fish 25.8427
# reptiles 18.92655


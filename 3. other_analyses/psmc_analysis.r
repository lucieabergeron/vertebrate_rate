patern = read.csv("patern_a.txt", sep="\t")
# Per generation rate
data_g = read.csv("../all_sp.txt", na.strings=c("","NA"), sep="\t")
ave_g=vector()
for(sp in unique(data_g$species_name)){
  g_ave=mean(as.numeric(as.character(data_g[which(data_g$species_name==sp),10])))
  ave_g=rbind(ave_g, c(sp, g_ave))
}
grate = as.numeric(ave_g[order(ave_g[,1]),2])
#

full_data = cbind(patern, grate)

# Mean over the whole period:
max_mean_hmean = vector()
for(i in seq(1,68)){
  sp = full_data$species[i]
  time = read.table(paste0(sp,"/estimate_t.txt"), sep=",")
  time_1 = replace(time[1,],1,sub(".", "", time[1,1]))
  time_2 = as.numeric(time_1)[-length(time_1)]
  idx = which(time_2<1000000 & time_2>30000)
  time_3 = time_2[idx]
  ne = read.table(paste0(sp,"/estimate.txt"), sep=",")
  ne_1 = replace(ne[1,],1,sub(".", "", ne[1,1]))
  ne_2 = as.numeric(ne_1)[-length(ne_1)]
  ne_3 = (as.numeric(ne_2[idx]))
  line = c(max(ne_3), mean(ne_3), length(ne_3)/sum(ne_3^(-1)))
  max_mean_hmean=rbind(max_mean_hmean,line)
}


cor.test(full_data$grate, max_mean_hmean[,1])
cor.test(full_data$grate, max_mean_hmean[,2])
cor.test(full_data$grate, max_mean_hmean[,3])

# With phylo
library(caper)
tree<-read.tree("../phylogeny/tree_calibrated_UCE.nwk")
tree_new = drop.tip(tree, c("Macaca_mulatta","Elephas_maximus", "Cercocebus_lunulatus", "Muntiacus_reevesi", "Oryzias_latipes", "Callithrix_jacchus"))

tree_mam = drop.tip(tree_new, c("Aptenodytes_forsteri","Ara_glaucogularis","Bubo_scandiacus","Chauna_torquata","Coturnix_japonica","Cyanistes_caeruleus","Gallus_gallus","Gyps_fulvus","Larus_argentatus","Larus_marinus","Pelecanus_crispus","Phoenicopterus_roseus","Platalea_ajaja","Pygoscelis_adeliae","Rhea_pennata","Saxicola_maurus","Taeniopygia_guttata","Turdus_merula","Amphiprion_ocellaris","Betta_splendens","Cynoglossus_semilaevis","Cyprinus_carpio","Larimichthys_crocea","Paralichthys_olivaceus","Salmo_salar","Syngnathus_scovelli","Chrysemys_picta","Coleonyx_brevis","Eublepharis_macularius","Pogona_vitticeps","Sphaerodactylus_macrolepis","Thamnophis_sirtalis"))


tree_bird = drop.tip(tree_new, c("Amphiprion_ocellaris","Betta_splendens","Cynoglossus_semilaevis","Cyprinus_carpio","Larimichthys_crocea","Paralichthys_olivaceus","Salmo_salar","Syngnathus_scovelli","Ailurus_fulgens","Arctocephalus_gazella","Callithrix_jacchus","Canis_lupus_familiaris","Capra_hircus","Cavia_aperea","Ceratotherium_simum_simum","Cervus_elaphus_yarkandensis","Cervus_nippon","Felis_catus","Fukomys_damarensis","Giraffa_camelopardalis","Hippopotamus_amphibius","Homo_sapiens","Hylobates_lar","Mandrillus_leucophaeus","Monodelphis_domestica","Moschus_berezovskii","Mus_musculus","Neovison_vison","Odobenus_rosmarus","Orcinus_orca","Pan_troglodytes","Panthera_pardus","Panthera_tigris","Pithecia_pithecia","Procavia_capensis","Rangifer_tarandus","Rousettus_aegyptiacus","Saimiri_boliviensis_boliviensis","Sarcophilus_harrisii","Sus_scrofa","Tapirus_indicus","Tupaia_belangeri","Tursiops_truncatus","Vicugna_pacos","Vulpes_vulpes","Chrysemys_picta","Coleonyx_brevis","Eublepharis_macularius","Pogona_vitticeps","Sphaerodactylus_macrolepis","Thamnophis_sirtalis"))

#####

matrix_phylo = cbind(full_data, max_mean_hmean[,1]/10000, max_mean_hmean[,2]/10000, max_mean_hmean[,3])
row.names(matrix_phylo) = full_data$species
colnames(matrix_phylo) = c(names(full_data),"max", "mean", "hmean")
sub_data = as.data.frame(matrix_phylo)
# na omit
comp.data<-comparative.data(tree_new, sub_data, names.col=species, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)

model<-pgls(grate~max, data=comp.data)
summary(model)
model<-pgls(grate~mean, data=comp.data)
summary(model)
model<-pgls(grate~hmean, data=comp.data)
summary(model)


# Mammals
matrix_phylo = cbind(full_data, max_mean_hmean[,1]/10000, max_mean_hmean[,2]/10000, max_mean_hmean[,3])
row.names(matrix_phylo) = full_data$species
colnames(matrix_phylo) = c(names(full_data),"max", "mean", "hmean")
sub_data = as.data.frame(matrix_phylo[which(as.character(matrix_phylo$col)=='firebrick3'),])
comp.data<-comparative.data(tree_mam, sub_data, names.col=species, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(grate~hmean, data=comp.data)
summary(model)
#Adjusted R-squared: 0.1239,  p-value: 0.02008, 8.2912e-09,-7.6338e-15
# Birds
sub_data = as.data.frame(matrix_phylo[which(as.character(matrix_phylo$col)=='goldenrod2'),])
comp.data<-comparative.data(tree_bird, sub_data, names.col=species, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(grate~hmean, data=comp.data)
summary(model)
# Not

png("grate_harmonic_mean.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(max_mean_hmean[,3], full_data$grate, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5)
points(max_mean_hmean[,3], full_data$grate, pch=19, col=as.character(full_data$col), cex=4)
text(c(1350000,1350000,1350000), c(3.8e-8,3.5e-8,3.2e-8), labels = c("Evolutionary correlation",expression('R'^2*' = 0.061'), "p.value = 0.024"), font=2, cex=3)
##abline(1.0500e-08, -4.8123e-15,lwd=4)
text(c(1350000,1350000,1350000), c(2.6e-8,2.3e-8,2e-8), labels = c("Evolutionary correlation mammals",expression('R'^2*' = 0.124'), "p.value = 0.020"), font=2, cex=3, col='firebrick3')
axis(2, at = c(0, 1e-8, 2e-8, 3e-8, 4e-8), labels = c(0,1,2,3,4), cex.axis=3.5, las=2)
mtext("Ne harmonic mean", 1, line=5, cex=4)
mtext(expression('Mutation rate per site'), 2, line=8, cex=4)
mtext(expression('per generation (10'^-8*')'), 2, line=4, cex=4)
dev.off()


# Harmonic mean over the recent period (130,000-30,000)
hmean_recent = vector()
for(i in seq(1,68)){
  sp = full_data$species[i]
  time = read.table(paste0(sp,"/estimate_t.txt"), sep=",")
  time_1 = replace(time[1,],1,sub(".", "", time[1,1]))
  time_2 = as.numeric(time_1)[-length(time_1)]
  idx = which(time_2<150000 & time_2>30000)
  time_3 = time_2[idx]
  ne = read.table(paste0(sp,"/estimate.txt"), sep=",")
  ne_1 = replace(ne[1,],1,sub(".", "", ne[1,1]))
  ne_2 = as.numeric(ne_1)[-length(ne_1)]
  ne_3 = (as.numeric(ne_2[idx]))
  line = length(ne_3)/sum(ne_3^(-1))
  hmean_recent=c(hmean_recent,line)
}


cor.test(full_data$grate, hmean_recent)

matrix_phylo = cbind(full_data, hmean_recent)
row.names(matrix_phylo) = full_data$species
colnames(matrix_phylo) = c(names(full_data),"hmean")
sub_data = as.data.frame(matrix_phylo)
# na omit
comp.data<-comparative.data(tree_new, sub_data, names.col=species, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)

model<-pgls(grate~hmean, data=comp.data)
summary(model)

#### Not significant

# Pente over the rencent time.
patern_b = read.csv("patern_b.txt", sep="\t")
full_data_b = cbind(full_data, patern_b[,2:4])

plot(full_data_b$slope, full_data_b$grate)
plot(full_data_b$shape, full_data_b$grate)
plot(full_data_b$Half_decreased, full_data_b$grate)

matrix_phylo = cbind(full_data_b)
row.names(matrix_phylo) = full_data_b$species
colnames(matrix_phylo) = names(full_data_b)
sub_data = as.data.frame(matrix_phylo)
# na omit
comp.data<-comparative.data(tree_new, sub_data, names.col=species, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)

model<-pgls(grate~slope, data=comp.data)
summary(model)
model<-pgls(grate~shape, data=comp.data)
summary(model)
model<-pgls(grate~Half_decreased, data=comp.data)
summary(model)
# this yes
# Mammals
sub_data = as.data.frame(matrix_phylo[which(as.character(matrix_phylo$col)=='firebrick3'),])
comp.data<-comparative.data(tree_mam, sub_data, names.col=species, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(grate~Half_decreased, data=comp.data)
summary(model)
#Adjusted R-squared: 0.1903,  p-value: 0.004565
# Birds
sub_data = as.data.frame(matrix_phylo[which(as.character(matrix_phylo$col)=='goldenrod2'),])
comp.data<-comparative.data(tree_bird, sub_data, names.col=species, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(grate~Half_decreased, data=comp.data)
summary(model)
# Not


x_half = factor(full_data_b$Half_decreased, levels=c("y", "n"))
#png("grate_half_decreased.png", width = 1300, height = 900)
#par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
#plot(x_half, full_data_b$grate, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5)
#points(x_half, full_data_b$grate, pch=19, col=as.character(full_data_b$col), cex=4)
#stat
#lines(c(1,2),c(3.20e-08,3.20e-08), lwd=3)
#lines(c(1,1),c(3.20e-08,3.15e-08), lwd=3)
#lines(c(2,2),c(3.20e-08,3.15e-08), lwd=3)
#text(1.5,3.30e-8, labels = c("*"), cex=4)
#text(1.5,3.45e-8, labels = expression('R'^2*' = 0.062, p.value = 0.022'), font=3, cex=2)
#axis(2, at = c(0, 1e-8, 2e-8, 3e-8, 4e-8), labels = c(0,1,2,3,4), cex.axis=3.5, las=2)
#mtext("Reduction of Ne over recent period (30,000 to 130,000 years)", 1, line=5, cex=4)
#mtext(expression('Mutation rate per site per generation (10'^-8*')'), 2, line=5, cex=4)
#dev.off()


png("grate_half_decreased.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
boxplot(full_data_b$grate[which(full_data_b$Half_decreased=="y")],
   full_data_b$grate[which(full_data_b$Half_decreased=="n")],
   full_data_b$grate[which(full_data_b$Half_decreased=="y" & as.character(full_data_b$col) =='firebrick3')],
   full_data_b$grate[which(full_data_b$Half_decreased=="n" & as.character(full_data_b$col) =='firebrick3')],
   names=rep(c("yes","no"), 2),
   yaxt="n", cex.axis=3)
points(x_half, full_data_b$grate, pch=19, col=as.character(full_data_b$col), cex=4)
points(rep(3,24), full_data_b$grate[which(full_data_b$Half_decreased=="y" & as.character(full_data_b$col) =='firebrick3')], col= 'firebrick3', cex=4, font=2, pch=19)
points(rep(4,12), full_data_b$grate[which(full_data_b$Half_decreased=="n" & as.character(full_data_b$col) =='firebrick3')], col= 'firebrick3', cex=4, font=2, pch=19)
abline(v=2.5)
#
lines(c(1,2),c(3.20e-08,3.20e-08), lwd=3)
lines(c(1,1),c(3.20e-08,3.15e-08), lwd=3)
lines(c(2,2),c(3.20e-08,3.15e-08), lwd=3)
text(1.5,3.30e-8, labels = c("*"), cex=4)
text(1.5,3.45e-8, labels = c("p.value = 0.022"), font=1, cex=3)
text(1.5,3.65e-8, labels = c("All species"), font=1, cex=3)
#
lines(c(3,4),c(1.66e-08,1.66e-08), lwd=3)
lines(c(3,3),c(1.66e-08,1.61e-08), lwd=3)
lines(c(4,4),c(1.66e-08,1.61e-08), lwd=3)
text(3.5,1.77e-08, labels = c("**"), cex=4)
text(3.5,1.92e-8, labels = c("p.value = 0.0046"), font=1, cex=3)
text(3.5,2.12e-8, labels = c("Mammals"), font=1, cex=3)
#
axis(2, at = c(0, 1e-8, 2e-8, 3e-8, 4e-8), labels = c(0,1,2,3,4), cex.axis=3.5, las=2)
mtext("Recent significant Ne decline", 1, line=5, cex=4)
mtext(expression('Mutation rate per site'), 2, line=8, cex=4)
mtext(expression('per generation (10'^-8*')'), 2, line=4, cex=4)
dev.off()



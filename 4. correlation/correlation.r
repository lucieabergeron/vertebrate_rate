# source activate caper
library(phytools)
library(nlme)
library(caper)
tree<-read.tree("../../phylogeny/tree_calibrated_UCE.nwk")
tree_new = drop.tip(tree, c("Macaca_mulatta","Elephas_maximus", "Cercocebus_lunulatus", "Muntiacus_reevesi", "Oryzias_latipes", "Callithrix_jacchus"))
data = read.csv("data_per_sp.txt", na.strings=c("","NA"), sep="\t")

tree_mam = drop.tip(tree_new, c("Aptenodytes_forsteri","Ara_glaucogularis","Bubo_scandiacus","Chauna_torquata","Coturnix_japonica","Cyanistes_caeruleus","Gallus_gallus","Gyps_fulvus","Larus_argentatus","Larus_marinus","Pelecanus_crispus","Phoenicopterus_roseus","Platalea_ajaja","Pygoscelis_adeliae","Rhea_pennata","Saxicola_maurus","Taeniopygia_guttata","Turdus_merula","Amphiprion_ocellaris","Betta_splendens","Cynoglossus_semilaevis","Cyprinus_carpio","Larimichthys_crocea","Paralichthys_olivaceus","Salmo_salar","Syngnathus_scovelli","Chrysemys_picta","Coleonyx_brevis","Eublepharis_macularius","Pogona_vitticeps","Sphaerodactylus_macrolepis","Thamnophis_sirtalis"))


tree_bird = drop.tip(tree_new, c("Amphiprion_ocellaris","Betta_splendens","Cynoglossus_semilaevis","Cyprinus_carpio","Larimichthys_crocea","Paralichthys_olivaceus","Salmo_salar","Syngnathus_scovelli","Ailurus_fulgens","Arctocephalus_gazella","Callithrix_jacchus","Canis_lupus_familiaris","Capra_hircus","Cavia_aperea","Ceratotherium_simum_simum","Cervus_elaphus_yarkandensis","Cervus_nippon","Felis_catus","Fukomys_damarensis","Giraffa_camelopardalis","Hippopotamus_amphibius","Homo_sapiens","Hylobates_lar","Mandrillus_leucophaeus","Monodelphis_domestica","Moschus_berezovskii","Mus_musculus","Neovison_vison","Odobenus_rosmarus","Orcinus_orca","Pan_troglodytes","Panthera_pardus","Panthera_tigris","Pithecia_pithecia","Procavia_capensis","Rangifer_tarandus","Rousettus_aegyptiacus","Saimiri_boliviensis_boliviensis","Sarcophilus_harrisii","Sus_scrofa","Tapirus_indicus","Tupaia_belangeri","Tursiops_truncatus","Vicugna_pacos","Vulpes_vulpes","Chrysemys_picta","Coleonyx_brevis","Eublepharis_macularius","Pogona_vitticeps","Sphaerodactylus_macrolepis","Thamnophis_sirtalis"))


matrix_phylo = data
row.names(matrix_phylo) = data$species_name
colnames(matrix_phylo) = names(data)
sub_data = as.data.frame(matrix_phylo)
# na omit
comp.data<-comparative.data(tree_new, sub_data, names.col=species_name, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)


#################################################################
# Correlation
type_col = vector()
type_col[which(data$Group=="Mammal")]="firebrick3"
type_col[which(data$Group=="Bird")]="goldenrod2"
type_col[which(data$Group=="Fish")]="royalblue4"
type_col[which(data$Group=="Reptile")]="limegreen"

## Genome size
##png("y_rate_genome_size.png", width = 1300, height = 900)
par(mar=c(8,11,4,2), mgp=c(3, 1.5, 0))
plot(data$genome_size, data$mu, type="n", xlab="", ylab="", cex=3, xaxt="n", yaxt="n")
text(data$genome_size, data$mu, labels=data$species_name, col=type_col, cex=1.5, font=2)
axis(1, cex.axis=2.5)
axis(2, cex.axis=2.5, las=2)
mtext("Genome size", 1, line=5, cex=3)
mtext("Mutation rate (per year)", 2, line=8, cex=3)
##dev.off()

cor.test(data$genome_size, data$mu)
## Phylogeny correction
matrix_phylo = data
row.names(matrix_phylo) = data$species_name
colnames(matrix_phylo) = names(data)
#
fit.yx<-gls(mu~genome_size,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)

gen_size = matrix_phylo$genome_size/1000000000
sub_genome = as.data.frame(cbind(matrix_phylo,gen_size))
# na omit
comp.data<-comparative.data(tree_new, sub_genome, names.col=species_name, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)

model<-pgls(mu~gen_size, data=comp.data)
summary(model)
# Error

## Group
##png("yrate_group.png", width = 1300, height = 1200)
png("yrate_group.png", width = 900, height = 1200)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
##plot(data$Group, data$mu, type="n", xlab="", ylab="", cex=1.5, yaxt="n", cex.axis=2.5, ylim=c(-1e-8,2.2e-8))
plot(data$Group, data$mu, type="n", xlab="", ylab="", cex=2, yaxt="n", cex.axis=3, pch=19, outcol=c("goldenrod2","royalblue4","firebrick3","limegreen"), cex=2)
points(data$Group, data$mu, pch=19, col=type_col, cex=6)
##data_label = data[c(53,18,62,6,61,60,11,17,22,49,59,64,16),]
##text(data_label$Group, data_label$mu+0.05e-08, labels=data_label$species_name, col=type_col[c(53,18,62,6,61,60,11,17,22,49,59,64,16)], cex=1.5, font=2)
axis(2, at = c(-1e-8, -5e-9, 0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(-1.0,-0.5, 0.0, 0.5,1.0,1.5,2.0), cex.axis=3, las=2)
mtext("Group", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()

x = factor(data$Group, levels=c("Mammal", "Bird", "Fish", "Reptile"))
png("yrate_group.png", width = 900, height = 1200)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x, data$mu, type="n", xlab="", ylab="", cex=2, yaxt="n", cex.axis=3, pch=19, outcol=c("firebrick3","goldenrod2","royalblue4","limegreen"), cex=2)
points(x, data$mu, pch=19, col=type_col, cex=6)
axis(2, at = c(-1e-8, -5e-9, 0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(-1.0,-0.5, 0.0, 0.5,1.0,1.5,2.0), cex.axis=3, las=2)
mtext("Group", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()


#
#But difference between mono and poly
x1 = factor(data$mating_b, levels=c("monogamous", "poly"))
t.test(data$mu[which(x1=="monogamous")],data$mu[which(x1=="poly")])
png("y_rate_mating_mono_poly.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x1, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5, ylim=c(-3.820248e-10,1.5e-8))
points(x1, data$mu, pch=19, col=type_col, cex=4)
lines(c(1,2),c(1.36e-08,1.36e-08), lwd=3)
lines(c(1,1),c(1.36e-08,1.31e-08), lwd=3)
lines(c(2,2),c(1.36e-08,1.31e-08), lwd=3)
text(1.5,1.4e-8, labels = c("ns"), cex=3)
text(1.5,1.45e-8, labels = c("p.value = 0.253"), font=2, cex=3)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
mtext("Mating system", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()

fit.yx<-gls(mu~mating_b,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)
sub_m=matrix_phylo[which(matrix_phylo$Group=="Mammal"),]
fit.yx<-gls(mu~mating_b,data=as.data.frame(sub_m),correlation=corBrownian(1,tree_mam))
anova(fit.yx)
sub_b=matrix_phylo[which(matrix_phylo$Group=="Bird"),]
fit.yx<-gls(mu~mating_b,data=as.data.frame(sub_b),correlation=corBrownian(1,tree_bird))
anova(fit.yx)

model<-pgls(mu~mating_b, data=comp.data)
summary(model)
sub_m_m = as.data.frame(sub_m)
comp.data_m<-comparative.data(tree_mam, sub_m_m, names.col=species_name, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(mu~mating_b, data=comp.data_m)
summary(model)
sub_b_b = as.data.frame(sub_b)
comp.data_b<-comparative.data(tree_bird, sub_b_b, names.col=species_name, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(mu~mating_b, data=comp.data_b)
summary(model)





####### K and r #########################################################
t.test(data$mu[which(data$r_K=="r")],data$mu[which(data$r_K=="K")])
t.test(data$mu[which(data$Group=="Mammal" & data$r_K=="r")],data$mu[which(data$Group=="Mammal" & data$r_K=="K")])
t.test(data$mu[which(data$Group=="Bird" & data$r_K=="r")],data$mu[which(data$Group=="Bird" & data$r_K=="K")])
##
only_b = data$r_K[which(data$Group=="Bird")]
only_m = data$r_K[which(data$Group=="Mammal")]
xaxis_b=vector()
xaxis_b[which(only_b =="K")] <-5
xaxis_b[which(only_b =="r")] <-6
xaxis_m=vector()
xaxis_m[which(only_m =="K")] <-3
xaxis_m[which(only_m =="r")] <-4
#
png("y_rate_mating_combine_r_K.png", width = 900, height = 900)
par(mar=c(8,9,4,2), mgp=c(3, 1.5, 0))
boxplot(data$mu[which(data$r_K=="K")], data$mu[which(data$r_K=="r")],
    data$mu[which(data$Group=="Mammal" & data$r_K=="K")],data$mu[which(data$Group=="Mammal" & data$r_K=="r")],
    data$mu[which(data$Group=="Bird" & data$r_K=="K")],data$mu[which(data$Group=="Bird" & data$r_K=="r")],
    names=rep(c("K","r"), 3),
    yaxt="n", cex.axis=2,
    ylim=c(-3.820248e-10,1.5e-8))
points(data$r_K, data$mu, col=type_col, cex=3, font=2, pch=20)
points(xaxis_m, data$mu[which(data$Group=="Mammal")], col=type_col[which(data$Group=="Mammal")], cex=3, font=2, pch=20)
points(xaxis_b, data$mu[which(data$Group=="Bird")], col=type_col[which(data$Group=="Bird")], cex=3, font=2, pch=20)
abline(v=2.5)
abline(v=4.5)
#
lines(c(1,2),c(1.36e-08,1.36e-08), lwd=3)
lines(c(1,1),c(1.36e-08,1.31e-08), lwd=3)
lines(c(2,2),c(1.36e-08,1.31e-08), lwd=3)
text(1.5,1.4e-8, labels = c("***"), cex=3)
text(1.5,1.45e-8, labels = c("p.value = 0.0002"), font=2, cex=2)
#
lines(c(3,4),c(0.96e-08,0.96e-08), lwd=3)
lines(c(3,3),c(0.96e-08,0.91e-08), lwd=3)
lines(c(4,4),c(0.96e-08,0.91e-08), lwd=3)
text(3.5,1e-08, labels = c("***"), cex=2)
text(3.5,1.05e-8, labels = c("p.value = 0.0002"), font=2, cex=2)
#
lines(c(5,6),c(1.16e-08,1.16e-08), lwd=3)
lines(c(5,5),c(1.16e-08,1.11e-08), lwd=3)
lines(c(6,6),c(1.16e-08,1.11e-08), lwd=3)
text(5.5,1.2e-8, labels = c("ns"), cex=2)
text(5.5,1.25e-8, labels = c("p.value = 0.285"), font=2, cex=2)
axis(2, at = c(0, 0.5e-8,1e-8,1.5e-8,2e-8), labels = c("0", "0.5","1","1.5","2"), cex.axis=2.5, las=2)
mtext("Strategy", 1, line=5, cex=3)
mtext(expression("Mutation rate per site per year (10"^-8*")"), 2, line=5, cex=3)
dev.off()

fit.yx<-gls(mu~r_K,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)
fit.yx<-gls(mu~r_K,data=as.data.frame(sub_m),correlation=corBrownian(1,tree_mam))
anova(fit.yx)
fit.yx<-gls(mu~r_K,data=as.data.frame(sub_b),correlation=corBrownian(1,tree_bird))
anova(fit.yx)



lines(c(1,2),c(9.46e-09,1.46e-08), lwd=3)
lines(c(1,1),c(9.46e-09,1.41e-08), lwd=3)
lines(c(2,2),c(9.46e-09,1.41e-08), lwd=3)
text(1.5,1.4e-8, labels = c("ns"), cex=2)
text(1.5,1.45e-8, labels = c("p.value = 0.378"), font=2, cex=2)
#
lines(c(1,2),c(1.16e-08,1.16e-08), lwd=3)
lines(c(1,1),c(1.16e-08,1.11e-08), lwd=3)
lines(c(2,2),c(1.16e-08,1.11e-08), lwd=3)
text(1.5,1.4e-8, labels = c("ns"), cex=2)
text(1.5,1.45e-8, labels = c("p.value = 0.202"), font=2, cex=2)



# Separate by all


## Reproduction
##t.test(data$mu[which(data$r_K=="r")],data$mu[which(data$r_K=="K")])
png("y_rate_repro.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$reproduction, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5, ylim=c(-3.820248e-10,1.5e-8))
points(data$reproduction, data$mu, pch=19, col=type_col, cex=4)
##text(1, 5.3e-8, labels = c("t.test_p.value = 1.23e-06"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Reproduction", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

fit.yx<-gls(mu~reproduction,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)

model<-pgls(mu~reproduction, data=comp.data)
summary(model)

# Latitude
x3 = factor(data$latitude, levels=c("high", "low"))
t.test(data$mu[which(data$latitude=="high")],data$mu[which(data$latitude=="low")])
png("y_rate_lat.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$latitude_b, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5, ylim=c(-3.820248e-10,1.5e-8))
points(data$latitude_b, data$mu, pch=19, col=type_col, cex=4)
lines(c(1,2),c(1.1e-08,1.1e-08), lwd=3)
lines(c(1,1),c(1.05e-08,1.1e-08), lwd=3)
lines(c(2,2),c(1.05e-08,1.1e-08), lwd=3)
text(1.5,1.14e-8, labels = c("ns"), cex=3)
text(1.5,1.19e-8, labels = c("p.value = 0.424"), font=2, cex=3)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
mtext("Latitude", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()

model<-pgls(mu~latitude_b, data=comp.data)
summary(model)



## Habitat
t.test(data$mu[which(data$habitat=="aquatic")],data$mu[which(data$habitat=="land")])
png("y_rate_habitat.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$habitat, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5, ylim=c(-3.820248e-10,1.5e-8))
points(data$habitat, data$mu, pch=19, col=type_col, cex=4)
lines(c(1,2),c(1.36e-08,1.36e-08), lwd=3)
lines(c(1,1),c(1.36e-08,1.31e-08), lwd=3)
lines(c(2,2),c(1.36e-08,1.31e-08), lwd=3)
text(1.5,1.4e-8, labels = c("ns"), cex=3)
text(1.5,1.45e-8, labels = c("p.value = 0.432"), font=2, cex=3)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
mtext("Habitat", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()

fit.yx<-gls(mu~habitat,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)
model<-pgls(mu~habitat, data=comp.data)
summary(model)


## Mass
cor.test(data$mass, data$mu)
png("y_rate_mass_log.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$mass, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5, log="x")
points(data$mass, data$mu, pch=19, col=type_col, cex=4)
text(c(310000,310000), c(1.18e-8,1.25e-8), labels = c("p.value = 0.490", "Evolutionary correlation"), font=1, cex=3)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
mtext("Mass (log scale)", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()

fit.yx<-gls(mu~mass,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)
summary(fit.yx)

obj<-phyl.vcv(cbind(matrix_phylo$mu,matrix_phylo$mass),vcv(tree_new),1)
obj$R
r.xz<-cov2cor(obj$R)

model<-pgls(mu~mass, data=comp.data)
summary(model)


## Longevity
cor.test(data$lifespan_wild, data$mu)
png("y_rate_longevity.png", width = 1000, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$lifespan_wild, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5)
points(data$lifespan_wild, data$mu, pch=19, col=type_col, cex=4)
text(c(42,42,42), c(1.25e-8,1.18e-8,1.09e-8), labels = c("Evolutionary correlation", expression('R'^2*' = 0.068'), "p.value = 0.031"), font=2, cex=3)
##abline(4.3098e-09, -3.7447e-11,lwd=4)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
mtext("Life span (in the wild)", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()

#
fit.yx<-gls(mu~lifespan_wild,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)
#i more strict
##lPic <- pic(data$lifespan_wild, tree_new)
##mPic <- pic(data$mu, tree_new)
##picModel <- lm(mPic ~ lPic)
##summary(picModel)
## our covariance matrix
obj<-phyl.vcv(cbind(matrix_phylo$mu,matrix_phylo$lifespan_wild),vcv(tree_new),1)
obj$R
r.xz<-cov2cor(obj$R)

model<-pgls(mu~lifespan_wild, data=comp.data)
summary(model)

model<-pgls(mu~lifespan_wild, data=comp.data)
summary(model)

# Generation time
cor.test(data$generation_time, data$mu)
png("y_rate_gt.png", width = 1000, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$generation_time, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5)
points(data$generation_time, data$mu, pch=19, col=type_col, cex=4)
text(c(20,20,20), c(1.25e-8,1.18e-8,1.09e-8), labels = c("Evolutionary correlation",expression('R'^2*' = 0.106'),"p.value = 0.007"), font=2, cex=3)
##abline(4.4494e-09,-1.1847e-10,lwd=4)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
mtext("Generation time", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()

fit.yx<-gls(mu~generation_time,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)
obj<-phyl.vcv(cbind(matrix_phylo$mu,matrix_phylo$generation_time),vcv(tree_new),1)
obj$R
r.xz<-cov2cor(obj$R)
summary(fit.yx)

model<-pgls(mu~generation_time, data=comp.data)
summary(model)


# Age used
#png("y_rate_gt.png", width = 1300, height = 900)
#par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
#plot(data$generation_time, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5)
#points(data$generation_time, data$mu, pch=19, col=type_col, cex=4)
#text(c(20,20,20), c(1.25e-8,1.18e-8,1.11e-8), labels = c("Evolutionary correlation",expression('R'^2*' = 0.106'),"p.value = 0.007"), font=2, cex=3)
###abline(4.4494e-09,-1.1847e-10,lwd=4)
#axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
#mtext("Generation time", 1, line=5, cex=4)
#mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
#dev.off()

model<-pgls(mu~age_used, data=comp.data)
summary(model)



# Maturation time
cor.test((data$maturation_time_mother+data$maturation_time_father)/2, data$mu)
png("y_rate_maturation_time.png", width = 1000, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot((data$maturation_time_mother+data$maturation_time_father)/2, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5)
points((data$maturation_time_mother+data$maturation_time_father)/2, data$mu, pch=19, col=type_col, cex=4)
text(c(8,8,8), c(1.25e-8,1.18e-8,1.09e-8), labels = c("Evolutionary correlation",expression('R'^2*' = 0.112'), "p.value = 0.005"), font=2, cex=3)
##abline(4.308373e-09, -2.512330e-10,lwd=4)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
mtext("Maturation time", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=4)
dev.off()

mat=(data$maturation_time_mother+data$maturation_time_father)/2
more = cbind(matrix_phylo,mat)
fit.yx<-gls(mu~mat,data=as.data.frame(more),correlation=corBrownian(1,tree_new))
anova(fit.yx)
obj<-phyl.vcv(cbind(more$mu,more$mat),vcv(tree_new),1)
obj$R
r.xz<-cov2cor(obj$R)
plot(fit.yx)    # residual plot
summary(fit.yx) # coefficients, SE's, fit statistics
confint(fit.yx) # confidence intervals for coefficients

more = as.data.frame(cbind(matrix_phylo,mat))
comp.data_more<-comparative.data(tree_new, more, names.col=species_name, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(mu~mat, data=comp.data_more)
summary(model)

# Nb offspring
t.test(data$mu[which(data$Offspring_per_batch=="few")],data$mu[which(data$Offspring_per_batch=="lot")])
png("y_rate_nboff.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$Offspring_per_batch, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5, ylim=c(-3.820248e-10,1.5e-8))
points(data$Offspring_per_batch, data$mu, pch=19, col=type_col, cex=4)
#stat
lines(c(1,2),c(1.36e-08,1.36e-08), lwd=3)
lines(c(1,1),c(1.36e-08,1.31e-08), lwd=3)
lines(c(2,2),c(1.36e-08,1.31e-08), lwd=3)
text(1.5,1.4e-8, labels = c("***"), cex=3)
text(1.5,1.45e-8, labels = expression('p.value = 2.5e-7, R'^2*' = 0.330'), font=2, cex=2)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=2.5, las=2)
mtext("Number of offspring per batch", 1, line=5, cex=3)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=5, cex=3)
dev.off()

##
only_b = data$offspring_perg[which(data$Group=="Bird")]
only_m = data$offspring_perg[which(data$Group=="Mammal")]
xaxis_b=vector()
xaxis_b[which(only_b =="few")] <-5
xaxis_b[which(only_b =="lot")] <-6
xaxis_m=vector()
xaxis_m[which(only_m =="few")] <-3
xaxis_m[which(only_m =="lot")] <-4

png("y_rate_offspring_combine.png", width = 1100, height = 900)
par(mar=c(8,9,4,2), mgp=c(3, 1.5, 0))
boxplot(data$mu[which(data$offspring_perg=="few")], data$mu[which(data$offspring_perg=="lot")],
    data$mu[which(data$Group=="Mammal" & data$offspring_perg=="few")],data$mu[which(data$Group=="Mammal" & data$offspring_perg=="lot")],
    data$mu[which(data$Group=="Bird" & data$offspring_perg=="few")],data$mu[which(data$Group=="Bird" & data$offspring_perg=="lot")],
    names=rep(c("Few","Lot"), 3),
    yaxt="n", cex.axis=3,
    ylim=c(-3.820248e-10,1.3e-8))
points(data$offspring_perg, data$mu, col=type_col, cex=3, font=2, pch=20)
points(xaxis_m, data$mu[which(data$Group=="Mammal")], col=type_col[which(data$Group=="Mammal")], cex=3, font=2, pch=20)
points(xaxis_b, data$mu[which(data$Group=="Bird")], col=type_col[which(data$Group=="Bird")], cex=3, font=2, pch=20)
abline(v=2.5)
abline(v=4.5)
#
lines(c(1,2),c(1.07e-08,1.07e-08), lwd=3)
lines(c(1,1),c(1.07e-08,1.02e-08), lwd=3)
lines(c(2,2),c(1.07e-08,1.02e-08), lwd=3)
text(1.5,1.11e-8, labels = c("***"), cex=4)
text(1.5,1.16e-8, labels = c("p.value = 4.3e-5"), font=1, cex=3)
text(1.5,1.23e-8, labels = c("All species"), font=1, cex=3)
#
lines(c(3,4),c(0.96e-08,0.96e-08), lwd=3)
lines(c(3,3),c(0.96e-08,0.91e-08), lwd=3)
lines(c(4,4),c(0.96e-08,0.91e-08), lwd=3)
text(3.5,1e-08, labels = c("***"), cex=4)
text(3.5,1.05e-8, labels = c("p.value = 0.0002"), font=1, cex=3)
text(3.5,1.12e-8, labels = c("Mammals"), font=1, cex=3)
#
lines(c(5,6),c(0.96e-08,0.96e-08), lwd=3)
lines(c(5,5),c(0.96e-08,0.91e-08), lwd=3)
lines(c(6,6),c(0.96e-08,0.91e-08), lwd=3)
text(5.5,1e-8, labels = c("ns"), cex=3)
text(5.5,1.05e-8, labels = c("p.value = 0.173"), font=1, cex=3)
text(5.5,1.12e-8, labels = c("Birds"), font=1, cex=3)
axis(2, at = c(0, 0.5e-8,1e-8,1.5e-8,2e-8), labels = c("0", "0.5","1","1.5","2"), cex.axis=3, las=2)
mtext("Number of offspring per generation", 1, line=5, cex=3.5)
mtext(expression("Mutation rate per site per year (10"^-8*")"), 2, line=5, cex=3.5)
dev.off()



# Dummy variable 0 or 1 work the same
fit.yx<-gls(mu~Offspring_per_batch,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)
summary(fit.yx)
model<-pgls(mu~Offspring_per_batch, data=comp.data)
summary(model)
model<-pgls(mu~Offspring_per_batch, data=comp.data_m)
summary(model)
model<-pgls(mu~Offspring_per_batch, data=comp.data_b)
summary(model)

model<-pgls(mu~offspring_perg, data=comp.data)
summary(model)
model<-pgls(mu~offspring_perg, data=comp.data_m)
summary(model)
model<-pgls(mu~offspring_perg, data=comp.data_b)
summary(model)


##
Danger
x4 = factor(data$status, levels=c("dd", "lc", "nt", "vu", "en", "cr"))
png("y_rate_status.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x4, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
points(x4, data$mu, col=type_col, cex=1.5, font=2, pch=19)
axis(2, cex.axis=2.5, las=2)
mtext("Status", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

fit.yx<-gls(mu~status,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)
model<-pgls(mu~status, data=comp.data)
summary(model)



## Domestication
x5 = factor(data$domesticated, levels=c("yes", "no"))
t.test(data$mu[which(data$domesticated=="yes")],data$mu[which(data$domesticated=="no")])
png("y_rate_domestication.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x5, data$mu, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=3.5, ylim=c(-4.128715e-10,1.55e-8))
points(x5, data$mu, pch=19, col=type_col, cex=4)
#stat
lines(c(1,2),c(1.36e-08,1.36e-08), lwd=3)
lines(c(1,1),c(1.36e-08,1.31e-08), lwd=3)
lines(c(2,2),c(1.36e-08,1.31e-08), lwd=3)
text(1.5,1.4e-8, labels = c("***"), cex=4)
text(1.5,1.5e-8, labels = 'p.value = 0.0004', font=1, cex=3)
axis(2, at = c(0, 5e-9, 1e-8, 1.5e-8, 2e-8), labels = c(0.0, 0.5,1.0,1.5,2.0), cex.axis=3.5, las=2)
mtext("Domesticated", 1, line=5, cex=4)
mtext(expression('Mutation rate per site per year (10'^-8*')'), 2, line=6, cex=4)
dev.off()

fit.yx<-gls(mu~domesticated,data=as.data.frame(matrix_phylo),correlation=corBrownian(1,tree_new))
anova(fit.yx)

model<-pgls(mu~domesticated, data=comp.data)
summary(model)


# Multiple regresion
## phylotools fucked change order
model<-pgls(mu~generation_time+lifespan_wild+offspring_perg+mat, data=comp.data_more)
summary(model)
# Only mammals
model<-pgls(mu~generation_time+lifespan_wild+offspring_perg+mat, data=comp.data_more)
summary(model)

more_m=more[which(data$Group=="Mammal"),]
comp.data_more_m<-comparative.data(tree_mam, more_m, names.col=species_name, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(mu~generation_time+lifespan_wild+offspring_perg+mat, data=comp.data_more_m)
summary(model)

#
more_b=more[which(data$Group=="Bird"),]
comp.data_more_b<-comparative.data(tree_bird, more_b, names.col=species_name, vcv.dim=2, warn.dropped=TRUE, na.omit = FALSE)
model<-pgls(mu~generation_time+lifespan_wild+offspring_perg+mat, data=comp.data_more_b)
summary(model)


#########################################


data = read.csv("../all_sp.txt", na.strings=c("","NA"), sep="\t")


#################################################################
# Correlation
type_col = vector()
type_col[which(data$Group=="Mammal")]="firebrick3"
type_col[which(data$Group=="Bird")]="goldenrod2"
type_col[which(data$Group=="Fish")]="royalblue4"
type_col[which(data$Group=="Reptile")]="limegreen"
type_col[which(data$Group=="Turtle")]="darkmagenta"

########################## Father age ###################################
# rate generation and Replace NA by generation time
age_or_gtime_f <- vector()
age_or_gtime_f[which(is.na(data$father))] <- data$generation_time[which(is.na(data$father))]
age_or_gtime_f[which(!is.na(data$father))] <- data$father[which(!is.na(data$father))]

#
age_or_gtime_m <- vector()
age_or_gtime_m[which(is.na(data$mother))] <- data$generation_time[which(is.na(data$mother))]
age_or_gtime_m[which(!is.na(data$mother))] <- data$mother[which(!is.na(data$mother))]

# Yearly rate for the rest of the analysis + average per species
for_yearly <- (data$father_age+data$mother_age)/2
for_yearly[which(is.na(for_yearly))] <- data$generation_time[which(is.na(for_yearly))]
rate_y <- data$mu/for_yearly
write.table(cbind(as.character(data$sample), rate_y), "y_rate.txt", sep="\t", quote = FALSE, col.names=FALSE, row.names=FALSE)

##for_yearly_b <- data$father_age
##for_yearly_b[which(is.na(for_yearly_b))] <- data$generation_time[which(is.na(for_yearly_b))]
##rate_y <- data$mu/for_yearly_b
##write.table(cbind(as.character(data$sample), rate_y), "y_rate_only_father.txt", sep="\t", quote = FALSE, col.names=FALSE, row.names=FALSE)

# Yearly rate with contribution
data_contrib = read.csv("for_yearly_rate.txt", na.strings=c("","NA"), sep="\t")
for_yearly <- ((data_contrib$father_age*data_contrib$father_contr)+data_contrib$mother_age)/(1+data_contrib$father_contr)
for_yearly[which(is.na(for_yearly))] <- data_contrib$generation_time[which(is.na(for_yearly))]
rate_y <- data_contrib$mu/for_yearly
write.table(cbind(as.character(data_contrib$sample), rate_y), "y_rate_contrib.txt", sep="\t", quote = FALSE, col.names=FALSE, row.names=FALSE)


## Genome size
png("y_rate_genome_size.png", width = 1300, height = 900)
par(mar=c(8,11,4,2), mgp=c(3, 1.5, 0))
plot(data$genome_size, rate_y, type="n", xlab="", ylab="", cex=3, xaxt="n", yaxt="n")
text(data$genome_size, rate_y, labels=data$sample, col=type_col, cex=1.5, font=2)
axis(1, cex.axis=2.5)
axis(2, cex.axis=2.5, las=2)
mtext("Genome size", 1, line=5, cex=3)
mtext("Mutation rate (per year)", 2, line=8, cex=3)
dev.off()

## Population size ( in generation because of lit)
png("rate_population_size.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$population_size, data$mu, type="n", xlab="", ylab="", cex=3, yaxt="n", cex.axis=2.5)
text(data$population_size, data$mu, labels=data$sample, col=type_col, cex=1.5, font=2)
axis(2, cex.axis=2.5, las=2)
mtext("Population size", 1, line=5, cex=3)
mtext("Mutation rate per site per generation", 2, line=9, cex=3)
dev.off()

## Population size Ne ( in generation because of lit)
Ne = read.csv("Ne.txt", sep="\t", header=FALSE)
Ne_list=vector()
for(sp in data$species_name){
  ne = Ne[which(Ne[,1]==sp),2]
  Ne_list=c(Ne_list,ne)
}

png("rate_population_size_Ne.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(Ne_list, data$mu, type="n", xlab="", ylab="", cex=3, yaxt="n", cex.axis=2.5)
points(Ne_list, data$mu, pch=19, col=type_col, cex=1.5, font=2)
##text(data$population_size, data$mu, labels=data$sample, col=type_col, cex=1.5, font=2)
axis(2, cex.axis=2.5, las=2)
mtext("Population size", 1, line=5, cex=3)
mtext("Mutation rate per site per generation", 2, line=9, cex=3)
dev.off()

## Group
png("yrate_group.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$Group, rate_y, type="n", xlab="", ylab="", cex=3, yaxt="n", cex.axis=2.5)
text(data$Group, rate_y, labels=data$sample, col=type_col, cex=1.5, font=2)
axis(2, cex.axis=2.5, las=2)
mtext("Group", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

# Mating #############################################################
# ON all with the various strategy hard to see
##x2  = factor(data$mating, levels=c("monogamous", "polyandrous", "polygynous", "polygynandrous"))
##t.test(rate_y[which(x1=="monogamous")],rate_y[which(x1=="poly")])
##png("y_rate_mating_mono_poly.png", width = 1300, height = 900)
##par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
##plot(x2, rate_y, type="n", ylab="", xlab="", cex=3, yaxt="n")
##text(x2, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
##text(1, 5.3e-8, labels = c("t.test_p.value = 3.647e-05"), font=2, cex=2)
##axis(2, cex.axis=2.5, las=2)
##mtext("Mating", 1, line=5, cex=3)
##mtext("Mutation rate per site per year", 2, line=9, cex=3)
##dev.off()

# But difference between mono and poly
data_mating = read.csv("data_mating.txt", na.strings=c("","NA"), sep="\t")
x1 = factor(data_mating$mating, levels=c("monogamous", "poly"))
t.test(rate_y[which(x1=="monogamous")],rate_y[which(x1=="poly")])
png("y_rate_mating_mono_poly.png", width = 800, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x1, rate_y, type="n", ylab="", xlab="", cex=2, yaxt="n", cex.axis=2.5,pch=20)
##text(x1, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
points(x1, rate_y, col=type_col, cex=2, font=2, pch=20)
text(1, 2e-8, labels = c("t.test_p.value = 4.2e-07"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Mating", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

# for mammals and birds separately
##x2  = factor(data$mating, levels=c("monogamous", "polyandrous", "polygynous", "polygynandrous"))
##x2_m=x2[which(data$Group=="Mammal")]
x2_m = x1[which(data$Group=="Mammal")]
t.test(rate_y[which(data$Group=="Mammal")][which(x2_m=="monogamous")],rate_y[which(data$Group=="Mammal")][which(x2_m=="poly")])
# Test between mono and polygynous
png("y_rate_mating_mammals_mono_poly.png", width = 800, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x2_m, rate_y[which(data$Group=="Mammal")], type="n", ylab="", xlab="", cex=2, yaxt="n", cex.axis=2.5, pch=20)
##text(x2_m, rate_y[which(data$Group=="Mammal")], labels=data$sample[which(data$Group=="Mammal")],col=type_col[which(data$Group=="Mammal")], cex=1.5, font=2)
points(x2_m, rate_y[which(data$Group=="Mammal")], col=type_col[which(data$Group=="Mammal")], cex=2, font=2, pch=20)
text(1, 1.2e-8, labels = c("t.test_p.value = 0.04"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Mating", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

##x2  = factor(data$mating, levels=c("monogamous", "polyandrous", "polygynous", "polygynandrous"))
##x2_b=x2[which(data$Group=="Bird")]
x2_b = x1[which(data$Group=="Bird")]
t.test(rate_y[which(data$Group=="Bird")][which(x2_b=="monogamous")],rate_y[which(data$Group=="Bird")][which(x2_b=="poly")])
# Test between mono and polygynous
png("y_rate_mating_bird_mono_poly.png", width = 800, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x2_b, rate_y[which(data$Group=="Bird")], type="n", ylab="", xlab="", cex=2, yaxt="n", cex.axis=2.5, pch=20)
##text(x2_b, rate_y[which(data$Group=="Bird")], labels=data$sample[which(data$Group=="Bird")],col=type_col[which(data$Group=="Bird")], cex=1.5, font=2)
points(x2_b, rate_y[which(data$Group=="Bird")], col=type_col[which(data$Group=="Bird")], cex=2, font=2, pch=20)
text(1, 9e-9, labels = c("t.test_p.value = 0.001"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Mating", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

# All together
# need some axis:
xaxis_b=vector()
xaxis_b[which(x2_b =="monogamous")] <-5
xaxis_b[which(x2_b =="poly")] <-6
xaxis_m=vector()
xaxis_m[which(x2_m =="monogamous")] <-3
xaxis_m[which(x2_m =="poly")] <-4
##
t.test(rate_y[which(x1=="monogamous")],rate_y[which(x1=="poly")])
t.test(rate_y[which(data$Group=="Mammal")][which(x2_m=="monogamous")],rate_y[which(data$Group=="Mammal")][which(x2_m=="poly")])
t.test(rate_y[which(data$Group=="Bird")][which(x2_b=="monogamous")],rate_y[which(data$Group=="Bird")][which(x2_b=="poly")])
## plot
png("y_rate_mating_combine_mono_poly.png", width = 900, height = 900)
par(mar=c(8,9,4,2), mgp=c(3, 1.5, 0))
boxplot(rate_y[which(x1=="monogamous")], rate_y[which(x1=="poly")],
  rate_y[which(data$Group=="Mammal")][which(x2_m=="monogamous")], rate_y[which(data$Group=="Mammal")][which(x2_m=="poly")],
  rate_y[which(data$Group=="Bird")][which(x2_b=="monogamous")], rate_y[which(data$Group=="Bird")][which(x2_b=="poly")],
  names=rep(c("mono","poly"), 3),
  yaxt="n", cex.axis=2)
points(x1, rate_y, col=type_col, cex=3, font=2, pch=20)
points(xaxis_m, rate_y[which(data$Group=="Mammal")], col=type_col[which(data$Group=="Mammal")], cex=3, font=2, pch=20)
points(xaxis_b, rate_y[which(data$Group=="Bird")], col=type_col[which(data$Group=="Bird")], cex=3, font=2, pch=20)
abline(v=2.5)
abline(v=4.5)
text(c(1.5,3.5,5.5), rep(1.7e-08, 3), labels = c("t-test all species", "t-test mammals", "t-test birds"), font=2, cex=1.5)
text(c(1.5,3.5,5.5), rep(1.6e-08, 3), labels = c("p-value = 4.2e-07", "p-value = 0.041", "p-value = 0.001"), font=2, cex=1.5)
axis(2, at = c(0, 0.5e-8,1e-8,1.5e-8,2e-8), labels = c("0", "0.5","1","1.5","2"), cex.axis=2.5, las=2)
mtext("Mating strategy", 1, line=5, cex=3)
mtext(expression("Mutation rate per site per year (10"^-8*")"), 2, line=5, cex=3)
dev.off()

## with sperm competition
##sperm = factor(data$sperm_comp, levels=c("yes", "no"))
##t.test(rate_y[which(sperm=="yes")][which(data$Group=="Bird")],rate_y[which(sperm=="no")][which(data$Group=="Bird")])
##png("y_rate_mating_mono_poly.png", width = 800, height = 900)
##par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
##plot(sperm[which(data$Group=="Bird")], rate_y[which(data$Group=="Bird")], type="n", ylab="", xlab="", cex=2, yaxt="n", cex.axis=2.5,pch=20)
##points(sperm[which(data$Group=="Bird")], rate_y[which(data$Group=="Bird")], col=type_col[which(data$Group=="Bird")], cex=2, font=2, pch=20)
##text(1, 2e-8, labels = c("t.test_p.value = 4.2e-07"), font=2, cex=2)
##axis(2, cex.axis=2.5, las=2)
##mtext("Mating", 1, line=5, cex=3)
##mtext("Mutation rate per site per year", 2, line=9, cex=3)
##dev.off()



####### K and r #########################################################
t.test(rate_y[which(data$r_K=="r")],rate_y[which(data$r_K=="K")])
png("y_rate_r_K.png", width = 800, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$r_K, rate_y, type="n", ylab="", xlab="", cex=2, pch=20,, yaxt="n", cex.axis=2.5)
##text(data$r_K, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
points(data$r_K, rate_y, col=type_col, cex=2,pch=20, font=2)
text(1, 2e-8, labels = c("t.test_p.value = 8.8e-12"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Strategy", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

# for mammals and birds separately
t.test(rate_y[which(data$Group=="Mammal" & data$r_K=="r")],rate_y[which(data$Group=="Mammal" & data$r_K=="K")])
png("y_rate_r_K_mammals.png", width = 800, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$r_K[which(data$Group=="Mammal")], rate_y[which(data$Group=="Mammal")], type="n", ylab="", xlab="", cex=2, pch=20, yaxt="n", cex.axis=2.5)
##text(data$r_K[which(data$Group=="Mammal")], rate_y[which(data$Group=="Mammal")], labels=data$sample[which(data$Group=="Mammal")],col=type_col[which(data$Group=="Mammal")], cex=1.5, font=2)
points(data$r_K[which(data$Group=="Mammal")], rate_y[which(data$Group=="Mammal")],col=type_col[which(data$Group=="Mammal")], cex=2, pch=20, font=2)
text(1, 1.2e-8, labels = c("t.test_p.value = 4.7e-07"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Strategy", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

x_off = factor(data$offspring, levels=c("low", "med", "lot", "huge"))
x_off_m = x_off[which(data$Group=="Mammal")]
png("y_rate_r_K_off_mammals.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x_off_m, rate_y[which(data$Group=="Mammal")], type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(x_off_m, rate_y[which(data$Group=="Mammal")], labels=data$sample[which(data$Group=="Mammal")],col=type_col[which(data$Group=="Mammal")], cex=1.5, font=2)
##text(1, 1.2e-8, labels = c("t.test_p.value = 3.779e-07"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Offspring", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

t.test(rate_y[which(data$Group=="Bird" & data$r_K=="r")],rate_y[which(data$Group=="Bird" & data$r_K=="K")])
png("y_rate_r_K_bird.png", width = 800, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$r_K[which(data$Group=="Bird")], rate_y[which(data$Group=="Bird")], type="n", ylab="", xlab="", cex=2, pch=20, yaxt="n", cex.axis=2.5)
##text(data$r_K[which(data$Group=="Bird")], rate_y[which(data$Group=="Bird")], labels=data$sample[which(data$Group=="Bird")],col=type_col[which(data$Group=="Bird")], cex=1.5, font=2)
points(data$r_K[which(data$Group=="Bird")], rate_y[which(data$Group=="Bird")], col=type_col[which(data$Group=="Bird")], cex=2, pch=20,font=2)
text(1, 9e-9, labels = c("t.test_p.value = 8.3e-05"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Strategy", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

x_off_b = x_off[which(data$Group=="Bird")]
png("y_rate_r_K_off_bird.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x_off_b, rate_y[which(data$Group=="Bird")], type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(x_off_b, rate_y[which(data$Group=="Bird")], labels=data$sample[which(data$Group=="Bird")],col=type_col[which(data$Group=="Bird")], cex=1.5, font=2)
##text(1, 1.5e-8, labels = c("t.test_p.value = 0.000553"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Offspring", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()


# All together
# need some axis:
only_b = data$r_K[which(data$Group=="Bird")]
only_m = data$r_K[which(data$Group=="Mammal")]
xaxis_b=vector()
xaxis_b[which(only_b =="K")] <-5
xaxis_b[which(only_b =="r")] <-6
xaxis_m=vector()
xaxis_m[which(only_m =="K")] <-3
xaxis_m[which(only_m =="r")] <-4
##
t.test(rate_y[which(data$r_K=="r")],rate_y[which(data$r_K=="K")])
t.test(rate_y[which(data$Group=="Mammal" & data$r_K=="r")],rate_y[which(data$Group=="Mammal" & data$r_K=="K")])
t.test(rate_y[which(data$Group=="Bird" & data$r_K=="r")],rate_y[which(data$Group=="Bird" & data$r_K=="K")])
## plot
png("y_rate_mating_combine_r_K.png", width = 900, height = 900)
par(mar=c(8,9,4,2), mgp=c(3, 1.5, 0))
boxplot(rate_y[which(data$r_K=="K")], rate_y[which(data$r_K=="r")],
    rate_y[which(data$Group=="Mammal" & data$r_K=="K")],rate_y[which(data$Group=="Mammal" & data$r_K=="r")],
    rate_y[which(data$Group=="Bird" & data$r_K=="K")],rate_y[which(data$Group=="Bird" & data$r_K=="r")],
    names=rep(c("K","r"), 3),
    yaxt="n", cex.axis=2)
points(data$r_K, rate_y, col=type_col, cex=3, font=2, pch=20)
points(xaxis_m, rate_y[which(data$Group=="Mammal")], col=type_col[which(data$Group=="Mammal")], cex=3, font=2, pch=20)
points(xaxis_b, rate_y[which(data$Group=="Bird")], col=type_col[which(data$Group=="Bird")], cex=3, font=2, pch=20)
abline(v=2.5)
abline(v=4.5)
text(c(1.5,3.5,5.5), rep(1.7e-08, 3), labels = c("t-test all species", "t-test mammals", "t-test birds"), font=2, cex=1.5)
text(c(1.5,3.5,5.5), rep(1.6e-08, 3), labels = c("p-value = 8.8e-12", "p-value = 4.7e-07", "p-value = 8.3e-05"), font=2, cex=1.5)
axis(2, at = c(0, 0.5e-8,1e-8,1.5e-8,2e-8), labels = c("0", "0.5","1","1.5","2"), cex.axis=2.5, las=2)
mtext("Strategy", 1, line=5, cex=3)
mtext(expression("Mutation rate per site per year (10"^-8*")"), 2, line=5, cex=3)
dev.off()




## Reproduction
##t.test(rate_y[which(data$r_K=="r")],rate_y[which(data$r_K=="K")])
png("y_rate_repro.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$reproduction, rate_y, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(data$reproduction, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
##text(1, 5.3e-8, labels = c("t.test_p.value = 1.23e-06"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Reproduction", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

# Latitude
x3 = factor(data$latitude, levels=c("high", "low"))
t.test(rate_y[which(data$latitude=="high")],rate_y[which(data$latitude=="low")])
png("y_rate_lat.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x3, rate_y, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(x3, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
text(1, 2e-8, labels = c("t.test_p.value = 2.2e-07"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Latitude", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

png("y_rate_lat_all.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$latitude, rate_y, type="n", ylab="", xlab="", cex=3, yaxt="n")
text(data$latitude, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
##text(1, 5.3e-8, labels = c("t.test_p.value = 1.13e-05"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Latitude", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

## Habitat
data_land = read.csv("data_land.txt", na.strings=c("","NA"), sep="\t")
t.test(rate_y[which(data_land$habitat=="aquatic")],rate_y[which(data_land$habitat=="land")])
png("y_rate_habitat.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data_land$habitat, rate_y, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(data_land$habitat, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
text(1, 2e-8, labels = c("t.test_p.value = 0.19"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Habitat", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

## Mass
##t.test(rate_y[which(data_land$habitat=="aquatic")],rate_y[which(data_land$habitat=="land")])
png("y_rate_mass.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$mass, rate_y, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(data$mass, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
##text(1, 5.3e-8, labels = c("t.test_p.value = 0.1066"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Mass", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

png("y_rate_mass_log.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$mass, rate_y, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5, log="x")
text(data$mass, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
##text(1, 5.3e-8, labels = c("t.test_p.value = 0.1066"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Mass (log)", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

## Longevity
png("y_rate_longevity.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$lifespan_wild, rate_y, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(data$lifespan_wild, rate_y, labels=data$sample,col=type_col, cex=1.5, font=2)
##text(1, 5.3e-8, labels = c("t.test_p.value = 0.1066"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("lifespan (wild)", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

png("g_rate_longevity.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data$lifespan_wild, data$mu, type="n", xlab="", ylab="", cex=3, yaxt="n", cex.axis=2.5)
text(data$lifespan_wild, data$mu, labels=data$sample, col=type_col, cex=1.5, font=2)
axis(2, cex.axis=2.5, las=2)
mtext("Longevity", 1, line=5, cex=3)
mtext("Mutation rate per site per generation", 2, line=9, cex=3)
dev.off()






#########################################
# Yearly rate only for the one with parental age!
for_yearly <- (data$father_age+data$mother_age)/2
rate_y_sub_all <- data$mu/for_yearly
rate_y_sub = rate_y_sub_all[-which(is.na(for_yearly))]
data_sub = data[-which(is.na(for_yearly)),]
write.table(cbind(as.character(data_sub$sample), rate_y_sub), "y_rate_sub.txt", sep="\t", quote = FALSE, col.names=FALSE, row.names=FALSE)

## Group
png("yrate_sub_group.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data_sub$Group, rate_y_sub, type="n", xlab="", ylab="", cex=3, yaxt="n", cex.axis=2.5)
text(data_sub$Group, rate_y_sub, labels=data_sub$sample, col=type_col[-which(is.na(for_yearly))], cex=1.5, font=2)
axis(2, cex.axis=2.5, las=2)
mtext("Group", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

# But difference between mono and poly
data_mating = read.csv("data_mating.txt", na.strings=c("","NA"), sep="\t")
data_mating_sub = data_mating[-which(is.na(for_yearly)),]
x1 = factor(data_mating_sub$mating, levels=c("monogamous", "poly"))
t.test(rate_y_sub[which(x1=="monogamous")],rate_y_sub[which(x1=="poly")])
png("y_rate_sub_mating_mono_poly.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x1, rate_y_sub, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(x1, rate_y_sub, labels=data_sub$sample,col=type_col[-which(is.na(for_yearly))], cex=1.5, font=2)
text(1, 1.2e-8, labels = c("t.test_p.value = 0.0003"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Mating", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

####### K and r #########################################################
t.test(rate_y_sub[which(data$r_K=="r")],rate_y_sub[which(data$r_K=="K")])
png("y_rate_sub_r_K.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data_sub$r_K, rate_y_sub, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(data_sub$r_K, rate_y_sub, labels=data_sub$sample,col=type_col[-which(is.na(for_yearly))], cex=1.5, font=2)
text(1, 1.2e-8, labels = c("t.test_p.value = 0.01"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Offspring", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

# Latitude
x3 = factor(data_sub$latitude, levels=c("high", "low"))
t.test(rate_y_sub[which(data_sub$latitude=="high")],rate_y_sub[which(data_sub$latitude=="low")])
png("y_rate_sub_lat.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(x3, rate_y_sub, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(x3, rate_y_sub, labels=data_sub$sample,col=type_col[-which(is.na(for_yearly))], cex=1.5, font=2)
text(1, 1.2e-8, labels = c("t.test_p.value = 0.0003"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Latitude", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

## Habitat
data_land = read.csv("data_land.txt", na.strings=c("","NA"), sep="\t")
data_land_sub = data_land[-which(is.na(for_yearly)),]
t.test(rate_y_sub[which(data_land_sub$habitat=="aquatic")],rate_y_sub[which(data_land_sub$habitat=="land")])
png("y_rate_sub_habitat.png", width = 1300, height = 900)
par(mar=c(8,12,4,2), mgp=c(3, 1.5, 0))
plot(data_land_sub$habitat, rate_y_sub, type="n", ylab="", xlab="", cex=3, yaxt="n", cex.axis=2.5)
text(data_land_sub$habitat, rate_y_sub, labels=data_sub$sample,col=type_col[-which(is.na(for_yearly))], cex=1.5, font=2)
text(1, 1.2e-8, labels = c("t.test_p.value = 0.24"), font=2, cex=2)
axis(2, cex.axis=2.5, las=2)
mtext("Habitat", 1, line=5, cex=3)
mtext("Mutation rate per site per year", 2, line=9, cex=3)
dev.off()

#########################################

data = read.csv("../all_sp.txt", na.strings=c("","NA"), sep="\t") # From supplementary table

#################################################################
type_col = vector()
type_col[which(data$Group=="Mammal")]="firebrick3"
type_col[which(data$Group=="Bird")]="goldenrod2"
type_col[which(data$Group=="Fish")]="royalblue4"
type_col[which(data$Group=="Reptile")]="limegreen"
type_col[which(data$Group=="Turtle")]="darkmagenta"

## Only when parental age available #############################
data_sub = data[-which(is.na(data$father)),]
type_col_sub = type_col[-which(is.na(data$father))]

# correlation ages
cor.test(data_sub$father_age, data_sub$mother_age)
fit = lm(data_sub$father_age~ data_sub$mother_age)
summary(fit)
png("correlation_age.png", width = 1300, height = 900)
par(mar=c(8,9,4,2), mgp=c(3, 1.5, 0))
plot(data_sub$mother_age, data_sub$father_age, type="n", xlab="", ylab="", cex=3, xaxt="n", yaxt="n")
points(data_sub$mother_age, data_sub$father_age, col=type_col_sub, cex=3, pch=19)
text(c(32,32), c(25, 21), labels = c(expression('adjusted R'^2*' = 0.77'), "p < 2.2e-16"), font=2, cex=3)
axis(1, cex.axis=2.5)
axis(2, cex.axis=2.5, las=2)
abline(0.9563,1.0063, lwd=4)
mtext("Maternal age", 1, line=5, cex=3)
mtext("Paternal age", 2, line=6, cex=3)
dev.off()

### Correlation with average age:
age_ave = (data_sub$father_age+data_sub$mother_age)/2

## cor.test(age_ave, data_sub$mu)
## fit = lm(data_sub$mu~age_ave)
## summary(fit)
## cor.test(age_ave[which(data_sub$Group=="Mammal")], data_sub$mu[which(data_sub$Group=="Mammal")])
## fit = lm(data_sub$mu[which(data_sub$Group=="Mammal")]~age_ave[which(data_sub$Group=="Mammal")])
## summary(fit)
##cor.test(age_ave[which(data_sub$Group=="Bird")][-23], data_sub$mu[which(data_sub$Group=="Bird")][-23])
## fit = lm(data_sub$mu[which(data_sub$Group=="Bird")][-23]~age_ave[which(data_sub$Group=="Bird")][-23])
## summary(fit)

png("g_rate_parental_age_all.png", width = 1300, height = 900)
par(mar=c(8,11,4,2), mgp=c(3, 1.5, 0))
plot(age_ave, data_sub$mu, type="n", xlab="", ylab="", cex=2, xaxt="n", yaxt="n", xlim=c(0,37), ylim=c(0,4e-8))
points(age_ave, data_sub$mu, col=type_col_sub, cex=4, pch=19)
#
text(c(33,33,33), c(2.2e-8, 2e-8, 1.8e-8), labels = c("105 samples", expression('adjusted R'^2*' = 0.14'), "p = 3.9e-5"), font=2, cex=2.5)
abline(6.327e-09,2.253e-10, lwd=4)
text(c(21,21, 21), c(2e-8, 1.8e-8, 1.6e-8), labels = c("Mammals", expression('adjusted R'^2*' = 0.37'), "p = 1.6e-7"), font=2, cex=2.5, col="firebrick3")
abline(6.117e-09,2.374e-10, col="firebrick3", lwd=4)
text(c(30,30,30), c(0.8e-8, 0.6e-8, 0.4e-8), labels = c("Birds without Rhea", expression('adjusted R'^2*' = 0.31'), "p = 0.0005"), font=2, cex=2.5, col="goldenrod2")
abline(5.389e-09,2.494e-10, col="goldenrod2", lwd=4)
axis(1, cex.axis=2.5)
axis(2, at = c(0,1e-8, 2e-8, 3e-8, 4e-8), labels = c(0,1,2,3,4), cex.axis=2.5, las=2)
mtext("Average parental age at the time of reproduction", 1, line=5, cex=3)
mtext(expression('Mutation rate per site per generation (10'^-8*')'), 2, line=4, cex=3)
dev.off()

fit = lm(data_sub$mu~data_sub$father_age+data_sub$mother_age)
summary(fit) # show results

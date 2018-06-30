agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")

setwd("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01")
##########################################################################
# # # # # # # # # #     Agave s.l.STEM    # # # # # # # # # # # #
# BlinDcons
load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/CIs/CI_cladeA1_BlinDcons_2016.10.28.RData")
nrow(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound) # 213 967

colFunc <- colorRampPalette(c("red","blue"))
lhCol <- colFunc(213967)
lhCol3 <- paste(lhCol, 30, sep="")
PARbounds <- CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound
x <- sample(1:213967, 10000)
pdf("CI_cladeA1_BlinDcons_2016.10.28_sppVSext_2018.05.08.pdf")
plot(x=NULL, main=c("Agave clade from stem best model (BlinDcons)", "Likelihood surface"), xlab="speciation rate", ylab="extinction rate", xlim=c(1,3.5), ylim=c(0,1.5))
for (i in x){
  points(PARbounds[i,"lamb0"],PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"lamb0"], CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"mu0"], col="black", pch=8, cex=1)
dev.off()

pdf("CI_cladeA1_BlinDcons_2016.10.28_divVSrelExt_2018.05.08.pdf")
plot(x=NULL, main=c("Agave clade from stem best model (BlinDcons)", "Likelihood surface"), xlab="diversification rate", ylab="relative extinction rate", xlim=c(-0.7,3), ylim=c(0,2))
abline(h=1, v=0, col="gray")
for (i in x){
  points(PARbounds[i,"r0"],PARbounds[i,"mu0"]/PARbounds[i,"lamb0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"r0"], CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"mu0"]/CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"lamb0"], col="black", pch=8, cex=1)
dev.off()

##########################################################################
# # # # # # # # # #     Agave s.l.STEM    # # # # # # # # # # # #
# AT STEM DIVERGENCE ONLY
# BlinDcons
agave_stem <- branching.times(agave)['57'] # 8.06559
pdf("CI_cladeA1_BlinDcons_2016.10.28_sppVSext_2018.05.08_at_stem.pdf")
plot(x=NULL, main=c("Agave clade best model (BlinDcons)", "Likelihood surface"), xlab="speciation rate", ylab="extinction rate", xlim=c(0,5), ylim=c(0,10), cex.axis = 2)
for (i in x){
  BlinAgave <- function(t){PARbounds[i, "lamb0"] + (abs(PARbounds[i, "alpha"]) * t)} 
  DconsAgave <- abs(PARbounds[i, "mu0"])
  # rAgave <- function(t){BlinAgave(t)-DconsAgave}
  points(BlinAgave(agave_stem),PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex=1.5)
}
BlinAgave <- function(t){CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"lamb0"] + (CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"alpha"] * t)} 
DconsAgave <- abs(CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"mu0"])
# rAgave <- function(t){BlinAgave(t)-DconsAgave}
points(BlinAgave(agave_stem), DconsAgave, col="black", pch=8, cex=2, lwd=3)
dev.off()

pdf("CI_cladeA1_BlinDcons_2016.10.28_divVSrelExt_2018.05.08_at_stem.pdf")
plot(x=NULL, main=c("woody Fouquieria clade best model (BexpDcons)", "Likelihood surface"), xlab="diversification rate at crown", ylab="relative extinction rate", xlim=c(-1,3), ylim=c(0,2), cex.axis = 2)
abline(h=1, v=0, col="gray")
for (i in x){
  BlinAgave <- function(t){PARbounds[i, "lamb0"] + (PARbounds[i, "alpha"] * t)} 
  DconsAgave <- abs(PARbounds[i, "mu0"])
  rAgave <- function(t){BlinAgave(t)-DconsAgave}
  points(rAgave(agave_stem),PARbounds[i,"mu0"]/BlinAgave(agave_stem), col=lhCol3[i], pch=16, cex=1.5)
}
BlinAgave <- function(t){CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"lamb0"] + (CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"alpha"] * t)} 
DconsAgave <- abs(CI_cladeA1_BlinDcons_2016.10.28$CIs[1,"mu0"])
rAgave <- function(t){BlinAgave(t)-DconsAgave}
points(rAgave(agave_stem), DconsAgave/BlinAgave(agave_stem), col="black", pch=8, cex=2, lwd=3)
dev.off()

##########################################################################
# # # # # # # # # #     Furcraea-Beschorneria STEM     # # # # # # # # # #
load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/CIs/CI_cladeA2_BconsDcons_2016.10.26.RData")

nrow(CI_cladeA2_BconsDcons$lik_vals_bound)  # 103 137
colFunc <- colorRampPalette(c("red","blue"))
lhCol <- colFunc(103137)
lhCol3 <- paste(lhCol, 30, sep="")
PARbounds <- CI_cladeA2_BconsDcons$lik_vals_bound
head(PARbounds)
x <- sample(1:103137, 10000)

pdf("CI_cladeA2_BconsDcons_sppVSext_2018.05.14.pdf")
plot(x=NULL, main=c("Furcraea-Beschorneria from stem BconsDcons model", 
                    "Likelihood surface"), xlab="speciation rate", ylab="extinction rate", 
     xlim=c(0,1), ylim=c(0,1))
for (i in x){
  points(PARbounds[i,"lamb0"], PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_cladeA2_BconsDcons$CIs[1,"lamb0"], abs(CI_cladeA2_BconsDcons$CIs[1,"mu0"]), col="black", pch=8, cex=1)
dev.off()

# more points: CI_cladeA2_BconsDcons_2
load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/CI_cladeA2_BconsDcons_2_2018.05.14.RData")
nrow(CI_cladeA2_BconsDcons_2$lik_vals_bound)  # 321364
colFunc <- colorRampPalette(c("red","blue"))
lhCol <- colFunc(321364)
lhCol3 <- paste(lhCol, 30, sep="")
PARbounds <- CI_cladeA2_BconsDcons_2$lik_vals_bound
head(PARbounds)
x <- sample(1:321364, 15000)

pdf("CI_cladeA2_BconsDcons_2_sppVSext_2018.05.14.pdf")
plot(x=NULL, main=c("Furcraea-Beschorneria from stem BconsDcons model", 
                    "Likelihood surface"), xlab="speciation rate", ylab="extinction rate", 
     xlim=c(0,1.5), ylim=c(0,1), labels = NA)
axis(side = 1, lwd = 0, line = 1, cex.axis = 3, at = c(0,0.5, 1,1.5), labels = as.character(c(0,0.5, 1,1.5)))
axis(side = 2, lwd = 0, line = -.2, cex.axis = 3, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = as.character(c(0, 0.2, 0.4, 0.6, 0.8, 1)))
for (i in x){
  points(PARbounds[i,"lamb0"], PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex= 1.5)
}
points(CI_cladeA2_BconsDcons_2$CIs[1,"lamb0"], abs(CI_cladeA2_BconsDcons_2$CIs[1,"mu0"]), col="black", pch=8, cex=2, lwd =3)
dev.off()

pdf("CI_cladeA2_BconsDcons_2_divVSrelExt_2018.05.14.pdf")
plot(x=NULL, main=c("Furcraea-Beschorneria from stem BconsDcons model", 
                    "Likelihood surface"), xlab="diversification rate", ylab="relative extinction rate", 
     xlim=c(-1,1), ylim=c(0,2), labels = NA)
axis(side = 1, lwd = 0, line = 1, cex.axis = 3, at = c(-1, -0.5, 0, 0.5, 1), labels = as.character(c(-1, -0.5, 0, 0.5, 1)))
axis(side = 2, lwd = 0, line = -.2, cex.axis = 3, at = c(0,0.5, 1,1.5,2), labels = as.character(c(0,0.5, 1,1.5,2)))
abline(h=1, v=0, col="gray")
for (i in x){
  points(PARbounds[i,"r0"], PARbounds[i,"mu0"]/PARbounds[i,"lamb0"], col=lhCol3[i], pch=16, cex= 1)
}
points(CI_cladeA2_BconsDcons_2$CIs[1,"r0"], CI_cladeA2_BconsDcons_2$CIs[1,"mu0"]/CI_cladeA2_BconsDcons_2$CIs[1,"lamb0"], col="black", pch=8, cex=2, lwd =3)
dev.off()

##########################################################################
# # # # # #     Backbone (Yucca + Hesperoyucca-Hesperaloe)     # # # # # #
load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/CI_BBcladeA1_A2_BconsDcons_1_2018.05.14.RData")
nrow(CI_BBcladeA1_A2_BconsDcons_1$lik_vals_bound)  #539
colFunc <- colorRampPalette(c("red","blue"))
lhCol <- colFunc(539)
lhCol3 <- paste(lhCol, 30, sep="")
PARbounds <- CI_BBcladeA1_A2_BconsDcons_1$lik_vals_bound
head(PARbounds)
x <- 1:539
pdf("CI_BBcladeA1_A2_BconsDcons_1_sppVSext_2018.05.14.pdf")
plot(x=NULL, main=c("Backbone (Yucca + Hesperoyucca-Hesperaloe) BconsDcons model", 
                    "Likelihood surface"), xlab="speciation rate", ylab="extinction rate", 
     xlim=c(0,1), ylim=c(0,1))
for (i in x){
  points(PARbounds[i,"lamb0"],PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_BBcladeA1_A2_BconsDcons_1$CIs[1,"lamb0"], abs(CI_BBcladeA1_A2_BconsDcons_1$CIs[1,"mu0"]), col="black", pch=8, cex=1)
dev.off()

pdf("CI_BBcladeA1_A2_BconsDcons_1_divVSrelExt_2018.05.14.pdf")
plot(x=NULL, main=c("Backbone (Yucca + Hesperoyucca-Hesperaloe) BconsDcons model", 
                    "Likelihood surface"), xlab="diversification rate", ylab="relative extinction rate", 
     xlim=c(-1,1), ylim=c(0,2))
abline(h=1, v=0, col="gray")
for (i in x){
  points(PARbounds[i,"r0"], PARbounds[i,"mu0"]/PARbounds[i,"lamb0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_BBcladeA1_A2_BconsDcons_1$CIs[1,"r0"], abs(CI_BBcladeA1_A2_BconsDcons_1$CIs[1,"mu0"])/CI_BBcladeA1_A2_BconsDcons_1$CIs[1,"lamb0"], col="black", pch=8, cex=1)
dev.off()


# WITH MORE POINTS:

load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/CI_BBcladeA1_A2_BconsDcons_2_2018.05.14.RData")
nrow(CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound)  #55296
colFunc <- colorRampPalette(c("red","blue"))
lhCol <- colFunc(55296)
lhCol3 <- paste(lhCol, 30, sep="")
PARbounds <- CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound
head(PARbounds)
x <- sample(1:nrow(CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound), 10000)

pdf("CI_BBcladeA1_A2_BconsDcons_2_sppVSext_2018.05.14.pdf")
plot(x=NULL, main=c("Backbone (Yucca + Hesperoyucca-Hesperaloe) BconsDcons model", 
                    "Likelihood surface"), xlab="speciation rate", ylab="extinction rate", 
     xlim=c(0,1), ylim=c(0,1), labels = NA)
axis(side = 1, lwd = 0, line = 1, cex.axis = 3, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = as.character(c(0, 0.2, 0.4, 0.6, 0.8, 1)))
axis(side = 2, lwd = 0, line = -.2, cex.axis = 3, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = as.character(c(0, 0.2, 0.4, 0.6, 0.8, 1)))
for (i in x){
  points(PARbounds[i,"lamb0"],PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_BBcladeA1_A2_BconsDcons_2$CIs[1,"lamb0"], abs(CI_BBcladeA1_A2_BconsDcons_2$CIs[1,"mu0"]), col="black", pch=8, cex=1)
dev.off()

pdf("CI_BBcladeA1_A2_BconsDcons_2_divVSrelExt_2018.05.14.pdf")
plot(x=NULL, main=c("Backbone (Yucca + Hesperoyucca-Hesperaloe) BconsDcons model", 
                    "Likelihood surface"), xlab="diversification rate", ylab="relative extinction rate", 
     xlim=c(-1,1), ylim=c(0,2), labels = NA)
axis(side = 1, lwd = 0, line = 1, cex.axis = 3, at = c(-1, -0.5, 0, 0.5, 1), labels = as.character(c(-1, -0.5, 0, 0.5, 1)))
axis(side = 2, lwd = 0, line = -.2, cex.axis = 3, at = c(0,0.5, 1,1.5,2), labels = as.character(c(0,0.5, 1,1.5,2)))
abline(h=1, v=0, col="gray")
for (i in x){
  points(PARbounds[i,"r0"], PARbounds[i,"mu0"]/PARbounds[i,"lamb0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_BBcladeA1_A2_BconsDcons_2$CIs[1,"r0"], abs(CI_BBcladeA1_A2_BconsDcons_2$CIs[1,"mu0"])/CI_BBcladeA1_A2_BconsDcons_2$CIs[1,"lamb0"], col="black", pch=8, cex=1)
dev.off()

##########################################################################
# # # # # # # # # #     Agave s.l.STEM    # # # # # # # # # # # #
# at crown
# BexpDcons
load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/CI_cladeA1_BexpDcons_2_2018.05.11.RData")
nrow(CI_cladeA1_BexpDcons_2$lik_vals_bound) # 90837

colFunc <- colorRampPalette(c("red","blue"))
lhCol <- colFunc(90837)
lhCol3 <- paste(lhCol, 30, sep="")
PARbounds <- CI_cladeA1_BexpDcons_2$lik_vals_bound
x <- sample(1:90837, 50000)
pdf("CI_cladeA1_BexpDcons_2_sppVSext_2018.05.14.pdf")
plot(x=NULL, main=c("Agave clade from stem BexpDcons model", "Likelihood surface"), 
     xlab="speciation rate", ylab="extinction rate", xlim=c(1,3), ylim=c(0,1))
for (i in x){
  points(PARbounds[i,"lamb0"],PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_cladeA1_BexpDcons_2$CIs[1,"lamb0"], CI_cladeA1_BexpDcons_2$CIs[1,"mu0"], col="black", pch=8, cex=1)
dev.off()

pdf("CI_cladeA1_BexpDcons_2_divVSrelExt_2018.05.14.pdf")
plot(x=NULL, main=c("Agave clade from stem BexpDcons model", "Likelihood surface"), 
     xlab="diversification rate", ylab="relative extinction rate", xlim=c(0,3), ylim=c(0,2))
abline(h=1, v=0, col="gray")
for (i in x){
  points(PARbounds[i,"r0"],PARbounds[i,"mu0"]/PARbounds[i,"lamb0"], col=lhCol3[i], pch=16, cex=0.7)
}
points(CI_cladeA1_BexpDcons_2$CIs[1,"r0"], CI_cladeA1_BexpDcons_2$CIs[1,"mu0"]/CI_cladeA1_BexpDcons_2$CIs[1,"lamb0"], col="black", pch=8, cex=1)
dev.off()

##########################################################################
# # # # # # # # # #     Agave s.l.STEM    # # # # # # # # # # # #
# AT STEM DIVERGENCE ONLY
# BexpDcons
colFunc <- colorRampPalette(c("red","blue"))
lhCol <- colFunc(90837)
lhCol3 <- paste(lhCol, 30, sep="")
PARbounds <- CI_cladeA1_BexpDcons_2$lik_vals_bound
x <- sample(1:90837, 10000)

agave_stem <- branching.times(agave)['57'] # 8.06559

pdf("CI_cladeA1_BexpDcons_2_sppVSext_2018.05.14_at_stem.pdf")
plot(x=NULL, main=c("Agave clade best BexpDcons model", "Likelihood surface"), 
     xlab="speciation rate", ylab="extinction rate", xlim=c(0,1), ylim=c(0,1), cex.axis = 2)
for (i in x){
  BlinAgave <- function(t){PARbounds[i, "lamb0"] * exp(PARbounds[i, "alpha"] * t)} 
  DconsAgave <- abs(PARbounds[i, "mu0"])
  # rAgave <- function(t){BlinAgave(t)-DconsAgave}
  points(BlinAgave(agave_stem),PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex=1.5)
}
BexpAgave <- function(t){CI_cladeA1_BexpDcons_2$CIs[1,"lamb0"] * exp(CI_cladeA1_BexpDcons_2$CIs[1,"alpha"] * t)} 
DconsAgave <- abs(CI_cladeA1_BexpDcons_2$CIs[1,"mu0"])
# rAgave <- function(t){BexpAgave(t)-DconsAgave}
points(BexpAgave(agave_stem), DconsAgave, col="black", pch=8, cex=2, lwd=3)
dev.off()

pdf("CI_cladeA1_BexpDcons_2_divVSrelExt_2018.05.14_at_stem.pdf")
plot(x=NULL, main=c("woody Fouquieria clade BexpDcons model", "Likelihood surface"), xlab="diversification rate at crown", ylab="relative extinction rate", xlim=c(-1,3), ylim=c(0,2), cex.axis = 2)
abline(h=1, v=0, col="gray")
for (i in x){
  BexpAgave <- function(t){PARbounds[i, "lamb0"] * exp(PARbounds[i, "alpha"] * t)} 
  DconsAgave <- abs(PARbounds[i, "mu0"])
  rAgave <- function(t){BexpAgave(t)-DconsAgave}
  points(rAgave(agave_stem), PARbounds[i,"mu0"]/BexpAgave(agave_stem), col=lhCol3[i], pch=16, cex=1.5)
}
BexpAgave <- function(t){CI_cladeA1_BexpDcons_2$CIs[1,"lamb0"] * exp(CI_cladeA1_BexpDcons_2$CIs[1,"alpha"] * t)} 
DconsAgave <- abs(CI_cladeA1_BexpDcons_2$CIs[1,"mu0"])
rAgave <- function(t){BexpAgave(t)-DconsAgave}
points(rAgave(agave_stem), DconsAgave/BexpAgave(agave_stem), col="black", pch=8, cex=2, lwd=3)
dev.off()

##########################################################################
# # # # # # # # # #     Agave s.l.STEM    # # # # # # # # # # # #
# AT PRESENT AND AT STEM DIVERGENCE
# BexpDcons
load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/CI_cladeA1_BexpDcons_2_2018.05.11.RData")
nrow(CI_cladeA1_BexpDcons_2$lik_vals_bound) # 90837
PARbounds <- CI_cladeA1_BexpDcons_2$lik_vals_bound

colFunc <- colorRampPalette(c("red","blue"))
lhCol <- colFunc(90837)
lhCol3 <- paste(lhCol, 30, sep="")

pdf("CI_cladeA1_BexpDcons_2_sppVSext_2018.05.14_present_and_stem.pdf")
plot(x=NULL, main=c("Agave clade from stem BexpDcons model", "Likelihood surface"), 
     xlab="speciation rate", ylab="extinction rate", xlim=c(0,3), ylim=c(0,1), labels = FALSE)
axis(side = 1, lwd = 0, line = 1, cex.axis = 3, at = c(0,1,2,3), labels = as.character(c(0,1,2,3)))
axis(side = 2, lwd = 0, line = -.2, cex.axis = 3, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = as.character(c(0, 0.2, 0.4, 0.6, 0.8, 1)))
#PRESENT:
x <- sample(1:90837, 90000)
for (i in x){
  points(PARbounds[i,"lamb0"],PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex= 1)
}
points(CI_cladeA1_BexpDcons_2$CIs[1,"lamb0"], CI_cladeA1_BexpDcons_2$CIs[1,"mu0"], 
       col="black", pch=8, cex=3, lwd =3)
#STEM:
x <- sample(1:90837, 50000)
for (i in x){
  BlinAgave <- function(t){PARbounds[i, "lamb0"] * exp(PARbounds[i, "alpha"] * t)} 
  DconsAgave <- abs(PARbounds[i, "mu0"])
  # rAgave <- function(t){BlinAgave(t)-DconsAgave}
  points(BlinAgave(agave_stem),PARbounds[i,"mu0"], col=lhCol3[i], pch=16, cex=1)
}
BexpAgave <- function(t){CI_cladeA1_BexpDcons_2$CIs[1,"lamb0"] * exp(CI_cladeA1_BexpDcons_2$CIs[1,"alpha"] * t)} 
DconsAgave <- abs(CI_cladeA1_BexpDcons_2$CIs[1,"mu0"])
# rAgave <- function(t){BexpAgave(t)-DconsAgave}
points(BexpAgave(agave_stem), DconsAgave, col="black", pch=8, cex=3, lwd=3)
dev.off()


pdf("CI_cladeA1_BexpDcons_2_divVSrelExt_2018.05.14_present_and_stem.pdf")
plot(x=NULL, main=c("Agave clade from stem BexpDcons model", "Likelihood surface"), 
     xlab="diversification rate", ylab="relative extinction rate", xlim=c(-1,3), ylim=c(0,3.5), labels = NA)
axis(side = 1, lwd = 0, line = 1, cex.axis = 3)
axis(side = 2, lwd = 0, line = -.2, cex.axis = 3, at = c(0,1,2,3), labels = as.character(c(0,1,2,3)))
abline(h=1, v=0, col="gray")
#PRESENT:
x <- sample(1:90837, 20000)
for (i in x){
  points(PARbounds[i,"r0"],PARbounds[i,"mu0"]/PARbounds[i,"lamb0"], col=lhCol3[i], pch=16, cex= 1)
}
points(CI_cladeA1_BexpDcons_2$CIs[1,"r0"], 
       CI_cladeA1_BexpDcons_2$CIs[1,"mu0"]/CI_cladeA1_BexpDcons_2$CIs[1,"lamb0"], 
       col="black", pch=8, cex=3, lwd =3)
#STEM:
x <- sample(1:90837, 20000)
for (i in x){
  BexpAgave <- function(t){PARbounds[i, "lamb0"] * exp(PARbounds[i, "alpha"] * t)} 
  DconsAgave <- abs(PARbounds[i, "mu0"])
  rAgave <- function(t){BexpAgave(t)-DconsAgave}
  points(rAgave(agave_stem), PARbounds[i,"mu0"]/BexpAgave(agave_stem), col=lhCol3[i], pch=16, cex= 1)
}
BexpAgave <- function(t){CI_cladeA1_BexpDcons_2$CIs[1,"lamb0"] * exp(CI_cladeA1_BexpDcons_2$CIs[1,"alpha"] * t)} 
DconsAgave <- abs(CI_cladeA1_BexpDcons_2$CIs[1,"mu0"])
rAgave <- function(t){BexpAgave(t)-DconsAgave}
points(rAgave(agave_stem), DconsAgave/BexpAgave(agave_stem), col="black", pch=8, cex=3, lwd=3)
dev.off()

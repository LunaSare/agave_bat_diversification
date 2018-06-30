# 7 dic 2016
# Gráficos mejor modelo de diversificación de Agavaceae: 3 dinámicas de diversificación

lapply(as.list(c("ape", "picante", "geiger","phyloch")), require, character.only=TRUE)
agave <- read.beast("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut_nombresPub")
agave.table <- read.beast.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut_nombresPub")
head(agave.table)

write.table(agave.table, file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/agave.nodes.xls", sep="\t", quote=FALSE, row.names=FALSE)


data(strat2012)
plot(agave)
edgelabels()

plot(ladderize(agave, right=FALSE), cex=0.5, show.node.label=FALSE, font = 4, label.offset = 0.2, x.lim=c(-4,24))
HPDbars(agave, col="#00000050")
nodelabels(round(agave$posterior[which(agave$posterior>=0.5)],2),(which(agave$posterior>=0.5))+49, frame="none", adj=c(1,-0.3), cex=0.5)

plot(ladderize(agave, right=FALSE), cex=0.5, show.node.label=FALSE, font = 4, label.offset = 0.2, x.lim=c(-4,24))
HPDbars(agave, col="#00000050")
nodelabels(round(agave$posterior[which(agave$posterior>=0.5)],2),(which(agave$posterior>=0.5))+49, frame="none", adj=c(-0.1,0.5), cex=0.5)

axisGeo(strat2012, ages=FALSE, col=c("grey80", "white"))
axis(side=1, cex.axis=0.75, at=c(15.7,10.7,5.7,0.7, 0), labels=c("0","5","10","1",""), line=1.5, tck=-0.01, mgp=c(0,0,0))
mtext("Time (MYA)", at = (max(get("last_plot.phylo",envir = .PlotPhyloEnv)$xx) * 0.5), side = 1, line = 2, cex = 0.75)

branch.colors<- rep("black", 96)
branch.colors[15:35]<- "orange"
branch.colors[36:96]<- "purple"
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/Tree_2shifts.pdf")
plot(ladderize(agave, right=FALSE), cex=0.5, show.node.label=FALSE, font = 4, label.offset = 0.2, x.lim=c(-4,24), edge.color=branch.colors, edge.width=1.5)
HPDbars(agave, col="#00000050")
nodelabels(round(agave$posterior[which(agave$posterior>=0.5)],2),(which(agave$posterior>=0.5))+49, frame="none", adj=c(-0.1,0.5), cex=0.5)
axis(side=1, cex.axis=0.75, at=c(15.7,10.7,5.7,0.7,0), labels=c("0","5","10","15",""), line=0.9, tck=-0.01, mgp=c(0,0,0), cex.axis=0.5)
axisGeo(strat2012, ages=TRUE, col=c("grey80", "white"), cex=0.6)
dev.off()

################# funciones de las tasas de especiación y extinción ############
################################################################################
agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")
stem.age1 <- branching.times(agave)['57']
stem.age1
stem.age2 <- stem.age1
stem.age2
     57 
8.06559 

whole.age <- branching.times(agave)['50']
whole.age
     # 50 
# 15.66983 

xA1 <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeA1_stem_2016.10.20.xls", header=TRUE)
xA1


BlinA1 <- function(t){xA1[1,"lambda"] + xA1[1,"alfa"] * t}

BexpA1 <- function(t){xA1[2,"lambda"] * exp( xA1[2,"alfa"] * t)}

DconsA1 <- xA1[1,"mu"]# puede ser la fila dos también, pero son prácticamente iguales
rlinA1<-function(t){BlinA1(t)-DconsA1}
rexpA1<-function(t){BexpA1(t)-DconsA1}

xA2 <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeA2_stem_2016.10.20.xls", header=TRUE)
xA2

BconsA2 <- xA2[1,"lambda"]# puede ser la fila dos también, pero son prácticamente iguales

DconsA2 <- xA2[1,"mu"]# puede ser la fila dos también, pero son prácticamente iguales

rA2 <- BconsA2-DconsA2

z <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves/2shift_BBcladeA1_A2_2016.10.22.xls", header=TRUE)
z
Blin_BB_A1A2 <- function(t){z[1,"lambda"] + z[1,"alfa"] * t}
Dcons_BB_A1A2 <- -z[1,"mu"]# puede ser la fila dos también, pero son prácticamente iguales
rlin_BB_A1A2<-function(t){Blin_BB_A1A2(t)-Dcons_BB_A1A2}

################################################################################
################################################################################

#############################       BEST MODEL      ############################

################################################################################
# Best Model: Two shifts, 3 dynamics
# cladeA1: Agaves s.l.
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/Agaves_exp_rpanda.pdf", height= 3.5, width=4.5)
plot(x=NULL, main=NULL, xlab="", ylab="", xlim=c(whole.age,0), ylim=c(-1,2), axes=FALSE)
axis(2, at=seq(-1,2,1), labels=T)
abline(h=0, col="gray")
# curve(expr=BlinA1, from=stem.age1, to=0, add=TRUE, col="purple", lwd=1.2, lty="dotted")
curve(expr=BexpA1, from=stem.age1, to=0, add=TRUE, col="purple", lwd=1.2, lty="dotted")
segments(x0=stem.age1, x1=0, y0=DconsA1, y1=DconsA1, col="purple", lwd=1.2, lty="longdash")
curve(expr=rexpA1, from=stem.age1, to=0, add=TRUE, col="purple", lwd=2)
dev.off()

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/Agaves_lin_rpanda.pdf", height= 3.5, width=4.5)
plot(x=NULL, main=NULL, xlab="", ylab="", xlim=c(whole.age,0), ylim=c(-1,2), axes=FALSE)
axis(2, at=seq(-1,2,1), labels=T)
abline(h=0, col="gray")
curve(expr=BlinA1, from=stem.age1, to=0, add=TRUE, col="purple", lwd=1.2, lty="dotted")
# curve(expr=BexpA1, from=stem.age1, to=0, add=TRUE, col="purple", lwd=1.2, lty="dotted")
segments(x0=stem.age1, x1=0, y0=DconsA1, y1=DconsA1, col="purple", lwd=1.2, lty="longdash")
curve(expr=rlinA1, from=stem.age1, to=0, add=TRUE, col="purple", lwd=2)
dev.off()

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/Agaves_rpanda_legend.pdf", height= 3.5, width=4.5)
plot(x=NULL, main=NULL, xlab="", ylab="", xlim=c(whole.age,0), ylim=c(-1,2), axes=FALSE)
legend(x=15, y=2, c("extinction rate", "diversification rate = speciation rate"), col="purple", bty='n', lty=c("dashed", "solid"), cex=0.5)
dev.off()

# cladeA2: Furcraea-Beschorneria
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/Furcrea_rpanda.pdf", height= 3.5, width=4.5)
plot(x=NULL, main=NULL, xlab="", ylab="", xlim=c(whole.age,0), ylim=c(-1,2), axes=FALSE)
axis(2, at=seq(-1,2,1), labels=T)
abline(h=0, col="gray")
segments(x0=stem.age2, x1=0, y0=BconsA2, y1=BconsA2, col="orange", lwd=1.2, lty="dotted")
segments(x0=stem.age2, x1=0, y0=rA2, y1=rA2, col="orange", lwd=2)
segments(x0=stem.age2, x1=0, y0=DconsA2, y1=DconsA2, col="orange", lwd=1.2, lty="longdash")
dev.off()

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/Furcrea_rpanda_legend.pdf", height= 3.5, width=4.5)
plot(x=NULL, main=NULL, xlab="", ylab="", xlim=c(whole.age,0), ylim=c(-1,2), axes=FALSE)
legend(x=15, y=2, c("extinction rate", "diversification rate = speciation rate"), col="orange", bty='n', lty=c("dashed", "solid"), cex=0.5)
dev.off()

# BB clade A1 and clade A2
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/Backbone_Agave_Furcrea_rpanda.pdf", height= 3.5, width=4.5)
plot(x=NULL, main=NULL, xlab="", ylab="", xlim=c(whole.age,0), ylim=c(-1,2), axes=FALSE)
axis(2, at=seq(-1,2,1), labels=T)
abline(h=0, col="gray")
curve(expr=Blin_BB_A1A2, from=whole.age, to=0, add=TRUE, col="black", lwd=1.2, lty="dotted")
segments(x0=whole.age, x1=0, y0=Dcons_BB_A1A2, y1=Dcons_BB_A1A2, col="black", lwd=1.2, lty="longdash")
curve(expr=rlin_BB_A1A2, from=whole.age, to=0, add=TRUE, col="black", lwd=2)
dev.off()

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/Backbone_Agave_Furcrea_rpanda_legend.pdf", height= 3.5, width=4.5)
plot(x=NULL, main=NULL, xlab="", ylab="", xlim=c(whole.age,0), ylim=c(-1,2), axes=FALSE)
legend(x=10, y=2, c("speciation rate", "extinction rate", "diversification rate"), col="black", bty='n', lty=c("dotted", "dashed", "solid"), cex=0.5)
dev.off()

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_NuriAgaves/time_axis.pdf", height= 3.5, width=4.5)
plot(x=NULL, main=NULL, xlab="", ylab="", xlim=c(whole.age,0), ylim=c(-1,2), axes=FALSE)
axis(side=1, at=c(0,5,10,15,15.7), labels=c("0","5","10","15",""), line=1, tck=-0.03, mgp=c(0,0.25,0), cex.axis=1)
mtext(c("15"),at=15, line=-10.6)
dev.off()
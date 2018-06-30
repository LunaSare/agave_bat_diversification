# 26 de octubre 2016
# Cambio en Agave s.l., en Furcraea-Beschorneria y en Hesperoyucca-Hesperaloe

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/3shifts_NuriAgaves")

lapply(c("ape", "picante", "geiger"), require, character.only=TRUE)

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/uncond.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC.R')


agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")
pdf("agavaceae.nodelabels.pdf", width=7, height=15)
plot(agave)
nodelabels()
dev.off()

cladeA1 <- extract.clade(agave, 58) #Agave s.l.
plot(cladeA1)


stem.age1 <- branching.times(agave)['57']
stem.age1
     # 57 
# 8.06559 

cladeA2 <- extract.clade(agave, 88) # Furcraea-Beschorneria
plot(cladeA2)


stem.age2 <- branching.times(agave)['57']
stem.age2
# [1] 8.06559


cladeB2 <- extract.clade(agave, 52) # Hesperoyucca-Hesperaloe
plot(cladeB2)


stem.age3 <- branching.times(agave)['51']
stem.age3
      # 51 
# 12.09561 

##################################################
BBcladeA1_A2_B2 <- drop.tip(agave, c(cladeA1$tip.label,cladeA2$tip.label, cladeB2$tip.label))
plot(BBcladeA1_A2_B2)
whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA1_A2_B2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2, stem.age3))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA1_A2_B2)/49) #hay 49 spp. de Yucca

name <- "3shift_BBcladeA1_A2_B2_2016.10.26"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')

DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

#  26 oct 2016: reestimaré todos los AICc de este clado:
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')
x <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/3shifts_NuriAgaves/3shift_BBcladeA1_A2_B2_2016.10.26.xls", header=TRUE)
x
pars <- c(4,4,2,3,3,3,3)
aicc <- c()
for(i in 1:nrow(x)){
	print(i)
	aicc <- c(aicc, jointAIC6(p=pars[i], jointLH=x[i,"LH"], nobs=2))# tomando como nobs un número que no anule nobs-p-1
}
aicc
[1] -49.03678  17.43337  14.28602 -46.23123  16.76671  16.28602  16.28602

daicc <- aicc-min(aicc)
xnew <- cbind(x,AICc_new=aicc, deltaAICc_new=daicc)
index <- order(aicc)
xnew <- xnew[index,]
xnew
file.name <- "3shift_BBcladeA1_A2_B2_2016.10.26_new.xls"
write.table(xnew, file=file.name, sep="\t", quote=FALSE, row.names=FALSE)

##################################################
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 4+3+2+2
jointLik <- 21.85172516+-46.4985469+-20.9989552+-6.942677992
jointLik
[1] -52.58845

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 134.312

# backbone BconsDcons, A1 BexpDcons, A2 same (BconsDcons), B2 same (BconsDcons) 
n.parameters <- 2+3+2+2
jointLik <- -11.14300877+-46.68371+-20.9989552+-6.942677992
jointLik
[1] -85.76835

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 194.1521

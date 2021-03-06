# 23 de octubre 2016
# Cambio en Agave s.l. y en Hesperoyucca-Hesperaloe 

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves")

lapply(c("ape", "picante", "geiger"), require, character.only=TRUE)

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/uncond.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC.R')


agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")

cladeA1 <- extract.clade(agave, 58) #Agave s.l.
plot(cladeA1)


stem.age1 <- branching.times(agave)['57']
stem.age1
     # 57 
# 8.06559 

cladeB2 <- extract.clade(agave, 52) # Hesperoyucca-Hesperaloe
plot(cladeB2)


stem.age2 <- branching.times(agave)['51']
stem.age2
# [1] 12.09561

##################################################
BBcladeA1_B2 <- drop.tip(agave, c(cladeA1$tip.label,cladeB2$tip.label))
plot(BBcladeA1_B2)
whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA1_B2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA1_B2)/(32+49)) #hay 81 spp. de Yucca y de Furcraea-Beschorneria descritas

name <- "2shift_BBcladeA1_B2_2016.10.23"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')

DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')
##################################################
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 3+3+2
jointLik <- -34.88019265+-46.4985469+-6.942677992
jointLik
[1] -88.32142

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 196.2428

BBcladeA1_B2 BexpDcons, A1 BexpDcons, B2 BconsDcons
n.parameters <- 3+3+2
jointLik <- -35.010162+-46.68371+-6.942677992
jointLik
[1] -88.63655

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 196.8731


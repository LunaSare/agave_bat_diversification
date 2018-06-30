# 26 de octubre 2016
# Cambio en Yucca y en Hesperoyucca-Hesperaloe 

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves")

lapply(c("ape", "picante", "geiger"), require, character.only=TRUE)

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/uncond.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC.R')


agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")

cladeB1 <- extract.clade(agave, 54) # Yucca
plot(cladeB1)


stem.age1 <- branching.times(agave)['51']
stem.age1
     # 57 
# 12.09561 

cladeB2 <- extract.clade(agave, 52) # Hesperoyucca-Hesperaloe
plot(cladeB2)


stem.age2 <- branching.times(agave)['51']
stem.age2
# [1] 12.09561

##################################################
BBcladeB1_B2 <- drop.tip(agave, c(cladeB1$tip.label,cladeB2$tip.label))
plot(BBcladeB1_B2)
whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeB1_B2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeB1_B2)/(208+32)) #hay 240 spp. de Agave s.l. y Furcraea-Beschorneria

name <- "2shift_BBcladeB1_B2_2016.10.26"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')

DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')
##################################################
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 3+3+2
jointLik <- -40.8398769+-8.256736602+-6.942677992
jointLik
[1] -56.03929

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 131.6786

# backbone and Yucca BexpDcons, HH BconsDcons
n.parameters <- 3+3+2
jointLik <- -72.05657689+-8.601220315+-6.942677992
jointLik
[1] -87.60048

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 194.801

# backbone BexpDcons, Yucca and HH BconsDcons
n.parameters <- 3+2+2
jointLik <- -72.05657689+-9.845259095+-6.942677992
jointLik
[1] -88.84451

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 194.4207

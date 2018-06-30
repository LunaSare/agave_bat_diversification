# 24 de octubre 2016
# Cambio en Furraea-Beschorneria y en Yucca 

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves")

lapply(c("ape", "picante", "geiger"), require, character.only=TRUE)

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/uncond.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC.R')


agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")

cladeA2 <- extract.clade(agave, 88) #Furcraea-Beschorneria
plot(cladeA2)


stem.age1 <- branching.times(agave)['57']
stem.age1
     # 57 
# 8.06559 

cladeB1 <- extract.clade(agave, 54) # Yucca
plot(cladeB1)


stem.age2 <- branching.times(agave)['51']
stem.age2
# [1] 12.09561

##################################################
BBcladeA2_B1 <- drop.tip(agave, c(cladeA2$tip.label,cladeB1$tip.label))
plot(BBcladeA2_B1)
whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA2_B1)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA2_B1)/(208+49)) #hay 257 spp. de Agave s.l. y de Yucca descritas

name <- "2shift_BBcladeA2_B1_2016.10.24"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')

DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')
##################################################
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 4+2+3
jointLik <- -63.39990509+-20.9989552+-8.256736602
jointLik
[1] -92.6556

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 207.9266

# backbone and F-B the same, Yucca B1 to BexpDcons
n.parameters <- 4+2+3
jointLik <- -63.39990509+-20.9989552+-8.601220315
jointLik
[1] -93.00008

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 208.6155

# backbone and F-B the same, Yucca B1 to BconsDcons
n.parameters <- 4+2+2
jointLik <- -63.39990509+-20.9989552+-9.845259095
jointLik
[1] -94.24412

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 208.0882

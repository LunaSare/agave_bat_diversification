# 23 de octubre 2016
# Cambio en Furcraea-Beschorneria + Yucca-Hesperoyucca-Hesperaloe (clade B)
# Tomé el flujo de trabajo de Pachy 1shift.cactaceae.crownR

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves")

lapply(c("ape", "picante", "geiger"), require, character.only=TRUE)

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/uncond.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC.R')


agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")

cladeA2 <- extract.clade(agave, 88) # Furcraea-Beschorneria
plot(cladeA2)

stem.age1 <- branching.times(agave)['57']
stem.age1
     # 57 
# 8.06559 

cladeB <- extract.clade(agave, 51) # Yucca-Hesperoyucca-Hesperaloe
plot(cladeB)


stem.age2 <- branching.times(agave)['50']
stem.age2
# [1] 15.66983

##################################################
BBcladeB_A2 <- drop.tip(agave, c(cladeB$tip.label,cladeA2$tip.label))
plot(BBcladeB_A2)


whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeB_A2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeB_A2)/208) #hay 208 spp. de Agave s.l.

name <- "2shift_BBcladeB_A2_2016.10.23"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################
# n.parameters <- 3+2+2
# jointLik <- -48.94710592+-18.69681914+-20.9989552
# jointLik
# [1] -88.64288

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 192.5585

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 3+2+2
jointLik <- -48.94710592+-18.69681914+-20.9989552
jointLik
[1] -88.64288

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 194.0175

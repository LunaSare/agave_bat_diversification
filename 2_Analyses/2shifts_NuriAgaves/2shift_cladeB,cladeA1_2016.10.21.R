# 21 de octubre 2016
# Cambio en Agave s.l. + Yucca-Hesperoyucca-Hesperaloe (clade B)
# Tomé el flujo de trabajo de Pachy 1shift.cactaceae.crownR

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves")

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

cladeB <- extract.clade(agave, 51)
plot(cladeB)


stem.age2 <- branching.times(agave)['50']
stem.age2
# [1] 15.66983

##################################################

# BBcladeB_A1 <- drop.tip(agave, c(cladeB$tip.label,cladeA1$tip.label))
# plot(BBcladeB_A1)
# whole.age <- branching.times(agave)['50']
# whole.age
# # [1] 15.66983
# PHYLO <- list(BBcladeB_A1)
# times <- list(whole.age)
# bb <- list("backbone1")
# spec.t <- list(stem.age1, stem.age2)
# branch.t <- list(NULL)
# condition <- list("crown")
# fraction <- list(Ntip(BBcladeB_A1)/32) #hay 32 spp. de Furcraea-Beschorneria descritas
# name <- "2shift_BBcladeB_A1_2016.10.21"
# source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')
# DATA <- get(name)
# file.name <- paste(name,".xls", sep="")
# source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')
# Made a mistake in assigning spec.t (did not concatenate), so had to repeat analyses, since 2shift_BBcladeB_A1_2016.10.21 are wrong
BBcladeB_A1 <- drop.tip(agave, c(cladeB$tip.label,cladeA1$tip.label))
plot(BBcladeB_A1)
whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeB_A1)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeB_A1)/32) #hay 32 spp. de Furcraea-Beschorneria descritas

name <- "2shift_BBcladeB_A1_2016.10.22"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

# da igual el orden de los spec times :)
# a continuación, la prueba:
spec.t <- list(c(stem.age2, stem.age1))
name <- "2shift_BBcladeB_A1_2016.10.23"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')
DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')
##################################################
# n.parameters <- 4+2+3
# jointLik <- 10.85420648+-18.69681914+-46.4985469
# jointLik
# [1] -54.34116

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 128.7753

# n.parameters <- 2+2+3
# jointLik <- -24.73611098+-18.69681914+-46.4985469
# jointLik
# [1] -89.93148

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 195.1357

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')
# original:
n.parameters <- 2+3+3
jointLik <- -24.73611098+-18.69681914+-46.4985469
jointLik
[1] -89.93148

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 199.463

# Backbone the same (BconsDcons, from 2shift_BBcladeB_A1_2016.10.23.xls), Clade B the same (BconsDcons), clade A1 to BexpDcons
n.parameters <- 2+2+3
jointLik <- -24.73611098+-18.69681914+-46.68371
jointLik
[1] -90.11664

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 196.965

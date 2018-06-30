# 19 de octubre 2016
# Cambio en Furcraea-Beschorneria (clade A2)
# Tomé el flujo de trabajo de Pachy 1shift.cactaceae.crownR

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10")

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

cladeA2 <- extract.clade(agave, 88) # Furcraea-Beschorneria
plot(cladeA2)
Ntip(cladeA2)

stem.age <- branching.times(agave)['57']
stem.age
     # 57 
# 8.06559 
PHYLO <- list(cladeA2)
times <- list(stem.age)
bb <- list(FALSE)
spec.t <- list(NULL)
branch.t <- list(NULL)
condition <- list(FALSE)
fraction <- list(Ntip(cladeA2)/32) #hay 32 spp. de Furcraea-Beschorneria

name <- "1shift_cladeA2_stem_2016.10.20"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################

BBcladeA2 <- drop.tip(agave, cladeA2$tip.label)
plot(BBcladeA2)
Ntip(BBcladeA2)

whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(stem.age)
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA2)/263) #hay 263 spp. de Agave s.l., Yucca + Hesperoyucca-Hesperalöe descritas

name <- "1shift_BBcladeA2_2016.10.20"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################
# n.parameters <- 3+2
# jointLik <- -74.52807124+-20.9989552
# jointLik
# [1] -95.52703

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 201.7207

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 3+2
jointLik <- -74.52807124+-20.9989552
jointLik
[1] -95.52703

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 202.4494 


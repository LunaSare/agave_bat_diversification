# 20 de octubre 2016
# Cambio en Agave s.l.  (clade A1)
# Tomé el flujo de trabajo de Pachy 1shift.cactaceae.crownR

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves")

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

cladeA1 <- extract.clade(agave, 58)  #Agave s.l.
plot(cladeA1)


stem.age <- branching.times(agave)['57']
stem.age
     # 57 
# 8.06559 
PHYLO <- list(cladeA1)
times <- list(stem.age)
bb <- list(FALSE)
spec.t <- list(NULL)
branch.t <- list(NULL)
condition <- list(FALSE)
fraction <- list(Ntip(cladeA1)/208) #hay 208 spp. de Agave s.l. 

name <- "1shift_cladeA1_stem_2016.10.20"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################

BBcladeA1 <- drop.tip(agave, cladeA1$tip.label)
plot(BBcladeA1)
Ntip(BBcladeA1)

whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA1)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(stem.age)
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA1)/87) #hay 87 spp. de Furcraea-Beschorneria, Yucca + Hesperoyucca-Hesperalöe descritas

name <- "1shift_BBcladeA1_2016.10.20"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################
# n.parameters <- 3+3
# jointLik <- -46.4985469+-43.33566519
# jointLik
# [1] -89.83421

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 192.6122

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')
# BexpDcons + BlinDcons
n.parameters <- 3+3
jointLik <- -43.33566519+-46.4985469
jointLik
# [1] -89.83421

jointAIC6(p = n.parameters, jointLH = jointLik, nobs = Ntip(agave))
# [1] 49
# [1] 193.6684


# BexpDcons + BexpDcons
PAR <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeA1_stem_2016.10.20.xls", header=TRUE)
n.parameters <- 3+3
jointLik <- -43.33566519+-46.68371
jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
# [1] 49
# [1] 194.0388
Ntip(cladeA1)
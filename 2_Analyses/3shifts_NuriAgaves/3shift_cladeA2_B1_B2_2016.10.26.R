# 26 de octubre 2016
# Cambio en Furcraea-Beschorneria, en Yucca y en Hesperoyucca-Hesperaloe

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

cladeA2 <- extract.clade(agave, 88) # Furcraea-Beschorneria
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


cladeB2 <- extract.clade(agave, 52) # Hesperoyucca-Hesperaloe
plot(cladeB2)


stem.age3 <- branching.times(agave)['51']
stem.age3
      # 51 
# 12.09561 

##################################################
BBcladeA2_B1_B2 <- drop.tip(agave, c(cladeA2$tip.label,cladeB1$tip.label, cladeB2$tip.label))
plot(BBcladeA2_B1_B2)
whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA2_B1_B2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2, stem.age3))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA2_B1_B2)/208) #hay 208 spp. de Agave s.l.

name <- "3shift_BBcladeA2_B1_B2_2016.10.26"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')

DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')
##################################################
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 3+2+3+2
jointLik <- -20.07689565+-20.9989552+-8.256736602+-6.942677992
jointLik
[1] -56.27527

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 138.34

# backbone BexpDcons, A2 same BconsDcons, B1 BconsDcons, B2 same (BconsDcons) 
n.parameters <- 3+2+2+2
jointLik <- -50.50751666+-20.9989552+-9.845259095+-6.942677992
jointLik
[1] -88.29441

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 199.2042

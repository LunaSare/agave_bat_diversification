# 23 de octubre 2016
# Cambio en Agave s.l.+ Furcraea-Beschorneria (cladeA) y en Yucca 

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves")

lapply(c("ape", "picante", "geiger"), require, character.only=TRUE)

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/uncond.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC.R')


agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")


cladeA <- extract.clade(agave, 57)
plot(cladeA)

stem.age1 <- branching.times(agave)['50']
stem.age1
     # 50 
# 15.66983  

cladeB1 <- extract.clade(agave, 54)
plot(cladeB1)


stem.age2 <- branching.times(agave)['51']
stem.age2
# [1] 12.09561

##################################################

BBcladeA_B1 <- drop.tip(agave, c(cladeA$tip.label,cladeB1$tip.label))
plot(BBcladeA_B1)


whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA_B1)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA_B1)/6) #hay 6 spp. de Hesperoyucca-Hesperaloe descritas

name <- "2shift_BBcladeA_B1_2016.10.23"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

#  24 oct 2016: reestimaré todos los AICc de este clado:
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')
x <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves/2shift_BBcladeA_B1_2016.10.23.xls", header=TRUE)
x
pars <- c(3,3,3,3,4,4,2)
aicc <- c()
for(i in 1:nrow(x)){
	print(i)
	aicc <- c(aicc, jointAIC6(p=pars[i], jointLH=x[i,"LH"], nobs=2))# tomando como nobs un número que no anule nobs-p-1
}
aicc
[1]  9.074581  9.074582 12.408942 12.553360 13.078636 13.220027  7.074581

daicc <- aicc-min(aicc)
xnew <- cbind(x,AICc_new=aicc, deltaAICc_new=daicc)
index <- order(aicc)
xnew <- xnew[index,]
xnew
file.name <- "2shift_BBcladeA_B1_2016.10.24.xls"
write.table(xnew, file=file.name, sep="\t", quote=FALSE, row.names=FALSE)

##################################################
# n.parameters <- 3+3+4
# jointLik <- -7.537290508+-70.44821155+-8.256960764
# jointLik
# [1] -86.24246

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 195.0732

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 2+3+3
jointLik <- -7.537290581+-70.44821155+-8.256736602
jointLik
[1] -86.24224

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 192.0845

Backbone BconsDcons, A BexpDcons, B1 BconsDcons
n.parameters <- 2+3+2
jointLik <- -7.537290581+-70.44821155+-9.845259095
jointLik
[1] -87.83076

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 192.3932

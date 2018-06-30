# 23 de octubre 2016
# Cambio en Agave s.l.+ Furcraea-Beschorneria (cladeA) y en Hesperoyucca-Hesperaloe 

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

cladeB2 <- extract.clade(agave, 52) # Hesperoyucca-Hesperaloe
plot(cladeB2)


stem.age2 <- branching.times(agave)['51']
stem.age2
# [1] 12.09561

##################################################

BBcladeA_B2 <- drop.tip(agave, c(cladeA$tip.label,cladeB2$tip.label))
plot(BBcladeA_B2)


whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA_B2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA_B2)/49) #hay 49 spp. de Yucca descritas

name <- "2shift_BBcladeA_B2_2016.10.23"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

#  24 oct 2016: reestimaré todos los AICc de este clado:
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')
x <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves/2shift_BBcladeA_B2_2016.10.23.xls", header=TRUE)
x
pars <- c(4,4,2,3,3,3,3)
aicc <- c()
for(i in 1:nrow(x)){
	print(i)
	aicc <- c(aicc, jointAIC6(p=pars[i], jointLH=x[i,"LH"], nobs=2))# tomando como nobs un número que no anule nobs-p-1
}
aicc
[1] 16.79119 16.88138 14.28602 16.12485 16.21486 16.28602 16.28602

daicc <- aicc-min(aicc)
xnew <- cbind(x,AICc_new=aicc, deltaAICc_new=daicc)
index <- order(aicc)
xnew <- xnew[index,]
xnew
file.name <- "2shift_BBcladeA_B2_2016.10.24.xls"
write.table(xnew, file=file.name, sep="\t", quote=FALSE, row.names=FALSE)

##################################################
# n.parameters <- 4+3+3
# jointLik <- -11.06226286+-70.44821155+-6.892976447
# jointLik
# [1] -88.40345

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 199.3951

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 2+3+2
jointLik <- -11.14300877+-70.44821155+-6.942677992
jointLik
[1] -88.5339

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 193.7995

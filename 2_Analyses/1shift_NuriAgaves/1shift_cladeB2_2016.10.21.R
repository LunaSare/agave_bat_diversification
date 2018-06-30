# 21 de octubre 2016
# Cambio en Hesperoyucca-Hesperaloe (clade B2)

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

cladeB2 <- extract.clade(agave, 52) # Hesperoyucca-Hesperaloe
plot(cladeB2)
Ntip(cladeB2)

stem.age <- branching.times(agave)['51']
     # 51 
#  12.09561
PHYLO <- list(cladeB2)
times <- list(stem.age)
bb <- list(FALSE)
spec.t <- list(NULL)
branch.t <- list(NULL)
condition <- list(FALSE)
fraction <- list(Ntip(cladeB2)/6) #hay 6 spp. de Hesperoyucca-Hesperaloe 

name <- "1shift_cladeB2_stem_2016.10.21"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V2.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V3.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V4.R')
jointAIC(p=3, jointLH=-6.892976447, phylo=list(cladeB2), tot_time= stem.age)
[1] 8
[1] 25.78595
jointAIC2(p=3, jointLH=-6.892976447, phylo=list(cladeB2), tot_time= stem.age)
[1] 8
[1] 25.78595
jointAIC3(p=3, jointLH=-6.892976447, phylo=list(cladeB2), tot_time= stem.age)
[1] 7
[1] 27.78595
jointAIC4(p=3, jointLH=-6.892976447, phylo=list(cladeB2), tot_time= stem.age)
[1] 4
[1] Inf
jointAIC5(p=3, jointLH=-6.892976447, phylo=list(cladeB2), tot_time= stem.age)
[1] 3
[1] -4.214047

jointAIC2(p=2, jointLH=-6.942677992, phylo=list(cladeB2), tot_time= stem.age)
[1] 8
[1] 20.28536
jointAIC3(p=2, jointLH=-6.942677992, phylo=list(cladeB2), tot_time= stem.age)
[1] 7
[1] 20.88536
jointAIC4(p=2, jointLH=-6.942677992, phylo=list(cladeB2), tot_time= stem.age)
[1] 4
[1] 29.88536
jointAIC5(p=2, jointLH=-6.942677992, phylo=list(cladeB2), tot_time= stem.age)
[1] 3
[1] Inf
#  23 oct 2016: reestimaré todos los AICc de este clado:
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')
x <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeB2_stem_2016.10.21.xls", header=TRUE)
n <- Ntip(cladeB2)*2-1 # esto es el número de ramas, incluyendo a la troncal
pars <- c(3,3,3,3,4,4,2)
aicc <- c()
for(i in 1:nrow(x)){
	print(i)
	aicc <- c(aicc, jointAIC6(p=pars[i], jointLH=x[i,"LH"], nobs=n))
}
aicc
[1] 43.78595 43.85972 43.88536 43.88536      Inf      Inf 23.88536
aicc <- c()
for(i in 1:nrow(x)){
	print(i)
	aicc <- c(aicc, jointAIC6(p=pars[i], jointLH=x[i,"LH"], nobs=6))# tomando con nobs el número de nodos (incluyendo las puntas y la raíz)
}
aicc
[1] 31.78595 31.85972 31.88536 31.88536 61.82313 61.85971 21.88536
aicc <- c()
for(i in 1:nrow(x)){
	print(i)
	aicc <- c(aicc, jointAIC6(p=pars[i], jointLH=x[i,"LH"], nobs=2))# tomando como nobs un número que no anule nobs-p-1
}
aicc
[1] 7.785953 7.859717 7.885356 7.885363 8.489801 8.526381 5.885356
daicc <- aicc-min(aicc)
xnew <- cbind(x,AICc_new=aicc, deltaAICc_new=daicc)
index <- order(aicc)
xnew <- xnew[index,]
xnew
file.name <- "1shift_cladeB2_stem_2016.10.23.xls"
write.table(xnew, file=file.name, sep="\t", quote=FALSE, row.names=FALSE)
##################################################

BBcladeB2 <- drop.tip(agave, cladeB2$tip.labels)
plot(BBcladeB2)

whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeB2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(stem.age)
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeB2)/289) #hay 289 spp. de Agave s.l. + Furcraea-Beschorneria, Yucca descritas

name <- "1shift_BBcladeB2_2016.10.21"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################
# n.parameters <- 3+3
# jointLik <- -6.892976447+-97.47416412
# jointLik
# [1] -104.3671

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 221.6781

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

# BexpDcons + BconsDcons
n.parameters <- 3+2
jointLik <- -97.47416412+-6.942677992
jointLik
[1] -104.4168

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 220.229

# 19 de octubre 2016
# Cambio en Agave s.l. + Furcraea-Beschorneria (clade A)
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

cladeA <- extract.clade(agave, 57)
plot(cladeA)


stem.age <- branching.times(agave)['50']
# [1] 15.66983

PHYLO <- list(cladeA)
times <- list(stem.age)
bb <- list(FALSE)
spec.t <- list(NULL)
branch.t <- list(NULL)
condition <- list(FALSE)
fraction <- list(Ntip(cladeA)/240) #hay 240 spp. de Agave s.l. + Furcraea-Beschorneria descritas

name <- "1shift_cladeA_stem_2016.10.19"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################

BBcladeA <- extract.clade(agave, 51)
plot(BBcladeA)


stem.age <- branching.times(agave)['50']
# [1] 15.66983

PHYLO <- list(BBcladeA)
times <- list(stem.age)
bb <- list("backbone1")
spec.t <- list(stem.age)
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA)/55) #hay 55 spp. de Yucca + Hesperoyucca-Hesperalöe descritas

name <- "1shift_BBcladeA_2016.10.19"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################
# n.parameters <- 2+3
# jointLik <- -18.84428883+-70.44821155
# [1] -89.2925

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 189.2517
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 2+3
jointLik <- -18.84428883+-70.44821155
jointLik
[1] -89.2925

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 189.9803


######### june 2018
# aiccs
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

# BACKBONE
PAR <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_BBcladeA_2016.10.19.xls", header=TRUE)
PAR
n.par <- c(2,3,3,3,3,4,4)
nuevos <- c()
# nobs es el tamaño de muestra, que en estudios comparativos corresponde al número de taxa en el árbol, aquí  puntas
for (i in 1:nrow(PAR)){
	print(i)
	aicci <- jointAIC6(p=n.par[i], jointLH=PAR[i, "LH"], nobs = Ntip(BBcladeA))
	print(aicci)
	nuevos <- c(nuevos, aicci)
}
nuevos
# nobs = 7: 44.68858 50.33662 50.50219 51.68858 51.68860 64.35005 64.50219

nuevos-min(nuevos)
# same aicc before and after :)

# CLADE A
PAR <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeA_stem_2016.10.19.xls", header=TRUE)
PAR
n.par <- c(3,3,4,4,3,2,3)
nuevos <- c()
# nobs es el tamaño de muestra, que en estudios comparativos corresponde al número de taxa en el árbol, aquí  puntas
for (i in 1:nrow(PAR)){
	print(i)
	aicci <- jointAIC6(p=n.par[i], jointLH=PAR[i, "LH"], nobs = Ntip(cladeA))
	print(aicci)
	nuevos <- c(nuevos, aicci)
}
nuevos
# nobs = 42: 147.5280 149.7522 149.9775 155.9237 156.9908 159.9633 167.1905

nuevos-min(nuevos)
# same aicc before and after :)

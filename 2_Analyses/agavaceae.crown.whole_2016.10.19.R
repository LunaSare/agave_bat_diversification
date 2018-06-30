# 19 de octubre 2016

# Volví a correr los análisis de Nuri Agaves para estar segura de los resultados.

# Tomé el flujo de trabajo de Pachy cactaceae.crown.whole.R

# modelos ajustados a toda la filogenia de Agavaceae

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10")

lapply(c("ape", "picante", "geiger"), require, character.only=TRUE)

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/uncond.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC.R')


agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")
plot(agave)
agave

crown.age <- max(branching.times(agave))
# [1] 15.66983

PHYLO <- list(agave)
times <- list(crown.age)
bb <- list(FALSE)
spec.t <- list(NULL)
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(agave)/295) #hay 295 spp. de Agavaceae descritas

name <- "agavaceae.crown.whole"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- "agavaceae.crown.whole_2016.10.19.xls"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

######### june 2018
# aiccs
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')
PAR <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/agavaceae.crown.whole_2016.10.19.xls", header=TRUE)
PAR
n.par <- c(3,4,3,4,2,3,3)
nuevos <- c()
# nobs es el tamaño de muestra, que en estudios comparativos corresponde al número de taxa en el árbol, aquí  puntas
for (i in 1:nrow(PAR)){
	print(i)
	aicci <- jointAIC6(p=n.par[i], jointLH=PAR[i, "LH"], nobs = Ntip(agave))
	print(aicci)
	nuevos <- c(nuevos, aicci)
}
nuevos
# nobs = 49: [1] 196.2890 198.6649 202.9773 205.0403 218.7530 219.9481 220.2807

nuevos-min(nuevos)
# same aicc before and after :)
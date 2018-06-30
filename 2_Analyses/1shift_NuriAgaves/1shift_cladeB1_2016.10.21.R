# 21 de octubre 2016
# Cambio en Yucca (clade B1)

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

cladeB1 <- extract.clade(agave, 54) # Yucca
plot(cladeB1)
Ntip(cladeB1)

stem.age <- branching.times(agave)['51']
stem.age
     # 51 
#  12.09561
PHYLO <- list(cladeB1)
times <- list(stem.age)
bb <- list(FALSE)
spec.t <- list(NULL)
branch.t <- list(NULL)
condition <- list(FALSE)
fraction <- list(Ntip(cladeB1)/49) #hay 49 spp. de Yucca 

name <- "1shift_cladeB1_stem_2016.10.21"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

# 23 oct 2016: reestimaré los AICc de este clado:
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')
x <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeB1_stem_2016.10.21.xls", header=TRUE)
n <- Ntip(cladeB1)*2-1 # esto es el número de ramas, incluyendo a la troncal, son 7
pars <- c(4,4,2,3,3,3,3)
aicc <- c()
for(i in 1:nrow(x)){
	print(i)
	aicc <- c(aicc, jointAIC6(p=pars[i], jointLH=x[i,"LH"], nobs=n))
}
aicc
[1] 44.51392 45.20753 26.69052 30.51347 31.20244 33.60907 33.69052
aicc <- c()
for(i in 1:nrow(x)){
	print(i)
	aicc <- c(aicc, jointAIC6(p=pars[i], jointLH=x[i,"LH"], nobs=2))
}
aicc
[1] 11.18059 11.87420 11.69052 10.51347 11.20244 13.60907 13.69052
daicc <- aicc-min(aicc)
xnew <- cbind(x,AICc_new=aicc, deltaAICc_new=daicc)
index <- order(aicc)
xnew <- xnew[index,]
xnew
file.name <- "1shift_cladeB1_stem_2016.10.23.xls"
write.table(xnew, file=file.name, sep="\t", quote=FALSE, row.names=FALSE)

##################################################

BBcladeB1 <- drop.tip(agave, cladeB1$tip.label)
plot(BBcladeB1)


whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeB1)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(stem.age)
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeB1)/246) #hay 246 spp. de Agave sl.l+Furcraea-Beschorneria, Hesperoyucca-Hesperalöe descritas

name <- "1shift_BBcladeB1_2016.10.21"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')


DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')

##################################################
# n.parameters <- 4+3
# jointLik <- -8.256960764+-84.19675471
# jointLik
# [1] -92.45372

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 200.1802

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 3+3
jointLik <- -84.19675471+-8.256736602
jointLik
[1] -92.45349

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 198.907

# BexpDcons B1
n.parameters <- 3+3
jointLik <- -84.19675471+-8.601220315
jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
# [1] 49
# [1] 199.596

# BconsDcons B1 
n.parameters <- 3+2
jointLik <- -84.19675471+-9.845259095
[1] -94.04201

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
# [1] 49
# [1] 199.4794

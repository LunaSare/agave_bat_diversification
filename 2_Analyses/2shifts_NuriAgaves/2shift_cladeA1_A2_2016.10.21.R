# 21 de octubre 2016
# Cambio en Agave s.l. y en Furcraea-Beschorneria

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

cladeA2 <- extract.clade(agave, 88) # Furcraea-Beschorneria
plot(cladeA2)


stem.age2 <- branching.times(agave)['57']
stem.age2
# [1] 8.06559

##################################################
# lo siguiente no lo usé, es incorrecto, porque el argumento spec.t está mal, tiene que ser un vector:
# BBcladeA1_A2 <- drop.tip(agave, c(cladeA1$tip.label,cladeA2$tip.label))
# plot(BBcladeA1_A2)
# whole.age <- branching.times(agave)['50']
# whole.age
# # [1] 15.66983
# PHYLO <- list(BBcladeA1_A2)
# times <- list(whole.age)
# bb <- list("backbone1")
# spec.t <- list(stem.age1, stem.age2)
# branch.t <- list(NULL)
# condition <- list("crown")
# fraction <- list(Ntip(BBcladeA1_A2)/55) #hay 55 spp. de Yucca-Hesperoyucca-HEsperaloe descritas
# name <- "2shift_BBcladeA1_A2_2016.10.21"
# source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')
# DATA <- get(name)
# file.name <- paste(name,".xls", sep="")
# source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')


BBcladeA1_A2 <- drop.tip(agave, c(cladeA1$tip.label,cladeA2$tip.label))
plot(BBcladeA1_A2)
whole.age <- branching.times(agave)['50']
whole.age
# [1] 15.66983

PHYLO <- list(BBcladeA1_A2)
times <- list(whole.age)
bb <- list("backbone1")
spec.t <- list(c(stem.age1, stem.age2))
branch.t <- list(NULL)
condition <- list("crown")
fraction <- list(Ntip(BBcladeA1_A2)/55) #hay 55 spp. de Yucca + Hesperoyucca-Hesperaloe descritas

name <- "2shift_BBcladeA1_A2_2016.10.22"
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/02.2015/mor&cond.Estimations.V4.R')
# [1] "BconsDcons"
# [1] "BlinDcons"
# [1] "BexpDcons"
# [1] "BconsDlin"
# [1] "BconsDexp"
# [1] "BlinDlin"
# [1] "BexpDexp"
# Error in integrate(Vectorize(r.int.0), x, y, stop.on.error = FALSE) : 
  # non-finite function value
  # sale este error, por lo que debo cambiar los parámetros iniciales de exploración de verosimilitud en las funciones exponenciales:
lamb.exp <- list(f.lamb=function(t,y){y[1] * exp( y[2] * t)}, lamb_par=c(0.1, 0.0001), cst.lamb=FALSE, expo.lamb=TRUE)
mu.exp <- list(f.mu=function(t,y){y[1] * exp( y[2] * t)}, mu_par=c(0.1, 0.0001), cst.mu=FALSE, expo.mu=TRUE)
# y volver a correr el modelo BexpDexp:
print("BexpDexp")
BexpDexp <- fitLikelihood(phylo=PHYLO, tot_time=times, backbone=bb, spec_times=spec.t, branch_times=branch.t, f.lamb=lamb.exp[[1]], f.mu=mu.exp[[1]], lamb_par=lamb.exp[[2]], mu_par=mu.exp[[2]], f=fraction, meth = "Nelder-Mead", cst.lamb=lamb.exp[[3]], cst.mu=mu.exp[[3]], expo.lamb=lamb.exp[[4]], expo.mu=mu.exp[[4]], fix.mu=F, cond=condition) 
#ahora sí a guardar los resultados:
results <- list(BconsDcons=BconsDcons, BlinDcons=BlinDcons, BexpDcons=BexpDcons, BconsDlin=BconsDlin, BconsDexp=BconsDexp, BlinDlin=BlinDlin, BexpDexp=BexpDexp)
assign(name, results)
save(list=name, file=paste(name, ".RData", sep=""))
print(name)

DATA <- get(name)
file.name <- paste(name,".xls", sep="")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/04.2015/mor&cond.Summary.V1.R')
##################################################
# resultados con 2shift_BBcladeA1_A2_2016.10.21:
# n.parameters <- 3+3+2
# jointLik <- 17.47520037+-46.4985469+-20.9989552
# jointLik
# [1] -50.0223

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 117.6998

# resultados con la función antigua de jointAIC
# n.parameters <- 3+3+2
# jointLik <- 15.4352473+-46.4985469+-20.9989552
# jointLik
# [1] -52.06225

# jointAIC(p=n.parameters, jointLH=jointLik, phylo=list(agave), tot_time=max(branching.times(agave)))
# [1] 96
# [1] 121.7797

source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/jointAIC_2016.10.23_V6.R')

n.parameters <- 3+3+2
jointLik <- 15.4352473+-46.4985469+-20.9989552
jointLik
[1] -52.06225

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 123.7245

# BconsCons Backbone, BexpDcons A1
n.parameters <- 2+3+2
jointLik <- -18.84428883+-46.68371+-20.9989552
jointLik
[1] -86.52695

jointAIC6(p=n.parameters, jointLH=jointLik, nobs=Ntip(agave))
[1] 49
[1] 189.7856

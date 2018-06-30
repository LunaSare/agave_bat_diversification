# May 11 2018

### CIs best model, two shifts: Agave s.l. and Furcraea-Beschorneria.


lapply(as.list(c("ape", "picante", "geiger")), require, character.only=TRUE)
#require("phyloch")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')


source("~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/maxlikbounds_BexpDcons.V8.R")

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01")
agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")
cladeA1 <- extract.clade(agave, 58)  # Agave s.l.
plot(cladeA1)
stem.age1 <- branching.times(agave)['57']
stem.age1
##########################################################################
# Agave s.l.STEM 

# Modelo BexpDcons
PAR <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeA1_stem_2016.10.20.xls", header=TRUE)
PAR
PHY <- list(cladeA1)
TIM <- list(stem.age1)
BB <- list(FALSE)
STIM <- list(NULL)
BTIM <- list(NULL)
COND <- list(FALSE)
FRA <- list(Ntip(cladeA1)/208) #hay 208 spp. de Agave s.l. 
RES <- list(LH = PAR[2,"LH"], lambda = PAR[2,"lambda"], alfa = PAR[2,"alfa"], mu = abs(PAR[2,"mu"]))
as.numeric(RES["LH"])
CI_cladeA1_BexpDcons_1 <- maxlikbounds_BexpDcons(phylo=PHY, tot_time=TIM, backbone=BB, 
                                                 spec_times=STIM, branch_times=BTIM, f=FRA, 
                                                 cond=COND, res=RES, chi=5.991, l.min=0.001, 
                                                 l.max=6, m.min=0, m.max=3, a.min=-2, 
                                                 a.max=2, l.by=0.1, a.by=0.1, m.by=0.1)
save(CI_cladeA1_BexpDcons_1, file="CI_cladeA1_BexpDcons_1_2018.05.11.RData")
CI_cladeA1_BexpDcons_1$CIs
# r0    lamb0      alpha          mu0        LH
# est       1.777452 1.777453 -0.2809011 3.088475e-07 -46.68371
# lowbound  1.101000 1.601000 -0.1000000 5.000000e-01 -48.99291
# highbound 2.401000 2.501000 -0.4000000 1.000000e-01 -49.11691

nrow(CI_cladeA1_BexpDcons_1$lik_vals_bound) # 162, few observations, we want around 1000
head(CI_cladeA1_BexpDcons_1$lik_vals_bound)
min(CI_cladeA1_BexpDcons_1$lik_vals_bound[,"alpha"]) # -0.5
max(CI_cladeA1_BexpDcons_1$lik_vals_bound[,"alpha"]) # -0.1
min(CI_cladeA1_BexpDcons_1$lik_vals_bound[,"mu0"]) # 0.1
max(CI_cladeA1_BexpDcons_1$lik_vals_bound[,"mu0"]) # 1.2
min(CI_cladeA1_BexpDcons_1$lik_vals_bound[,"lamb0"]) # 1.201
max(CI_cladeA1_BexpDcons_1$lik_vals_bound[,"lamb0"]) # 2.501

CI_cladeA1_BexpDcons_2 <- maxlikbounds_BexpDcons(phylo=PHY, tot_time=TIM, backbone=BB, 
                                                 spec_times=STIM, branch_times=BTIM, f=FRA, 
                                                 cond=COND, res=RES, chi=5.991, l.min=1.2, 
                                                 l.max=2.6, m.min=0, m.max=1.3, a.min=-0.6, 
                                                 a.max=-0.2, l.by=0.01, a.by=0.01, m.by=0.01)
save(CI_cladeA1_BexpDcons_2, file="CI_cladeA1_BexpDcons_2_2018.05.11.RData")
CI_cladeA1_BexpDcons_2$CIs
# r0    lamb0      alpha          mu0        LH
# est       1.777452 1.777453 -0.2809011 3.088475e-07 -46.68371
# lowbound  1.210000 1.320000 -0.2000000 1.100000e-01 -49.66682
# highbound 2.590000 2.600000 -0.5000000 1.000000e-02 -49.32512

nrow(CI_cladeA1_BexpDcons_2$lik_vals_bound) # 90837
head(CI_cladeA1_BexpDcons_2$lik_vals_bound)
# lamb0 alpha  mu0   r0        LH
# [1,]  1.78 -0.28 0.01 1.77 -46.71679
# [2,]  1.79 -0.28 0.01 1.78 -46.71692
# [3,]  1.76 -0.27 0.01 1.75 -46.71847
# [4,]  1.77 -0.27 0.01 1.76 -46.71912
# [5,]  1.77 -0.28 0.01 1.76 -46.71935
# [6,]  1.80 -0.28 0.01 1.79 -46.71973

min(CI_cladeA1_BexpDcons_2$lik_vals_bound[,"alpha"]) # -0.58
max(CI_cladeA1_BexpDcons_2$lik_vals_bound[,"alpha"]) # -0.2
min(CI_cladeA1_BexpDcons_2$lik_vals_bound[,"mu0"]) # 0.01
max(CI_cladeA1_BexpDcons_2$lik_vals_bound[,"mu0"]) # 0.82
min(CI_cladeA1_BexpDcons_2$lik_vals_bound[,"lamb0"]) # 1.23
max(CI_cladeA1_BexpDcons_2$lik_vals_bound[,"lamb0"]) # 2.6

##########################################################################
# BBcladeA1_A2
cladeA2 <- extract.clade(agave, 88) # Furcraea-Beschorneria
plot(cladeA2)
stem.age2 <- branching.times(agave)['57']  # [1] 8.06559

BBcladeA1_A2 <- drop.tip(agave, c(cladeA1$tip.label,cladeA2$tip.label))
plot(BBcladeA1_A2)
BBcladeA1_A2$tip.label
whole.age <- branching.times(agave)['50']  # 15.66983

# Modelo BlinDcons
source("~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/maxlikbounds_BlinDcons.V8.R")

PAR <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves/2shift_BBcladeA1_A2_2016.10.22.xls", header=TRUE)
PAR
PHY <- list(BBcladeA1_A2)
TIM <- list(15.66983)
BB <- list("backbone1")
STIM <- list(c(stem.age1, stem.age2))
BTIM <- list(NULL)
FRA <- list(Ntip(BBcladeA1_A2)/55)  # Yucca + Hesperoyucca-Hesperaloe diversity count
COND <- list("crown")
RES <- list(LH = PAR[1,"LH"], lambda = PAR[1,"lambda"], alfa = PAR[1,"alfa"], mu = abs(PAR[1,"mu"]))

# l.min=0.001, 
# l.max=2, m.min=0, m.max=1, a.min=-1, 
# a.max=1, l.by=0.01, a.by=0.01, m.by=0.01
# "set another range of initial parameter values"

# screen -r May 13, 22hrs 
CI_BBcladeA1_A2_BlinDcons_1 <- maxlikbounds_BlinDconsV8(phylo=PHY, tot_time=TIM, backbone=BB, 
                                                        spec_times=STIM, branch_times=BTIM, f=FRA, 
                                                        cond=COND, res=RES, chi=5.991, l.min=0.001, 
                                                        l.max=1, m.min=0, m.max=0.5, a.min=-0.5, 
                                                        a.max=0.5, l.by=0.01, a.by=0.001, m.by=0.01)

save(CI_BBcladeA1_A2_BlinDcons_1, file="CI_BBcladeA1_A2_BlinDcons_1_2018.05.11.RData")


# Modelo BconsDcons
source("~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/maxlikbounds_BconsDcons.V8.R")

PAR <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves/2shift_BBcladeA1_A2_2016.10.22.xls", header=TRUE)
PAR
PHY <- list(BBcladeA1_A2)
TIM <- list(15.66983)
BB <- list("backbone1")
STIM <- list(c(stem.age1, stem.age2))
BTIM <- list(NULL)
FRA <- list(Ntip(BBcladeA1_A2)/55)  # Yucca + Hesperoyucca-Hesperaloe diversity count
COND <- list("crown")
RES <- list(LH = PAR[3,"LH"], lambda = PAR[3,"lambda"], mu = abs(PAR[3,"mu"]))
CHI <- 3.841  # for a 2 params model
CI_BBcladeA1_A2_BconsDcons_1 <- maxlikbounds_BconsDconsV8(phylo=PHY, tot_time=TIM, backbone=BB, 
                                                        spec_times=STIM, branch_times=BTIM, f=FRA, 
                                                        cond=COND, res=RES, chi=CHI, l.min=0, 
                                                        l.max=1, m.min=0, m.max=1, l.by=0.01, m.by=0.01)
# l.min=0, 
# l.max=1, m.min=0, m.max=0.1, l.by=0.01, m.by=0.000000001
save(CI_BBcladeA1_A2_BconsDcons_1, file="CI_BBcladeA1_A2_BconsDcons_1_2018.05.14.RData")
CI_BBcladeA1_A2_BconsDcons_1$CIs
# r0     lamb0          mu0        LH
# est       0.3224192 0.3224193 9.123916e-08 -18.84429
# lowbound  0.2000000 0.2400000 4.000000e-02 -20.67638
# highbound 0.4800000 0.4900000 1.000000e-02 -20.62479

nrow(CI_BBcladeA1_A2_BconsDcons_1$lik_vals_bound) # 539
head(CI_BBcladeA1_A2_BconsDcons_1$lik_vals_bound)
# lamb0  mu0   r0        LH
# [1,]  0.33 0.01 0.32 -18.92123
# [2,]  0.34 0.01 0.33 -18.92868
# [3,]  0.32 0.01 0.31 -18.93173
# [4,]  0.35 0.01 0.34 -18.95298
# [5,]  0.31 0.01 0.30 -18.96144
# [6,]  0.36 0.01 0.35 -18.99310

min(CI_BBcladeA1_A2_BconsDcons_1$lik_vals_bound[,"mu0"]) # 0.01
max(CI_BBcladeA1_A2_BconsDcons_1$lik_vals_bound[,"mu0"]) # 0.27
min(CI_BBcladeA1_A2_BconsDcons_1$lik_vals_bound[,"lamb0"]) # 0.21
max(CI_BBcladeA1_A2_BconsDcons_1$lik_vals_bound[,"lamb0"]) # 0.58

CI_BBcladeA1_A2_BconsDcons_2 <- maxlikbounds_BconsDconsV8(phylo=PHY, tot_time=TIM, backbone=BB, 
                                                          spec_times=STIM, branch_times=BTIM, f=FRA, 
                                                          cond=COND, res=RES, chi=CHI, l.min=0.16, 
                                                          l.max=0.63, m.min=0, m.max=0.32, l.by=0.001, m.by=0.001)
save(CI_BBcladeA1_A2_BconsDcons_2, file="CI_BBcladeA1_A2_BconsDcons_2_2018.05.14.RData")
CI_BBcladeA1_A2_BconsDcons_2$CIs
# r0     lamb0          mu0        LH
# est       0.3224192 0.3224193 9.123916e-08 -18.84429
# lowbound  0.1980000 0.2330000 3.500000e-02 -20.73628
# highbound 0.4910000 0.4920000 1.000000e-03 -20.74866

nrow(CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound) # 55296
head(CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound)
# lamb0   mu0    r0        LH
# [1,] 0.323 0.001 0.322 -18.85201
# [2,] 0.324 0.001 0.323 -18.85206
# [3,] 0.322 0.001 0.321 -18.85215
# [4,] 0.325 0.001 0.324 -18.85228
# [5,] 0.321 0.001 0.320 -18.85247
# [6,] 0.326 0.001 0.325 -18.85268

min(CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound[,"mu0"]) # 0.001
max(CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound[,"mu0"]) # 0.277
min(CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound[,"lamb0"]) # 0.199
max(CI_BBcladeA1_A2_BconsDcons_2$lik_vals_bound[,"lamb0"]) # 0.588


##########################################################################
# cladeA2
lapply(as.list(c("ape", "picante", "geiger")), require, character.only=TRUE)
#require("phyloch")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')
source("~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/maxlikbounds_BconsDcons.V8.R")

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01")
agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")
cladeA2 <- extract.clade(agave, 88) # Furcraea-Beschorneria
plot(cladeA2)
stem.age2 <- branching.times(agave)['57']  # [1] 8.06559

PAR <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeA2_stem_2016.10.20.xls", header=TRUE)
PAR
PHY <- list(cladeA2)
TIM <- list(stem.age2)
BB <- list(FALSE)
STIM <- list(NULL)
BTIM <- list(NULL)
COND <- list(FALSE)
FRA <- list(Ntip(cladeA2)/32) # 32 total spp 
RES <- list(LH = PAR[1,"LH"], lambda = PAR[1,"lambda"], mu = abs(PAR[1,"mu"]))
CHI <- 3.841  # for a 2 params model

CI_cladeA2_BconsDcons_2 <- maxlikbounds_BconsDconsV8(phylo=list(cladeA2), tot_time=list(stem.age2), 
                                        backbone=list(F), spec_times=list(NULL), 
                                        branch_times=list(NULL), f=list(Ntip(cladeA2)/32), 
                                        cond=list(FALSE), 
                                        res=RES, chi=CHI, l.min=0, 
                                        l.max=2,m.min=0,m.max=1,l.by=0.001,m.by=0.001)
save(CI_cladeA2_BconsDcons_2, file="CI_cladeA2_BconsDcons_2_2018.05.14.RData")


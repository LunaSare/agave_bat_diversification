### 26 octubre 2016
### CIs mejor modelo de dos cambios: uno en Agave s.l. y otro en Furcraea-Beschorneria Nuri Agaves.

setwd("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/CIs")

lapply(as.list(c("ape", "picante", "geiger")), require, character.only=TRUE)
#require("phyloch")
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/fitLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/Condamine.3.12.2013/getLikelihood.R')
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011Framework/MorlonPnas2011codes/Psi_Phi_timevar.R')


source("~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/maxlikbounds_BlinDcons.V8.R")


agave <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut")
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

BBcladeA1_A2 <- drop.tip(agave, c(cladeA1$tip.label,cladeA2$tip.label))
plot(BBcladeA1_A2)

whole.age <- branching.times(agave)['50']
whole.age
      # 50 
# 15.66983 

##########################################################################
# Agave s.l.STEM 

# Modelo BlinDcons
x <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeA1_stem_2016.10.20.xls", header=TRUE)
x
PHYLO <- list(cladeA1)
times <- list(stem.age1)
bb <- list(FALSE)
spec.t <- list(NULL)
branch.t <- list(NULL)
condition <- list(FALSE)
fraction <- list(Ntip(cladeA1)/208) #hay 208 spp. de Agave s.l. 
RES <- list(LH=x[1,"LH"], lambda=x[1,"lambda"], alfa=x[1,"alfa"], mu=-x[1,"mu"])
as.numeric(RES["LH"])

# 26 oct 2016 17:40 MExico
# compu susana screen -r 21124
CI_cladeA1_BlinDcons <- maxlikbounds_BlinDconsV8(phylo=PHYLO, tot_time=times, backbone=bb, spec_times=spec.t, branch_times=branch.t, f=fraction, cond=condition, res=RES, chi=5.991,l.min=0,l.max=4,m.min=0.01,m.max=4,a.min=-3,a.max=3,l.by=0.1,m.by=0.1,a.by=0.1)

save(CI_cladeA1_BlinDcons, file="CI_cladeA1_BlinDcons_2016.10.26.RData")
load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/CIs/CI_cladeA1_BlinDcons_2016.10.26.RData")
ls(CI_cladeA1_BlinDcons)
CI_cladeA1_BlinDcons$probas
                                          # negative positive null non null
# net diversification rate, bound                 0       NA   NA       NA
# rate of change in speciation rate, bound        1        0   NA       NA
# extinction rate, bound                         NA       NA    0        1

CI_cladeA1_BlinDcons$CIs
                # r0    lamb0      alpha          mu0        LH
# est       1.477838 1.477838 -0.2108164 5.421874e-08 -46.49855
# lowbound  0.990000 1.800000 -0.1000000 8.100000e-01 -49.19375
# highbound 2.490000 3.200000 -0.5000000 7.100000e-01 -49.37190

nrow(CI_cladeA1_BlinDcons$lik_vals_bound) # 226, few observations, we want around 1000
head(CI_cladeA1_BlinDcons$lik_vals_bound)
min(CI_cladeA1_BlinDcons$lik_vals_bound[,"alpha"]) #-0.5
max(CI_cladeA1_BlinDcons$lik_vals_bound[,"alpha"]) #-0.1
min(CI_cladeA1_BlinDcons$lik_vals_bound[,"mu0"]) #0.01
max(CI_cladeA1_BlinDcons$lik_vals_bound[,"mu0"]) #1.31
min(CI_cladeA1_BlinDcons$lik_vals_bound[,"lamb0"]) #1.2
max(CI_cladeA1_BlinDcons$lik_vals_bound[,"lamb0"]) #3.2

# 27 oct 2016 16:12 MExico
# compu susana screen -r 21124
CI_cladeA1_BlinDcons_fine <- maxlikbounds_BlinDconsV8(phylo=PHYLO, tot_time=times, backbone=bb, spec_times=spec.t, branch_times=branch.t, f=fraction, cond=condition, res=RES, chi=5.991,l.min=1.1,l.max=3.3,m.min=0,m.max=1.4,a.min=-0.6,a.max=0,l.by=0.001,m.by=0.001,a.by=0.001)
print(Sys.time)
# lo paré porque iba super lento, un día después iba en lambda = 1.106
# lo volveré a hacer con by´s de 0.01

# 28 oct 2016 14:02 MExico
# compu susana screen -r 21124
CI_cladeA1_BlinDcons_fine <- maxlikbounds_BlinDconsV8(phylo=PHYLO, tot_time=times, backbone=bb, spec_times=spec.t, branch_times=branch.t, f=fraction, cond=condition, res=RES, chi=5.991,l.min=1.1,l.max=3.3,m.min=0,m.max=1.4,a.min=-0.6,a.max=0,l.by=0.01,m.by=0.01,a.by=0.01)
CI_cladeA1_BlinDcons_2016.10.28 <- CI_cladeA1_BlinDcons_fine
save(CI_cladeA1_BlinDcons_2016.10.28, file="CI_cladeA1_BlinDcons_2016.10.28.RData")
load("/Users/luna/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/CIs/CI_cladeA1_BlinDcons_2016.10.28.RData")

CI_cladeA1_BlinDcons_2016.10.28$CIs
#                r0    lamb0      alpha          mu0        LH
# est       1.477838 1.477838 -0.2108164 5.421874e-08 -46.49855
# lowbound  0.930000 1.940000 -0.0800000 1.010000e+00 -49.48789
# highbound 2.930000 3.230000 -0.4600000 3.000000e-01 -48.04073
CI_cladeA1_BlinDcons_2016.10.28$probas
# negative positive null non null
# net diversification rate, bound                 0        1    0        1
# rate of change in speciation rate, bound        1        0   NA       NA
# extinction rate, bound                         NA       NA    1        0

nrow(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound) # 213 967
head(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound)
min(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound[,"alpha"]) #-0.6
max(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound[,"alpha"]) #-0.04
min(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound[,"mu0"]) #0.01
max(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound[,"mu0"]) #1.38
min(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound[,"lamb0"]) #1.1
max(CI_cladeA1_BlinDcons_2016.10.28$lik_vals_bound[,"lamb0"]) #3.3

##########################################################################
# Furcraea-Beschorneria STEM 
source('~/Google Drive/GoogleDrive.7nov2014/Analysis/Diversification/Morlon2011&Condamine/maxlikbounds_BconsDcons.V8.R', chdir = TRUE)
# Modelo BconsDcons
y <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/1shift_NuriAgaves/1shift_cladeA2_stem_2016.10.20.xls", header=TRUE)

# 26 oct 2016 17:45 MExico
# compu susana screen -r 21039
CI_cladeA2 <- maxlikbounds_BconsDconsV8(phylo=list(cladeA2), tot_time=list(stem.age2), backbone=list(F), spec_times=list(NULL), branch_times=list(NULL), f=list(Ntip(cladeA2)/32), cond=lsit(FALSE), res=list(LH=y[1,"LH"], lambda=y[1,"lambda"], mu=y[1,"mu"]), chi=3.841,l.min=0,l.max=2,m.min=y[1,"mu"],m.max=0.1,l.by=0.01,m.by=0.01)
# CIs cannot be estimated for the range of parameters, try new ones:
CI_cladeA2 <- maxlikbounds_BconsDconsV8(phylo=list(cladeA2), tot_time=list(stem.age2), backbone=list(F), spec_times=list(NULL), branch_times=list(NULL), f=list(Ntip(cladeA2)/32), cond=lsit(FALSE), res=list(LH=y[1,"LH"], lambda=y[1,"lambda"], mu=y[1,"mu"]), chi=3.841,l.min=0,l.max=2,m.min=y[1,"mu"],m.max=1,l.by=0.001,m.by=0.001)
[1] "LHs can't be calculated from initial parameter range; change initial parameters"
# try new ones: NOOOO, estaba mal el argumento en cond, puse lsit en vez de list!!! Volveré a intentar con los parámetros del primer intento:

# 27 oct 2016 15:48 MExico
# compu susana screen -r 21039
CI_cladeA2 <- maxlikbounds_BconsDconsV8(phylo=list(cladeA2), tot_time=list(stem.age2), 
                                        backbone=list(F), spec_times=list(NULL), 
                                        branch_times=list(NULL), f=list(Ntip(cladeA2)/32), 
                                        cond=list(FALSE), 
                                        res=list(LH=y[1,"LH"], lambda=y[1,"lambda"], 
                                                 mu=y[1,"mu"]), chi=3.841, l.min=0, 
                                        l.max=2,m.min=0,m.max=0.2,l.by=0.001,m.by=0.001)
# este vez sí salió
CI_cladeA2_BconsDcons <- CI_cladeA2
save(CI_cladeA2_BconsDcons, file="CI_cladeA2_BconsDcons_2016.10.27.RData")
CI_cladeA2_BconsDcons$CIs
#                 r0     lamb0          mu0        LH
# est       0.563784 0.5637846 5.645275e-07 -20.99896
# lowbound  0.272000 0.4720000 2.000000e-01 -22.91935
# highbound 0.861000 0.8620000 1.000000e-03 -22.91684

CI_cladeA2_BconsDcons$probas
#                                 negative       null  non null
# net diversification rate, bound        0         NA        NA
# extinction rate, bound                NA 0.04562135 0.9543786

min(CI_cladeA2_BconsDcons$lik_vals_bound[,"mu0"])  # 0.001
max(CI_cladeA2_BconsDcons$lik_vals_bound[,"mu0"])  # 0.2
min(CI_cladeA2_BconsDcons$lik_vals_bound[,"lamb0"])  # 0.337
max(CI_cladeA2_BconsDcons$lik_vals_bound[,"lamb0"])  # 0.971


nrow(CI_cladeA2_BconsDcons$lik_vals_bound)  # 103137

# 28 oct 2016 13:59 MExico
# compu susana screen -r 21039
CI_cladeA2_BconsDcons <- maxlikbounds_BconsDconsV8(phylo=list(cladeA2), tot_time=list(stem.age2), backbone=list(F), spec_times=list(NULL), branch_times=list(NULL), f=list(Ntip(cladeA2)/32), cond=list(FALSE), res=list(LH=y[1,"LH"], lambda=y[1,"lambda"], mu=y[1,"mu"]), chi=3.841,l.min=0,l.max=2,m.min=0,m.max=2,l.by=0.001,m.by=0.001)
# Error in `[<-`(`*tmp*`, "net diversification rate, bound", "positive",  : 
  # subscript out of bounds
##########################################################################
# BBcladeA1_A2

# Modelo BlinDcons

z <- read.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/2shifts_NuriAgaves/2shift_BBcladeA1_A2_2016.10.22.xls", header=TRUE)
z

# 26 oct 2016 17:50 MExico
# compu susana screen -r 21159
CI_BBcladeA1_A2_BlinDcons <- maxlikbounds_BlinDconsV8(phylo=list(cladeA1), tot_time=list(whole.age), backbone=list(F), spec_times=list(NULL), branch_times=list(NULL), f=list(Ntip(cladeA1)/208), cond=list(FALSE), res=list(LH=z[1,"LH"], lambda=z[1,"lambda"], alfa=z[1,"alfa"], mu=-z[1,"mu"]), chi=5.991,l.min=0,l.max=2,m.min=0.00000001,m.max=0.1,a.min=-1,a.max=1,l.by=0.01,m.by=0.001,a.by=0.01)
[1] "cannot estimate LH values with the proposed range of parameters"
[1] "set another range of initial parameter values"
# me equivoqué de clados :o
# lo volveré a lanzar

# 27 oct 2016 16:22 MExico
# compu susana screen -r 21159
CI_BBcladeA1_A2_BlinDcons 
<- maxlikbounds_BlinDconsV8(phylo=list(BBcladeA1_A2), tot_time=list(whole.age), 
                            backbone=list("backbone1"), spec_times=list(c(stem.age1, stem.age2)), branch_times=list(NULL), f=list(Ntip(BBcladeA1_A2)/55), cond=list("crown"), res=list(LH=z[1,"LH"], lambda=z[1,"lambda"], alfa=z[1,"alfa"], mu=-z[1,"mu"]), chi=5.991,l.min=0,l.max=4,m.min=0,m.max=4,a.min=-3,a.max=3,l.by=0.01,m.by=0.01,a.by=0.01)
print(Sys.time())
save(CI_BBcladeA1_A2_BlinDcons, file="CI_BBcladeA1_A2_BlinDcons_2016.10.26.RData")
# [1] "cannot estimate LH values with the proposed range of parameters"
# [1] "set another range of initial parameter values"
# There were 50 or more warnings (use warnings() to see the first 50)

# no salió; volver a intentar
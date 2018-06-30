# 6 de mayo 2018

#figuras filogenia Phyllostomidae de Roberto
lapply(as.list(c("ape", "picante", "geiger","phyloch")), require, character.only=TRUE)
data(strat2012)
######################################### filogenia ####################################
branch.colors<- rep("black", 244)
branch.colors[25:29]<- "orange"
branch.colors[122:244]<- "green"
branch.colors[53:105]<- "purple"
plot(ladderize(phyllos, right=FALSE), cex=0.5, show.node.label=FALSE, font = 3, label.offset = 0.3, edge.color=branch.colors, x.lim=c(-5,50))
HPDbars(ladderize(phyllos), col="#00000050")


# 2013 tree
phyllos <- read.beast("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/ConcatenadoTimes.TreeAnnOut")
branch.colors<- rep("gray", 246)
# branch.colors[51:103]<- "orange"
branch.colors[54:105]<- "black" # clado de nectarívoros (también incluye a dos especies que son frugívoras)
tip.colors <- rep("gray", length(phyllos$tip.label))
tip.colors[c(76:77, 44, 82, 36:37, 51, 54:59, 85)]<- "black"
tip.colors[c(45:48, 52, 53, 38, 2:4, 79)]<- "purple"
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/bats_tree_ConcatenadoTimes.TreeAnnOut_2018.05.06.pdf", height=11)
plot(ladderize(phyllos, right=FALSE), cex=0.4, show.node.label=FALSE, font = 4,
     label.offset = 0.5, edge.color=branch.colors, edge.width=1.3, x.lim=c(-3,55),
     tip.color = tip.colors)
HPDbars(phyllos, col="#00000050")
nodelabels(round(phyllos$posterior[which(phyllos$posterior>=0.5)],2),(which(phyllos$posterior>=0.5))+123, frame="none", adj=c(-0.15,0.5), cex=0.3)
axis(side=1, cex.axis=0.75, at=c(36.2,31.2,26.2,21.2,16.2,11.2,6.2,1.2,0), labels=c("0","5","10","15","20","25","30","35",""), line=0, tck=-0.01, mgp=c(0,0,0), cex.axis=0.5)
axisGeo(strat2012, ages=TRUE, col=c("grey80", "white"), cex=0.6)
mtext(c("35"), at=1.2, line=-46.8, cex=0.53)
dev.off()


tip.colors[c(67, 66, 68, 13, 12, 11, 8, 10, 7, 9, 14)] <- "red"  # problematic names
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/bats_tree_ConcatenadoTimes.TreeAnnOut_2018.05.06_original_names.pdf", height=11)
plot(ladderize(phyllos, right=FALSE), cex=0.4, show.node.label=FALSE, font = 4,
label.offset = 0.5, edge.color=branch.colors, edge.width=1.3, x.lim=c(-3,55),
tip.color = tip.colors)
HPDbars(phyllos, col="#00000050")
nodelabels(round(phyllos$posterior[which(phyllos$posterior>=0.5)],2),(which(phyllos$posterior>=0.5))+123, frame="none", adj=c(-0.15,0.5), cex=0.3)
axis(side=1, cex.axis=0.75, at=c(36.2,31.2,26.2,21.2,16.2,11.2,6.2,1.2,0), labels=c("0","5","10","15","20","25","30","35",""), line=0, tck=-0.01, mgp=c(0,0,0), cex.axis=0.5)
axisGeo(strat2012, ages=TRUE, col=c("grey80", "white"), cex=0.6)
mtext(c("35"), at=1.2, line=-46.8, cex=0.53)
dev.off()


problematic_names <- phyllos$tip.label[c(67, 66, 68, 13, 12, 11, 8, 10, 7, 9, 14)]
gsub(".*_", "", problematic_names[1])
problematic_names <- c("Lampronycteris_brachyotis", "Micronycteris_Glyphonycteris_daviesi",
"Trinycteris_nicefori", "Enchisthenes_hartii", "Dermanura_tolteca", "Dermanura_phaeotis",
"Dermanura_aztecus", "Artibeus_glaucus", "Dermanura_anderseni",
"Dermanura_cinereus", "Artibeus_concolor")
# phyllos <- drop.tip(phyllos, 66)
phyllos$tip.label[c(67, 66, 68, 13, 12, 11, 8, 10, 7, 9, 14)] <- problematic_names

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/bats_tree_ConcatenadoTimes.TreeAnnOut_2018.05.06_corrected_names.pdf", height=11)
plot(ladderize(phyllos, right=FALSE), cex=0.4, show.node.label=FALSE, font = 4,
label.offset = 0.5, edge.color=branch.colors, edge.width=1.3, x.lim=c(-3,55),
tip.color = tip.colors)
HPDbars(phyllos, col="#00000050")
nodelabels(round(phyllos$posterior[which(phyllos$posterior>=0.5)],2),(which(phyllos$posterior>=0.5))+123, frame="none", adj=c(-0.15,0.5), cex=0.3)
axis(side=1, cex.axis=0.75, at=c(36.2,31.2,26.2,21.2,16.2,11.2,6.2,1.2,0), labels=c("0","5","10","15","20","25","30","35",""), line=0, tck=-0.01, mgp=c(0,0,0), cex.axis=0.5)
axisGeo(strat2012, ages=TRUE, col=c("grey80", "white"), cex=0.6)
mtext(c("35"), at=1.2, line=-46.8, cex=0.53)
dev.off()


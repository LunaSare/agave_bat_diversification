# 15 ago 2018

#figuras agave
lapply(as.list(c("ape", "picante", "geiger","phyloch")), require, character.only=TRUE)
agave <- read.beast("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut_nombresPub")

# tree
branch.colors<- rep("black", 96)
branch.colors[15:35]<- "orange"
branch.colors[36:96]<- "purple"

agave$tip.label[which(agave$tip.label == "Beschorneria_rigidae")] <- "Beschorneria_rigida"
agave$tip.label[which(agave$tip.label == "Beschorneria_calcicula")] <- "Beschorneria_calcicola"
agave$tip.label[which(agave$tip.label == "Furcraea_longeva")] <- "Furcraea_longaeva"
agave$tip.label[which(agave$tip.label == "Hesperalöe_funifera")] <- "Hesperaloe_funifera"
agave$tip.label[which(agave$tip.label == "Hesperalöe_nocturna")] <- "Hesperaloe_nocturna"
agave$tip.label[which(agave$tip.label == "Furcraea_vendrhausin")] <- "Furcraea_bedinghausii"


pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/agave_tree_2shifts.pdf")
plot(ladderize(agave, right=FALSE), cex=0.5, show.node.label=FALSE, font = 4, label.offset = 0.2, x.lim=c(-4,24), edge.color=branch.colors, edge.width=1.5)
HPDbars(agave, col="#00000050")
nodelabels(round(agave$posterior[which(agave$posterior>=0.5)],2),(which(agave$posterior>=0.5))+49, frame="none", adj=c(-0.1,0.5), cex=0.5)
axis(side=1, cex.axis=0.75, at=c(15.7,10.7,5.7,0.7,0), labels=c("0","5","10","15",""), line=0.9, tck=-0.01, mgp=c(0,0,0), cex.axis=0.5)
axisGeo(strat2012, ages=TRUE, col=c("grey80", "white"), cex=0.6)
dev.off()
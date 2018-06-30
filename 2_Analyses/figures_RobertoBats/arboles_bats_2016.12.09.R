# 9 dic 2016

#figuras filogenia Phyllostomidae de Roberto
lapply(as.list(c("ape", "picante", "geiger","phyloch")), require, character.only=TRUE)
######################################### filogenia ####################################
branch.colors<- rep("black", 244)
branch.colors[25:29]<- "orange"
branch.colors[122:244]<- "green"
branch.colors[53:105]<- "purple"
plot(ladderize(phyllos, right=FALSE), cex=0.5, show.node.label=FALSE, font = 3, label.offset = 0.3, edge.color=branch.colors, x.lim=c(-5,50))
HPDbars(ladderize(phyllos), col="#00000050")
nodelabels(round(phyllos$posterior[which(phyllos$posterior>=0.5)],2),(which(phyllos$posterior>=0.5))+49, frame="none", adj=c(1,-0.3), cex=0.5)
#nodelabels(round(agaveB$posterior,2), frame="none", adj=c(1,-0.3), cex=0.5)
#nodelabels(as.character(c(50:98)), 50:98)
#edgelabels(as.character(c(1:244)), 1:244)
axisPhylo(cex.axis=0.75)
mtext("Time (MYA)", at = (max(get("last_plot.phylo",envir = .PlotPhyloEnv)$xx) * 0.5), side = 1, line = 2, cex = 0.75)


# este es un árbol del 2013 en formato beast
phyllos <- read.beast("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/ConcatenadoTimes.TreeAnnOut")
phyllos.table <- read.beast.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/ConcatenadoTimes.TreeAnnOut")
write.table(phyllos.table, file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_RobertoBats/phyllostomidae.nodes.xls", sep="\t", quote=FALSE, row.names=FALSE)
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_RobertoBats/phyllos_nodelabels.pdf", height=30)
plot(ladderize(phyllos, right=FALSE), cex=0.5)
nodelabels(cex=0.5)
dev.off()

# Este es un árbol del 2016 en formato beast pero desconfigurado porque no puedo leerlo con read.beast
# De alguna manera lo convertí a newick, pero no sé como. El árbol en newick es wholeOK.new
wholeOK.nex <- read.nexus("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/datos.2marzo2016/Archivos adjuntos_201632/ConcatenadoTimes")
wholeOK.table <- read.beast.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/datos.2marzo2016/ConcatenadoTimes con colores")

TREE <- read.tree("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/datos.2marzo2016/wholeOK.new")

plot(ladderize(phyllos, right=FALSE))
plot(ladderize(TREE, right=FALSE))

plot(ladderize(phyllos, right=FALSE), show.tip.label=FALSE)
plot(ladderize(TREE, right=FALSE), show.tip.label=FALSE)

phyllos # tiene 123 puntas
TREE # tiene 120 puntas

phyllos$tip.label
TREE$tip.label

x <- match(phyllos$tip.label, TREE$tip.label)
y <- match(TREE$tip.label, phyllos$tip.label)
# Parece ser que el correcto es wholeOK, porque es el árbol en el que quitamos las especies con nombres raros y mal identificadas

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_RobertoBats/wholeOK_edgelabels.pdf", height=30)
plot(ladderize(TREE, right=FALSE), cex=0.5)
edgelabels(cex=0.5)
dev.off()
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_RobertoBats/ConcatenadoTimes_edgelabels.pdf", height=30)
plot(ladderize(phyllos, right=FALSE), cex=0.5)
edgelabels(cex=0.5)
dev.off()
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_RobertoBats/ConcatenadoTimes_tiplabels.pdf", height=30)
plot(ladderize(phyllos, right=FALSE), cex=0.5, label.offset= 3)
tiplabels(cex=0.5)
dev.off()

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_RobertoBats/wholeOK_nodelabels.pdf", height=30)
plot(ladderize(TREE, right=FALSE), cex=0.5)
nodelabels(cex=0.5)
dev.off()

branch.colors<- rep("black", 240)
# branch.colors[51:103]<- "orange"
branch.colors[51:103]<- "purple" # clado de nectarívoros (también incluye a dos especies que son frugívoras)
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_RobertoBats/wholeOK_nectarívoros_2016.12.09.pdf", height=11)
plot(ladderize(TREE, right=FALSE), cex=0.5, show.node.label=FALSE, font = 4, label.offset = 0.2, edge.color=branch.colors, edge.width=1.5)
# HPDbars(TREE, col="#00000050")
# nodelabels(round(TREE$posterior[which(TREE$posterior>=0.5)],2),(which(TREE$posterior>=0.5))+120, frame="none", adj=c(-0.1,0.5), cex=0.5)
axis(side=1, cex.axis=0.75, at=c(36.2,31.2,26.2,21.2,16.2,11.2,6.2,1.2,0), labels=c("0","5","10","15","20","25","30","35",""), line=0, tck=-0.01, mgp=c(0,0,0), cex.axis=0.5)
axisGeo(strat2012, ages=TRUE, col=c("grey80", "white"), cex=0.6)
dev.off()

branch.colors<- rep("black", 246)
# branch.colors[51:103]<- "orange"
branch.colors[53:105]<- "purple" # clado de nectarívoros (también incluye a dos especies que son frugívoras)
pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2016.10/Revision_Diversificacion_NuriAgaves_2016.10/figures_RobertoBats/concatenadoTimes_nectarívoros_2016.12.09.pdf", height=11.5)
plot(ladderize(phyllos, right=FALSE), cex=0.5, show.node.label=FALSE, font = 4, label.offset = 0.2, edge.color=branch.colors, edge.width=1.3, x.lim=c(-3,55))
HPDbars(phyllos, col="#00000050")
nodelabels(round(phyllos$posterior[which(phyllos$posterior>=0.5)],2),(which(phyllos$posterior>=0.5))+123, frame="none", adj=c(-0.15,0.5), cex=0.3)
axis(side=1, cex.axis=0.75, at=c(36.2,31.2,26.2,21.2,16.2,11.2,6.2,1.2,0), labels=c("0","5","10","15","20","25","30","35",""), line=0, tck=-0.01, mgp=c(0,0,0), cex.axis=0.5)
axisGeo(strat2012, ages=TRUE, col=c("grey80", "white"), cex=0.6)
mtext(c("35"), at=1.2, line=-49.3, cex=0.5)
dev.off()

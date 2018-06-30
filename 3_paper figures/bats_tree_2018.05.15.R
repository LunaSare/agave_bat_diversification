# 15 de mayo 2018

#figuras filogenia Phyllostomidae de Roberto\
devtools::load_all("~/Desktop/datelife")
lapply(as.list(c("ape", "picante", "geiger","phyloch")), require, character.only=TRUE)
data(strat2012)
#
# 2013 tree
phyllos <- read.beast("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/ConcatenadoTimes.TreeAnnOut")

problematic_names_index <- c(67, 66, 68, 13, 12, 11, 8, 10, 7, 9, 14)
new_names <- c("Lampronycteris_brachyotis", "Micronycteris_Glyphonycteris_daviesi",
                       "Trinycteris_nicefori", "Enchisthenes_hartii", "Dermanura_tolteca", "Dermanura_phaeotis",
                       "Dermanura_aztecus", "Artibeus_glaucus", "Dermanura_anderseni",
                       "Dermanura_cinereus", "Artibeus_concolor")
phyllos$tip.label[problematic_names_index] <- new_names
phyllos_final <- drop.tip(phyllos, 66)  # "Micronycteris_Glyphonycteris_daviesi" probable missidentification
phyllos_final <- ladderize(phyllos_final, right=FALSE)
names(phyllos_final)
phyllos_nodes <- tree_get_node_data(phyllos, c("node_number", "descendant_tips_label"))
phyllos_nodes2 <- phyllos_nodes 
phyllos_nodes2$descendant_tips_label <- sapply(phyllos_nodes$descendant_tips_label, function(x) {
  rm <- which(x == "Micronycteris_Glyphonycteris_daviesi")
  if(length(rm >0)) x <- x[-rm]
  return(x)
  })
# interesting behaviour of rapply:
# remove <- rapply(phyllos_nodes$descendant_tips_label, function(x) which(x == "Micronycteris_Glyphonycteris_daviesi"))
# phyllos_nodes2$descendant_tips_label[names(remove)]
                                                              
phyllos_final_nodes <- tree_get_node_data(phyllos_final, c("node_number", "descendant_tips_label"))
match(phyllos_nodes2$descendant_tips_label, phyllos_final_nodes$descendant_tips_label)
phyllos_nodes2$descendant_tips_label[105]
rapply(phyllos_final_nodes$descendant_tips_label, function(x) which(x== "Choeroniscus_godmani"))
length(phyllos_final_nodes$descendant_tips_label) + x = 243
243-121
226-122
match(phyllos_final_nodes$descendant_tips_label[104], phyllos_nodes2$descendant_tips_label[105])
match(phyllos_final_nodes$descendant_tips_label[104:105], phyllos_nodes2$descendant_tips_label[105:106])
match(phyllos_final_nodes$descendant_tips_label[105], phyllos_nodes2$descendant_tips_label[106])
# the order of names matter!!!!!
# we need to sorth them:
phyllos_nodes2$descendant_tips_label[106]
lapply(phyllos_nodes2$descendant_tips_label, sort)[106]
phyllos_nodes2$descendant_tips_label <- lapply(phyllos_nodes2$descendant_tips_label, sort)
phyllos_final_nodes$descendant_tips_label[105]
lapply(phyllos_final_nodes$descendant_tips_label, sort)[105]
phyllos_final_nodes$descendant_tips_label <- lapply(phyllos_final_nodes$descendant_tips_label, sort)

phyllos_nodes2$descendant_tips_label[31]
new_nodes_pos <- match(phyllos_final_nodes$descendant_tips_label, phyllos_nodes2$descendant_tips_label)  # nodes from original study in the order of new nodes
# yay this worked!!!
# basically, we have to identify the nodes that disappear when dropping a branch, and then 
phyllos_final$`height_95%_HPD_MIN` <- phyllos$`height_95%_HPD_MIN`[new_nodes_pos]
phyllos_final$`height_95%_HPD_MAX` <- phyllos$`height_95%_HPD_MAX`[new_nodes_pos]
phyllos_final$`posterior` <- phyllos$`posterior`[new_nodes_pos]
length(phyllos_final$`posterior`)
length(phyllos$`posterior`)
branch.colors<- rep("black", 245)
branch.colors[65:105]<- "purple" # nectarivorous clade Glossophaginae (includes two frugivore species)

pdf(file="~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/2018.01/bats_tree_ConcatenadoTimes.TreeAnnOut_2018.05.14_black_tips_final_tips.pdf", height=11)
plot(ladderize(phyllos_final, right=FALSE), cex=0.4, show.node.label=FALSE, font = 4,
     label.offset = 0.5, edge.color=branch.colors, edge.width=1.3, x.lim=c(-3,55))
HPDbars(ladderize(phyllos_final, right=FALSE), col="#00000050")
nodelabels(round(phyllos_final$posterior[which(phyllos_final$posterior>=0.5)],2),
           (which(phyllos_final$posterior>=0.5)) + length(phyllos_final$tip.label), # you need to sum the total number of tip labels to get the actual node positions
           frame="none", adj=c(-0.15,0.5), cex=0.3)
axis(side=1, cex.axis=0.75, at=c(36.2,31.2,26.2,21.2,16.2,11.2,6.2,1.2,0), labels=c("0","5","10","15","20","25","30","35",""), line=0, tck=-0.01, mgp=c(0,0,0), cex.axis=0.5)
axisGeo(strat2012, ages=TRUE, col=c("grey80", "white"), cex=0.6)
mtext(c("35"), at=1.2, line=-46.8, cex=0.53)
dev.off()


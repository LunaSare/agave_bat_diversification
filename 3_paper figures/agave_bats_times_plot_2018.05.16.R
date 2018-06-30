# may 16 2018

#

lapply(as.list(c("ape", "picante", "geiger","phyloch")), require, character.only=TRUE)

phyllos <- read.beast("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/ConcatenadoTimes.TreeAnnOut")
names(phyllos)
beast_phylo <- phyllos
phyllos_table <- read.beast.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Roberto/ConcatenadoTimes.TreeAnnOut")
beast_table <- phyllos_table
names(beast_table)
class(beast_table)
head(beast_table)
node_name <- c("Glossophagini", "Choeronycterini", "Anurina", "Glossophaga + Leptonycteris", "Glossophaginae")
is.character(node_name)
is.null(node_name)
is.null(agave_table)
# visually identify nodes of each clade on figure phyllos_nodelabels.pdf
crown_node <- c(236, 224, 231, 237, 223)
stem_node <- c(229, 223, 230, 236, 217)

plot_beast_ages <- function(beast_file = NULL, beast_object = NULL, crown_node, stem_node, node_name, node_color) {
  if(is.character(beast_file)){
    beast_table <- phyloch::read.beast.table(file = beast_file)
  } 
  # if(!is.table(beast_object)){
  # check that it has required age data
  #   try to convert to table if it has beast information
  #   otherwise stop()
  # }
  if(is.null(beast_file) & is.null(beast_object)){
    stop("no data to plot")
  }
  crown_node_index <- match(crown_node, beast_table[,"node"])
  stem_node_index <- match(stem_node, beast_table[,"node"])
  x_lim_max <- round(max(beast_table[,"height_95%_HPD_MAX"], na.rm = TRUE) + 5, digits = -1)  #rounding to the nearest 10
  # plot.phylo(beast_phylo, plot = FALSE, show.tip.label = FALSE, x.lim = c(0, x_lim_max))
  # phyloch::axisGeo(GTS = strat2012, cex = 0.5)
  # par(new = TRUE)  # necessary if you want to overlay two plots together
  node_order <- order(beast_table[crown_node_index, "height_median"], decreasing = T)
  par (xpd = T)  # necessary to show node names
  plot(x = beast_table[crown_node_index[node_order], "height_median"],  y = 1:length(node_order), 
       xlim = c(x_lim_max, 0), axes = F, xlab = "", ylab = "nodes", col= "#CC0000")
  points(x = beast_table[stem_node_index[node_order], "height_median"],  y = 1:length(node_order), col= "#660000")
  axis(side = 1, line = 2, col = "gray", col.ticks = "gray")
  mtext(side = 1, line = 4, text = "Time (myr)")
  text(x = (beast_table[crown_node_index[node_order], "height_median"] + 
         beast_table[stem_node_index[node_order], "height_median"]) / 2,
       y = 1:length(node_order), labels = node_name[node_order], col = "black", pos = 3)
  index <- 0
  # for(i in node_order){
  #   index <- index + 1
  #   segments(x0 = beast_table[crown_node_index[i], "height_95%_HPD_MAX"], x1 = beast_table[crown_node_index[i], "height_95%_HPD_MIN"], y0 = index, y1 = index, col = "#CC000050", lwd = 3)
  #   segments(x0 = beast_table[stem_node_index[i], "height_95%_HPD_MAX"], x1 = beast_table[stem_node_index[i], "height_95%_HPD_MIN"], y0 = index, y1 = index, col = "#66000050", lwd = 3)
  # }
  # this is way better:
  
}           
plot_beast_ages(beast_table = agave_table, crown_node = crown_node, stem_node = stem_node, node_name = node_name)

# plot of bats + agave:

agave <- read.beast("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut_nombresPub")
plot(agave)
nodelabels()

pdf("agave_bats_main_times_plot_2018.05.16.pdf", height = 8, width = 8)
par (xpd = T, oma = c(3,1,2,4))
CEX <- 2 
bats_y_pos <- c(1,2,3,5,6)
agave_y_pos <- 4
text_y_pos <- 0.5
CI_lwd <- 4
plot(x = beast_table[crown_node_index[node_order], "height_median"],  y = bats_y_pos, 
     xlim = c(40, 0), axes = F, xlab = "", ylab = "", col= "#CC0000", cex = CEX, lwd =CI_lwd)
points(pch = 2, x = beast_table[stem_node_index[node_order], "height_median"],  y = bats_y_pos, col= "#CC0000", cex = CEX, lwd =CI_lwd)  #"#660000"
axis(side = 1, line = 0.5, col = "gray", col.ticks = "gray", cex.axis = CEX, labels = NA)
axis(side = 1, line = 1, cex.axis = CEX, lwd = 0)
mtext(side = 1, line = 4, text = "Time (myr)", cex = CEX)
text(x = (beast_table[crown_node_index[node_order], "height_median"] + 
            beast_table[stem_node_index[node_order], "height_median"]) / 2,  
     y = c(1,2,3,5,6) + text_y_pos, labels = node_name[node_order], col = "black", font = c(1,1,1,3,1), cex = CEX)  #to automatically set them above the points, use pos = 3, 
# segments can be applied to a vector
segments(x0 = beast_table[crown_node_index[node_order], "height_95%_HPD_MAX"], 
         x1 = beast_table[crown_node_index[node_order], "height_95%_HPD_MIN"], 
         y0 = bats_y_pos, y1 = bats_y_pos, col = "#CC000060", lwd = CI_lwd)
segments(x0 = beast_table[stem_node_index[node_order], "height_95%_HPD_MAX"], 
         x1 = beast_table[stem_node_index[node_order], "height_95%_HPD_MIN"], 
         y0 = bats_y_pos, y1 = bats_y_pos, col = "#CC000060", lwd = CI_lwd)
# agave estimated crown node is 58
# agave stem node is 57
agave_table <- read.beast.table("~/Google Drive/GoogleDrive.7nov2014/Analysis/Nuria/Archivos11.2013/Agave5_14Nov13.TreeAnnOut_nombresPub")
head(agave_table)
match(58, agave_table[,"node"]) # 9
match(57, agave_table[,"node"]) # 8
segments(x0 = agave_table[9:8, "height_95%_HPD_MAX"], x1 = agave_table[9:8, "height_95%_HPD_MIN"], y0 = rep(agave_y_pos ,2), y1 = rep(agave_y_pos ,2), col = "#7f00ff60", lwd = CI_lwd)
points(pch = 1:2, x = agave_table[9:8, "height_median"],  y = rep(agave_y_pos ,2), col= "#7f00ff", cex = CEX, lwd =CI_lwd)
text(x = mean(agave_table[9:8, "height_median"]), y = agave_y_pos + text_y_pos, col = "black", labels = "Agave sensu lato", font = 3, cex = CEX)
legend(pch = 1:2, lty = 1, lwd = CI_lwd, x = 48, y = 4.8, 
       legend = c("Node median age", "+ 95% HPD CI", "crown node", "stem node"), # expression(italic(Agave~sensu~lato)*~crown~and~stem~node~age*","~95*"%"*~HPD)
       bty = "n", pt.cex = 2, col = c(rep("white",2), rep("#CC000050", 2)), y.intersp = 1.1, cex = 1.5)  
dev.off()

# Anurina node ages:
beast_table[crown_node_index[3], c("height_median", "height_95%_HPD_MIN", "height_95%_HPD_MAX")]
beast_table[stem_node_index[3], c("height_median", "height_95%_HPD_MIN", "height_95%_HPD_MAX")]

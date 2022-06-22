############################################################################
#                    ECOBIO - PhD
#
#       Classification and PCA on recombination maps
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com

# Supervisor: Sylvain GLEMIN
#             ECOBIO Lab


#============================================================================
# LOADING ENVIRONMENT
#============================================================================

# clear global environment: remove all variables
rm(list=ls(all=TRUE))

#----------------------
# Loading packages
library(rstudioapi)
library(ggplot2)
library(MareyMap)
library(stringr)
library(ade4)
library(factoextra)
library(ggdendro)
library(pvclust) # Clustering with bootstrap
# library(adegraphics)
library(egg)

#============================================================================
# Loading variables & objects
#============================================================================

# Get the directory of the file & set working directory
wd = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

# Custom functions
source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
source("sources/MareyMap.R") # Consolidate a MareyMap input file
source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)
source("sources/stats.R") # Statistical functions implemented

# Loading data
metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
chromosomes = read.table("tables/chromosome.stats.csv", header = TRUE, sep = ";")

#############################################################################
#============================================================================
#       ANALYSES OF SUMMARY STATISTICS
#============================================================================
#############################################################################
# Make a data frame with all quantitative summary statistics relevant
data.pca = chromosomes[,-c(1:7, 24:27)]
colnames(data.pca)

# Check if there is correlation between some variables - redundant information
matcor = cor(data.pca, method = "spearman")
# Remove diagonal values with correlation == 1
diag(matcor) = NA
# Change to tidy format
library(reshape2)
matcor = melt(matcor)
colnames(matcor) = c('var1', 'var2', 'correlation')
# Which variables are higly correlated (>0.8)?
matcor[which(abs(matcor$correlation) > 0.8),]
# density.markers.cM & nb.markers & density.markers.cM are highyl correlated
# Keep only one statistic of map quality (hence quality checking is important but it is not what we are looking for)
data.pca = chromosomes[,c(17, 20, 21, 22, 25, 26, 27)]
colnames(data.pca)
rownames(data.pca) = paste(chromosomes$set, chromosomes$chromosome, sep ="_")
# pretty names
colnames(data.pca) = c("Mean recombination rate", "Max recombination rate",  "Coefficient of variation", "Gini", "Centromeric index", "Periphery-bias ratio", "Brokenstick variance")

# Transform the centromeric index in a distance to the center
data.pca$`Centromeric index` = abs(data.pca$`Centromeric index`-0.5)
# No need to log-transform, because all variables are centered-reduced

#============================================================================
# PCA
#============================================================================
# PCA uses reduced-centered variables: same dimension allows the comparison of maps of different lengths
# pca = dudi.pca(data.pca, scannf = TRUE)
# Put NA values to 0
data.pca[is.na(data.pca)] = 0
pca = dudi.pca(df = data.pca, scannf = FALSE, nf = 3)
save(data.pca, file = "output/classification/data.pca.chromosomes.Rda")
save(pca, file = "output/classification/pca.chromosomes.Rda")

# Informations contained by each axis
pca$eig
barplot(pca$eig)
# If we kept axes 1&2, we had 73% of variance explained
# By keeping axes one to three, we had 86% of variance explained
pca$eig[1]/sum(pca$eig)
pca$eig[2]/sum(pca$eig)
pca$eig[3]/sum(pca$eig)
sum(pca$eig[1:2])/sum(pca$eig)
sum(pca$eig[1:3])/sum(pca$eig)

screeplot(pca, main = "Screeplot - Eigenvalues")
# scatter(pca)
s.label(pca$li, xax = 1, yax = 2, label = (1:nrow(data.pca)), boxes = FALSE)
pal = c(hcl(0,100, seq(20,100, length.out=30)), hcl(240,100, seq(100,20, length.out=30)))
ade4::s.class(pca$li, xax = 1, yax = 2, chromosomes$species, add.plot = FALSE, label = "", col = pal[1:nlevels(chromosomes$species)])

s.label(pca$li, xax = 1, yax = 3, label = (1:nrow(data.pca)), boxes = FALSE)
ade4::s.class(pca$li, xax = 1, yax = 3, chromosomes$species, add.plot = FALSE, label = "", col = pal[1:nlevels(chromosomes$species)])

# Wich variables were correlated to axes?
pca$co
s.corcircle(pca$co, xax = 1, yax = 2)
s.corcircle(pca$co, xax = 1, yax = 3)

thresh = 1/nrow(pca$co)
rownames(pca$co)[which(abs(pca$co$Comp1) > thresh)]
# "peripherybias_ratio" "mean.recrate"  "cv.recrate" and "max" were the most correlated to first axis
rownames(pca$co)[which(abs(pca$co$Comp2) > thresh)]
# "linkage.map.length.correctedHW" "max" "cv.recrate" were the most correlated to second axis
rownames(pca$co)[which(abs(pca$co$Comp3) > thresh)]
# "linkage.map.length.correctedHW" "phys.map.length" "centromeric_index" were the most correlated to third axis


# Blurry scatterplot with classes, yet, species seem globally to cluster together

# Axis one: chromosomes with higher recombination rate (mean and max)
# are also those with the lower heterogeneity (cv.recrate) and spatial distribution toward the periphery
# two groups of variables separated on axis one:
# mean rec rate, max rec rate, linkage map length vs phys map length, peripherybias_ratio and cv rec rate

# Axis two opposes centromeric index to all other variables

# Axis three opposes centromeric index to linkage map length
# Are longer genetic maps more asymetrical? (centromeric index is here an absolute distance to the center of the chromosome)

# ARRANGE PLOTS
# var.pca = fviz_pca_var(pca,
#                        col.var = "contrib", 
#                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                        repel = TRUE     
# )
# 
# fviz_pca_ind(pca, axes = c(1, 2), geom.ind = "point")
# 
# fviz_pca_ind(pca, axes = c(1, 2),
#              geom.ind = "point", # show points only (but not "text")
#              col.ind = chromosomes$species, # color by groups
#              repel = TRUE,
#              addEllipses = TRUE
# ) +
#   scale_shape_manual(values=c(rep(19, nlevels(chromosomes$species)))) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=14),
#         axis.title.y = element_text(color="black", size=14),
#         axis.text=element_text(size=14, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(5,"line"),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=14),
#         legend.position='none')

# detach(package:adegraphics)
# # AXES 1&2
# png(filename = "figures/article_one/PCA_Chromosomes_Axes1-2.png", units = "cm", width = 30, height = 10,
#     res = 150)
# # Set plot layout
# layout(mat = matrix(c(1, 2),
#                     nrow = 1,
#                     ncol = 2),
#        heights = c(4),    # Heights of the two rows
#        widths = c(6, 4))     # Widths of the two columns
# par(mar =  c(5.1, 4.1, 4.1, 2.1))
# layout.show(2)

# par(mfrow = c(1, 2))
# Bottom left - scatter(pca) with labels
# s.label(pca$li, xax = 1, yax = 2, label = (1:nrow(data.pca)), boxes = FALSE)
# Color gradient for chromosome size
color.gradient = function(x, colors = c("lightgrey", "black"), colsteps = 100) {
  return(colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
# list_sp = levels(chromosomes$species)
# list_chrsize = NA
# for (i in 1:length(list_sp)) {
#   list_chrsize[i] = data.pca.pooled$phys.map.length[rownames(data.pca.pooled) == list_sp[i]]
# }
# list_color = color.gradient(log(list_chrsize), colors = c("Blue", "Green", "Red"))
# # pal = c(hcl(0,100, seq(20,100, length.out=30)), hcl(240,100, seq(100,20, length.out=30)))
# s.class(pca$li, xax = 1, yax = 2, chromosomes$species, label = "", add.plot = FALSE, col = list_color)
# # Top left - Screeplot
# # screeplot(pca, main = "Screeplot - Eigenvalues")
# # Bottom right - Correlation circles
# par(mar =  c(2, 2, 4.1, 2.1))
# s.corcircle(pca$co, xax = 1, yax = 2, box = TRUE, clabel = 0.8)
# # par(mfrow = c(1, 1))
# dev.off()
# 
# # AXES 1&3
# png(filename = "figures/article_one/PCA_Chromosomes_Axes1-3.png", units = "cm", width = 30, height = 10,
#     res = 150)
# # Set plot layout
# layout(mat = matrix(c(1, 2),
#                     nrow = 1,
#                     ncol = 2),
#        heights = c(4),    # Heights of the two rows
#        widths = c(6, 4))     # Widths of the two columns
# par(mar =  c(5.1, 4.1, 4.1, 2.1))
# layout.show(2)
# # par(mfrow = c(1, 2))
# # Bottom left - scatter(pca) with labels
# # s.label(pca$li, xax = 1, yax = 2, label = (1:nrow(data.pca)), boxes = FALSE)
# # pal = c(hcl(0,100, seq(20,100, length.out=30)), hcl(240,100, seq(100,20, length.out=30)))
# list_color = color.gradient(log(list_chrsize), colors = c("Blue", "Green", "Red"))
# s.class(pca$li, xax = 1, yax = 3, chromosomes$species, add.plot = FALSE, label = "", col = list_color)
# # Top left - Screeplot
# # screeplot(pca, main = "Screeplot - Eigenvalues")
# # Bottom right - Correlation circles
# par(mar =  c(2, 2, 4.1, 2.1))
# s.corcircle(pca$co, xax = 1, yax = 3, box = TRUE, clabel = 0.8)
# # par(mfrow = c(1, 1))
# dev.off()
# # library(adegraphics)



#####################################
# # PLOTS
# # AXES 1&2
# Color gradient for chromosome size
list_chrsize = chromosomes$phys.map.length
# Plot of individuals 
plot.ind = fviz_pca_ind(pca, axes = c(1, 2), label = "none", # hide individual labels
                        col.ind = log(list_chrsize/1000000),
                        legend.title = "Mb") +
  scale_colour_viridis_c(option = "viridis", breaks = c(log(20), log(50), log(200), log(400)), labels = c(20, 50, 200, 400)) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
plot.ind
# Control variable colors using their contributions
plot.cor = fviz_pca_var(pca, axes = c(1, 2), col.var="contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE, legend.title = "Contrib.", labelsize = 4) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
plot.cor

# png("figures/article_one/PCA_Chromosomes_Axes1-2.png", units = "cm", width = 30, height = 10, res = 150)
# ggarrange(plot.ind, plot.cor, ncol = 2, nrow = 1, widths = c(2,1),
#           labels = c("(a)", "(b)"))
# dev.off()


# AXES 1&3
list_chrsize = chromosomes$phys.map.length
# Plot of individuals 
plot.ind = fviz_pca_ind(pca, axes = c(1, 3), label = "none", # hide individual labels
                        col.ind = log(list_chrsize/1000000),
                        legend.title = "Mb") +
  scale_colour_viridis_c(option = "viridis", breaks = c(log(20), log(50), log(200), log(400)), labels = c(20, 50, 200, 400)) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
plot.ind
# Control variable colors using their contributions
plot.cor = fviz_pca_var(pca, axes = c(1, 3), col.var="contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE, legend.title = "Contrib.", labelsize = 4) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
plot.cor
# png("figures/article_one/PCA_Chromosomes_Axes1-3.png", units = "cm", width = 30, height = 10, res = 150)
# ggarrange(plot.ind, plot.cor, ncol = 2, nrow = 1, widths = c(2,1),
#           labels = c("(a)", "(b)"))
# dev.off()



#============================================================================
# Hierarchical Classification
#============================================================================
# Same method as clustering with DAPC in 'adegenet'
load(file = "output/classification/data.pca.species.Rda")
### Clustering based on PCA individuals
d1 = dist(pca$li)
h1 = hclust(d1, method="ward.D")
plot(h1)
ggdendrogram(h1, rotate = FALSE, size = 2)
# order of species with the classification
h1$order
# We observed at least four distinct groups
ngroupes = 4
rect.hclust(h1,k=ngroupes)
# Clustering in four groups
groupes.h1 <- as.factor(cutree(h1,k=ngroupes))
# Individuals in groups
names(groupes.h1[which(groupes.h1 == 1)])
names(groupes.h1[which(groupes.h1 == 2)])
names(groupes.h1[which(groupes.h1 == 3)])
names(groupes.h1[which(groupes.h1 == 4)])

# The representation of clustering and pca were very blurry, with two much individuals
# Hence we performed the same analysis on pooled chromosomes per species
# With averaged summary statistics



#============================================================================
# POOL CHROMOSOMES AND AVERAGE SUMMARY STATISTICS PER SPECIES
#============================================================================
# Make a data frame with all quantitative summary statistics relevant
data.pca = chromosomes[,-c(1:5)]
colnames(data.pca)
# POOL PER SPECIES - MEAN()
data.pca.pooled = aggregate(data.pca, by = list(species = data.pca$species), mean)
rownames(data.pca.pooled) = data.pca.pooled$species
data.pca.pooled = data.pca.pooled[-c(1, 20)]

# Check if there is correlation between some variables - redundant information
matcor = cor(data.pca.pooled)
# Remove diagonal values with correlation == 1
diag(matcor) = NA
# Change to tidy format
library(reshape2)
matcor = melt(matcor)
colnames(matcor) = c('var1', 'var2', 'correlation')
# Which variables are higly correlated (>0.8)?
matcor[which(abs(matcor$correlation) > 0.8),]
# density.markers.cM & nb.markers & density.markers.cM are highyl correlated
# Keep only one statistic of map quality (hence quality checking is important but it is not what we are looking for)
data.pca.pooled = data.pca.pooled[,c(12, 15, 16, 17, 19, 20, 21)]
colnames(data.pca.pooled)
# pretty names
colnames(data.pca.pooled) = c("Mean recombination rate", "Max recombination rate",  "Coefficient of variation", "Gini", "Centromeric index", "Periphery-bias ratio", "Brokenstick variance")

# Transform the centromeric index in a distance to the center
data.pca.pooled$`Centromeric index` = abs(data.pca.pooled$`Centromeric index`-0.5)
# No need to log-transform, because all variables are centered-reduced

#============================================================================
# PCA
#============================================================================
# PCA uses reduced-centered variables: same dimension allows the comparison of maps of different lengths
# pca = dudi.pca(data.pca, scannf = TRUE)
# Put NA values to 0
data.pca.pooled[is.na(data.pca.pooled)] = 0
pca = dudi.pca(df = data.pca.pooled, scannf = FALSE, nf = 3)
save(data.pca.pooled, file = "output/classification/data.pca.species.Rda")
save(pca, file = "output/classification/pca.species.Rda")

# Informations contained by each axis
pca$eig
barplot(pca$eig)
# If we kept axes 1&2, we had 74% of variance explained
# By keeping axes one to three, we had 88% of variance explained
pca$eig[1]/sum(pca$eig)
pca$eig[2]/sum(pca$eig)
pca$eig[3]/sum(pca$eig)
sum(pca$eig[1:2])/sum(pca$eig)
sum(pca$eig[1:3])/sum(pca$eig)

# scatter(pca)
s.label(pca$li, xax = 1, yax = 2, label = (rownames(data.pca.pooled)), boxes = FALSE)
# pal = c(hcl(0,100, seq(20,100, length.out=30)), hcl(240,100, seq(100,20, length.out=30)))
# s.class(pca$li, xax = 1, yax = 2, chromosomes$species, add.plot = FALSE, label = "", col = pal[1:nlevels(chromosomes$species)])

s.label(pca$li, xax = 1, yax = 3, label = (rownames(data.pca.pooled)), boxes = FALSE)
# s.class(pca$li, xax = 1, yax = 3, chromosomes$species, add.plot = FALSE, label = "", col = pal[1:nlevels(chromosomes$species)])

# Wich variables were correlated to axes?
pca$co
s.corcircle(pca$co, xax = 1, yax = 2)
s.corcircle(pca$co, xax = 1, yax = 3)

thresh = 1/nrow(pca$co)
rownames(pca$co)[which(abs(pca$co$Comp1) > thresh)]
# "peripherybias_ratio" "mean.recrate"  "cv.recrate" and "max" were the most correlated to first axis
rownames(pca$co)[which(abs(pca$co$Comp2) > thresh)]
# "linkage.map.length.correctedHW" "max" "cv.recrate" were the most correlated to second axis
rownames(pca$co)[which(abs(pca$co$Comp3) > thresh)]
# "linkage.map.length.correctedHW" "phys.map.length" "centromeric_index" were the most correlated to third axis

save(data.pca.pooled, file = "output/classification/data.pca.species.Rda")

# Blurry scatterplot with classes, yet, species seem globally to cluster together

# Axis one: chromosomes with higher recombination rate (mean and max)
# are also those with the lower heterogeneity (cv.recrate) and spatial distribution toward the periphery
# two groups of variables separated on axis one:
# mean rec rate, max rec rate, linkage map length vs phys map length, peripherybias_ratio and cv rec rate

# Axis two opposes centromeric index to all other variables

# Axis three opposes centromeric index to linkage map length
# Are longer genetic maps more asymetrical? (centromeric index is here an absolute distance to the center of the chromosome)

# # PLOTS
# # AXES 1&2
# Color gradient for chromosome size
list_sp = levels(chromosomes$species)
list_chrsize = NA
for (i in 1:length(list_sp)) {
  list_chrsize[i] = mean(chromosomes$phys.map.length[chromosomes$species == list_sp[i]], na.rm = TRUE)
}
# Plot of individuals 
plot.ind = fviz_pca_ind(pca, axes = c(1, 2), repel = FALSE, geom.ind = c("point", "text"), col.ind = log(list_chrsize/1000000),
             legend.title = "Mb", labelsize = 4) +
  scale_colour_viridis_c(breaks = c(log(20), log(50), log(200), log(400)), labels = c(20, 50, 200, 400)) +
  xlim(-6, 6) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
plot.ind
# Control variable colors using their contributions
plot.cor = fviz_pca_var(pca, axes = c(1, 2), col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, legend.title = "Contrib.", labelsize = 4) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
plot.cor

# png("figures/article_one/PCA_Species_Axes1-2.png", units = "cm", width = 30, height = 10, res = 150)
# ggarrange(plot.ind, plot.cor, ncol = 2, nrow = 1, widths = c(2,1),
#           labels = c("(a)", "(b)"))
# dev.off()


# # AXES 1&3
# Plot of individuals 
plot.ind = fviz_pca_ind(pca, axes = c(1, 3), repel = FALSE, geom.ind = c("point", "text"), col.ind = log(list_chrsize/1000000),
                        legend.title = "Mb", labelsize = 4) +
  scale_colour_viridis_c(breaks = c(log(20), log(50), log(200), log(400)), labels = c(20, 50, 200, 400)) +
  xlim(-6, 6) +
  ylim(-2, 2.5) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=18, colour="black"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
plot.ind
# Control variable colors using their contributions
plot.cor = fviz_pca_var(pca, axes = c(1, 3), col.var="contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE, legend.title = "Contrib.", labelsize = 4) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
plot.cor

# png("figures/article_one/PCA_Species_Axes1-3.png", units = "cm", width = 30, height = 10, res = 150)
# ggarrange(plot.ind, plot.cor, ncol = 2, nrow = 1, widths = c(2,1),
#           labels = c("(a)", "(b)"))
# dev.off()





#============================================================================
# Hierarchical Classification
#============================================================================
load(file = "output/classification/data.pca.species.Rda")

#----------------------------------------------------------------------------
### Clustering based on PCA individuals
d1 = dist(pca$li)
h1 = hclust(d1, method="ward.D2")
plot(h1)
save(h1, file = "output/classification/h1.clust.Rda")

hclustering = as.dendrogram(h1)
list_sp = labels(hclustering)
list_chrsize = NA
for (i in 1:length(list_sp)) {
  list_chrsize[i] = mean(chromosomes$phys.map.length[chromosomes$species == list_sp[i]], na.rm = TRUE)
}
color.gradient = function(x, colors = c("lightgrey", "black"), colsteps = 100) {
  return(colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
list_color = color.gradient(log(list_chrsize))
# ggdendrogram(hclustering, rotate = FALSE, size = 2, tip.color = list_chrsize) + geom_text(colour = list_color)
library(dendextend)
dend = as.dendrogram(h1)
dend %>% dendextend::set("labels_col", list_color) %>% plot
dend %>% rect.dendrogram(k=5, border = 8, lty = 5, lwd = 2)

# order of species with the classification
h1$order
h1$height
# We observed at least five distinct groups
ngroupes = 5
plot(h1)
rect.hclust(h1,k=ngroupes)
# Clustering in groups
groupes.h1 = as.factor(cutree(h1,k=ngroupes))
# Individuals in groups
names(groupes.h1[which(groupes.h1 == 1)])
names(groupes.h1[which(groupes.h1 == 2)])
names(groupes.h1[which(groupes.h1 == 3)])
names(groupes.h1[which(groupes.h1 == 4)])
names(groupes.h1[which(groupes.h1 == 5)])

png(filename = "figures/article_one/HierarchicalClustering_Species.png", units = "cm", width = 20, height = 10,
    res = 150)
load(file = "output/classification/h1.clust.Rda")
hclustering = as.dendrogram(h1)
list_sp = labels(hclustering)
list_chrsize = NA
for (i in 1:length(list_sp)) {
  list_chrsize[i] = mean(chromosomes$phys.map.length[chromosomes$species == list_sp[i]], na.rm = TRUE)
}
color.gradient = function(x, colors = c("lightgrey", "black"), colsteps = 100) {
  return(colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
list_color = color.gradient(log(list_chrsize))
list_color = color.gradient(log(list_chrsize), colors = viridis(3))
# ggdendrogram(hclustering, rotate = FALSE, size = 2, tip.color = list_chrsize) + geom_text(colour = list_color)
library(dendextend)
par(mar=c(12,3,3,1))
dend = as.dendrogram(h1)
dend %>% dendextend::set("labels_col", value = list_color) %>% plot
dend %>% rect.dendrogram(k=6, border = "Darkgrey", lty = 5, lwd = 2)
# dend %>% rect.dendrogram(k=5, lty = 5, lwd = 0, which = c(1, 3, 5), col = rgb(0.1, 0.2, 0.4, 0.1) ) 
par(mar =  c(5.1, 4.1, 4.1, 2.1))
dev.off()















#############################################################################
#============================================================================
#       ANALYSES OF RECOMBINATION MAPS (BINS)
#============================================================================
#############################################################################


#============================================================================
# Loading all maps (1000 bins) in a concatenated dataset
#============================================================================
# Individuals are maps (species)
# Variables are 1000 bins of recombination rate along the genome

# Get the list of recombination maps in the directory
path = "output/recombination_maps/loess/bins/"
list = system(paste("ls ", path, sep = ""), intern = TRUE)
# Concatenate all single recombination maps files in a single 'data' data frame
# Make paths
list = paste(path, list, sep = "")
# Read all chromosomes at once
data = lapply(list, function(x){tmp = read.table(file = x,header = TRUE,  sep = "\t")
  tmp = c(gsub("output/recombination_maps/loess/bins/", "",gsub(".txt", "", x)), tmp$rec.rate)
  return(tmp)
  })
data = do.call("rbind", data)
# data = as.data.frame(data)
rownames(data) = as.character(data[,1])
data = data[,-1]
# rm(list)
# NA transformed to 0 values
# 0 values does not count in PCA (centered-reduced)
data[is.na(data)] = 0
data = apply(data, c(1,2), as.numeric)
is.numeric(data[1])

#----------------------------------------------------------------------------
# Trimming
# PCA shows an extreme outlier with Glycine_max_SoyBase_chromosomeB1
# so, remove it from dataset
data = data[!(rownames(data) == "Glycine_max_SoyBase_chromosomeB1"),]

# Aegilops is not in the barplot
data = data[!(rownames(data) == "Aegilops_speltoides_Zhang2019_chromosome2S"),]
# Some Populus chromosomes have errors

# Populus_deltoides_Mousavi2016_chromosome10
# Populus_deltoides_Mousavi2016_chromosome19
# Populus_deltoides_Mousavi2016_chromosome6
# Populus_deltoides_Mousavi2016_chromosome1
# Populus_deltoides_Mousavi2016_chromosome8
# Populus_deltoides_Mousavi2016_chromosome9
# Populus_deltoides_Mousavi2016_chromosome11
data = data[!(rownames(data) %in% c("Populus_deltoides_Mousavi2016_chromosome10",
              "Populus_deltoides_Mousavi2016_chromosome19",
              "Populus_deltoides_Mousavi2016_chromosome6",
              "Populus_deltoides_Mousavi2016_chromosome1",
              "Populus_deltoides_Mousavi2016_chromosome8",
              "Populus_deltoides_Mousavi2016_chromosome9",
              "Populus_deltoides_Mousavi2016_chromosome11")),]

# Populus_simonii_Mousavi2016_chromosome1
# Populus_simonii_Mousavi2016_chromosome11
# Populus_simonii_Mousavi2016_chromosome13
# Populus_simonii_Mousavi2016_chromosome14
# Populus_simonii_Mousavi2016_chromosome15
# Populus_simonii_Mousavi2016_chromosome2
# Populus_simonii_Mousavi2016_chromosome3
# Populus_simonii_Mousavi2016_chromosome5
# Populus_simonii_Mousavi2016_chromosome6
# Populus_simonii_Mousavi2016_chromosome9
data = data[!(rownames(data) %in% c("Populus_simonii_Mousavi2016_chromosome1",
                                    "Populus_simonii_Mousavi2016_chromosome11",
                                    "Populus_simonii_Mousavi2016_chromosome13",
                                    "Populus_simonii_Mousavi2016_chromosome14",
                                    "Populus_simonii_Mousavi2016_chromosome15",
                                    "Populus_simonii_Mousavi2016_chromosome2",
                                    "Populus_simonii_Mousavi2016_chromosome3",
                                    "Populus_simonii_Mousavi2016_chromosome5",
                                    "Populus_simonii_Mousavi2016_chromosome6",
                                    "Populus_simonii_Mousavi2016_chromosome9")),]



#----------------------------------------------------------------------------
# Pool chromosomes per species
# Mean of all bins for a species

# Get the list of species
list_species = unique(gsub("_chromosome[A-Za-z0-9]*", "", row.names(data)))
for (i in 1:length(list_species)) {
  list_species[i] = as.character(metadata$species[which(metadata$id == list_species[i])])
}
list_species = unique(list_species)
list_species

# Compute the matrix of binned recombination rates
data.pooled = matrix(NA, length(list_species), 999)
row.names(data.pooled) = list_species
for (i in 1:nrow(data.pooled)) {
  for (j in 1:ncol(data.pooled)) {
    data.pooled[i,j] = mean(data[grep(row.names(data.pooled)[i], row.names(data)),j], na.rm = TRUE)
  }
}

# Now do classification on pooled chromosomes

#============================================================================
# PCA
#============================================================================
# PCA uses reduced-centered variables: same dimension allows the comparison of maps of different lengths
#réalisation ACP sur données normalisées
acp = dudi.pca(data.pooled, scannf = FALSE, nf = 3)
s.label(acp$li)


#informations de chaque axe de l'ACP
acp$eig
barplot(acp$eig)#information contenue dans les axes
# Si on garde les deux premiers axes, on garde 49% de la variance : 
sum(acp$eig[1:2])/sum(acp$eig)
sum(acp$eig[1:3])/sum(acp$eig)

# Pour afficher toutes les informations en même temps
# scatter(acp)

#============================================================================
# Hierarchical Classification on recombination maps
#============================================================================
# Same method as clustering with DAPC in 'adegenet'

### Clustering based on PCA individuals
d1 = dist(acp$li)
h1 = hclust(d1, method="ward.D")
plot(h1)
ggdendrogram(h1, rotate = FALSE, size = 2)
# order of species with the classification
h1$order
# Save the results of the hierarchical classification
save(h1, file = "output/classification/hclustPooledChromosomes.Rda")
load(file = "output/classification/hclustPooledChromosomes.Rda")

# We observed at least four distinct groups
ngroupes = 4
rect.hclust(h1,k=ngroupes)
# Clustering in four groups
groupes.h1 <- as.factor(cutree(h1,k=ngroupes))
# Individuals in groups
names(groupes.h1[which(groupes.h1 == 1)])
names(groupes.h1[which(groupes.h1 == 2)])
names(groupes.h1[which(groupes.h1 == 3)])
names(groupes.h1[which(groupes.h1 == 4)])

# # Species clustered together
# unique(gsub("_chromosome[A-Za-z0-9]*", "", names(groupes.h1[which(groupes.h1 == 1)])))
# unique(gsub("_chromosome[A-Za-z0-9]*", "", names(groupes.h1[which(groupes.h1 == 2)])))
# unique(gsub("_chromosome[A-Za-z0-9]*", "", names(groupes.h1[which(groupes.h1 == 3)])))
# unique(gsub("_chromosome[A-Za-z0-9]*", "", names(groupes.h1[which(groupes.h1 == 4)])))
# 
# 
# # Some species are in multiple groups
# # Percentage of assignation to each group for each species
# cluster = groupes.h1
# names(cluster) = gsub("_chromosome[A-Za-z0-9]+", "", names(cluster))
# table(names(cluster))
# table(names(cluster[which(cluster == 1)]))
# table(names(cluster[which(cluster == 2)]))
# table(names(cluster[which(cluster == 3)]))
# table(names(cluster[which(cluster == 4)]))
# 
# prop = data.frame(species = unique(names(cluster)), cluster1 = numeric(length(unique(names(cluster)))),
#                   cluster2 = numeric(length(unique(names(cluster)))), cluster3 = numeric(length(unique(names(cluster)))),
#                   cluster4 = numeric(length(unique(names(cluster)))))
# for (i in 1:nrow(prop)) {
#   prop$cluster1[i] = sum(names(cluster[which(cluster == 1)]) == prop$species[i])/sum(names(cluster) == prop$species[i])
#   prop$cluster2[i] = sum(names(cluster[which(cluster == 2)]) == prop$species[i])/sum(names(cluster) == prop$species[i])
#   prop$cluster3[i] = sum(names(cluster[which(cluster == 3)]) == prop$species[i])/sum(names(cluster) == prop$species[i])
#   prop$cluster4[i] = sum(names(cluster[which(cluster == 4)]) == prop$species[i])/sum(names(cluster) == prop$species[i])
# }
# 
# df1 = data.frame(species = prop$species, cluster = rep(1, nrow(prop)), prop = prop$cluster1)
# df2 = data.frame(species = prop$species, cluster = rep(2, nrow(prop)), prop = prop$cluster2)
# df3 = data.frame(species = prop$species, cluster = rep(3, nrow(prop)), prop = prop$cluster3)
# df4 = data.frame(species = prop$species, cluster = rep(4, nrow(prop)), prop = prop$cluster4)
# df = rbind(df1, df2, df3, df4)
# # Identify the major cluster (cluster with the max of individuals of a given species)
# # More than 50% of proportion
# df$major = NA
# for (i in 1:nrow(df)) {
#   if (df$prop[i] > 0.5) {df$major[i] = df$cluster[i]} else {df$major[i] = NA}
# }
# # Order species by clusters
# # Levels ordered by clusters
# df$species = factor(df$species, levels = c(levels(df$species)[which(table(df$species, df$major)[,1] == 1)],
#                                            levels(df$species)[which(table(df$species, df$major)[,2] == 1)],
#                                            levels(df$species)[which(table(df$species, df$major)[,3] == 1)],
#                                            levels(df$species)[which(table(df$species, df$major)[,4] == 1)]))
# 
# # Stacked + percent
# ggplot(df, aes(fill = cluster, y = prop, x = species)) + 
#   geom_bar(position="fill", stat="identity") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.text.x=element_text(angle=-90),
#         legend.position = "none")
# 
# 
# #### CAH ####
# # Classifier les cartes
# # Classifier les espèces (summary statistics)
# # Heatmap des deux classifications
# donnees1=read.table("data1.dat",header=T)
# head(donnees1)
# dim(donnees1)
# # 23 patients et 33409 sondes
# 
# # ACP sur le tableau de données
# acp1=dudi.pca(donnees1, scale=T, scannf = FALSE, nf = 2)
# 
# #affichage des patients
# s.label(acp1$co)
# #on distingue 2 groupes
# 
# typecancer=read.table('typecancer.dat',header=T)
# s.class(acp1$co, typecancer$type, add.plot = TRUE, col = rainbow(2))
# #les deux groupes qu'on voit correspondent aux types de cancers
# 
# #jeu de données réduit
# donnees1_red=read.table("data1_red.dat",header=T)
# head(donnees1_red)
# 
# #classification sur les sondes
# d1=dist(donnees1_red)
# h1=hclust(d1,method="ward.D")
# plot(h1)
# 
# #classification sur les patients
# d2=dist(t(donnees1_red))
# h2=hclust(d2,method="ward.D")
# plot(h2)
# 
# #visualisation
# H1 = as.dendrogram(h1)
# H2 = as.dendrogram(h2)
# 
# # Sondes en lignes, individus en colonnes
# gplots::heatmap.2(as.matrix(donnees1_red),Rowv = H1,Colv = H2,trace = "none")



#############################################################################
#============================================================================
#       ANALYSES OF BROCKEN STICKS
#============================================================================
#############################################################################


#============================================================================
# Loading all maps in a concatenated dataset
#============================================================================
# Individuals are maps (species)
# Variables are 10 segments of the brocken stick model

data = read.table(paste(wd, "/tables/brockenstick/brockenstick_k10.txt", sep = ""), header = TRUE, sep = "\t")

#----------------------------------------------------------------------------
# Pool chromosomes per species
# Mean of all bins for a species

# Get the list of species
list_species = unique(gsub("_chromosome[A-Za-z0-9]*$", "", data$set))
for (i in 1:length(list_species)) {
  list_species[i] = as.character(metadata$species[which(metadata$id == list_species[i])])
}
list_species = unique(list_species)
list_species

# Compute the matrix of binned recombination rates
data.pooled = matrix(NA, length(list_species), 10)
row.names(data.pooled) = list_species
for (i in 1:nrow(data.pooled)) {
  for (j in 1:ncol(data.pooled)) {
    data.pooled[i,j] = mean(data[grep(row.names(data.pooled)[i], data$set),j+2], na.rm = TRUE)
  }
}

data.pooled[data.pooled == "NaN"] = 0
# data.pooled = data.pooled[-which(row.names(data.pooled) == c("Panicum_hallii", "Triticum_urartu")),]
  
# Now do classification on pooled chromosomes


#============================================================================
# PCA
#============================================================================
# PCA uses reduced-centered variables: same dimension allows the comparison of maps of different lengths
#réalisation ACP sur données normalisées
acp = dudi.pca(data.pooled, scannf = FALSE, nf = 3)
s.label(acp$li)

#informations de chaque axe de l'ACP
acp$eig
barplot(acp$eig)#information contenue dans les axes
# Si on garde les deux premiers axes, on garde 49% de la variance : 
sum(acp$eig[1:2])/sum(acp$eig)
sum(acp$eig[1:3])/sum(acp$eig)

# Pour afficher toutes les informations en même temps
# scatter(acp)

#============================================================================
# Hierarchical Classification
#============================================================================
# Same method as clustering with DAPC in 'adegenet'

### Clustering based on PCA individuals
d1 = dist(acp$li)
# d1 = as.matrix(dist(data.pooled))
h1 = hclust(d1, method = "ward.D")
# pvclust for bootstrapping
# h1 = clust(d1, method.hclust = "ward.D")
plot(h1)

# order of species with the classification
h1$order
# Save the results of the hierarchical classification
save(h1, file = "output/classification/hclustPooledChromosomes.Rda")
load(file = "output/classification/hclustPooledChromosomes.Rda")

# We observed at least three distinct groups
ngroupes = 4
rect.hclust(h1, k = ngroupes)
# Clustering in three groups
groupes.h1 = as.factor(cutree(h1, k = ngroupes))
# Individuals in groups
names(groupes.h1[which(groupes.h1 == 1)])
names(groupes.h1[which(groupes.h1 == 2)])
names(groupes.h1[which(groupes.h1 == 3)])
names(groupes.h1[which(groupes.h1 == 4)])

# Species clustered together
unique(gsub("_[A-Za-z0-9]*$", "", names(groupes.h1[which(groupes.h1 == 1)])))
unique(gsub("_[A-Za-z0-9]*$", "", names(groupes.h1[which(groupes.h1 == 2)])))
unique(gsub("_[A-Za-z0-9]*$", "", names(groupes.h1[which(groupes.h1 == 3)])))
unique(gsub("_[A-Za-z0-9]*$", "", names(groupes.h1[which(groupes.h1 == 4)])))

# Some species are in multiple groups
# Percentage of assignation to each group for each species
cluster = groupes.h1
names(cluster) = gsub("_[A-Za-z0-9]+$", "", names(cluster))
table(names(cluster))
table(names(cluster[which(cluster == 1)]))
table(names(cluster[which(cluster == 2)]))
table(names(cluster[which(cluster == 3)]))
table(names(cluster[which(cluster == 4)]))

prop = data.frame(species = unique(names(cluster)), cluster1 = numeric(length(unique(names(cluster)))),
                  cluster2 = numeric(length(unique(names(cluster)))), cluster3 = numeric(length(unique(names(cluster)))),
                  cluster4 = numeric(length(unique(names(cluster)))))
for (i in 1:nrow(prop)) {
  prop$cluster1[i] = sum(names(cluster[which(cluster == 1)]) == prop$species[i])/sum(names(cluster) == prop$species[i])
  prop$cluster2[i] = sum(names(cluster[which(cluster == 2)]) == prop$species[i])/sum(names(cluster) == prop$species[i])
  prop$cluster3[i] = sum(names(cluster[which(cluster == 3)]) == prop$species[i])/sum(names(cluster) == prop$species[i])
  prop$cluster4[i] = sum(names(cluster[which(cluster == 4)]) == prop$species[i])/sum(names(cluster) == prop$species[i])
}

df1 = data.frame(species = prop$species, cluster = rep(1, nrow(prop)), prop = prop$cluster1)
df2 = data.frame(species = prop$species, cluster = rep(2, nrow(prop)), prop = prop$cluster2)
df3 = data.frame(species = prop$species, cluster = rep(3, nrow(prop)), prop = prop$cluster3)
df4 = data.frame(species = prop$species, cluster = rep(4, nrow(prop)), prop = prop$cluster4)
df = rbind(df1, df2, df3, df4)
# Identify the major cluster (cluster with the max of individuals of a given species)
# More than 50% of proportion
df$major = NA
for (i in 1:nrow(df)) {
  if (df$prop[i] > 0.5) {df$major[i] = df$cluster[i]} else {df$major[i] = NA}
}
# Order species by clusters
# Levels ordered by clusters
df$species = factor(df$species, levels = c(levels(df$species)[which(table(df$species, df$major)[,1] == 1)],
                                           levels(df$species)[which(table(df$species, df$major)[,2] == 1)],
                                           levels(df$species)[which(table(df$species, df$major)[,3] == 1)],
                                           levels(df$species)[which(table(df$species, df$major)[,4] == 1)]))

# Stacked + percent
ggplot(df, aes(fill = cluster, y = prop, x = species)) + 
  geom_bar(position="fill", stat="identity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle=-90),
        legend.position = "none")


#============================================================================
# END
#============================================================================
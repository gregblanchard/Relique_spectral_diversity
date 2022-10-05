
####################################################################################
####################################################################################
####################################################################################
############## Validation of spectral diversity from field plots  ############## 
####################################################################################
####################################################################################
####################################################################################

####################################################################################
#### load files and data ####
####################################################################################
#### spatial data forest ####
library(rgdal)
library(sf)
library(rgeos)
# forest 
forest <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.shp")
zone <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/zone_plaine_des_lacs/zone.shp")
# forest in the zone
CP <- as(extent(zone), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(forest))
forest_zone <-  gIntersection(forest, CP, byid=TRUE)
forest_zone$id_patch <- 1:length(forest_zone)

#### # CHOOSE and get spectral diversity indices from plots ####

# 100 spectral sp. from PC : 1,6,7,8 (imported from server)
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL"
Biodiv_Indicators <- readRDS(paste0(file_path,"/VALIDATIONBiodiv_Indicators.rds"))
# this is the plot file of reference (from which the spectral indices are extracted)
Path_Vector <- paste0(file_path, "/plot_data_ok")

# 20 spectral sp. from PC : 1,2,6,7,8 (generated on local)
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL"
Biodiv_Indicators <- readRDS(paste0(file_path,"/VALIDATIONBiodiv_Indicators.rds"))
# this is the plot file of reference (from which the spectral indices are extracted)
Path_Vector <- paste0(file_path, "/plot_data_ok")

#### get objects and data for field plots ####

# get field diversity indices from plots 
library(rgdal)
# plots adjusted to match with sentinel reolution (3*3 cells)
plot_data_ok <- readOGR(paste0(Path_Vector, "/geom_plots_PDL_div_indices_utm_ok.shp"))
plot_data_true <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data/geom_plots_PDL_div_indices_utm.shp")
# keep same plot set than for S2 images                
plot_data_true <- plot_data_true[plot_data_true$localty %in% plot_data_ok$localty,]

# get data table from field plots 
plots_PDL_div_indices <- readRDS('/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_PDL_div_indices.rds')

# get data distance matrices from field plots 
plot_data_matrices <- readRDS(file = '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plot_data_matrices.rds')

# filter matrices with selected plots (during the biodivmapR pipeline)
plots_PDL_div_indices <- plots_PDL_div_indices[plots_PDL_div_indices$locality %in% plot_data_ok$localty,]
library(usedist)
plotok <- rownames(plot_data_matrices$community_matrix)[rownames(plot_data_matrices$community_matrix) %in% plot_data_ok$localty]
plot_data_matrices$community_matrix_plotok <- plot_data_matrices$community_matrix[plotok,]
plot_data_matrices$dist_mat_taxo_plotok <- dist_subset(plot_data_matrices$dist_mat_taxo, plotok)
plot_data_matrices$dist_mat_taxo_BA_plotok <- dist_subset(plot_data_matrices$dist_mat_taxo_BA, plotok)

plot_data_matrices$dist_mat_functio_plotok <- dist_subset(plot_data_matrices$dist_mat_functio, plotok)
plot_data_matrices$dist_mat_functio_trans_BA_plotok <- dist_subset(plot_data_matrices$dist_mat_functio_trans_BA, plotok)
plot_data_matrices$dist_mat_functio_trans_synth_BA_plotok <- dist_subset(plot_data_matrices$dist_mat_functio_trans_synth_BA, plotok)

dist_mat_taxo <- plot_data_matrices$dist_mat_taxo_plotok
dist_mat_taxo_BA <- plot_data_matrices$dist_mat_taxo_BA_plotok

dist_mat_functio <- plot_data_matrices$dist_mat_functio_plotok
dist_mat_functio_trans_BA <- plot_data_matrices$dist_mat_functio_trans_BA_plotok
dist_mat_functio_trans_synth_BA <- plot_data_matrices$dist_mat_functio_trans_synth_BA_plotok

#### make sure that we get same order of plots for all data ####
# reference data from field 
plot_order <- plots_PDL_div_indices$locality
# spatial data used for extracting spectral data
plot_data_ok <- plot_data_ok[match(plot_order, plot_data_ok$localty),]
plot_data_true <- plot_data_true[match(plot_order, plot_data_true$localty),]
# biodiv indicators from spectral data
Biodiv_Indicators$Richness <- data.frame(Biodiv_Indicators$Richness[match(plot_order, Biodiv_Indicators$localty),])
Biodiv_Indicators$Fisher <- data.frame(Biodiv_Indicators$Fisher[match(plot_order, Biodiv_Indicators$localty),])
Biodiv_Indicators$Shannon <- data.frame(Biodiv_Indicators$Shannon[match(plot_order, Biodiv_Indicators$localty),])
Biodiv_Indicators$Simpson <- data.frame(Biodiv_Indicators$Simpson[match(plot_order, Biodiv_Indicators$localty),])
Biodiv_Indicators$fisher.All <- Biodiv_Indicators$fisher.All[match(plot_order, Biodiv_Indicators$localty),]
Biodiv_Indicators$Shannon.All <- Biodiv_Indicators$Shannon.All[match(plot_order, Biodiv_Indicators$localty),]
Biodiv_Indicators$Simpson.All <- Biodiv_Indicators$Simpson.All[match(plot_order, Biodiv_Indicators$localty),]
Biodiv_Indicators$BCdiss <- Biodiv_Indicators$BCdiss[match(plot_order, Biodiv_Indicators$localty),match(plot_order, Biodiv_Indicators$localty)]
for (i in 1:length(Biodiv_Indicators$BCdiss.All))
  Biodiv_Indicators$BCdiss.All[[i]] <- as.matrix(Biodiv_Indicators$BCdiss.All[[i]])[match(plot_order, Biodiv_Indicators$localty),match(plot_order, Biodiv_Indicators$localty)]

Biodiv_Indicators$localty <- Biodiv_Indicators$localty[match(plot_order, Biodiv_Indicators$localty)]
#### make geographic groups from plots localities ####
plot_grps <- plots_PDL_div_indices$localty_grp
# merge the only plot from Kuebini (donwside) with Wadjana
Biodiv_Indicators$localty_grp <- plot_grps
plot_data_ok$localty_grp <- plot_grps
plot_data_true$localty_grp <- plot_grps

#### extract diversity indices ####
BC_spectral <- Biodiv_Indicators$BCdiss
# beta div : apply ordination using PCoA (same as done for map_beta_div)
library(labdsv)
dist_mat_spectral <- as.dist(BC_spectral, diag = FALSE, upper = FALSE)
names(dist_mat_spectral) <- Biodiv_Indicators$localty
spectral_betadiv <- labdsv::pco(dist_mat_spectral, k = 3)
spectral_betadiv_ok <- data.frame(spectral_betadiv$points)
names(spectral_betadiv_ok) <- c("PC1_spectral", "PC2_spectral", "PC3_spectral")

# alpha div
Shannon_RS <- c(Biodiv_Indicators$Shannon)[[1]]
FRic <- c(Biodiv_Indicators$FunctionalDiversity$FRic)
FEve <- c(Biodiv_Indicators$FunctionalDiversity$FEve)
FDiv <- c(Biodiv_Indicators$FunctionalDiversity$FDiv)
# spectral_div_indices <- Shannon_RS, FRic, 
#### arrange data  ####
taxo_betadiv <- data.frame(PC1_taxo = plots_PDL_div_indices$PCoA_PC1_taxo,
                           PC2_taxo = plots_PDL_div_indices$PCoA_PC2_taxo,
                           PC3_taxo = plots_PDL_div_indices$PCoA_PC3_taxo)

taxo_betadiv_BA <- data.frame(PC1_taxo_BA = plots_PDL_div_indices$PCoA_PC1_taxo_BA,
                           PC2_taxo_BA = plots_PDL_div_indices$PCoA_PC2_taxo_BA,
                           PC3_taxo_BA = plots_PDL_div_indices$PCoA_PC3_taxo_BA)

functio_betadiv <- data.frame(PC1_functio = plots_PDL_div_indices$PCoA_PC1_functio,
                              PC2_functio = plots_PDL_div_indices$PCoA_PC2_functio,
                              PC3_functio = plots_PDL_div_indices$PCoA_PC3_functio)

functio_betadiv_trans_BA <- data.frame(PC1_functio_trans_BA = plots_PDL_div_indices$PCoA_PC1_functio_trans_BA,
                              PC2_functio_trans_BA = plots_PDL_div_indices$PCoA_PC2_functio_trans_BA,
                              PC3_functio_trans_BA = plots_PDL_div_indices$PCoA_PC3_functio_trans_BA)


###########################################################################
#####               plot maps             #####  
###########################################################################
library(RColorBrewer)
library(viridis)     
library(dichromat)     
library(rasterVis)     
library(colorspace)
library(viridisLite)
library(colorRamps)
library(raster)
library(rgeos)
library(sf)

# rgb image
aerial_photo <- brick("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/rapport/figures/map_grandsud.tif")
# transform multiband image to raster

aerial_photo_ok_r <- raster(aerial_photo)
cols_bg_image <- factor(rgb(aerial_photo[], maxColorValue=255))
Sys.sleep(5)
aerial_photo_ok_r[] <- cols_bg_image
levelplot(aerial_photo_ok_r, col.regions=as.character(levels(cols_bg_image)), colorkey=FALSE)

# get same CRS for plots 
plot_data_ok
plot_data_ok_wgs <- spTransform(plot_data_ok, CRS("+proj=longlat +datum=WGS84 +no_defs"))

spdf <- SpatialPointsDataFrame(coords = cbind(plot_data_ok$longitd, plot_data_ok$latitud), data = data.frame(plot_data_ok_wgs),
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

#### plot ####
bg_image <- levelplot(aerial_photo_ok_r, col.regions=as.character(levels(cols_bg_image)), colorkey=FALSE)
forest_layer <- latticeExtra::layer(sp.polygons(forest , lwd=2, fill = "forestgreen", alpha = 0.3)) 
plots_KB <- latticeExtra::layer(sp.points(spdf[spdf$localty_grp == "Kuebini",], pch= 15, cex=1, col= "dodgerblue1", alpha = .8)) 
plots_CORIFOR <- latticeExtra::layer(sp.points(spdf[spdf$localty_grp == "CORIFOR",], pch= 15, cex=1, col= "coral2", alpha = .8)) 
plots_ForetNord <- latticeExtra::layer(sp.points(spdf[spdf$localty_grp == "Forêt Nord",], pch= 15, cex=1, col= "darkolivegreen3", alpha = .8)) 
plots_Wadjana <- latticeExtra::layer(sp.points(spdf[spdf$localty_grp == "Wadjana",], pch= 15, cex=1, col= "magenta", alpha = .8)) 

plots_map <- bg_image + forest_layer + plots_KB + plots_CORIFOR  + plots_ForetNord + plots_Wadjana 
plots_map

####################################################################################
#### Based on spectral beta div: produce figures in order to locate the different types of vegetation in the PCoA space ####
####################################################################################
# assign a type of vegetation to each plot, assuming that the type of vegetation 
# is defined by the name of the shapefile
library(rgdal)
library(tools)
library(ggplot2)
library(gridExtra)
library(grid)
nbSamples <- plotName <- c()
for (i in 1:length(Path_Vector)){
  shp <- Path_Vector[i]
  nbSamples[i] <- length(rgdal::readOGR(shp,verbose = FALSE))
  plotName[i] <- "PDL"
}

Type_Vegetation = c()
for (i in 1: length(nbSamples)){
  for (j in 1:nbSamples[i]){
    Type_Vegetation = c(Type_Vegetation,Biodiv_Indicators$localty_grp[j])
  }
}

# create data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('vgtype'=Type_Vegetation,'pco1'= spectral_betadiv_ok[,1],'pco2'= spectral_betadiv_ok[,2],'pco3' = spectral_betadiv_ok[,3],
                      'shannon'=Shannon_RS,'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)

# plot field data in the PCoA space, with size corresponding to shannon index
g1 <-ggplot (Results, aes (x=pco1, y=pco2, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

g2 <-ggplot (Results, aes (x=pco1, y=pco3, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

g3 <-ggplot (Results, aes (x=pco2, y=pco3, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

#extract legend
get_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(g3)
gAll <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                                 g2 + theme(legend.position="none"),
                                 g3 + theme(legend.position="none"),
                                 nrow=1),legend,nrow=1, widths =c(5, 1))

####################################################################################
#### Based on taxo beta div: produce figures in order to locate the different types of vegetation in the PCoA space ####
####################################################################################
# create data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('vgtype'=Type_Vegetation,'pco1'= taxo_betadiv_BA[,1],'pco2'= taxo_betadiv_BA[,2],'pco3' = taxo_betadiv_BA[,3],
                      'shannon'=Shannon_RS,'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)

# plot field data in the PCoA space, with size corresponding to shannon index
g1 <-ggplot (Results, aes (x=pco1, y=pco2, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

g2 <-ggplot (Results, aes (x=pco1, y=pco3, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

g3 <-ggplot (Results, aes (x=pco2, y=pco3, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

#extract legend
get_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(g3)
gAll <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                                 g2 + theme(legend.position="none"),
                                 g3 + theme(legend.position="none"),
                                 nrow=1),legend,nrow=1, widths =c(5, 1))

####################################################################################
#### Based on functio beta div: produce figures in order to locate the different types of vegetation in the PCoA space ####
####################################################################################
# create data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('vgtype'=Type_Vegetation,'pco1'= functio_betadiv_trans_BA[,1],'pco2'= functio_betadiv_trans_BA[,2],'pco3' = functio_betadiv_trans_BA[,3],
                      'shannon'=Shannon_RS,'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)

# plot field data in the PCoA space, with size corresponding to shannon index
g1 <-ggplot (Results, aes (x=pco1, y=pco2, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

g2 <-ggplot (Results, aes (x=pco1, y=pco3, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

g3 <-ggplot (Results, aes (x=pco2, y=pco3, color=vgtype,size=shannon)) + 
  geom_point(alpha=0.6) 

#extract legend
get_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(g3)
gAll <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                                 g2 + theme(legend.position="none"),
                                 g3 + theme(legend.position="none"),
                                 nrow=1),legend,nrow=1, widths =c(5, 1))


####################################################################################
#### beta div: compare field data with spectral data ####
####################################################################################

##### CANONICAL CORRELATION ANALYSIS : CCA package #####
# see https://stats.oarc.ucla.edu/r/dae/canonical-correlation-analysis/#:~:text=Canonical%20correlation%20analysis%20is%20used,are%20multiple%20intercorrelated%20outcome%20variables.
library(CCA)
library(GGally)


cc_tax_spec_beta <- cc(taxo_betadiv_BA, spectral_betadiv_ok)
summary(cc_tax_spec_beta)
# display the canonical correlations
cc_tax_spec_beta$cor 
barplot(cc_tax_spec_beta$cor, xlab = "Dimension")
# tests of canonical dimensions
rho <- cc_tax_spec_beta$cor
## Define taxo_betadiv of observations, number of variables in first set, and number of variables in the second set.
n <- dim(taxo_betadiv)[1]
p <- dim(taxo_betadiv)[2]
q <- dim(spectral_betadiv_ok)[2]

## Calculate p-values using the F-approximations of different test statistics:
library(CCP)
p.asym(rho, n, p, q, tstat = "Wilks")                          
p.perm(taxo_betadiv, spectral_betadiv_ok, nboot = 999, rhostart = 1, type = "Wilks")                        
p.asym(rho, n, p, q, tstat = "Roy")

plt.cc(cc_tax_spec_beta, var.label = TRUE)

##### CANONICAL CORRELATION ANALYSIS : ade4 package #####
library(ade4)
dudipco_spectral <- dudi.pco(quasieuclid(dist_mat_spectral), scannf = FALSE, nf = 3)
dudipco_taxo <- dudi.pco(quasieuclid(dist_mat_taxo_BA), scannf = FALSE, nf = 3)

plot(dudipco_taxo$tab$A1 ~ dudipco_spectral$tab$A1)
plot(dudipco_taxo$tab$A2 ~ dudipco_spectral$tab$A2)
plot(dudipco_taxo$tab$A3 ~ dudipco_spectral$tab$A3)

coinert_spect_tax <- coinertia(dudipco_spectral, dudipco_taxo, scannf = FALSE, nf = 2)

summary(coinert_spect_tax)
plot(coinert_spect_tax)

round(coinert_spect_tax$eig/sum(coinert_spect_tax$eig)*100, digits = 1)[1:10]

rv1 <- RV.rtest(dudipco_spectral$tab, dudipco_taxo$tab, 999)
plot(rv1)

#### linear model ####

summary(lm(spectral_betadiv_ok$PC1_spectral ~ taxo_betadiv$PC1_taxo + taxo_betadiv$PC2_taxo + taxo_betadiv$PC3_taxo))

summary(lm(taxo_betadiv$PC1_taxo ~ spectral_betadiv_ok$PC1_spectral + spectral_betadiv_ok$PC2_spectral + spectral_betadiv_ok$PC3_spectral))
summary(lm(taxo_betadiv$PC2_taxo ~ spectral_betadiv_ok$PC1_spectral + spectral_betadiv_ok$PC2_spectral + spectral_betadiv_ok$PC3_spectral))
summary(lm(taxo_betadiv$PC3_taxo ~ spectral_betadiv_ok$PC1_spectral + spectral_betadiv_ok$PC2_spectral + spectral_betadiv_ok$PC3_spectral))

summary(lm(taxo_betadiv$PC1_taxo~ spectral_betadiv_ok$PC1_spectral))
plot(taxo_betadiv$PC1_taxo ~ spectral_betadiv_ok$PC1_spectral)

summary(lm(functio_betadiv$PC1_functio ~ spectral_betadiv_ok$PC1_spectral + spectral_betadiv_ok$PC2_spectral + spectral_betadiv_ok$PC3_spectral))
summary(lm(functio_betadiv$PC2_functio ~ spectral_betadiv_ok$PC1_spectral + spectral_betadiv_ok$PC2_spectral + spectral_betadiv_ok$PC3_spectral))
summary(lm(functio_betadiv$PC3_functio ~ spectral_betadiv_ok$PC1_spectral + spectral_betadiv_ok$PC2_spectral + spectral_betadiv_ok$PC3_spectral))

plot(functio_betadiv$PC1_functio ~ spectral_betadiv_ok$PC1_spectral)

##### Multiple Regression on distance Matrices ##### 
library(MDMR)
dist_mat_taxo_test <- plot_data_matrices$dist_mat_taxo_BA_plotok
dist_mat_functio_test <- plot_data_matrices$dist_mat_functio_trans_synth_BA_plotok
mdmr1 <- mdmr(spectral_betadiv_ok, dist_mat_taxo_test)
mdmr1
summary(mdmr1)

mdmr2 <- mdmr( taxo_betadiv_BA, dist_mat_spectral)
mdmr2
summary(mdmr2)

mdmr3 <- mdmr( functio_betadiv_trans_BA, dist_mat_spectral)
mdmr3
summary(mdmr3)

t_f_beta_div <- cbind(taxo_betadiv_BA, functio_betadiv_trans_BA)
names(t_f_beta_div) <- c("PC1t","PC2t","PC3t","PC1f","PC2f","PC3f")
mdmr2 <- mdmr(t_f_beta_div, dist_mat_spectral)
summary(mdmr2)

##### Mantel test ##### 
library(ade4)
mantel.rtest(dist_mat_taxo, dist_mat_functio, nrepet = 999)
plot(dist_mat_taxo ~ dist_mat_functio)

mantel.rtest(dist_mat_taxo_BA, dist_mat_functio_trans_synth_BA, nrepet = 999)
plot(dist_mat_taxo_BA ~ dist_mat_functio_trans_synth_BA)

mantel.rtest(dist_mat_taxo, dist_mat_spectral, nrepet = 999)
plot(dist_mat_spectral ~ dist_mat_taxo)

mantel.rtest(dist_mat_taxo_BA, dist_mat_spectral, nrepet = 999)
plot(dist_mat_spectral ~ dist_mat_taxo_BA)

mantel.rtest(dist_mat_functio, dist_mat_spectral, nrepet = 999)
plot(dist_mat_spectral ~ dist_mat_functio)

mantel.rtest(dist_mat_functio_trans_synth_BA, dist_mat_spectral, nrepet = 999)
plot(dist_mat_spectral ~ dist_mat_functio_trans_synth_BA)

library(ecodist)
mantel(dist_mat_taxo, dist_mat_functio)
mantel(dist_mat_taxo_BA, dist_mat_spectral)
mantel(dist_mat_spectral, dist_mat_functio_trans_synth_BA)

# mantel residuals
library(daee)
mantres <- mantel.residuals(dist_mat_spectral, dist_mat_functio_trans_synth_BA )
mantres$statistic
mantres$residuals
mantel(mantres$residuals, dist_mat_taxo)

####                  geographic distance decay                        ####
# spectral distance matrix to pairs of distance 
dist_mat_spectral_mat <- as.matrix(dist_mat_spectral)
dist_mat_spectral_mat[upper.tri(dist_mat_spectral_mat)] <- NA
df_dist_spectral <- reshape2::melt(as.matrix(dist_mat_spectral_mat), varnames = c("row", "col"))
df_dist_spectral <- df_dist_spectral[!is.na(df_dist_spectral$value) ,]

# taxo distance matrix to pairs of distance 
dist_mat_taxo_mat <- as.matrix(dist_mat_taxo)
dist_mat_taxo_mat[upper.tri(dist_mat_taxo_mat)] <- NA
df_dist_taxo <- reshape2::melt(as.matrix(dist_mat_taxo_mat), varnames = c("row", "col"))
df_dist_taxo <- df_dist_taxo[!is.na(df_dist_taxo$value) ,]

# functio distance matrix to pairs of distance 
dist_mat_functio_mat <- as.matrix(dist_mat_functio)
dist_mat_functio_mat[upper.tri(dist_mat_functio_mat)] <- NA
df_dist_functio <- reshape2::melt(as.matrix(dist_mat_functio_mat), varnames = c("row", "col"))
df_dist_functio <- df_dist_functio[!is.na(df_dist_functio$value) ,]

# goegraphic distance matrix
library(sf)
dist_geo <- st_distance(st_as_sf(plot_data_ok))
dist_mat_geo <- as.matrix(dist_geo)
dist_mat_geo[upper.tri(dist_mat_geo)] <- NA
rownames(dist_mat_geo) <- colnames(dist_mat_geo) <- plot_data_ok$localty
units(dist_mat_geo) <- NULL
df_dist <- reshape2::melt(dist_mat_geo, varnames = c("row", "col"))
df_dist <- df_dist[!is.na(df_dist$value) ,]

#### merge selected dataframes of pairs ###
sel_mat <- df_dist_taxo
sel_mat <- df_dist_functio
sel_mat <- df_dist_spectral

df_dist_bc <- cbind(df_dist, sel_mat$value)
names(df_dist_bc) <- c("row", "col", "m", "bc")
df_dist_bc <- df_dist_bc[complete.cases(df_dist_bc),]
# no zero distance (same plot)
df_dist_bc <- df_dist_bc[df_dist_bc$m>0,]
####  linear model #### 
plot(df_dist_bc$bc ~ df_dist_bc$m)
lm1 <- lm(bc ~ m, data = df_dist_bc)
smry1 <- summary(lm1)
smry1
plot(df_dist_spectral$value ~ df_dist_taxo$value)
plot(df_dist_spectral$value ~ df_dist_functio$value)

####  mantel test #### 
mantel(dist_mat_spectral, dist_geo)

####################################################################################
#### alpha div: compare field data with spectral data ####
####################################################################################

summary(lm(Shannon_RS ~ plot_data_ok$sp_shnn))
summary(lm(Shannon_RS ~ plot_data_ok$rr20_txsh))
summary(lm(Shannon_RS ~ plot_data_ok$FRic))
summary(lm(Shannon_RS ~ plot_data_ok$FDiv))
summary(lm(FDiv ~ plot_data_ok$FDiv))
summary(lm(FRic ~ plot_data_ok$FRic))

###########################################################################
###                                                                     ###
####           Influence of fragmentation on spectral diversity        ####
###                                                                     ###
###########################################################################
####  distance to edge  #### 

# use metric distance 
forest_zone_utm <- spTransform(forest_zone, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# ID patch for plots
patchs_over_point <- over(plot_data_ok,forest_zone_utm) 
plot_data_ok$id_patch <- patchs_over_point$id_patch
# get distance to edge for the centroid of each plot
dist_edge <- c()
for (i in 1:nrow(plot_data_ok)) {
  pts_tmp <- plot_data_ok[i,]
  ## NA if points are out of forest (reference stations "cagnar")
  if(!is.na(pts_tmp$id_patch)){
    patch_tmp <- forest_zone_utm[forest_zone_utm$id_patch == pts_tmp$id_patch,]
    forest_edge <- as(patch_tmp, "SpatialLines") 
    dist_edge <- c(dist_edge, gDistance(rgeos::gCentroid(pts_tmp,byid=TRUE), forest_edge, byid=TRUE)) 
  }else{
    dist_edge <- c(dist_edge, NA) 
  }
}
plot_data_ok$dist_edge <- dist_edge
hist(dist_edge)
####   ####   ####   ####   ####   ####   ####   #### 
#### spatral diversity VS. distance to edge  #### 
####   ####   ####   ####   ####   ####   ####   #### 
##### Linear model ##### 
summary(lm(spectral_betadiv_ok$PC1_spectral ~ log(plot_data_ok$dist_edge)))
plot(spectral_betadiv_ok$PC1_spectral ~ plot_data_ok$dist_edge)
plot(spectral_betadiv_ok$PC1_spectral ~ log(plot_data_ok$dist_edge))

summary(lm(Shannon_RS ~ log(plot_data_ok$dist_edge)))
plot(Shannon_RS ~ plot_data_ok$dist_edge)
summary(lm(FRic ~ log(plot_data_ok$dist_edge)))
plot(FRic ~ plot_data_ok$dist_edge)
summary(lm(FEve ~ log(plot_data_ok$dist_edge)))
plot(FEve ~ plot_data_ok$dist_edge)
summary(lm(FDiv ~ log(plot_data_ok$dist_edge)))
plot(FDiv ~ plot_data_ok$dist_edge)

##############################################################################
###                                                                        ###
#### Relationship between biological components and spectral components   ####
###                                                                        ###
##############################################################################
library(exactextractr)
library(raster)
#### get raw spectral cbands ####
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021"
raw_spectral <- stack(paste0(file_path,"/T58KFA_raw_int_forest_zone_PDL.tif"))
# rename bands 
names(raw_spectral) <- c("B2",     "B3",     "B4",     "B5",     "B6",     "B7",     "B8",    "B8A",    "B11",    "B12")
#### extract PCA values as functional traits for plots #### 
## use spatial plots with true area ##
# Community mean 
plot_CWM_raw_spectral <-  exact_extract(raw_spectral, plot_data_true, "mean")
colnames(plot_CWM_raw_spectral) <- paste0("S2.",colnames(plot_CWM_raw_spectral))
# Community variance
plot_CWV_raw_spectral <-  exact_extract(raw_spectral, plot_data_true, "variance")
colnames(plot_CWV_raw_spectral) <- paste0("S2.",colnames(plot_CWV_raw_spectral))

#### get PCA spectral  ####
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL"
pca_spectral <- stack(paste0(file_path,"/SPCA/PCA/OutputPCA_8_PCs"))
#### extract PCA values as functional traits for plots #### 
## use spatial plots with true area ##
# Community mean 
plot_CWM_PCA_spectral <-  exact_extract(pca_spectral, plot_data_true, "mean")
colnames(plot_CWM_PCA_spectral) <- paste0("PCAS2.",colnames(plot_CWM_PCA_spectral))
# Community variance
plot_CWV_PCA_spectral <-  exact_extract(pca_spectral, plot_data_true, "variance")
colnames(plot_CWV_PCA_spectral) <- paste0("PCAS2.",colnames(plot_CWV_PCA_spectral))

#### compil data ####
spectral_data_plots <- cbind(Shannon_RS_S2 = Shannon_RS,
                             FRic_S2 = FRic,
                             FEve_S2 = FEve,
                             FDiv_S2 = FDiv,
                             betadiv_PC1_S2 = spectral_betadiv_ok$PC1,
                             betadiv_PC2_S2 = spectral_betadiv_ok$PC2,
                             betadiv_PC3_S2 = spectral_betadiv_ok$PC3,
                             plot_CWM_raw_spectral, 
                             plot_CWV_raw_spectral,
                             plot_CWM_PCA_spectral,
                             plot_CWV_PCA_spectral)
 
plots_field_spectral <- cbind(plots_PDL_div_indices, spectral_data_plots)

#### relationship CWM traits ~ CWM PCA Axis for plots #### 
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_BA_SLA)
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_BA_LA)
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_BA_LDMC)
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_BA_WD)
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_synth_BA_Axis1)

#### multivariate model selection  ####


library(MuMIn)
library(parallel)
#### select variables ####
# choose set of explanatory variable
spectral_div_indices_alpha <- c("Shannon_RS_S2", "FRic_S2", "FEve_S2", "FDiv_S2")
spectral_div_indices_beta <-c("betadiv_PC1_S2", "betadiv_PC2_S2", "betadiv_PC3_S2")

expl_v_list <-  list(c(spectral_div_indices_alpha))
expl_v_list <-  list(c(spectral_div_indices_beta))

expl_v_list <- list(colnames(plot_CWM_raw_spectral))
expl_v_list <- list(colnames(plot_CWV_raw_spectral))
expl_v_list <- list(colnames(plot_CWM_PCA_spectral))
expl_v_list <- list(colnames(plot_CWV_PCA_spectral))

# set of response variables
resp_v_list <- colnames(plots_PDL_div_indices)[8:ncol(plots_PDL_div_indices)]

# total set of variable
data_for_mod_all <- plots_field_spectral
#### model selection ####
list_model_tab <- list()
list_models_best <- list()
list_first_best_model <- list()
list_model_avg <- list()
list_model_full <- list()
for(i in 1:length(resp_v_list)){
  data_for_mod <- plots_field_spectral
    print(paste("model selection for model", i, "on", length(resp_v_list)))
    resp_tmp <- resp_v_list[i]
    expl_tmp <- expl_v_list[[1]]
    
    #### all models based on plot data exept for microclimate ? ####
    # if(!resp_tmp == "max_VPD"){data_for_mod <- data_for_mod[!is.na(data_for_mod$sp_richness),]}
    data_for_mod <- data_for_mod_all
    #### select variables ####
    data_for_mod <- data_for_mod[,c(resp_tmp, expl_tmp)]
    # no NA values 
    complet_cases <- complete.cases(data_for_mod) 
    data_for_mod <- data_for_mod[complet_cases,]
    
    #### final number of points ####
    dim(data_for_mod)
    #### standardize variables ?  ####
    data_for_mod_toscale <- data_for_mod
    data_for_mod_scaled <- data_for_mod
    data_for_mod_scaled <- as.data.frame(apply(data_for_mod_toscale, 2, scale))
    
    #### model selection for each variable ####
    # with gaussiian family distribution 
    # do not keep Variance Inflation Factors > 5
    glm1 <- glm( formula(paste0(resp_tmp,  " ~ .")), data = data_for_mod_scaled, na.action = "na.fail")
    
    list_model_full[[i]] <- glm1
    
    vif <- car::vif(glm1)
    data_for_mod_scaled_ok <- data_for_mod_scaled
    while(max(vif)>10){
      data_for_mod_scaled_ok <- data_for_mod_scaled_ok[, !colnames(data_for_mod_scaled_ok) %in% names(vif)[vif == max(vif)]]
      glm1 <- glm(formula(paste0(resp_tmp,  " ~ .")), data = data_for_mod_scaled_ok, na.action = "na.fail")
      vif <- car::vif(glm1)
    }
    vif
    glm1 <- glm(formula(paste0(resp_tmp,  " ~ .")), data = data_for_mod_scaled_ok, na.action = "na.fail")
    
    rm(dd)
    # cl <- makeCluster(3)
    # clusterExport(cl, "data_for_mod_scaled_ok")
    dd <- pdredge(glm1,rank = "AICc", beta = "partial.sd",
                  extra = list(
                    "R^2", "*" = function(x) {
                      s <- summary(x)
                      c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
                        F = s$fstatistic[[1]])
                    })
    )
    # stopCluster(cl)
    dd
    list_model_tab[[i]] <- dd 
    # selection few best models based on delta AIC
    rm(models.list)
    models.list <- get.models(dd,subset =  delta < 1.5) 
    print(paste(i, length(models.list)))
    list_models_best[[i]] <- models.list
    list_first_best_model[[i]] <-  models.list[[1]]
    # average models and coeficiants based on weigth
    if(length(models.list) > 1){   muminavg <- model.avg(models.list, revised.var = TRUE, fit = TRUE, beta = "partial.sd")
    list_model_avg[[i]] <- muminavg
    }else{list_model_avg[[i]] <- NA}
}

names(list_model_tab) =  names(list_models_best) =  names(list_model_avg) = names(list_first_best_model) <- resp_v_list
# list_model_tab
# list_first_best_model

# saveRDS(list_first_best_model,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_first_best_model_div_indices_S2.rds") 

# saveRDS(list_first_best_model,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_first_best_model_CWM_S2.rds") 
# saveRDS(list_first_best_model,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_first_best_model_CVM_S2.rds") 
# saveRDS(list_first_best_model,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_first_best_model_CWM_PCAS2.rds") 

# noice! 
plot(plots_field_spectral$CWM_SLA~plots_field_spectral$betadiv_PC1_S2)

##### export results of best models #####
coef_best <- c()
for(i in 1: length(list_first_best_model)){
  coef_tmp <- rep(NA, length( expl_v_list[[1]]))
  names(coef_tmp) <- expl_v_list[[1]]
  sm_tmp <- summary(list_first_best_model[[i]])
  if(nrow(sm_tmp$coefficients)>1){
    for(j in 2: nrow(sm_tmp$coefficients)){
      coef_tmp[names(coef_tmp) == rownames(sm_tmp$coefficients)[j]] <- 
        paste0(round(sm_tmp$coefficients[j,1], digits = 2))
    }
  }
  # add R2
  R2 <- round(list_model_tab[[i]]$`R^2`[1], digits = 2)
  coef_tmp <- c(coef_tmp, R2 = R2)
  coef_tmp <- as.numeric(coef_tmp)
  names(coef_tmp) <- c(expl_v_list[[1]],"R2")
  coef_best <- rbind(coef_best, coef_tmp)
  rownames(coef_best)[i] <- names(list_first_best_model)[i]
}
coef_best <- data.frame(coef_best)
# explore results 
coef_best[rev(order(coef_best$R2)),]

#### partial least square regression ####
library(pls)
# arrange data
CWM_S2 <- as.matrix(plot_CWM_raw_spectral)
CWM_PCAS2 <- as.matrix(plot_CWM_PCA_spectral)
CWM_field <- as.matrix(cbind(WD = plots_field_spectral$CWM_trans_BA_WD,
                                  SLA = plots_field_spectral$CWM_trans_BA_SLA,
                                  LA = plots_field_spectral$CWM_trans_BA_LA,
                                  LDMC = plots_field_spectral$CWM_trans_BA_LDMC))

CWM_synth_field <- as.matrix(cbind(Ax1 = plots_field_spectral$CWM_trans_synth_BA_Axis1,
                                   Ax2 = plots_field_spectral$CWM_trans_synth_BA_Axis2))
#### transform data? ####
# CWM_field[,"SLA"] <- log(CWM_field[,"SLA"])
# raw spectral bands vs. field data
model <- plsr(CWM_field ~ CWM_S2, data=plots_field_spectral, scale=TRUE, validation="CV")
model <- plsr(CWM_field[,"SLA"] ~ CWM_S2, data=plots_field_spectral, scale=TRUE, validation="CV")
model <- plsr(CWM_synth_field ~ CWM_S2, data=plots_field_spectral, scale=TRUE, validation="CV")
model <- plsr(CWM_synth_field[,"Ax1"] ~ CWM_S2, data=plots_field_spectral, scale=TRUE, validation="CV")

# PCA spectral bands vs. field data
model <- plsr(CWM_synth_field[,"Ax1"] ~ CWM_PCAS2, data=plots_field_spectral, scale=TRUE, validation="CV")


summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

plot(RMSEP(model), legendpos = "topright")

#### subset plots? ####
plots_field_spectral$localty_grp
subset <- plots_field_spectral$localty_grp %in% c("CORIFOR", "Kuebini")
subset <- plots_field_spectral$localty_grp %in% c("Kuebini")

CWM_field[subset,]

model <- plsr(CWM_synth_field[subset,] ~ CWM_S2[subset,], data=plots_field_spectral[subset,], scale=TRUE, validation="CV")
summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")










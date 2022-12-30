
###########################################################################
###                                                                     ###
####                             Load data                             ####
###                                                                     ###
###########################################################################
library(rgdal)
library(sf)
library(dplyr)

data_ncpipn_taxon_sf <- readRDS( "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/data_ncpipn_taxon_sf_all.rds")
plots_field_spectral_ok_env <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_field_spectral_ok_env.rds")

# get dataframe
data_ncpipn_taxon_df_PDL <- data_ncpipn_taxon_sf_PDL %>% st_drop_geometry()

###########################################################################
###                                                                     ###
####         UPDATE PLOT SELECTION (same dataset as for PLSR)          ####
###                                                                     ###
###########################################################################

data_ncpipn_taxon_df_PDL_update <- data_ncpipn_taxon_df_PDL[data_ncpipn_taxon_df_PDL$locality %in% plots_field_spectral_ok_env$locality,]
data_ncpipn_plots_sf_PDL_update <- unique(data_ncpipn_taxon_df_PDL_update[,c("longitude", "latitude", "source","locality", "plot")])

community_matrix_update <- as.data.frame.matrix(t(table(data.frame(data_ncpipn_taxon_df_PDL_update[,c("nom_taxon", "locality")]))))

#### get same order than goemetry object ####
community_matrix_update <- community_matrix_update[data_ncpipn_plots_sf_PDL_update$locality,]

community_matrix_update <- community_matrix_update[, colSums(community_matrix_update)>0]

#### community matrix with relative crown area ####
# from the log-linear relationship in Blanchard et al. 2016
# y = b * x^a  <=> log(y) = log(b) + a * log(x)
data_ncpipn_taxon_df_PDL$CA <- (0.169 * (data_ncpipn_taxon_df_PDL$dbh^1.354 ) )
# OR : data_ncpipn_taxon_df_PDL$CA <- exp(log(0.169) + 1.354 * log(data_ncpipn_taxon_df_PDL$dbh) ) # (same)
community_matrix_CA_update <- community_matrix_update
community_matrix_CA_update[ ] <- NA
for(i in 1:length(rownames(community_matrix_update))){
  data_loc_tmp <-  data_ncpipn_taxon_df_PDL[data_ncpipn_taxon_df_PDL$locality == rownames(community_matrix_update)[i],]
  for(j in 1:length(colnames(community_matrix_update))){
    tx_tmp <- data_loc_tmp[data_loc_tmp$nom_taxon == colnames(community_matrix_update)[j],]
    # sum individual crown area for each species
    rel_CA <- sum(tx_tmp$CA, na.rm = T)
    community_matrix_CA_update[i,j] <- rel_CA
  }
}
# round CA to interger
community_matrix_CA_update <- round(community_matrix_CA_update)
# verif order plots
cbind(rownames(community_matrix_CA_update),data_ncpipn_plots_sf_PDL_update$locality)


###########################################################################
####                         taxo  betadiv                             ####        
###########################################################################

#### beta div taxo based on a similar method as in BiodivmapR ####
library(ecodist)
library(vegan)
# distance matrix
dist_mat_taxo <- vegdist(community_matrix_update, method = "bray")
dist_mat_taxo <- vegdist(community_matrix_update, method = "bray")

# verif order plots
cbind(names(dist_mat_taxo),data_ncpipn_plots_sf_PDL_update$locality)

# PCoA
pcoa_taxo <- pco(dist_mat_taxo, negvals = "zero", dround = 0)

# proportion of variance explained 
plot(pcoa_taxo$values/sum(pcoa_taxo$values))
(pcoa_taxo$values/sum(pcoa_taxo$values))[1]

# keep the 3 first PC
PCoA_PCs_taxo_update <- pcoa_taxo$vectors[,1:3]

# saveRDS(PCoA_PCs_taxo_update, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data/PCoA_PCs_taxo_update.rds")



#### SAME USING relative CANOPY AREA : beta div taxo based on a similar method as in BiodivmapR ####
library(ecodist)
library(vegan)
# distance matrix
# dist_mat_taxo <- vegdist(community_matrix_CA_update, method = "bray")
community_matrix_CA_update_relat <- community_matrix_CA_update/rowSums(community_matrix_CA_update)     
dist_mat_taxo_CA <- vegdist(community_matrix_CA_update_relat, method = "bray")
# dist_mat_taxo_CA <- vegdist(community_matrix_CA_update, method = "bray")

# verif order plots
cbind(names(dist_mat_taxo_CA),data_ncpipn_plots_sf_PDL_update$locality)

# PCoA
pcoa_taxo_CA <- pco(dist_mat_taxo_CA, negvals = "zero", dround = 0)

# proportion of variance explained 
plot(pcoa_taxo_CA$values/sum(pcoa_taxo_CA$values))
(pcoa_taxo_CA$values/sum(pcoa_taxo_CA$values))[1]
sum((pcoa_taxo_CA$values/sum(pcoa_taxo_CA$values))[1:3])
percent_var_pc <-  (pcoa_taxo_CA$values/sum(pcoa_taxo_CA$values))[1:3]
percent_var_pc <- round(percent_var_pc*100, digits = 1)
# keep the 3 first PC
PCoA_PCs_taxo_CA_update <- pcoa_taxo_CA$vectors[,1:3]

# saveRDS(PCoA_PCs_taxo_CA_update, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data/PCoA_PCs_taxo_CA_update.rds")

# using the ape package
pcoa_t <- ape::pcoa(dist_mat_taxo_CA)

plot(pcoa_t$values$Relative_eig)
plot(pcoa_t)
####################################################################################
#### nice plot ####
####################################################################################

plot(pcoa_taxo_CA$vectors[,1], pcoa_taxo_CA$vectors[,2], type = "n", xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA plots")
text(pcoa_taxo_CA$vectors[,1], pcoa_taxo_CA$vectors[,2], labels(dist_mat_taxo_CA), 
     cex = 0.9, xpd = TRUE)
####################################################################################
#### Based on taxo beta div: produce figures in order to locate the different types of vegetation in the PCoA space ####
####################################################################################
# change name for soil: Iron crust vs. eroded ferritic soil (see mcCoy 1999)
plots_field_spectral_ok_env$geology <- ifelse(plots_field_spectral_ok_env$geology == "Péridotites", "Eroded ferritic soils", "Iron crust" )

# create data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('pco1'= PCoA_PCs_taxo_CA_update[,1],'pco2'= PCoA_PCs_taxo_CA_update[,2],'pco3' = PCoA_PCs_taxo_CA_update[,3],
                      'Shannon'=plots_field_spectral_ok_env$sp_shannon, "Soil" = plots_field_spectral_ok_env$geology2 )

# plot field data in the PCoA space, with size corresponding to shannon index
my_colors <- c( "springgreen4", "orangered4")

g1 <- ggplot (Results, aes (x=pco1, y=pco2, color= Soil)) + #, size=Shanno
  geom_point(alpha=0.8, size = 1.8)+ 
  scale_color_manual(values = my_colors) +
  xlab(paste0("Dim1 (", percent_var_pc[1], "%)")) +
  ylab(paste0("Dim2 (", percent_var_pc[2], "%)")) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title=element_blank()) 

g1

# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/paper/PCoA_taxo.png", width=130, height=100, units = 'mm', res = 300) 
g1
dev.off()

####################################################################################
#### NMDS ####
####################################################################################
nmds <- metaMDS(community_matrix_CA_update, k=3,trymax=100)
plot(nmds)
stressplot(nmds)

# saveRDS(nmds, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data/nmds_taxo_CA_update.rds")



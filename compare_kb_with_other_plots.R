
###########################################################################
###                                                                     ###
####                             Load data                             ####
###                                                                     ###
###########################################################################
library(rgdal)
library(sf)
library(dplyr)

#### WITH ORIGINAL DATASET ####
data_ncpipn_taxon_sf <- readRDS( "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/data_ncpipn_taxon_sf_all.rds")

colnames(data_ncpipn_taxon_sf)

data_ncpipn_taxon_sf$geo_pt
data_ncpipn_taxon_sf$geo_pt..89
st_geometry(data_ncpipn_taxon_sf) <- "geo_pt..89"
plot(st_geometry(data_ncpipn_taxon_sf))

# keep only geometry of individuals 
data_ncpipn_taxon_sf <- data_ncpipn_taxon_sf[,!(colnames(data_ncpipn_taxon_sf) == "geo_pt")]
# get dataframe
data_ncpipn_taxon_df <- data_ncpipn_taxon_sf %>% st_drop_geometry()

# export simple shape file for NCPIPPN
plot_geo_ncpipn <- unique(data.frame(locality = data_ncpipn_taxon_sf$locality, data_ncpipn_taxon_sf$geo_pt..89))
st_write(plot_geo_ncpipn, paste0("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/", "plot_geo_ncpipn.shp"))

###########################################################################
###                                                                     ###
####                          filter data                              ####
###                                                                     ###
###########################################################################
library(rgdal)
library(sf)
#### only plaine des lacs ####
zone <- readOGR("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/GIS/forest_map/zone_ncpippn.shp")
zone <- st_as_sf(zone)

data_ncpipn_taxon_sf_PDL <- st_intersection(data_ncpipn_taxon_sf, zone)
# saveRDS(data_ncpipn_taxon_sf_PDL, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/data_ncpipn_taxon_sf_PDL.rds")
# data_ncpipn_taxon_sf_PDL <- readRDS( "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/data_ncpipn_taxon_sf_PDL.rds")

#### only NCPIPPN ####
data_ncpipn_taxon_sf_PDL <- data_ncpipn_taxon_sf_PDL %>% filter(source == "nc-pippn")

#### only dbh >= 10 cm ####
colnames(data_ncpipn_taxon_sf_PDL)
data_ncpipn_taxon_sf_PDL <- data_ncpipn_taxon_sf_PDL %>% filter(dbh >= 10)
data_ncpipn_taxon_sf_PDL$strate
#### only living trees ####
data_ncpipn_taxon_sf_PDL <- data_ncpipn_taxon_sf_PDL %>% filter(isliving == T)

#### only living trees ####
data_ncpipn_taxon_sf_PDL <- data_ncpipn_taxon_sf_PDL %>% filter(isliving == T)

#### ALL COUNT FOR THE PAPER: total number of plots, of species, genus, family ####

data_ncpipn_taxon_sf_PDL_tree_count <- data_ncpipn_taxon_sf_PDL[data_ncpipn_taxon_sf_PDL$is_tree == T,]
dim(data_ncpipn_taxon_sf_PDL_tree_count)

# total number of plot in the landscape (before removing plots)

data_ncpipn_taxon_sf_PDL_tree_count <- data_ncpipn_taxon_sf_PDL_tree_count[data_ncpipn_taxon_sf_PDL_tree_count$plot %in% c("circle_R10" ,"circle_R11.28"  ,"plot_20x20") ,]
length(unique(data_ncpipn_taxon_sf_PDL_tree_count$locality))

####  removed plots ####
# these plots are too much outside forest (~open area = 1/3 of the plot) 
plot_remove <- c("10-10-BIS","Forêt Nord 35","Forêt Nord 18")
data_ncpipn_taxon_sf_PDL_tree_count <- data_ncpipn_taxon_sf_PDL_tree_count[!data_ncpipn_taxon_sf_PDL_tree_count$locality %in% plot_remove  ,]

# load the final plot list to get the same number of plots (from script plsr_uncertainties.R)
plots_field_spectral_ok_final <-readRDS( '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_field_spectral_ok.rds')
data_ncpipn_taxon_sf_PDL_tree_count <- data_ncpipn_taxon_sf_PDL_tree_count[data_ncpipn_taxon_sf_PDL_tree_count$locality %in% plots_field_spectral_ok_final$locality,]
# nb of individuals and species in total and of plots 
dim(data_ncpipn_taxon_sf_PDL_tree_count)
length(unique(data_ncpipn_taxon_sf_PDL_tree_count$locality))
length(unique(data_ncpipn_taxon_sf_PDL_tree_count$nom_taxon))
# number of genus and families
length(unique(data_ncpipn_taxon_sf_PDL_tree_count$genus))
length(unique(data_ncpipn_taxon_sf_PDL_tree_count$family))

# nb of individual and species in plots
sort(table(data_ncpipn_taxon_sf_PDL_tree_count$locality))
mean(sort(table(data_ncpipn_taxon_sf_PDL_tree_count$locality)))
data_ncpipn_taxon_sf_PDL_tree_count_species <- unique(data_ncpipn_taxon_sf_PDL_tree_count[c("locality", "plot", "date_det", "nom_taxon")])
sort(table(data_ncpipn_taxon_sf_PDL_tree_count_species$locality))
mean(sort(table(data_ncpipn_taxon_sf_PDL_tree_count_species$locality)))
# type of plot
data_ncpipn_taxon_sf_PDL_tree_count_plot <- unique(data_ncpipn_taxon_sf_PDL_tree_count[c("locality", "plot", "date_det")])
table(data_ncpipn_taxon_sf_PDL_tree_count_plot$plot)
table(data_ncpipn_taxon_sf_PDL_tree_count_plot$date_det)




###########################################################################
###                                                                     ###
####            create geometry for plot area                          ####
###                                                                     ###
###########################################################################

# get only plot data
names(data_ncpipn_taxon_sf_PDL)
data_ncpipn_plots_sf_PDL <- unique(data_ncpipn_taxon_sf_PDL[,c("longitude", "latitude", "source","locality", "plot","geo_pt..89")])
data_ncpipn_plots_sf_PDL_utm <- st_transform(data_ncpipn_plots_sf_PDL, "+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# plot type
plot_type <-  unique(data_ncpipn_plots_sf_PDL_utm$plot)
# create buffer around each plot
geom_plots <- c()
for (i in 1:length(plot_type)){
  plot_type_tmp <- plot_type[i]
  data_ncpipn_plots_sf_PDL_utm_tmp <- data_ncpipn_plots_sf_PDL_utm %>% filter(data_ncpipn_plots_sf_PDL_utm$plot == plot_type_tmp)
  # create buffers around plots
  if(plot_type_tmp == "circle_R10")
    geom_plots_tmp <- st_buffer(data_ncpipn_plots_sf_PDL_utm_tmp, dist = 10)
  if(plot_type_tmp == "circle_R11.28")
    geom_plots_tmp <- st_buffer(data_ncpipn_plots_sf_PDL_utm_tmp, dist = 11.28)
  if(plot_type_tmp == "plot_100x100")
    geom_plots_tmp <- st_buffer(data_ncpipn_plots_sf_PDL_utm_tmp, dist = 50, endCapStyle = "SQUARE")
  if(plot_type_tmp == "plot_20x20")
    geom_plots_tmp <- st_buffer(data_ncpipn_plots_sf_PDL_utm_tmp, dist = 10, endCapStyle = "SQUARE")
  if(plot_type_tmp == "plot_70x70")
    geom_plots_tmp <- st_buffer(data_ncpipn_plots_sf_PDL_utm_tmp, dist = 35, endCapStyle = "SQUARE")
  if(plot_type_tmp == "plot_50x50")
    geom_plots_tmp <- st_buffer(data_ncpipn_plots_sf_PDL_utm_tmp, dist = 25, endCapStyle = "SQUARE")
  # merge resulting buffers
  geom_plots <- rbind(geom_plots, geom_plots_tmp)
}
plot(geom_plots$geo_pt..89)

#### keep same order of plot ####
geom_plots <- geom_plots[match(data_ncpipn_plots_sf_PDL_utm$locality, geom_plots$locality),]

###########################################################################
###                                                                     ###
####                          community matrix                         ####
###                                                                     ###
###########################################################################

data_ncpipn_taxon_sf <- readRDS( "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/data_ncpipn_taxon_sf_all.rds")

# get dataframe
data_ncpipn_taxon_df_PDL <- data_ncpipn_taxon_sf_PDL %>% st_drop_geometry()

community_matrix <- as.data.frame.matrix(t(table(data.frame(data_ncpipn_taxon_df_PDL[,c("nom_taxon", "locality")]))))

#### get same order than goemetry object ####
community_matrix <- community_matrix[data_ncpipn_plots_sf_PDL_utm$locality,]

rownames(community_matrix)
cbind(rownames(community_matrix),data_ncpipn_plots_sf_PDL_utm$locality)
colnames(community_matrix)

dim(community_matrix)
colSums(community_matrix)
rowSums(community_matrix)

#### community matrix with relative basal area ####
data_ncpipn_taxon_df_PDL$BA <- 2*pi*data_ncpipn_taxon_df_PDL$dbh
community_matrix_BA <- community_matrix
community_matrix_BA[ ] <- NA
for(i in 1:length(rownames(community_matrix))){
  data_loc_tmp <-  data_ncpipn_taxon_df_PDL[data_ncpipn_taxon_df_PDL$locality == rownames(community_matrix)[i],]
  for(j in 1:length(colnames(community_matrix))){
    tx_tmp <- data_loc_tmp[data_loc_tmp$nom_taxon == colnames(community_matrix)[j],]
    rel_BA <- sum(tx_tmp$BA, na.rm = T)
    community_matrix_BA[i,j] <- rel_BA
  }
}
# round BA to interger
community_matrix_BA <- round(community_matrix_BA)
# verif order plots
cbind(rownames(community_matrix_BA),data_ncpipn_plots_sf_PDL_utm$locality)

#### community matrix with relative crown area ####
# from the log-linear relationship in Blanchard et al. 2016
# y = b * x^a  <=> log(y) = log(b) + a * log(x)
data_ncpipn_taxon_df_PDL$CA <- (0.169 * (data_ncpipn_taxon_df_PDL$dbh^1.354 ) )
# OR : data_ncpipn_taxon_df_PDL$CA <- exp(log(0.169) + 1.354 * log(data_ncpipn_taxon_df_PDL$dbh) ) # (same)
community_matrix_CA <- community_matrix
community_matrix_CA[ ] <- NA
for(i in 1:length(rownames(community_matrix))){
  data_loc_tmp <-  data_ncpipn_taxon_df_PDL[data_ncpipn_taxon_df_PDL$locality == rownames(community_matrix)[i],]
  for(j in 1:length(colnames(community_matrix))){
    tx_tmp <- data_loc_tmp[data_loc_tmp$nom_taxon == colnames(community_matrix)[j],]
    # sum individual crown area for each species
    rel_CA <- sum(tx_tmp$CA, na.rm = T)
    community_matrix_CA[i,j] <- rel_CA
  }
}
# round CA to interger
community_matrix_CA <- round(community_matrix_CA)
# verif order plots
cbind(rownames(community_matrix_CA),data_ncpipn_plots_sf_PDL_utm$locality)

#### export png ####
# how many plot per species?
hist(colSums(community_matrix>0), 
     breaks = seq(min(colSums(community_matrix>0)), max(colSums(community_matrix>0)), length.out = 100),
     main = "",
     xlab = "Number of plots",
     ylab = "Number of species")

table(colSums(community_matrix))
table(colSums(community_matrix>0))

###########################################################################
###                                                                     ###
####        compute diversity indices for field data from plots        ####
###                                                                     ###
###########################################################################

colnames(data_ncpipn_taxon_sf_PDL)
###########################################################################
####         taxo div (R, shanon, betadiv) ####        
###########################################################################
#### hill numbers ####
library(hillR)
taxo_div <- data.frame(sp_richness = hill_taxa(community_matrix,  q = 0), 
                       sp_shannon = hill_taxa(community_matrix,  q = 1),
                       sp_simpson = hill_taxa(community_matrix,  q = 2))

taxo_div
#### rarefied hill numbers ####
library(iNEXT)
sample_size_stand_2 <- 20
raref <- iNEXT(t(community_matrix),  size = c(30,sample_size_stand_2), q = c(0,1,2) )
# verif order plots
cbind(names(raref$iNextEst),data_ncpipn_plots_sf_PDL_utm$locality)

rar_div2 <- c()
for(i in 1:length(raref$iNextEst)){
  plot_tmp <- raref$iNextEst[[i]]
  plot_tmp <- unique(plot_tmp)
  div_tmp <- plot_tmp[plot_tmp$m == sample_size_stand_2,]$qD
  rar_div2 <- rbind(rar_div2, div_tmp)
}
rar_div2 <- data.frame(rar_div2)
colnames(rar_div2) <- c(paste0("rar",sample_size_stand_2,"_taxrich"),
                        paste0("rar",sample_size_stand_2,"_taxsha"),
                        paste0("rar",sample_size_stand_2,"_taxsim"))

#### beta div taxo based on a similar method as in BiodivmapR ####
library(ecodist)
library(vegan)
# distance matrix
dist_mat_taxo <- vegdist(community_matrix, method = "bray")
# verif order plots
cbind(names(dist_mat_taxo),data_ncpipn_plots_sf_PDL_utm$locality)

# PCoA
pcoa_taxo <- pco(dist_mat_taxo, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_taxo$values/sum(pcoa_taxo$values))
sum((pcoa_taxo$values/sum(pcoa_taxo$values))[1:3])
PCoA_PCs_taxo <- pcoa_taxo$vectors[,1:3]

#### SAME WITH BA: beta div taxo based on a similar method as in BiodivmapR ####
# distance matrix
dist_mat_taxo_BA <- vegdist(community_matrix_BA, method = "bray")
# verif order plots
cbind(names(dist_mat_taxo_BA),data_ncpipn_plots_sf_PDL_utm$locality)

# PCoA
pcoa_taxo_BA <- pco(dist_mat_taxo_BA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_taxo_BA$values/sum(pcoa_taxo_BA$values))
sum((pcoa_taxo_BA$values/sum(pcoa_taxo_BA$values))[1:3])
PCoA_PCs_taxo_BA <- pcoa_taxo_BA$vectors[,1:3]


#### SAME WITH CA: beta div taxo based on a similar method as in BiodivmapR ####
# distance matrix
dist_mat_taxo_CA <- vegdist(community_matrix_CA, method = "bray")
# verif order plots
cbind(names(dist_mat_taxo_CA),data_ncpipn_plots_sf_PDL_utm$locality)

# PCoA
pcoa_taxo_CA <- pco(dist_mat_taxo_CA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_taxo_CA$values/sum(pcoa_taxo_CA$values))
sum((pcoa_taxo_CA$values/sum(pcoa_taxo_CA$values))[1:3])
PCoA_PCs_taxo_CA <- pcoa_taxo_CA$vectors[,1:3]

###########################################################################
#### Functional diversity ####
###########################################################################
# FT_all_compil <- read.csv( '/home/thesardfou/Documents/data/FT_all_122021.csv')
FT_all_compil <- readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/data_taxa.rds")
names(FT_all_compil)
FT_all_compil_ok <- FT_all_compil[order(FT_all_compil$taxaname),]
list_sp <- unique(data.frame(nom_taxon = data_ncpipn_taxon_sf_PDL$nom_taxon, taxaname = data_ncpipn_taxon_sf_PDL$taxaname))
sum(!list_sp$nom_taxon %in% FT_all_compil_ok$taxaname)
tax_not_in_db <- list_sp[!list_sp$taxaname %in% FT_all_compil_ok$taxonref,]
tax_not_in_db
dim(tax_not_in_db)

#### get trait data at the species level ####
names(FT_all_compil_ok)
FT_all_compil_ok_traits <- FT_all_compil_ok[,c("family","genus","taxonref","taxaname",
                                               "wood_density_avg","leaf_sla_avg", "leaf_area_avg", "leaf_ldmc_avg")]

#### complete with traits at the genus level (if available) ####
colnames(FT_all_compil_ok_traits)
col_traits <- c("wood_density_avg", "leaf_sla_avg", "leaf_area_avg", "leaf_ldmc_avg")

FT_all_compil_ok_traits <- FT_all_compil_ok_traits[!is.na(FT_all_compil_ok_traits$genus),]
FT_all_traits_complet_genus <- FT_all_compil_ok_traits

for(i in 1:length(col_traits)){
  trait_tmp <- col_traits[i]
  for(j in 1:nrow(FT_all_traits_complet_genus)){
    if(is.na(FT_all_traits_complet_genus[j,trait_tmp])){
      FT_all_traits_complet_genus[j,trait_tmp] <- 
        FT_all_traits_complet_genus[FT_all_traits_complet_genus["taxaname"] == FT_all_traits_complet_genus[j,"genus"],trait_tmp]
    }
  }
}

####  get same number of species for trait and site-species matrices #### 
FT_select_sp <- FT_all_compil_ok_traits[,c("wood_density_avg", "leaf_sla_avg", "leaf_area_avg", "leaf_ldmc_avg")]
FT_select <- FT_all_traits_complet_genus[,c("wood_density_avg", "leaf_sla_avg", "leaf_area_avg", "leaf_ldmc_avg")]

rownames(FT_select) <- rownames(FT_select_sp) <-  FT_all_traits_complet_genus$taxaname
colnames(FT_select) <- colnames(FT_select_sp) <-  c("WD", "SLA", "LA", "LDMC" )

FT_select_sp <- FT_select_sp[complete.cases(FT_select_sp),]
FT_select <- FT_select[complete.cases(FT_select),]

community_matrix_for_traits_sp <- community_matrix[,colnames(community_matrix) %in% rownames(FT_select_sp)]
FT_select_sp <- FT_select_sp[rownames(FT_select_sp) %in% colnames(community_matrix_for_traits_sp), ]

community_matrix_for_traits <- community_matrix[,colnames(community_matrix) %in% rownames(FT_select)]
FT_select <- FT_select[rownames(FT_select) %in% colnames(community_matrix_for_traits), ]

dim(FT_select)
# verif order plots
cbind(rownames(community_matrix_for_traits),data_ncpipn_plots_sf_PDL_utm$locality)

# get trait matrice for only SLA
FT_select_SLA <- data.frame(SLA = FT_select[ , "SLA"])
rownames(FT_select_SLA) <- rownames(FT_select)

#### percentage complete trait data for plots at the species level ####
# in plots
hist(rowSums(community_matrix_for_traits_sp) / rowSums(community_matrix))
# total
sum(community_matrix_for_traits_sp)/sum(community_matrix)
# species 
ncol(community_matrix_for_traits_sp)/ncol(community_matrix)

#### percentage complete trait data for plots with trait completed at genus level ####
# in plots
hist(rowSums(community_matrix_for_traits) / rowSums(community_matrix))
# total
sum(community_matrix_for_traits)/sum(community_matrix)
# species 
ncol(community_matrix_for_traits)/ncol(community_matrix)


#### transform traits values ####
FT_select_trans <- FT_select
FT_select_trans$SLA <- log(FT_select$SLA)
FT_select_trans$LA <- log(FT_select$LA)
FT_select_SLA_trans <- FT_select_SLA
FT_select_SLA_trans$SLA <- log(FT_select_SLA$SLA)
# save 
# saveRDS(FT_select_trans, '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/FT_select_trans.rds')
# saveRDS(community_matrix_for_traits, '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/community_matrix_for_traits.rds')

#### synthetic traits from PCA on traits ####
library(ade4)
pca_traits <- dudi.pca(FT_select_trans, scannf = F, nf = 2)
library(factoextra)
fviz_pca_ind(pca_traits)
fviz_pca_var(pca_traits)
synthetic_trait <- data.frame(pca_traits$li)
#### functional diversity indices ####
#### get functional alpha div indices and get functional beta div based on a similar method as in BiodivmapR ####
library(FD)
# get functional diversity
FD_PDL <- dbFD(FT_select, community_matrix_for_traits, stand.FRic = T, calc.CWM = T)
names(FD_PDL$CWM) <- paste0("CWM_",names(FD_PDL$CWM))
FD_PDL_trans <- dbFD(FT_select_trans, community_matrix_for_traits, stand.FRic = T, calc.CWM = T)
names(FD_PDL_trans$CWM) <- paste0("CWM_trans_",names(FD_PDL_trans$CWM))

# only for SLA
FD_PDL_SLA <- dbFD(FT_select_SLA, community_matrix_for_traits , calc.FRic = F, calc.CWM = F)
FD_PDL_trans_SLA <- dbFD(FT_select_SLA_trans, community_matrix_for_traits , calc.FRic = F, calc.CWM = F)

#### get functional beta div based on a similar method as in BiodivmapR ####
# distance matrix
library(BAT)
dist_mat_functio <- beta(community_matrix_for_traits,FT_select)
dist_mat_functio <- dist_mat_functio$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio <- pco(dist_mat_functio, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio$values/sum(pcoa_functio$values))
sum((pcoa_functio$values/sum(pcoa_functio$values))[1:3])
PCoA_PCs_functio <- pcoa_functio$vectors[,1:3]

#### get functional beta div based on a similar method as in BiodivmapR with transformed traits ####
# distance matrix
library(BAT)
dist_mat_functio_trans <- beta(community_matrix_for_traits,FT_select_trans)
dist_mat_functio_trans <- dist_mat_functio_trans$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_trans <- pco(dist_mat_functio_trans, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_trans$values/sum(pcoa_functio_trans$values))
sum((pcoa_functio_trans$values/sum(pcoa_functio_trans$values))[1:3])
PCoA_PCs_functio_trans <- pcoa_functio_trans$vectors[,1:3]

## same for only SLA ##
# distance matrix
library(BAT)
dist_mat_functio_sla <- beta(community_matrix_for_traits,FT_select_SLA)
dist_mat_functio_sla <- dist_mat_functio_sla$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_sla <- pco(dist_mat_functio_sla, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_sla$values/sum(pcoa_functio_sla$values))
sum((pcoa_functio_sla$values/sum(pcoa_functio_sla$values))[1:3])
PCoA_PCs_functio_SLA <- pcoa_functio_sla$vectors[,1:3]

# distance matrix
library(BAT)
dist_mat_functio_trans_sla <- beta(community_matrix_for_traits,FT_select_SLA_trans)
dist_mat_functio_trans_sla <- dist_mat_functio_trans_sla$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_trans_sla <- pco(dist_mat_functio_trans_sla, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_trans_sla$values/sum(pcoa_functio_trans_sla$values))
sum((pcoa_functio_trans_sla$values/sum(pcoa_functio_trans_sla$values))[1:3])
PCoA_PCs_functio_trans_SLA <- pcoa_functio_trans_sla$vectors[,1:3]

#### use basal area to compute functional metrics ####
community_matrix_for_traits_BA <- 
  community_matrix_BA[,colnames(community_matrix_BA) %in% colnames(community_matrix_for_traits)]
## functional diversity indices ##
library(FD)
# get functional diversity
FD_PDL_BA <- dbFD(FT_select, community_matrix_for_traits_BA, stand.FRic = T, calc.CWM = T)
names(FD_PDL_BA$CWM) <- paste0("CWM_BA_",names(FD_PDL_BA$CWM))

FD_PDL_trans_BA <- dbFD(FT_select_trans, community_matrix_for_traits_BA, stand.FRic = T, calc.CWM = T)
names(FD_PDL_trans_BA$CWM) <- paste0("CWM_trans_BA_",names(FD_PDL_trans_BA$CWM))

FD_PDL_trans_synthetic_BA  <- dbFD(synthetic_trait, community_matrix_for_traits_BA, stand.FRic = T, calc.CWM = T)
names(FD_PDL_trans_synthetic_BA$CWM) <- paste0("CWM_trans_synth_BA_",names(FD_PDL_trans_synthetic_BA$CWM))

# only for SLA
FD_PDL_BA_SLA <- dbFD(FT_select_SLA, community_matrix_for_traits_BA , calc.FRic = F, calc.CWM = F)
FD_PDL_trans_BA_SLA <- dbFD(FT_select_SLA_trans, community_matrix_for_traits_BA , calc.FRic = F, calc.CWM = F)

## get functional beta div based on a similar method as in BiodivmapR ##
# distance matrix
library(BAT)
dist_mat_functio_BA <- beta(community_matrix_for_traits_BA,FT_select)
dist_mat_functio_BA <- dist_mat_functio_BA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_BA <- pco(dist_mat_functio_BA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_BA$values/sum(pcoa_functio_BA$values))
sum((pcoa_functio_BA$values/sum(pcoa_functio_BA$values))[1:3])
PCoA_PCs_functio_BA <- pcoa_functio_BA$vectors[,1:3]

## get functional beta div based on a similar method as in BiodivmapR with transformed traits ##
# distance matrix
library(BAT)
dist_mat_functio_trans_BA <- beta(community_matrix_for_traits_BA,FT_select_trans)
dist_mat_functio_trans_BA <- dist_mat_functio_trans_BA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_trans_BA <- pco(dist_mat_functio_trans_BA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_trans_BA$values/sum(pcoa_functio_trans_BA$values))
sum((pcoa_functio_trans_BA$values/sum(pcoa_functio_trans_BA$values))[1:3])
PCoA_PCs_functio_trans_BA <- pcoa_functio_trans_BA$vectors[,1:3]

## get functional beta div based on a similar method as in BiodivmapR with synthetic traits ##
# distance matrix
library(BAT)
dist_mat_functio_trans_synth_BA <- beta(community_matrix_for_traits_BA,synthetic_trait)
dist_mat_functio_trans_synth_BA <- dist_mat_functio_trans_synth_BA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_trans_synth_BA <- pco(dist_mat_functio_trans_synth_BA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_trans_synth_BA$values/sum(pcoa_functio_trans_synth_BA$values))
sum((pcoa_functio_trans_synth_BA$values/sum(pcoa_functio_trans_synth_BA$values))[1:3])
PCoA_PCs_functio_trans_synt_BA <- pcoa_functio_trans_synth_BA$vectors[,1:3]

## same for only SLA ##
# distance matrix
library(BAT)
dist_mat_functio_sla_BA <- beta(community_matrix_for_traits_BA,FT_select_SLA)
dist_mat_functio_sla_BA <- dist_mat_functio_sla_BA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_sla_BA <- pco(dist_mat_functio_sla_BA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_sla_BA$values/sum(pcoa_functio_sla_BA$values))
sum((pcoa_functio_sla_BA$values/sum(pcoa_functio_sla_BA$values))[1:3])
PCoA_PCs_functio_SLA_BA <- pcoa_functio_sla_BA$vectors[,1:3]

# distance matrix
library(BAT)
dist_mat_functio_trans_sla_BA <- beta(community_matrix_for_traits_BA,FT_select_SLA_trans)
dist_mat_functio_trans_sla_BA <- dist_mat_functio_trans_sla_BA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_trans_sla_BA <- pco(dist_mat_functio_trans_sla_BA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_trans_sla_BA$values/sum(pcoa_functio_trans_sla_BA$values))
sum((pcoa_functio_trans_sla_BA$values/sum(pcoa_functio_trans_sla_BA$values))[1:3])
PCoA_PCs_functio_trans_SLA_BA <- pcoa_functio_trans_sla_BA$vectors[,1:3]


#### use crown area to compute functional metrics ####
community_matrix_for_traits_CA <- 
  community_matrix_CA[,colnames(community_matrix_CA) %in% colnames(community_matrix_for_traits)]
## functional diversity indices ##
library(FD)
# get functional diversity
FD_PDL_CA <- dbFD(FT_select, community_matrix_for_traits_CA, stand.FRic = T, calc.CWM = T)
names(FD_PDL_CA$CWM) <- paste0("CWM_CA_",names(FD_PDL_CA$CWM))

FD_PDL_trans_CA <- dbFD(FT_select_trans, community_matrix_for_traits_CA, stand.FRic = T, calc.CWM = T)
names(FD_PDL_trans_CA$CWM) <- paste0("CWM_trans_CA_",names(FD_PDL_trans_CA$CWM))

FD_PDL_trans_synthetic_CA  <- dbFD(synthetic_trait, community_matrix_for_traits_CA, stand.FRic = T, calc.CWM = T)
names(FD_PDL_trans_synthetic_CA$CWM) <- paste0("CWM_trans_synth_CA_",names(FD_PDL_trans_synthetic_CA$CWM))

# only for SLA
FD_PDL_CA_SLA <- dbFD(FT_select_SLA, community_matrix_for_traits_CA , calc.FRic = F, calc.CWM = F)
FD_PDL_trans_CA_SLA <- dbFD(FT_select_SLA_trans, community_matrix_for_traits_CA , calc.FRic = F, calc.CWM = F)

## get functional beta div based on a similar method as in BiodivmapR ##
# distance matrix
library(BAT)
dist_mat_functio_CA <- beta(community_matrix_for_traits_CA,FT_select)
dist_mat_functio_CA <- dist_mat_functio_CA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_CA <- pco(dist_mat_functio_CA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_CA$values/sum(pcoa_functio_CA$values))
sum((pcoa_functio_CA$values/sum(pcoa_functio_CA$values))[1:3])
PCoA_PCs_functio_CA <- pcoa_functio_CA$vectors[,1:3]

## get functional beta div based on a similar method as in BiodivmapR with transformed traits ##
# distance matrix
library(BAT)
dist_mat_functio_trans_CA <- beta(community_matrix_for_traits_CA,FT_select_trans)
dist_mat_functio_trans_CA <- dist_mat_functio_trans_CA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_trans_CA <- pco(dist_mat_functio_trans_CA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_trans_CA$values/sum(pcoa_functio_trans_CA$values))
sum((pcoa_functio_trans_CA$values/sum(pcoa_functio_trans_CA$values))[1:3])
PCoA_PCs_functio_trans_CA <- pcoa_functio_trans_CA$vectors[,1:3]

## get functional beta div based on a similar method as in BiodivmapR with synthetic traits ##
# distance matrix
library(BAT)
dist_mat_functio_trans_synth_CA <- beta(community_matrix_for_traits_CA,synthetic_trait)
dist_mat_functio_trans_synth_CA <- dist_mat_functio_trans_synth_CA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_trans_synth_CA <- pco(dist_mat_functio_trans_synth_CA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_trans_synth_CA$values/sum(pcoa_functio_trans_synth_CA$values))
sum((pcoa_functio_trans_synth_CA$values/sum(pcoa_functio_trans_synth_CA$values))[1:3])
PCoA_PCs_functio_trans_synt_CA <- pcoa_functio_trans_synth_CA$vectors[,1:3]


## same for only SLA ##
# distance matrix
library(BAT)
dist_mat_functio_sla_CA <- beta(community_matrix_for_traits_CA,FT_select_SLA)
dist_mat_functio_sla_CA <- dist_mat_functio_sla_CA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_sla_CA <- pco(dist_mat_functio_sla_CA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_sla_CA$values/sum(pcoa_functio_sla_CA$values))
sum((pcoa_functio_sla_CA$values/sum(pcoa_functio_sla_CA$values))[1:3])
PCoA_PCs_functio_SLA_CA <- pcoa_functio_sla_CA$vectors[,1:3]

# distance matrix
library(BAT)
dist_mat_functio_trans_sla_CA <- BAT::beta(community_matrix_for_traits_CA,FT_select_SLA_trans)
dist_mat_functio_trans_sla_CA <- dist_mat_functio_trans_sla_CA$Btotal
library(ecodist)
# distance matrix
# PCoA
pcoa_functio_trans_sla_CA <- pco(dist_mat_functio_trans_sla_CA, negvals = "zero", dround = 0)
# keep the first 3 PCs
plot(pcoa_functio_trans_sla_CA$values/sum(pcoa_functio_trans_sla_CA$values))
sum((pcoa_functio_trans_sla_CA$values/sum(pcoa_functio_trans_sla_CA$values))[1:3])
PCoA_PCs_functio_trans_SLA_CA <- pcoa_functio_trans_sla_CA$vectors[,1:3]




###########################################################################
#### make groups ####
###########################################################################
pattern <- c("-", "COR", "Forêt Nord", "PGK", "Grand Lac", "Kuebini", "Wadjana")

locality_names  <- geom_plots$locality
plot_grps <- c(rep(NA,length(locality_names)))
for(i in 1:length(locality_names)){
  for(j in pattern){
    if( grepl(j,locality_names[i])) plot_grps[i] <- j
  }
}
# merge the only plot from Kuebini (donwside) with Wadjana
plot_grps[plot_grps == "Kuebini"] <- "Wadjana"
plot_grps[plot_grps == "-"] <- "Kuebini"
plot_grps[plot_grps == "COR"] <- "CORIFOR"

localty_grp <- plot_grps

###########################################################################
###                                                                     ###
####            Compile diversity indices with spatial data            ####
###                                                                     ###
###########################################################################
# verif order plots
cbind(geom_plots$locality,data_ncpipn_plots_sf_PDL_utm$locality)

# export plot geometry (caanot include diversity indices bbeacause theare are to much and colnames do not fit for a shapefile)
geom_plots_PDL_div_indices <- cbind(geom_plots,
                                    localty_grp
      )
# export
library(sf)
st_write(geom_plots_PDL_div_indices, dsn = '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data',
         layer = 'geom_plots_PDL_div_indices_utm', driver = "ESRI Shapefile", append=F)

#### complil indices ####
plots_PDL_div_indices <- cbind(data.frame(geom_plots),
                                    localty_grp,
                                    taxo_div,
                                    rar_div2,
                                    PCoA_PC1_taxo = PCoA_PCs_taxo[,1],
                                    PCoA_PC2_taxo = PCoA_PCs_taxo[,2],
                                    PCoA_PC3_taxo = PCoA_PCs_taxo[,3],
                                    
                                    PCoA_PC1_taxo_BA = PCoA_PCs_taxo_BA[,1],
                                    PCoA_PC2_taxo_BA = PCoA_PCs_taxo_BA[,2],
                                    PCoA_PC3_taxo_BA = PCoA_PCs_taxo_BA[,3],
                                    
                                    PCoA_PC1_taxo_CA = PCoA_PCs_taxo_CA[,1],
                                    PCoA_PC2_taxo_CA = PCoA_PCs_taxo_CA[,2],
                                    PCoA_PC3_taxo_CA = PCoA_PCs_taxo_CA[,3],
                                    
                                    FRic = FD_PDL$FRic,
                                    FEve = FD_PDL$FEve,
                                    FDiv = FD_PDL$FDiv,
                                    FDis = FD_PDL$FDis,
                                    RaoQ = FD_PDL$RaoQ,
                                    FD_PDL$CWM,
                               
                                     FEve_SLA = FD_PDL_SLA$FEve,
                                     FDis_SLA = FD_PDL_SLA$FDis,
                                     RaoQ_SLA = FD_PDL_SLA$RaoQ,
                                    
                                    FRic_BA = FD_PDL_BA$FRic,
                                    FEve_BA = FD_PDL_BA$FEve,
                                    FDiv_BA = FD_PDL_BA$FDiv,
                                    FDis_BA = FD_PDL_BA$FDis,
                                    RaoQ_BA = FD_PDL_BA$RaoQ,
                                    FD_PDL_BA$CWM,
                               
                                     FEve_BA_SLA = FD_PDL_BA_SLA$FEve,
                                     FDis_BA_SLA = FD_PDL_BA_SLA$FDis,
                                     RaoQ_BA_SLA = FD_PDL_BA_SLA$RaoQ,
                                           
                                     FRic_CA = FD_PDL_CA$FRic,
                                     FEve_CA = FD_PDL_CA$FEve,
                                     FDiv_CA = FD_PDL_CA$FDiv,
                                     FDis_CA = FD_PDL_CA$FDis,
                                     RaoQ_CA = FD_PDL_CA$RaoQ,
                                     FD_PDL_CA$CWM,
                               
                                     FEve_CA_SLA = FD_PDL_CA_SLA$FEve,
                                     FDis_CA_SLA = FD_PDL_CA_SLA$FDis,
                                     RaoQ_CA_SLA = FD_PDL_CA_SLA$RaoQ,
                                    
                                    FRic_trans = FD_PDL_trans$FRic,
                                    FEve_trans = FD_PDL_trans$FEve,
                                    FDiv_trans = FD_PDL_trans$FDiv,
                                    FDis_trans = FD_PDL_trans$FDis,
                                    RaoQ_trans = FD_PDL_trans$RaoQ,
                                    FD_PDL_trans$CWM,
                               
                                     FEve_trans_SLA = FD_PDL_trans_SLA$FEve,
                                     FDis_trans_SLA = FD_PDL_trans_SLA$FDis,
                                     RaoQ_trans_SLA = FD_PDL_trans_SLA$RaoQ,
                                    
                                    FRic_trans_BA = FD_PDL_trans_BA$FRic,
                                    FEve_trans_BA = FD_PDL_trans_BA$FEve,
                                    FDiv_trans_BA = FD_PDL_trans_BA$FDiv,
                                    FDis_trans_BA = FD_PDL_trans_BA$FDis,
                                    RaoQ_trans_BA = FD_PDL_trans_BA$RaoQ,
                                    FD_PDL_trans_BA$CWM,
                               
                                     FEve_trans_BA_SLA = FD_PDL_trans_BA_SLA$FEve,
                                     FDis_trans_BA_SLA = FD_PDL_trans_BA_SLA$FDis,
                                     RaoQ_trans_BA_SLA = FD_PDL_trans_BA_SLA$RaoQ,
                                     
                                     FRic_trans_CA = FD_PDL_trans_CA$FRic,
                                     FEve_trans_CA = FD_PDL_trans_CA$FEve,
                                     FDiv_trans_CA = FD_PDL_trans_CA$FDiv,
                                     FDis_trans_CA = FD_PDL_trans_CA$FDis,
                                     RaoQ_trans_CA = FD_PDL_trans_CA$RaoQ,
                                     FD_PDL_trans_CA$CWM,
                               
                                     FEve_trans_CA_SLA = FD_PDL_trans_CA_SLA$FEve,
                                     FDis_trans_CA_SLA = FD_PDL_trans_CA_SLA$FDis,
                                     RaoQ_trans_CA_SLA = FD_PDL_trans_CA_SLA$RaoQ,
                                    
                                    FRic_trans_synth_BA = FD_PDL_trans_synthetic_BA$FRic,
                                    FEve_trans_synth_BA = FD_PDL_trans_synthetic_BA$FEve,
                                    FDiv_trans_synth_BA = FD_PDL_trans_synthetic_BA$FDiv,
                                    FDis_trans_synth_BA = FD_PDL_trans_synthetic_BA$FDis,
                                    RaoQ_trans_synth_BA = FD_PDL_trans_synthetic_BA$RaoQ,
                                    FD_PDL_trans_synthetic_BA$CWM,
                               
                                     FRic_trans_synth_CA = FD_PDL_trans_synthetic_CA$FRic,
                                     FEve_trans_synth_CA = FD_PDL_trans_synthetic_CA$FEve,
                                     FDiv_trans_synth_CA = FD_PDL_trans_synthetic_CA$FDiv,
                                     FDis_trans_synth_CA = FD_PDL_trans_synthetic_CA$FDis,
                                     RaoQ_trans_synth_CA = FD_PDL_trans_synthetic_CA$RaoQ,
                                     FD_PDL_trans_synthetic_CA$CWM,
                                    
                                    PCoA_PC1_functio = PCoA_PCs_functio[,1],
                                    PCoA_PC2_functio = PCoA_PCs_functio[,2],
                                    PCoA_PC3_functio = PCoA_PCs_functio[,3],
                               
                                     PCoA_PC1_functio_SLA = PCoA_PCs_functio_SLA[,1],
                                     PCoA_PC2_functio_SLA = PCoA_PCs_functio_SLA[,2],
                                     PCoA_PC3_functio_SLA = PCoA_PCs_functio_SLA[,3],
                                          
                                    PCoA_PC1_functio_BA = PCoA_PCs_functio_BA[,1],
                                    PCoA_PC2_functio_BA = PCoA_PCs_functio_BA[,2],
                                    PCoA_PC3_functio_BA = PCoA_PCs_functio_BA[,3],
                               
                                     PCoA_PC1_functio_SLA_BA = PCoA_PCs_functio_SLA_BA[,1],
                                     PCoA_PC2_functio_SLA_BA = PCoA_PCs_functio_SLA_BA[,2],
                                     PCoA_PC3_functio_SLA_BA = PCoA_PCs_functio_SLA_BA[,3],
                                     
                                     PCoA_PC1_functio_CA = PCoA_PCs_functio_CA[,1],
                                     PCoA_PC2_functio_CA = PCoA_PCs_functio_CA[,2],
                                     PCoA_PC3_functio_CA = PCoA_PCs_functio_CA[,3],
                               
                                     PCoA_PC1_functio_SLA_CA = PCoA_PCs_functio_SLA_CA[,1],
                                     PCoA_PC2_functio_SLA_CA = PCoA_PCs_functio_SLA_CA[,2],
                                     PCoA_PC3_functio_SLA_CA = PCoA_PCs_functio_SLA_CA[,3],
                                          
                                    PCoA_PC1_functio_trans = PCoA_PCs_functio_trans[,1],
                                    PCoA_PC2_functio_trans = PCoA_PCs_functio_trans[,2],
                                    PCoA_PC3_functio_trans = PCoA_PCs_functio_trans[,3],
                                    
                                    PCoA_PC1_functio_trans_BA = PCoA_PCs_functio_trans_BA[,1],
                                    PCoA_PC2_functio_trans_BA = PCoA_PCs_functio_trans_BA[,2],
                                    PCoA_PC3_functio_trans_BA = PCoA_PCs_functio_trans_BA[,3],
                               
                                     PCoA_PC1_functio_trans_SLA = PCoA_PCs_functio_trans_SLA[,1],
                                     PCoA_PC2_functio_trans_SLA = PCoA_PCs_functio_trans_SLA[,2],
                                     PCoA_PC3_functio_trans_SLA = PCoA_PCs_functio_trans_SLA[,3],
                                     
                                     PCoA_PC1_functio_trans_CA = PCoA_PCs_functio_trans_CA[,1],
                                     PCoA_PC2_functio_trans_CA = PCoA_PCs_functio_trans_CA[,2],
                                     PCoA_PC3_functio_trans_CA = PCoA_PCs_functio_trans_CA[,3],
                                    
                                    PCoA_PC1_functio_trans_synt_BA = PCoA_PCs_functio_trans_synt_BA[,1],
                                    PCoA_PC2_functio_trans_synt_BA = PCoA_PCs_functio_trans_synt_BA[,2],
                                    PCoA_PC3_functio_trans_synt_BA = PCoA_PCs_functio_trans_synt_BA[,3],
                              
                                      PCoA_PC1_functio_trans_SLA_BA = PCoA_PCs_functio_trans_SLA_BA[,1],
                                     PCoA_PC2_functio_trans_SLA_BA = PCoA_PCs_functio_trans_SLA_BA[,2],
                                     PCoA_PC3_functio_trans_SLA_BA = PCoA_PCs_functio_trans_SLA_BA[,3],
                               
                                     PCoA_PC1_functio_trans_synt_CA = PCoA_PCs_functio_trans_synt_CA[,1],
                                     PCoA_PC2_functio_trans_synt_CA = PCoA_PCs_functio_trans_synt_CA[,2],
                                     PCoA_PC3_functio_trans_synt_CA = PCoA_PCs_functio_trans_synt_CA[,3],
                               
                                     PCoA_PC1_functio_trans_SLA_CA = PCoA_PCs_functio_trans_SLA_CA[,1],
                                     PCoA_PC2_functio_trans_SLA_CA = PCoA_PCs_functio_trans_SLA_CA[,2],
                                     PCoA_PC3_functio_trans_SLA_CA = PCoA_PCs_functio_trans_SLA_CA[,3]
                               
                               
)
# export
# saveRDS(plots_PDL_div_indices,'/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_PDL_div_indices.rds')

###########################################################################
###                                                                     ###
####               Save RDS object with matrix type data               ####
###                                                                     ###
###########################################################################

plot_data_matrices <- list(locality = geom_plots$locality,
                           community_matrix = community_matrix,
                           community_matrix_BA = community_matrix_BA,
                           dist_mat_taxo = dist_mat_taxo,
                           dist_mat_taxo_BA = dist_mat_taxo_BA,
                           dist_mat_functio = dist_mat_functio,
                           dist_mat_functio_BA = dist_mat_functio_BA,
                           dist_mat_functio_trans = dist_mat_functio_trans,
                           dist_mat_functio_trans_BA = dist_mat_functio_trans_BA,
                           dist_mat_functio_trans_CA = dist_mat_functio_trans_CA,
                           dist_mat_functio_trans_synth_BA = dist_mat_functio_trans_synth_BA,
                           dist_mat_functio_trans_synth_CA = dist_mat_functio_trans_synth_CA,
                           dist_mat_functio_sla = dist_mat_functio_sla,
                           dist_mat_functio_trans_sla = dist_mat_functio_trans_sla,
                           dist_mat_functio_sla_BA = dist_mat_functio_sla_BA,
                           dist_mat_functio_trans_sla_BA = dist_mat_functio_trans_sla_BA,
                           dist_mat_functio_sla_CA = dist_mat_functio_sla_CA,
                           dist_mat_functio_trans_sla_CA = dist_mat_functio_trans_sla_CA,
                           FT_select = FT_select,
                           FT_select_trans = FT_select_trans,
                           community_matrix_for_traits = community_matrix_for_traits,
                           community_matrix_for_traits_BA = community_matrix_for_traits_BA)
dir.create('/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices')
# saveRDS(plot_data_matrices, '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plot_data_matrices.rds')








######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
##### Same but for filtered dataset  ####
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

#### WITH DATASET FILTERED : ONLY PLOTS THAT ARE IN SPECTRAL SPECTIES RASTER (see) ####


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
##### OLD SCRIPT FOR VALE KUEBINI REPORT  AND OTHER STUFF ####
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

###########################################################################
#### Functional diversity with other packages (only for pres/abs data) ####
###########################################################################

library(betapart)
dist_mat_functio <- functional.beta.pair(community_matrix_for_traits,FT_select)

Library(mDF)
# Computing distances between species based on functional traits
trat_cat <- data.frame(trait_name = colnames(FT_select), trait_type = (rep("Q", length(colnames(FT_select)))))
sp_dist_funct <- mFD::funct.dist(
  sp_tr         = FT_select,
  tr_cat        = trat_cat,
  metric        = "euclidean",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)
# Compute multimensional functional spaces and assess their quality
fspaces_quality_funct <- mFD::quality.fspaces(
  sp_dist             = sp_dist_funct,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

sp_faxes_coord_funct <- fspaces_quality_funct$"details_fspaces"$"sp_pc_coord"

beta_fd_indices_funct <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_funct[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_occ       = community_matrix_for_traits,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)


###########################################################################
###                                                                     ###
####                  Make geographic groups                           ####
###                                                                     ###
###########################################################################

#### change plot names ####

library(stringr)

rownames(community_matrix)
remove <- c("-BIS")

locality_names <- str_remove_all(rownames(community_matrix), paste(remove, collapse = "|"))
locality_names <- strsplit(locality_names,"-")
kbn=1
for(i in 1:length(locality_names)){
  if(length(locality_names[[i]])>1){
    locality_names[[i]] <- paste0("KB", kbn,"-", locality_names[[i]][2])
    kbn = kbn+1
  } 
}
locality_names <- unlist(locality_names)
rownames(community_matrix) <- locality_names

#### make groups ####

pattern <- c("KB", "COR", "Forêt Nord", "PGK", "Grand Lac", "Kuebini", "Wadjana")

plot_grps <- c(rep(NA,length(locality_names)))
for(i in 1:length(locality_names)){
  for(j in pattern){
     if( grepl(j,locality_names[i])) plot_grps[i] <- j
  }
}
# merge the only plot from Kuebini (donwside) with Wadjana

plot_grps[plot_grps == "Kuebini"] <- "Wadjana"
plot_grps[plot_grps == "KB"] <- "Kuebini"
plot_grps[plot_grps == "COR"] <- "CORIFOR"

###########################################################################
###                                                                     ###
####                             NMDS                                  ####
###                                                                     ###
###########################################################################
library(vegan)
df <- community_matrix
# log transform abundances? 
df <- sqrt(community_matrix)


# number of dimensions
n = 5
stress <- vector(length = n)
for (i in 1:n) {
  stress[i] <- metaMDS(df, distance = "bray", k = i)$stress
}
names(stress) <- paste0(1:n, "Dim")
# x11(width = 10/2.54, height = 7/2.54)
par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")

# choice for NMDS
NMDS=metaMDS(df,
             k=3,
             trymax=100)

plot(NMDS)
stressplot(NMDS)
# plot(NMDS, choices = c(1,3))
# ordihull(NMDS,groups=plot_grps,draw="polygon",col="grey90",label=F, choices =c(1,3))

ordiplot(NMDS,type="n")
orditorp(NMDS,display="species",col="red",air=0.01)
orditorp(NMDS,display="sites",cex=1,air=0.01)

ordihull(NMDS,groups=plot_grps,draw="polygon",col="grey90",label=F)

############################################################################
#### Nice plot NMDS ####
############################################################################

#### extract NMDS parameters ####
data.scores <- as.data.frame(scores(NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
# distance edge class
data.scores$plot_grps <- plot_grps   #  add the grp variable created earlier
data.scores$plot_grps <- factor(data.scores$plot_grps, levels = unique(data.scores$plot_grps))

head(data.scores)  #look at the data
species.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

####  display only some species #### 
library(tidyverse)
plot_sp_all <- scores(NMDS, display = "species") %>% 
  as.data.frame() %>% 
  rownames_to_column("species")
# most abundant species in all the dataset
most_abund_sp <- sort(names(sort(colSums(df), decreasing = T)[1:60]))
most_abund_sp
plot_sp <- plot_sp_all[plot_sp_all$species %in% most_abund_sp,]
# most abundant species in each sub-dataset
pattern <- c("KB", "COR", "Forêt Nord", "PGK", "Grand Lac", "Wadjana")
most_abund_sp <- c()
for(j in pattern){
  df_tmp <- df[grepl(j, rownames(df)),]
  most_abund_sp_tmp <- sort(names(sort(colSums(df_tmp), decreasing = T)[1:8]))
  if(j %in% c("Kuebini") )  most_abund_sp_tmp <- sort(names(sort(colSums(df_tmp), decreasing = T)[1:15]))
  if(j %in% c("PGK", "Grand Lac") )  most_abund_sp_tmp <- sort(names(sort(colSums(df_tmp), decreasing = T)[1:7]))
  if(j %in% c("Forêt Nord") ){
    most_abund_sp_tmp <- sort(names(sort(colSums(df_tmp), decreasing = T)[1:8]))
    other_abund_sp_tmp <- sort(names(sort(colSums(df_tmp), decreasing = T)[8:52]))
    other_abund_sp_tmp <- other_abund_sp_tmp[other_abund_sp_tmp %in% plot_sp_all$species[plot_sp_all$NMDS1>.75 & plot_sp_all$NMDS2>.5]]
    most_abund_sp_tmp <- c(most_abund_sp_tmp,
                           other_abund_sp_tmp)
  }
  most_abund_sp <- c(most_abund_sp, most_abund_sp_tmp)
}

sp_keep <- unique(plot_sp_all$species[plot_sp_all$species %in% most_abund_sp])
# significant sp in NMDS among most abundants
# df_abundant <- df[,colnames(df) %in% most_abund_sp]
# fit_sp <- envfit(NMDS, df_abundant)
# fit_sp_axe1 <- envfit(NMDS, df_abundant, choices = 1)
# fit_sp_axe2 <- envfit(NMDS, df_abundant, choices = 2)
# keep only significant sp 
# sp_keep <- names(fit_sp$vectors$pvals[fit_sp$vectors$pvals<0.05])
sp_keep[sp_keep == "Plerandra {gordonii}"] <- "Plerandra gordonii"
sp_keep
species_to_plot <- species.scores[species.scores$species %in% sp_keep,]
# only species identified
species_to_plot <- 
  species_to_plot[unlist(lapply(species_to_plot$species,
                                function(x) length(unlist(strsplit(x, " "))))) > 1 ,]
# use centroids of plot groups
scrs <-
  scores(NMDS, display = "sites", "species")
centroid_plot_grps <-
  aggregate(scrs ~plot_grps, FUN = "mean")
colnames(centroid_plot_grps)[1] <- "plot_grps"
centroid_plot_grps$plot_grps <- factor(centroid_plot_grps$plot_grps, levels = centroid_plot_grps$plot_grps)
#### create hulls for edge class ####
library(data.table)
hulls_plot <- data.table(data.scores)[, .SD[chull(NMDS1, NMDS2)], by = plot_grps]

#### nice plot ####
sp_size <- log10(colSums(df)+10)*5
species.scores$size <- sp_size
species.scores$size <- as.numeric(sp_size)

species.scores[species.scores$NMDS2>0.8,]
##### plot all  #####
library(ggplot2)
library(ggrepel)

plot_NMDS_all <- ggplot() + 
  geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2), alpha=0)  +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, fill = plot_grps), pch = 22, alpha = .7, size = 6.3) +
  geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2, size = size), color = "grey",alpha=0.3)   +
   guides(size= "none")  +
  theme_bw() +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position=c(.15,.82),
        legend.box.background = element_rect(color="darkgrey", size=1)
  ) +
  geom_label_repel(data = species_to_plot,
                   aes(x = NMDS1, y = NMDS2, label = species,
                       hjust=.5), alpha=0.7, 
                   size = 3.3, fontface = "italic", lineheight = .75,
                   show.legend  = FALSE, max.overlaps = 20,
                   colour = "black") +
  annotate("text", x=0.8, y= -1.1, label= paste0("Stress = ", round(NMDS$stress, digits = 2))) +
  xlim(c(-1.6,1.8)) +
  ylim(c(-1.05,1.2)) +
  guides(fill=guide_legend(title="Plot groups")) 
  
plot_NMDS_all


#### export png ####
png("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/figures/plot_NMDS_PDL2", width=265, height=265, units = 'mm', res = 300) 

plot_NMDS_all
dev.off()

###########################################################################
###                                                                     ###
####                           map                                     ####
###                                                                     ###
###########################################################################


# forest 
forest <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/GIS/forest_map/forest_present_PDL.shp")

aerial_photo <- brick("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/rapport/figures/map_grandsud.tif")

#### groups of plot for map ####
data_ncpipn_taxon_sf_PDL_for_map <- data_ncpipn_taxon_sf_PDL
data_ncpipn_taxon_sf_PDL_for_map <- unique(data_ncpipn_taxon_sf_PDL_for_map[,c("locality","geo_pt..89")])

#### change plot names ####
library(stringr)
remove <- c("-BIS")
data_ncpipn_taxon_sf_PDL_for_map$locality <- str_remove_all(data_ncpipn_taxon_sf_PDL_for_map$locality, paste(remove, collapse = "|"))
data_ncpipn_taxon_sf_PDL_for_map$locality <- strsplit(data_ncpipn_taxon_sf_PDL_for_map$locality,"-")
kbn=1
for(i in 1:length(data_ncpipn_taxon_sf_PDL_for_map$locality)){
  if(length(data_ncpipn_taxon_sf_PDL_for_map$locality[[i]])>1){
    data_ncpipn_taxon_sf_PDL_for_map$locality[[i]] <- paste0("KB", kbn,"-", data_ncpipn_taxon_sf_PDL_for_map$locality[[i]][2])
    kbn = kbn+1
  } 
}
data_ncpipn_taxon_sf_PDL_for_map$locality <- unlist(data_ncpipn_taxon_sf_PDL_for_map$locality)

#### make groups ####
plots_PDL_for_map <- data_ncpipn_taxon_sf_PDL_for_map
pattern <- c("KB", "COR", "Forêt Nord", "PGK", "Grand Lac", "Kuebini", "Wadjana")

for(i in 1:length(plots_PDL_for_map$locality)){
  for(j in pattern){
    if( grepl(j,plots_PDL_for_map$locality[i])) plots_PDL_for_map$locality [i] <- j
  }
}
# merge the only plot from Kuebini (donwside) with Wadjana

plots_PDL_for_map$locality[plots_PDL_for_map$locality == "Kuebini"] <- "Wadjana"
plots_PDL_for_map$locality[plots_PDL_for_map$locality == "KB"] <- "Kuebini"
plots_PDL_for_map$locality[plots_PDL_for_map$locality == "COR"] <- "CORIFOR"
plots_PDL_for_map <- as_Spatial(plots_PDL_for_map)
##### transform multiband image to raster ##### 
library(raster)
library(rgeos)
library(sf)
library(rasterVis)
  
aerial_photo_ok_r <- raster(aerial_photo)
cols_bg_image <- factor(rgb(aerial_photo[], maxColorValue=255))
Sys.sleep(5)
aerial_photo_ok_r[] <- cols_bg_image
levelplot(aerial_photo_ok_r, col.regions=as.character(levels(cols_bg_image)), colorkey=FALSE)

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

bg_image <- levelplot(aerial_photo_ok_r, col.regions=as.character(levels(cols_bg_image)), colorkey=FALSE)
forest_layer <- latticeExtra::layer(sp.polygons(forest , lwd=2, fill = "forestgreen", alpha = 0.3)) 
plots_KB <- latticeExtra::layer(sp.points(plots_PDL_for_map[plots_PDL_for_map$locality == "Kuebini",], pch= 15, cex=.8, col= "Brown2", alpha = .8)) 
plots_CORIFOR <- latticeExtra::layer(sp.points(plots_PDL_for_map[plots_PDL_for_map$locality == "CORIFOR",], pch= 15, cex=.8, col= "gold3", alpha = .8)) 
plots_GrandLac <- latticeExtra::layer(sp.points(plots_PDL_for_map[plots_PDL_for_map$locality == "Grand Lac",], pch= 15, cex=.8, col= "turquoise2", alpha = .8)) 
plots_ForetNord <- latticeExtra::layer(sp.points(plots_PDL_for_map[plots_PDL_for_map$locality == "Forêt Nord",], pch= 15, cex=.8, col= "springgreen3", alpha = .8)) 
plots_Wadjana <- latticeExtra::layer(sp.points(plots_PDL_for_map[plots_PDL_for_map$locality == "Wadjana",], pch= 15, cex=.8, col= "dodgerblue1", alpha = .8)) 
plots_PGK <- latticeExtra::layer(sp.points(plots_PDL_for_map[plots_PDL_for_map$locality == "PGK",], pch= 15, cex=.8, col= "magenta", alpha = .8)) 


plots_map <- bg_image + forest_layer + plots_KB + plots_CORIFOR + plots_GrandLac + plots_ForetNord + plots_Wadjana + plots_PGK
plots_map

#### export png ####
png("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/figures/plots_map", width=200, height=200, units = 'mm', res = 300) 

plots_map
dev.off()

###########################################################################
###                                                                     ###
####                     distance to forest edge                       ####
###                                                                     ###
###########################################################################
data_ncpipn_taxon_sf_PDL_for_map
data_PDL_for_map <- as_Spatial(data_ncpipn_taxon_sf_PDL_for_map)

####  distance to edge  #### 
# ID patch
patchs_over_point <- over(data_PDL_for_map,forest) 
data_PDL_for_map$id <- patchs_over_point$id
# use metric distance 
forest_utm <- spTransform(forest, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
data_PDL_for_map_utm <- spTransform(data_PDL_for_map, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

dist_edge <- c()
for (i in 1:nrow(data_PDL_for_map)) {
  pts_tmp <- data_PDL_for_map_utm[i,]
  ## NA if points are out of forest (reference stations "cagnar")
  if(!is.na(pts_tmp$id)){
    patch_tmp <- forest_utm[forest_utm$id == pts_tmp$id,]
    forest_edge <- as(patch_tmp, "SpatialLines") 
    dist_edge <- c(dist_edge, gDistance(pts_tmp, forest_edge, byid=TRUE)) 
  }else{
    dist_edge <- c(dist_edge, NA) 
  }
}
data_PDL_for_map$dist_edge <- dist_edge

###########################################################################
###                                                                     ###
####                  geographic distance decay                        ####
###                                                                     ###
###########################################################################
library(betapart)
# same order 
community_matrix <- community_matrix[order(rownames(community_matrix)),]
data_ncpipn_taxon_sf_PDL_for_map <- data_ncpipn_taxon_sf_PDL_for_map[order(data_ncpipn_taxon_sf_PDL_for_map$locality),]

# dissimilarity matrix
dist_mat_bray <- vegdist(community_matrix)
dist_mat_bray
bray_part <- bray.part(community_matrix)

dist_mat_bray <- as.matrix(bray_part$bray)
dist_mat_bray[upper.tri(dist_mat_bray)] <- NA
# distance matrix
dist_mat_geo <- st_distance(data_ncpipn_taxon_sf_PDL_for_map)
dist_mat_geo <- as.matrix(dist_mat_geo)
dist_mat_geo[upper.tri(dist_mat_geo)] <- NA
rownames(dist_mat_geo) <- colnames(dist_mat_geo) <- data_ncpipn_taxon_sf_PDL_for_map$locality
units(dist_mat_geo) <- NULL

# merge matrices
df_bc <- reshape2::melt(dist_mat_bray, varnames = c("row", "col"))
df_bc <- df_bc[!is.na(df_bc$value) ,]

df_dist <- reshape2::melt(dist_mat_geo, varnames = c("row", "col"))
df_dist <- df_dist[!is.na(df_dist$value) ,]

df_dist_bc <- cbind(df_dist, df_bc$value)
names(df_dist_bc) <- c("row", "col", "m", "bc")
df_dist_bc <- df_dist_bc[complete.cases(df_dist_bc),]
# no zero distance (same plot)
df_dist_bc <- df_dist_bc[df_dist_bc$m>0,]

# linear model
lm1 <- lm(bc ~ m, data = df_dist_bc)
smry1 <- summary(lm1)

# mant_test1 <- mantel(dist_mat_bray, dist_mat_geo, method = "spearman", permutations = 9999, na.rm = TRUE)
# saveRDS(mant_test1, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/mantel_test_PDL.rds")

# saveRDS(mant_test1_log, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/mant_test_log_PDL.rds")


#### plot ####

#with confidence intervals
ts = 3.8
plot_geo <- ggplot(data = df_dist_bc, aes(x = m, y = bc)) +
  geom_point(color='darkolivegreen', size = 1, alpha =.1) +
  geom_smooth( formula = y ~ log(x),
              method = "lm", color='darkred', se=T, fill = "darkred" , alpha = 0.3) +
  annotate("text", x=12000, y= .3,
           label= paste0(" Mantel r = ", round(mant_test1$statistic, digits = 2), gtools::stars.pval(mant_test1$signif)),
           color='darkred',  size = ts) +
  
  xlab("Geographic distance (m)") + ylab("Bray-Curtis dissimilarity") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
  ) +
  theme_bw() 

plot_geo

#### export png ####
png("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/figures/distance_decay_PDL.png", width=120, height=90, units = 'mm', res = 300) 

plot_geo
dev.off()

###########################################################################
###                                                                     ###
####                        distance decay from edge                  ####
###                                                                     ###
###########################################################################
library(betapart)
# same order 
community_matrix <- community_matrix[order(rownames(community_matrix)),]
data_ncpipn_taxon_sf_PDL_for_map <- data_ncpipn_taxon_sf_PDL_for_map[order(data_ncpipn_taxon_sf_PDL_for_map$locality),]

# dissimilarity matrix
dist_mat_bray <- vegdist(community_matrix)
dist_mat_bray
bray_part <- bray.part(community_matrix)

dist_mat_bray <- as.matrix(bray_part$bray)
dist_mat_bray[upper.tri(dist_mat_bray)] <- NA
# distance matrix
dist_mat_edge <- dist(data_PDL_for_map$dist_edge)
dist_mat_edge <- as.matrix(dist_mat_edge)
dist_mat_edge[upper.tri(dist_mat_edge)] <- NA
rownames(dist_mat_edge) <- colnames(dist_mat_edge) <- data_ncpipn_taxon_sf_PDL_for_map$locality

# merge matrices
df_bc <- reshape2::melt(dist_mat_bray, varnames = c("row", "col"))
df_bc <- df_bc[!is.na(df_bc$value) ,]
# no zero dissimilarity (same plot)
df_bc <- df_bc[df_bc$value > 0 ,]

df_dist <- reshape2::melt(dist_mat_edge, varnames = c("row", "col"))
df_dist <- df_dist[!is.na(df_dist$value) ,]

df_dist_bc <- df_dist %>% dplyr::right_join( df_bc, by=c("row","col"))

names(df_dist_bc) <- c("row", "col", "m", "bc")
df_dist_bc <- df_dist_bc[complete.cases(df_dist_bc),]
df_dist_bc <- df_dist_bc[df_dist_bc$m>0,]

# linear model
lm1 <- lm(bc ~ m, data = df_dist_bc)
smry1 <- summary(lm1)

 # mant_test2 <- mantel(dist_mat_bray, dist_mat_geo, method = "spearman", permutations = 9999, na.rm = TRUE)
# saveRDS(mant_test2, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/mantel_test_edge_PDL.rds")


#### plot ####

#with confidence intervals
ts = 3.8
plot_edge <- ggplot(data = df_dist_bc, aes(x = m, y = bc)) +
  geom_point(color='darkolivegreen', size = 1, alpha =.1) +
  geom_smooth( formula = y ~ log(x+1),
               method = "lm", color='darkred', se=T, fill = "darkred" , alpha = 0.3) +
  annotate("text", x=300, y= .3,
           label= paste0(" Mantel r = ", round(mant_test2$statistic, digits = 2), gtools::stars.pval(mant_test1$signif)),
           color='darkred',  size = ts) +
  
  xlab("Difference in distance from edge (m)") + ylab("Bray-Curtis dissimilarity") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
  ) +
  theme_bw() 

plot_edge

#### export png ####
png("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/figures/distance_decay_edge_PDL.png", width=120, height=90, units = 'mm', res = 300) 

plot_edge
dev.off()


#### export png combined ####

png("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/figures/distance_decay_geo_edge_PDL.png", width=220, height=90, units = 'mm', res = 300) 

do.call(gridExtra::grid.arrange,c(list(plot_geo, plot_edge), ncol = 2))

dev.off()


###########################################################################
###                                                                     ###
####       dissimilarity by groups of distance from edge               ####
###                                                                     ###
###########################################################################

library(betapart)
# same order 
community_matrix <- community_matrix[order(rownames(community_matrix)),]
data_ncpipn_taxon_sf_PDL_for_map <- data_ncpipn_taxon_sf_PDL_for_map[order(data_ncpipn_taxon_sf_PDL_for_map$locality),]
data_PDL_for_map <- data_PDL_for_map[order(data_PDL_for_map$locality),]
# groups of distance from edge

data_PDL_for_map$dist_edge_class <- NA
# data_PDL_for_map$dist_edge_class[data_PDL_for_map$dist_edge<25] <-"0-25"
# data_PDL_for_map$dist_edge_class[data_PDL_for_map$dist_edge>25 & data_PDL_for_map$dist_edge<50] <-"25-50"
# data_PDL_for_map$dist_edge_class[data_PDL_for_map$dist_edge>50 & data_PDL_for_map$dist_edge<100] <-"50-100"
# data_PDL_for_map$dist_edge_class[data_PDL_for_map$dist_edge>100] <-">100"

data_PDL_for_map$dist_edge_class[data_PDL_for_map$dist_edge<50] <-"0-50"
data_PDL_for_map$dist_edge_class[data_PDL_for_map$dist_edge>50 & data_PDL_for_map$dist_edge<100] <-"50-100"
data_PDL_for_map$dist_edge_class[data_PDL_for_map$dist_edge>100] <-">100"

# no NA
data_PDL_for_map <- data_PDL_for_map[!is.na(data_PDL_for_map$dist_edge_class),]
# dissimilarity matrice for each group
df_bc_class_edge <- c()
for (i in 1:length(unique(data_PDL_for_map$dist_edge_class))){
  class_tmp <- unique(data_PDL_for_map$dist_edge_class)[i]
  comat_tmp <- community_matrix[rownames(community_matrix) %in% data_PDL_for_map$locality[data_PDL_for_map$dist_edge_class == class_tmp],]
  # no empty species
  comat_tmp <- comat_tmp[,colSums(comat_tmp)>0]
  # presabs
  comat_tmp[comat_tmp>0] <- 1
  bray_part_tmp <- bray.part(comat_tmp)
  dist_mat_bray_tmp <- as.matrix(bray_part_tmp$bray)
  dist_mat_bray_tmp[upper.tri(dist_mat_bray_tmp)] <- NA
  df_bc_tmp <- reshape2::melt(dist_mat_bray_tmp, varnames = c("row", "col"))
  df_bc_tmp <- df_bc_tmp[!is.na(df_bc_tmp$value) ,]
  df_bc_tmp <- df_bc_tmp[df_bc_tmp$value > 0 ,]
  df_bc_class_edge_tmp <- cbind(df_bc_tmp, grp = rep(class_tmp, nrow(df_bc_tmp)))
  df_bc_class_edge <- rbind(df_bc_class_edge, df_bc_class_edge_tmp)
}

df_bc_class_edge <- df_bc_class_edge[order(df_bc_class_edge$grp),]

colnames(df_bc_class_edge)

#### plot ####
boxplot_BC <- ggplot(df_bc_class_edge, aes(x=grp, y=value, fill=grp)) +
  # geom_violin() +
  geom_boxplot(width=0.5, color="black", alpha=.3,outlier.shape=NA, lwd=.5, notch = TRUE, coef = .25) + 
  geom_point(color='black', size = 1, alpha =.1) +
  # scale_x_discrete(limits=c("Hydro_and_other", "Anemo", "Zoo"), labels = NULL) +
  # labels = c("Anemochorous\n& zoochorous", "Hydrochorous", "Other species") ) +
  scale_fill_manual(values= c( "darkolivegreen","dodgerblue3", "darkorange" ,"darkred")) +
  #  theme_ipsum() +
  theme_bw()
boxplot_BC
  
theme(
    legend.position="none",
    plot.title = element_text(size=14),
    plot.subtitle = element_text(size=13),
    axis.text.x = element_text(size=10, colour = "white"),
    axis.title.y = element_text(size=14, hjust = .6),
    axis.ticks.x=element_blank()
  ) +
  xlab("") +
  ylab("") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", inherit.aes = FALSE, tip.length = 0.005, size = 5)+
  coord_cartesian(ylim=c(0.12,.65)) 

boxplot_BC
###########################################################################
###                                                                     ###
####                        compare diversity                          ####
###                                                                     ###
###########################################################################

#### But plots does not have the same size !!!!! ####

community_matrix01 <- community_matrix
community_matrix01[community_matrix01>1] <-1
RS <- rowSums(community_matrix01)
RS_df <- data.frame(cbind(RS, plot_grps))
RS_df$RS <- as.numeric(RS_df$RS)
# library
library(ggplot2)
library(dplyr)
library(hrbrthemes)
p <- RS_df %>%
  ggplot( aes(x=RS, fill=plot_grps)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',  stat="count") 
p

p <- RS_df %>%
  ggplot( aes(x=RS, fill=plot_grps)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',  stat="count") +
  # scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")
p

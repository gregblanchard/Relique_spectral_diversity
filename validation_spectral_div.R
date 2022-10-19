
####################################################################################
####################################################################################
####################################################################################
############## Validation of spectral diversity from field plots  ############## 
####################################################################################
####################################################################################
####################################################################################

####################################################################################
#### Choose and load files and data ####
####################################################################################
#### spatial data forest ####
library(rgdal)
library(sf)
library(rgeos)
library(raster)
# forest 
forest <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.shp")
zone <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/zone_plaine_des_lacs/zone.shp")
# forest in the zone
CP <- as(extent(zone), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(forest))
forest_zone <-  gIntersection(forest, CP, byid=TRUE)
forest_zone$id_patch <- 1:length(forest_zone)

#### get raw spectral bands from S2 image  ####
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021"
raw_spectral <- stack(paste0(file_path,"/T58KFA_raw_int_forest_zone_PDL.tif"))
# rename bands 
names(raw_spectral) <- c("B2",     "B3",     "B4",     "B5",     "B6",     "B7",     "B8",    "B8A",    "B11",    "B12")

#### # CHOOSE and get PCA spectral and spectral diversity indices from plots otained with different parameters (biodivmapR) ####
# see scropts : "sriptR_BiodivMap_58FA_local" or "sriptR_BiodivMap_58FA_server"
## Products from biodimapR pipeline use on server with 100sp and selected PCs: 1,6,7,8, and 5*5 cells spectral communities ##
# PCA from server
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL"
# PCA from same image from local 
pca_spectral <- stack(paste0(file_path,"/SPCA/PCA/OutputPCA_8_PCs"))
# 100 spectral sp. from PC : 1,6,7,8 (imported from server)
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL"
Biodiv_Indicators <- readRDS(paste0(file_path,"/VALIDATIONBiodiv_Indicators.rds"))
# this is the plot file of reference (from which the spectral indices are extracted)
Path_Vector <- paste0(file_path, "/plot_data_ok")

## Products from biodimapR pipeline use on local with 20sp and selected PCs: 1,2,6,7,8, and 5*5 cells spectral communities ##
# # PCA from same image from local
# file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL"
# pca_spectral <- stack(paste0(file_path,"/SPCA/PCA/OutputPCA_8_PCs"))
# # 20 spectral sp. from PC : 1,2,6,7,8 (generated on local)
# file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL"
# Biodiv_Indicators <- readRDS(paste0(file_path,"/SPCA/VALIDATIONBiodiv_Indicators.rds"))
# # this is the plot file of reference (from which the spectral indices are extracted)
# Path_Vector <- paste0(file_path, "/plot_data_ok")

#### get objects and data for field plots ####

# get field diversity indices from plots 
library(rgdal)
# polygons of plots adjusted to match with sentinel reolution (3*3 cells)
plot_data_ok <- readOGR(paste0(Path_Vector, "/geom_plots_PDL_div_indices_utm_ok.shp"))
#### no outlier (feedback from PLS models, this plot is out of forest at the time of S2 image 2021)
plot_data_ok <- readOGR(paste0(Path_Vector, "/geom_plots_PDL_div_indices_utm_ok.shp"))
plot_data_ok <- plot_data_ok[!plot_data_ok$localty == "Forêt Nord 35",]
# foret nord 18 is 1/2 out of forest... 
plot_data_ok <- plot_data_ok[!plot_data_ok$localty == "Forêt Nord 18",]

# polygons of plots with true shape and area
plot_data_true <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data/geom_plots_PDL_div_indices_utm.shp")
# keep same plot set than for S2 images                
plot_data_true <- plot_data_true[plot_data_true$localty %in% plot_data_ok$localty,]

# get data table with all diversity indices from field plots 
plots_PDL_div_indices <- readRDS('/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_PDL_div_indices.rds')

# get data with all diversity distance matrices from field plots 
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

for (i in 1:length(Biodiv_Indicators$BCdiss.All)){
  Biodiv_Indicators$BCdiss.All[[i]] <- as.matrix(Biodiv_Indicators$BCdiss.All[[i]])[match(plot_order, Biodiv_Indicators$localty),match(plot_order, Biodiv_Indicators$localty)]
}

Biodiv_Indicators$FunctionalDiversity <- Biodiv_Indicators$FunctionalDiversity[match(plot_order, Biodiv_Indicators$localty),]

# finaly get the same ID order as all is OK now 
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
Richness <-Biodiv_Indicators$Richness[[1]]
Shannon_RS <- c(Biodiv_Indicators$Shannon)[[1]]
Simpson_RS <- Biodiv_Indicators$Simpson[[1]]
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

#### noice! #### 
plot(plots_PDL_div_indices$CWM_SLA~spectral_betadiv_ok$PC1)
summary(lm(plots_PDL_div_indices$CWM_SLA~spectral_betadiv_ok$PC1))
plot(plots_PDL_div_indices$CWM_trans_SLA~spectral_betadiv_ok$PC1)
summary(lm(plots_PDL_div_indices$CWM_trans_SLA~spectral_betadiv_ok$PC1))

summary(lm(plots_PDL_div_indices$CWM_trans_BA_SLA~spectral_betadiv_ok$PC1))


##############################################################################
###                                                                        ###
#### Relationship between biological components and spectral components   ####
###                                                                        ###
##############################################################################
library(exactextractr)
library(raster)

# #### extract raw and PCA values as functional traits for plots #### 
# ####  use spatial plots with true area #### 
# # Community mean raw values
# plot_CWM_raw_spectral <-  exact_extract(raw_spectral, plot_data_true, "mean")
# colnames(plot_CWM_raw_spectral) <- paste0("S2.",colnames(plot_CWM_raw_spectral))
# # Community variance raw values
# plot_CWV_raw_spectral <-  exact_extract(raw_spectral, plot_data_true, "variance")
# colnames(plot_CWV_raw_spectral) <- paste0("S2.",colnames(plot_CWV_raw_spectral))
# 
# # Community mean PCs values
# plot_CWM_PCA_spectral <-  exact_extract(pca_spectral, plot_data_true, "mean")
# colnames(plot_CWM_PCA_spectral) <- paste0("PCAS2.",colnames(plot_CWM_PCA_spectral))
# # Community variance PCs values
# plot_CWV_PCA_spectral <-  exact_extract(pca_spectral, plot_data_true, "variance")
# colnames(plot_CWV_PCA_spectral) <- paste0("PCAS2.",colnames(plot_CWV_PCA_spectral))
# 

#### use spatial plots addapted to get 3*3 S2 cells #### 
# Community mean raw values
plot_CWM_raw_spectral <-  exact_extract(raw_spectral, plot_data_ok, "mean")
colnames(plot_CWM_raw_spectral) <- paste0("S2.",colnames(plot_CWM_raw_spectral))
# Community mean PCs values
plot_CWV_raw_spectral <-  exact_extract(raw_spectral, plot_data_ok, "variance")
colnames(plot_CWV_raw_spectral) <- paste0("S2.",colnames(plot_CWV_raw_spectral))

# Community mean PCs values
plot_CWM_PCA_spectral <-  exact_extract(pca_spectral, plot_data_ok, "mean")
colnames(plot_CWM_PCA_spectral) <- paste0("PCAS2.",colnames(plot_CWM_PCA_spectral))
# Community variance PCs values
plot_CWV_PCA_spectral <-  exact_extract(pca_spectral, plot_data_ok, "variance")
colnames(plot_CWV_PCA_spectral) <- paste0("PCAS2.",colnames(plot_CWV_PCA_spectral))


#### compil data ####
spectral_data_plots <- cbind(Ric_S2 = Richness,
                             Shannon_RS_S2 = Shannon_RS,
                             Simpson_RS_S2 = Simpson_RS,
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

#### save big table for analysis ####
# plots_field_spectral_trueplots_100sp_formPC1678 <- plots_field_spectral
# saveRDS(plots_field_spectral_trueplots_100sp_formPC1678,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_trueplots_100sp_formPC1678.rds")
# plots_field_spectral_30x30plots_100sp_formPC1678 <- plots_field_spectral
# saveRDS(plots_field_spectral_30x30plots_100sp_formPC1678,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_100sp_formPC1678.rds")
# plots_field_spectral_trueplots_25sp_formPC12678 <- plots_field_spectral
# saveRDS(plots_field_spectral_trueplots_25sp_formPC12678,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_trueplots_25sp_formPC12678.rds")
# plots_field_spectral_30x30plots_25sp_formPC12678 <- plots_field_spectral
# saveRDS(plots_field_spectral_30x30plots_25sp_formPC12678,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_25sp_formPC12678.rds")

#### choose big table for analysis ####
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_trueplots_100sp_formPC1678.rds")
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_100sp_formPC1678.rds")
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_trueplots_25sp_formPC12678.rds")
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_25sp_formPC12678.rds")

#### relationship CWM traits ~ CWM PCA Axis for plots #### 
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_BA_SLA)
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_BA_LA)
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_BA_LDMC)
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_BA_WD)
plot(plot_CWM_PCA_spectral$PCAS2.mean.PC.1 ~ plots_PDL_div_indices$CWM_trans_synth_BA_Axis1)
# noice! 
plot(plots_field_spectral$CWM_SLA~plots_field_spectral$betadiv_PC1_S2)
summary(lm(plots_field_spectral$CWM_SLA ~ plots_field_spectral$betadiv_PC1_S2))

#### multivariate model selection  ####

library(MuMIn)
library(parallel)
#### select variables ####
# choose set of explanatory variable
spectral_div_indices_alpha <- cbind("Ric_S2", "Shannon_RS_S2", "Simpson_RS_S2", "FEve_S2", "FDiv_S2")
spectral_div_indices_beta <-c("betadiv_PC1_S2", "betadiv_PC2_S2", "betadiv_PC3_S2")

expl_v_list <-  list(c(spectral_div_indices_alpha))
expl_v_list <-  list(c(spectral_div_indices_beta))
expl_v_list <- list(colnames(plot_CWM_raw_spectral))
expl_v_list <- list(colnames(plot_CWV_raw_spectral))
expl_v_list <- list(colnames(plot_CWM_PCA_spectral))
expl_v_list <- list(colnames(plot_CWV_PCA_spectral))

expl_v_list_list <- list(  c(spectral_div_indices_alpha),
                           c(spectral_div_indices_beta),
                           colnames(plot_CWM_raw_spectral),
                           colnames(plot_CWV_raw_spectral),
                           colnames(plot_CWM_PCA_spectral),
                           colnames(plot_CWV_PCA_spectral))

# set of response variables
resp_v_list <- colnames(plots_PDL_div_indices)[8:ncol(plots_PDL_div_indices)]

# total set of variable
data_for_mod_all <- plots_field_spectral
#### model selection ####
list_list_model_tab <- list()
list_list_models_best <- list()
list_list_first_best_model <- list()
list_list_model_avg <- list()
list_list_model_full <- list()
for(l in 1:length(expl_v_list_list)){
  expl_v_list <- list(expl_v_list_list[[l]])
  
  list_model_tab <- list()
  list_models_best <- list()
  list_first_best_model <- list()
  list_model_avg <- list()
  list_model_full <- list()
  for(i in 1:length(resp_v_list)){
    data_for_mod <- plots_field_spectral
    print(paste("model selection for model", i, "on", length(resp_v_list), "response variables, for", l, "on", length(expl_v_list_list), "groups of predictors"))
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
  list_list_model_tab[[l]] <- list_model_tab
  list_list_models_best[[l]] <- list_models_best
  list_list_first_best_model[[l]] <- list_first_best_model
  list_list_model_avg[[l]] <- list_model_avg
  list_list_model_full[[l]] <- list_model_full
}

# list_model_tab
# list_first_best_model

##### export results of best models #####

list_coef_best <- list()
for(l in 1:length(list_list_first_best_model)){
  list_first_best_model <- list_list_first_best_model[[l]]
  expl_v_list <- list(expl_v_list_list[[l]])
  list_model_tab <-  list_list_model_tab[[l]]
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
  # order results
  coef_best <- data.frame(coef_best)
  list_coef_best[[l]] <- coef_best[rev(order(coef_best$R2)),]
}

# explore results 
list_coef_best

# export
list_coef_best_trueplots_100sp_formPC1678 <- list_coef_best
# saveRDS(list_coef_best_trueplots_100sp_formPC1678,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_trueplots_100sp_formPC1678.rds") 
list_coef_best_30x30plots_100sp_formPC1678 <- list_coef_best
# saveRDS(list_coef_best_30x30plots_100sp_formPC1678,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_30x30plots_100sp_formPC1678.rds") 
list_coef_best_trueplots_25sp_formPC12678 <- list_coef_best
# saveRDS(list_coef_best_trueplots_25sp_formPC12678,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_trueplots_25sp_formPC12678.rds") 
list_coef_best_30x30plots_25sp_formPC12678 <- list_coef_best
# saveRDS(list_coef_best_30x30plots_25sp_formPC12678,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_30x30plots_25sp_formPC12678.rds") 

# list_coef_best <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_trueplots_100sp_formPC1678.rds")
# list_coef_best <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_30x30plots_100sp_formPC1678.rds")
# list_coef_best <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_trueplots_25sp_formPC12678.rds")
# list_coef_best <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_30x30plots_25sp_formPC12678.rds")

#### this is the best dataset according to linear modle selection ####
# list_coef_best <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/list_coef_best_30x30plots_100sp_formPC1678.rds")
#### IN CONSEQUENCE: choose big table for analysis ####
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_trueplots_100sp_formPC1678.rds")
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_100sp_formPC1678.rds")
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_trueplots_25sp_formPC12678.rds")
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_25sp_formPC12678.rds")


plot(plots_field_spectral$Ric_S2)
#### this is the best dataset according to linear modle selection ####
# plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_100sp_formPC1678.rds")
rownames(plots_field_spectral) <- 1:nrow(plots_field_spectral)
# foret nord 18 is 1/2 out of forest... 
# plots_field_spectral <- plots_field_spectral[!plots_field_spectral$locality == "Forêt Nord 18",]

#### partial least square regression ####
library(pls)
colnames(plots_field_spectral)

# use species diversity VS. bands
model <- plsr(sp_shannon ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

# use species diversity VS. spectral diversity
model <- plsr(sp_shannon ~ Ric_S2 + Shannon_RS_S2 + Simpson_RS_S2 +
                FRic_S2 + FEve_S2 ,
              data=plots_field_spectral, scale=TRUE, validation="CV")

# use one trait or synthertic trait 
model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
              S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

model <- plsr(CWM_SLA ~ PCAS2.mean.PC.1 + PCAS2.mean.PC.2+  PCAS2.mean.PC.3 +
                PCAS2.mean.PC.4 + PCAS2.mean.PC.5 + PCAS2.mean.PC.6 +
                PCAS2.mean.PC.7 + PCAS2.mean.PC.8,
              data=plots_field_spectral, scale=TRUE, validation="CV")

model <- plsr(CWM_CA_SLA ~ PCAS2.mean.PC.1 + PCAS2.mean.PC.2+  PCAS2.mean.PC.3 +
                PCAS2.mean.PC.4 + PCAS2.mean.PC.5 + PCAS2.mean.PC.6 +
                PCAS2.mean.PC.7 + PCAS2.mean.PC.8,
              data=plots_field_spectral, scale=TRUE, validation="CV")


summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

plot(RMSEP(model), legendpos = "topright")

### choose the number of components
plot(model, ncomp = 2, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 3, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 4, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 5, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 6, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 7, asp = 1, line = TRUE) # SLA

#### look for outliers in the model (plot that are not in forest dur to recent landscape changes?) ####
# negative SLA values in model prediction? 
model$fitted.values[model$fitted.values<0]

model$fitted.values[,1,2][model$fitted.values[,1,2]<0]
plots_field_spectral[rownames(plots_field_spectral) == 151,"locality"]

#### update the dataset (no outliers) ####
plots_field_spectral_no_outlier <- plots_field_spectral[!plots_field_spectral$locality == "Forêt Nord 35",]

model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_no_outlier, scale=TRUE, validation="CV")

summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

plot(RMSEP(model), legendpos = "topright")

### choose the number of components
plot(model, ncomp = 2, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 3, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 4, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 5, asp = 1, line = TRUE) # SLA
plot(model, ncomp = 6, asp = 1, line = TRUE) # SLA

R2(model, estimate = "all") 
R2(model, estimate = "train") 

#### subset plots? ####
plots_field_spectral$localty_grp
subset <- plots_field_spectral$localty_grp %in% c("CORIFOR", "Kuebini")
subset <- plots_field_spectral$localty_grp %in% c("Kuebini")

plots_field_spectral_subset <- plots_field_spectral[subset,]
CWM_PCAS2_subset <- CWM_PCAS2[subset,]

model <- plsr(CWM_SLA ~ CWM_PCAS2_subset, data= plots_field_spectral_subset, scale=TRUE, validation="CV")
summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")
plot(RMSEP(model), legendpos = "topright")

### choose the number of components
plot(model, ncomp = 2, asp = 1, line = TRUE) # SLA
R2(model, estimate = "all") 
R2(model, estimate = "train") 

##############################################################################
#### Use the Final PLSR Model to Make Predictions ####
##############################################################################
# get fragmentation metrics and mean environmental variables for 30*30m cells with mean S2 bands and mean spectral PCA values
poly_cells_30m_ok_landscape_env_S2 <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/poly_cells_30m_ok_landscape_env_S2.rds")
poly_cells_30m_ok_landscape_env_S2$log_dist_edge <- log(poly_cells_30m_ok_landscape_env_S2$dist_edge)

colnames(plots_field_spectral)
# reload package 
library(pls)
# taxo div
model <- plsr(sp_shannon ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

model <- plsr(sp_shannon ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12 +
              S2.variance.B2+S2.variance.B3 +            
              S2.variance.B4+S2.variance.B5+S2.variance.B6  +              
            S2.variance.B7+S2.variance.B8+S2.variance.B8A +              
              S2.variance.B11+S2.variance.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

model <- plsr(sp_shannon ~  S2.variance.B2+S2.variance.B3 +            
                S2.variance.B4+S2.variance.B5+S2.variance.B6  +              
                S2.variance.B7+S2.variance.B8+S2.variance.B8A +              
                S2.variance.B11+S2.variance.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

model <- plsr(sp_richness ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

model <- plsr(rar20_taxsha ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")


model <- plsr(PCoA_PC1_taxo ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

# functio div
model <- plsr(FDis ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

model <- plsr(FDis ~  S2.variance.B2+S2.variance.B3 +            
                S2.variance.B4+S2.variance.B5+S2.variance.B6  +              
                S2.variance.B7+S2.variance.B8+S2.variance.B8A +              
                S2.variance.B11+S2.variance.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

model <- plsr(PCoA_PC1_functio ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")


# WD
model <- plsr(CWM_WD ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")


# SLA
model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, validation="CV")

# summary model
summary(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")
plot(RMSEP(model), legendpos = "topright")

names(poly_cells_30m_ok_landscape_env_S2)
landscape_data_utm <- poly_cells_30m_ok_landscape_env_S2[,c("mean.B2","mean.B3","mean.B4",                              
                                                          "mean.B5","mean.B6","mean.B7",                              
                                                           "mean.B8","mean.B8A","mean.B11","mean.B12")]

landscape_data <- data.frame(landscape_data_utm)[,c("mean.B2","mean.B3","mean.B4",                              
                                                    "mean.B5","mean.B6","mean.B7",                              
                                                    "mean.B8","mean.B8A","mean.B11","mean.B12")]
colnames(landscape_data) <- paste0("S2.", colnames(landscape_data))

pcr_pred <- predict(model, landscape_data, ncomp=6)

dim(pcr_pred)

#### get spatial object from predictions ####
landscape_pred_utm <- landscape_data_utm
landscape_pred_utm <- cbind(landscape_pred_utm, pcr_pred)
library(dplyr)
landscape_pred_utm  <- landscape_pred_utm %>% select(CWM_SLA.6.comps)
plot(landscape_pred_utm)

# st_write(landscape_pred_utm, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/landscape_pred_utm.shp")
####  get saved preditions #### 
landscape_pred_utm <- st_read("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/landscape_pred_utm.shp")

######################################################################################################
#### Analyze predictions (extrapolated biological attributes from S2) ~ fragmentation indices ####
######################################################################################################
preds_CWM_SLA <- landscape_pred_utm$CWM_SLA.6.comps
preds_CWM_SLA <- landscape_pred_utm$CWM_SLA
poly_cells_30m_ok_landscape_env_S2$preds_CWM_SLA <- preds_CWM_SLA
poly_cells_30m_ok_landscape_env_S2_df <- data.frame(poly_cells_30m_ok_landscape_env_S2)
#### multivariate model selection  ####

library(MuMIn)
library(parallel)
#### select variables ####
# choose set of explanatory variable
colnames(poly_cells_30m_ok_landscape_env_S2_df)

expl_v_list_list <- list(c( "log_dist_edge", "prop.landscape_100_centroid", "prop.landscape_250_centroid", "prop.landscape_500_centroid",
                            "elevation","slope","aspect",                        
                            "curvature", "twi", "prec" ))

expl_v_list_list <- list(c( "log_dist_edge", "prop.landscape_500_centroid",
                            "elevation","slope",                        
                            "curvature", "twi", "prec" ))

expl_v_list_list <- list(c( "dist_edge", "prop.landscape_500_centroid",
                            "elevation","slope",                        
                            "curvature", "twi", "prec" ))


expl_v_list_list <- list(c( "dist_edge", "edge.density_500_centroid",
                            "elevation","slope",                        
                            "curvature", "twi", "prec" ))


expl_v_list_list <- list(c( "dist_edge", "edge.density_500_centroid", "prop.landscape_500_centroid","effective.mesh.size_250_centroid",
                            "effective.mesh.size_500_centroid", "total.edge_500_centroid", "n.patches_500_centroid", 
                            "elevation","slope",                        
                            "curvature", "twi", "prec" ))

expl_v_list_list <- list(c( "dist_edge", "total.edge_500_centroid", 
                            "elevation","slope",                        
                            "curvature", "twi" ))
expl_v_list_list <- list(c( "log_dist_edge", "total.edge_500_centroid", 
                            "elevation","slope",                        
                            "curvature", "twi" ))

expl_v_list_list <- list(c( "log_dist_edge", "total.edge_500_centroid", "effective.mesh.size_500_centroid",
                            "elevation","slope",                        
                            "curvature", "twi" ))


# set of response variables
resp_v_list <- c("preds_CWM_SLA")


#### remove outliers from SLA prediction (visual chack shows that negative values are on non forest zones) ####
poly_cells_30m_ok_landscape_env_S2_df <- poly_cells_30m_ok_landscape_env_S2_df[poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA>0,]
#### remove cells centroids a less than 15m from the edge ####
poly_cells_30m_ok_landscape_env_S2_df <- poly_cells_30m_ok_landscape_env_S2_df[poly_cells_30m_ok_landscape_env_S2_df$dist_edge>15,]

#### compute perimeter / area ratio for local landscapes ####
poly_cells_30m_ok_landscape_env_S2_df$perimeter.area.ratio_500_centroid <- 
poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid/poly_cells_30m_ok_landscape_env_S2_df$prop.landscape.core_500_centroid
# total set of variable
data_for_mod_all <- poly_cells_30m_ok_landscape_env_S2_df
#### model selection ####
list_list_model_tab <- list()
list_list_models_best <- list()
list_list_first_best_model <- list()
list_list_model_avg <- list()
list_list_model_full <- list()
for(l in 1:length(expl_v_list_list)){
  expl_v_list <- list(expl_v_list_list[[l]])
  
  list_model_tab <- list()
  list_models_best <- list()
  list_first_best_model <- list()
  list_model_avg <- list()
  list_model_full <- list()
  for(i in 1:length(resp_v_list)){
    data_for_mod <- data_for_mod_all
    print(paste("model selection for model", i, "on", length(resp_v_list), "response variables, for", l, "on", length(expl_v_list_list), "groups of predictors"))
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
    # do not keep Variance Inflation Factors > 10 or 5
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
  list_list_model_tab[[l]] <- list_model_tab
  list_list_models_best[[l]] <- list_models_best
  list_list_first_best_model[[l]] <- list_first_best_model
  list_list_model_avg[[l]] <- list_model_avg
  list_list_model_full[[l]] <- list_model_full
}

# list_model_tab
# list_first_best_model

##### export results of best models #####

list_coef_best <- list()
for(l in 1:length(list_list_first_best_model)){
  list_first_best_model <- list_list_first_best_model[[l]]
  expl_v_list <- list(expl_v_list_list[[l]])
  list_model_tab <-  list_list_model_tab[[l]]
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
  # order results
  coef_best <- data.frame(coef_best)
  list_coef_best[[l]] <- coef_best[rev(order(coef_best$R2)),]
}

# explore results 
list_coef_best

# plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$log_dist_edge)
# summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$log_dist_edge))

# plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge)
# summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge))
# 
# summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge))
# plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$elevation)
# summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$elevation))
# summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid))
# plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid)
# 
# plot(poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge)
# 
# summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$total.edge_500_centroid))

#### nice plot coeficient ####
library(GGally)
# choose wich set of restponse and explanatory variable to keep 
list_model_avg <- list_list_model_avg[[1]]
list_model_best <- list_list_first_best_model[[1]]
list_model_tab <-  list_list_model_tab[[1]]
plot_coef_list <- list()
for(i in 1:length(list_model_tab)){
  
  muminavg <- list_model_avg[[i]]
  model_best <- list_model_best[[i]]
  if(!is.na(muminavg)){
    weighted_confint <- confint(muminavg)
    weighted_coef <- coef(muminavg)
  }else{
    weighted_confint <- confint(model_best)
    weighted_coef <- coef(model_best)
  }
  
  
  df_coef <- data.frame(term = names(weighted_coef[2:length(weighted_coef)]),
                        estimate = weighted_coef[2:length(weighted_coef)],
                        conf.low = weighted_confint[2:length(weighted_coef),1],
                        conf.high = weighted_confint[2:length(weighted_coef),2])
  
  grp_term <- ifelse((grepl("edge", df_coef$term, fixed = TRUE) + grepl("centroid", df_coef$term, fixed = TRUE))>0, "frag", "topo")
  df_coef$grp_term <- grp_term
  
  plot_coef_list[[i]] <- ggcoef(df_coef, sort = "ascending", vline_linetype =  "dotted", mapping = aes(x = estimate, y = term, colour = grp_term),
                                errorbar_height = .2, cex = 1) + 
    scale_color_manual('grp_term', labels=c('Positive','Negative'),
                       values=c('forestgreen','darkorange')) +
    theme_bw() + 
    ggtitle(paste(names(list_model_tab[i]), 'from Plaine des Lacs')) +
    theme(legend.position="none",
          plot.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10)) +
    xlab('standardized coefficients') 
}

# change variable names on the plot
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "dist_edge"] <- "Distance to forest edge"
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "log_dist_edge"] <- "Distance to forest edge"
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "curvature"] <- "Curvature"
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "elevation"] <- "Elevation"
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "slope"] <- "Slope"
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "total.edge_500_centroid"] <- "Edge density"
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "twi"] <- "Topographic wetness index"
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "effective.mesh.size_500_centroid"] <- "Effective mesh size"

plot_coef_list[[1]]$labels$title <- "Influence of edge and topography on canopy SLA"
plot_coef_list
## export png ##
# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/plot_coef_pred_SLA.png", width=150, height=70, units = 'mm', res = 300) 
plot_coef_list
dev.off()

################################################################
#### variance partitionning  ####
################################################################
#### use vegan ####
library(vegan)

#### remove outliers from SLA prediction (visual chack shows that negative values are on non forest zones) ####
poly_cells_30m_ok_landscape_env_S2_df <- poly_cells_30m_ok_landscape_env_S2_df[poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA>0,]
#### remove cells centroids a less than 15m from the edge ####
poly_cells_30m_ok_landscape_env_S2_df <- poly_cells_30m_ok_landscape_env_S2_df[poly_cells_30m_ok_landscape_env_S2_df$dist_edge>15,]
#### no na  ####

data_for_varpart <- poly_cells_30m_ok_landscape_env_S2_df[complete.cases(poly_cells_30m_ok_landscape_env_S2_df[,c("elevation","slope",                        
                                         "curvature", "twi",
                                         "log_dist_edge", "total.edge_500_centroid", 
                                         "prop.landscape_500_centroid")]),]

colnames(data_for_varpart)

Topography <- data.frame(data_for_varpart[,c("elevation","slope",                        
                                                                  "curvature", "twi")])

Edge_influence <- data.frame(data_for_varpart[,c("log_dist_edge", "total.edge_500_centroid", "effective.mesh.size_500_centroid")])

# relative importance
relaimpo::calc.relimp(data.frame(data_for_varpart[c("preds_CWM_SLA",
                                                                         "elevation","slope",                        
                                                                         "curvature", "twi",
                                                                         "log_dist_edge", "total.edge_500_centroid",
                                                    "effective.mesh.size_500_centroid")]))
# variance partitioning 
library("rdacca.hp")
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13800
pred_SLA <- data.frame(pred_SLA <- data_for_varpart$preds_CWM_SLA)
hierpar_pred_SLA <- rdacca.hp(pred_SLA, list(Edge_influence, Topography), var.part = T)
hierpar_pred_SLA_all <- rdacca.hp(pred_SLA, cbind(Edge_influence, Topography ), var.part = T)

#### plot varpart ####
col_frag <- grDevices::adjustcolor( "forestgreen", alpha.f = 0.2)
col_topo <- grDevices::adjustcolor( "orange", alpha.f = 0.2)
# get custom venn diagram function
source("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/kuebini/showvarparts_custom.R")
Xnames <- c("Edge\ninfluence", "Topography")

hierpar_pred_SLA$Var.part <- round(hierpar_pred_SLA$Var.part, digits = 2)
labels <- c(paste0("R²=", hierpar_pred_SLA$Var.part[1,1], "         "),
            paste0("\nR²=", ifelse(hierpar_pred_SLA$Var.part[3,1]>0, hierpar_pred_SLA$Var.part[3,1], 0)),
            paste0("          R²=", hierpar_pred_SLA$Var.part[2,1]),
            1-hierpar_pred_SLA$Var.part[4,1])

showvarparts_custom(2, labels = labels, Xnames = Xnames, bg = c(col_frag, col_topo), lty = 0)
## export png ##
# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/vaprpart_pred_SLA.png", width=130, height=100, units = 'mm', res = 300) 
showvarparts_custom(2, labels = labels, Xnames = Xnames, bg = c(col_frag, col_topo), lty = 0)
title("Influence of edge and\ntopography on canopy SLA")
dev.off()

################################################################
#### Random forest  ####
################################################################
library(randomForest)
# select variables 
data_for_mod <- poly_cells_30m_ok_landscape_env_S2_df[,c("preds_CWM_SLA", "dist_edge","total.edge_500_centroid","elevation","slope","curvature","twi")]
# no outliers 
data_for_mod <- data_for_mod[data_for_mod$preds_CWM_SLA>0,]
# no na
data_for_mod <- data_for_mod[complete.cases(data_for_mod),]

data_for_mod_scaled <- as.data.frame(apply(data_for_mod, 2, scale))

rf.fit <- randomForest(preds_CWM_SLA ~ ., data=data_for_mod_scaled, ntree=200,
                     importance=TRUE)

# saveRDS(rf.fit,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/rf_fit.rds")

rf.fit
rf.fit$mse
rf.fit$mse/sqrt(rf.fit$mse)

plot(rf.fit)
### Visualize variable importance ----------------------------------------------

# Get variable importance from the model fit
ImpData <- as.data.frame(importance(rf.fit))
ImpData$Var.Names <- row.names(ImpData)

ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

#find number of trees that produce lowest test MSE
which.min(rf.fit$mse)
#find RMSE of best model
sqrt(rf.fit$mse[which.min(rf.fit$mse)]) 
#### other method to get variable importance ####
i_scores <- caret::varImp(rf.fit, conditional=TRUE, scale = T)
varImpPlot(rf.fit)
# standardize to 100
i_scores_stand <- i_scores/sum(i_scores)*100
i_scores_stand

################################################################
#### BRT ####
################################################################
library(gbm)
library(dismo)
#### change variable names ####
colnames(data_for_mod_scaled)

data_for_mod_brt <- data_for_mod_scaled

#### fit model ####
gbm_mod1 <- gbm.step(data= data_for_mod_brt, gbm.x = 2:ncol(data_for_mod_brt), gbm.y = 1,
                     family = "gaussian", tree.complexity = 5,
                     learning.rate = 0.01, bag.fraction = 0.5)

# saveRDS(gbm_mod1,"/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/gbm_mod1.rds")

summary(gbm_mod1)

gbm.plot(gbm_mod1)
gbm.plot.fits(gbm_mod1)

gbm_mod1_simp <- gbm.simplify(gbm_mod1, n.drops = 5)

find.int <- gbm.interactions(gbm_mod1)
find.int$rank.list

gbm.perspec(gbm_mod1, 5, 4)
gbm.perspec(gbm_mod1, 2, 1)
gbm.perspec(gbm_mod1, 3, 4)
gbm.perspec(gbm_mod1, 3, 5)

#### variance explained ####
#  2 = 1 – (residual deviance/total deviance)

################################################################
#### Moran's I test for residual spatial autocorrelation ####
################################################################
library(spdep)
nb <- knn2nb(knearneigh(st_coordinates(poly_cells_30m_ok_landscape_env_S2)[,1:2], k=2)) 
lstw <- nb2listw(nb, style="S")
moran_tes_res <- lm.morantest(best_model, lstw) #check for raw data
moran_tes_res
moran.mc(best_model$residuals, lstw, 999)

plot(coords_spl_pts)		
plot(nb,coords_spl_pts,add=TRUE)		

# other solution 
library(tidyverse)
library(gridExtra)
library(NLMR)
library(DHARMa)
# formal test
sims <- simulateResiduals(best_model)
testSpatialAutocorrelation(sims, x = coords_spl_pts$lat , y = coords_spl_pts$long, plot = FALSE)


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

#### old plot ####
# bg_image <- levelplot(aerial_photo_ok_r, col.regions=as.character(levels(cols_bg_image)), colorkey=FALSE)
# forest_layer <- latticeExtra::layer(sp.polygons(forest , lwd=2, fill = "forestgreen", alpha = 0.3)) 
# plots_KB <- latticeExtra::layer(sp.points(spdf[spdf$localty_grp == "Kuebini",], pch= 15, cex=1, col= "dodgerblue1", alpha = .8)) 
# plots_CORIFOR <- latticeExtra::layer(sp.points(spdf[spdf$localty_grp == "CORIFOR",], pch= 15, cex=1, col= "coral2", alpha = .8)) 
# plots_ForetNord <- latticeExtra::layer(sp.points(spdf[spdf$localty_grp == "Forêt Nord",], pch= 15, cex=1, col= "darkolivegreen3", alpha = .8)) 
# plots_Wadjana <- latticeExtra::layer(sp.points(spdf[spdf$localty_grp == "Wadjana",], pch= 15, cex=1, col= "magenta", alpha = .8)) 
# 
# plots_map <- bg_image + forest_layer + plots_KB + plots_CORIFOR  + plots_ForetNord + plots_Wadjana 
# plots_map

#### new plot ! ####

library(RStoolbox)
e <- extent(aerial_photo)
plot_image <- ggplot() + ggRGB(img = aerial_photo,
                               r = 1,
                               g = 2,
                               b = 3,
                               stretch = 'none',
                               ggLayer = T,
                               coord_equal = TRUE,
                               alpha =.7) + 
  xlim(e[1:2]) + ylim(e[3:4])
plot_image

#### with legend  ####
cols_fill <- c("Forest"="forestgreen")
cols <- c("Plots"="red4", "Microclimate\nloggers" = "gold2")

#####  add forest map ##### 
library(sf)
# forest_sf <- st_as_sf(forest)
# forest_sf <- st_transform(forest_sf, "+proj=longlat +datum=WGS84 +no_defs" )
# forest_sf <- st_make_valid(forest_sf)
# forest_ok <- st_crop( forest_sf$geometry, e)
# saveRDS(forest_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/forest_ok.rds" )
forest_ok_map <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/forest_ok.rds" )

plot_imag_map <- plot_image  + 
  geom_sf(aes( fill = "Forest"),  data = forest_ok_map,  colour = "grey25" , alpha = 0.1) 

## legend fill ##
plot_imag_map <- plot_imag_map +
  scale_fill_manual(name = " ", values=cols_fill )

##  add plots ##
my_colors <- c( "coral2", "darkolivegreen3","dodgerblue1", "magenta")
plots_ok <-  st_as_sf(spdf)
plot_imag_map <- plot_imag_map  + 
  geom_sf(aes( color = localty_grp), data = plots_ok,  alpha = 1, pch = 15, cex = 1.2) +
  scale_color_manual(values = my_colors, name="") 

##  add theme ###
plot_imag_map  <- plot_imag_map + theme_bw() 
plot_imag_map

# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/plots_map_PDL.png", width=160, height=140, units = 'mm', res = 300) 
plot_imag_map
dev.off()
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
  nbSamples[i] <- length(plot_data_ok)
  plotName[i] <- "PDL"
}

Type_Vegetation = c()
for (i in 1: length(nbSamples)){
  for (j in 1:nbSamples[i]){
    Type_Vegetation = c(Type_Vegetation,Biodiv_Indicators$localty_grp[j])
  }
}

# create data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('Inventory'=Type_Vegetation,'pco1'= spectral_betadiv_ok[,1],'pco2'= spectral_betadiv_ok[,2],'pco3' = spectral_betadiv_ok[,3],
                      'Spectral_Shannon'=Shannon_RS,'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)

# plot field data in the PCoA space, with size corresponding to shannon index
g1 <-ggplot (Results, aes (x=pco1, y=pco2, color=Inventory,size=Spectral_Shannon)) + 
  geom_point(alpha=0.6) + 
  theme_bw() 

g2 <-ggplot (Results, aes (x=pco1, y=pco3, color=Inventory,size=Spectral_Shannon)) + 
  geom_point(alpha=0.6) + 
  theme_bw() 

g3 <-ggplot (Results, aes (x=pco2, y=pco3, color=Inventory,size=Spectral_Shannon)) + 
  geom_point(alpha=0.6) + 
  theme_bw() 

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

# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/PCoA_spectral.png", width=280, height=75, units = 'mm', res = 300) 
grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                         g2 + theme(legend.position="none"),
                         g3 + theme(legend.position="none"),
                         nrow=1),legend,nrow=1, widths =c(5, 1)) 
dev.off()
####################################################################################
#### Based on taxo beta div: produce figures in order to locate the different types of vegetation in the PCoA space ####
####################################################################################
# create data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('Inventory'=Type_Vegetation,'pco1'= taxo_betadiv_BA[,1],'pco2'= taxo_betadiv_BA[,2],'pco3' = taxo_betadiv_BA[,3],
                      'Taxo_Shannon'=plots_PDL_div_indices$sp_shannon,'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)

# plot field data in the PCoA space, with size corresponding to shannon index
g1 <-ggplot (Results, aes (x=pco1, y=pco2, color=Inventory,size=Taxo_Shannon)) + 
  geom_point(alpha=0.6)+ 
  theme_bw()  

g2 <-ggplot (Results, aes (x=pco1, y=pco3, color=Inventory,size=Taxo_Shannon)) + 
  geom_point(alpha=0.6) + 
  theme_bw() 

g3 <-ggplot (Results, aes (x=pco2, y=pco3, color=Inventory,size=Taxo_Shannon)) + 
  geom_point(alpha=0.6) + 
  theme_bw() 

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
# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/PCoA_taxo.png", width=280, height=75, units = 'mm', res = 300) 
grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                         g2 + theme(legend.position="none"),
                         g3 + theme(legend.position="none"),
                         nrow=1),legend,nrow=1, widths =c(5, 1)) 
dev.off()
####################################################################################
#### Based on functio beta div: produce figures in order to locate the different types of vegetation in the PCoA space ####
####################################################################################
# create data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('Inventory'=Type_Vegetation,'pco1'= functio_betadiv_trans_BA[,1],'pco2'= functio_betadiv_trans_BA[,2],'pco3' = functio_betadiv_trans_BA[,3],
                      'Functio_Div'=plot_data_ok$FDiv,'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)

# plot field data in the PCoA space, with size corresponding to shannon index
g1 <-ggplot (Results, aes (x=pco1, y=pco2, color=Inventory,size=Functio_Div)) + 
  geom_point(alpha=0.6) + 
  theme_bw() 

g2 <-ggplot (Results, aes (x=pco1, y=pco3, color=Inventory,size=Functio_Div)) + 
  geom_point(alpha=0.6) + 
  theme_bw() 

g3 <-ggplot (Results, aes (x=pco2, y=pco3, color=Inventory,size=Functio_Div)) + 
  geom_point(alpha=0.6) + 
  theme_bw() 

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
# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/PCoA_functio.png", width=280, height=75, units = 'mm', res = 300) 
grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                         g2 + theme(legend.position="none"),
                         g3 + theme(legend.position="none"),
                         nrow=1),legend,nrow=1, widths =c(5, 1)) 
dev.off()

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
summary(lm(Richness ~ plots_PDL_div_indices$sp_richness))
summary(lm(Richness ~ plots_PDL_div_indices$rar20_taxrich))
summary(lm(Shannon_RS ~ plots_PDL_div_indices$sp_shannon))
summary(lm(Shannon_RS ~ plots_PDL_div_indices$FDis))
summary(lm(Shannon_RS ~ plots_PDL_div_indices$FRic))
summary(lm(Shannon_RS ~ plots_PDL_div_indices$FDiv))
summary(lm(FDiv ~ plots_PDL_div_indices$FDiv))
summary(lm(FRic ~ plots_PDL_div_indices$FRic))
summary(lm(FRic ~ plots_PDL_div_indices$FDis))

###########################################################################
###########################################################################
###########################################################################
###                                                                     ###
####           Influence of fragmentation on spectral diversity        ####
###                                                                     ###
###########################################################################
###########################################################################
###########################################################################


####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  
####  distance to edge from plot position only  #### 
####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  
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
# minimum 10m
dist_edge[dist_edge<10] <- 15
plots_field_spectral$dist_edge <- dist_edge
plots_field_spectral$log_dist_edge <- log(dist_edge)

hist(dist_edge)



####   ####   ####   ####   ####   ####   ####   #### 
#### taxo, funct comp / diversity VS. distance to edge  #### 
####   ####   ####   ####   ####   ####   ####   #### 

summary(lm(plots_field_spectral$sp_richness ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$sp_richness ~ plots_field_spectral$log_dist_edge)

summary(lm(plots_field_spectral$sp_shannon ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$sp_shannon ~ plots_field_spectral$log_dist_edge)

summary(lm(plots_field_spectral$FDis ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$FDis ~ plots_field_spectral$log_dist_edge)

summary(lm(plots_field_spectral$rar20_taxrich ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$rar20_taxrich ~ plots_field_spectral$log_dist_edge)


summary(lm(plots_field_spectral$rar20_taxsha ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$rar20_taxsha ~ plots_field_spectral$log_dist_edge)


summary(lm(plots_field_spectral$CWM_SLA ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$CWM_SLA ~ plots_field_spectral$log_dist_edge)

plot(plots_field_spectral$CWM_SLA ~ plots_field_spectral$dist_edge)

plot(plots_field_spectral$CWM_SLA[plots_field_spectral$localty_grp == "Kuebini"] ~ plots_field_spectral$dist_edge[plots_field_spectral$localty_grp == "Kuebini"])
summary(lm(plots_field_spectral$CWM_SLA[plots_field_spectral$localty_grp == "Kuebini"] ~ plots_field_spectral$log_dist_edge[plots_field_spectral$localty_grp == "Kuebini"]))

#### nice plot ####
# rar sp rich
library(ggplot2)
lm1 <- lm(rar20_taxrich ~ log_dist_edge, data = plots_field_spectral)
smry1 = summary(lm1)
plot1 <- ggplot(data = plots_field_spectral, aes(x = dist_edge, y = rar20_taxrich)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = plots_field_spectral, aes(x = dist_edge, y = rar20_taxrich), formula = y ~ log(x), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x=300, y= 4,
           label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("Rarefied sp. richness (20 trees)") + xlab("Distance to edge (log)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot1

# rar sp rich
lm1 <- lm(FDis ~ log_dist_edge, data = plots_field_spectral)
smry1 = summary(lm1)
plot2 <- ggplot(data = plots_field_spectral, aes(x = dist_edge, y = FDis)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = plots_field_spectral, aes(x = dist_edge, y = FDis), formula = y ~ log(x), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x=300, y= .5,
           label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("Functional dispersion") + xlab("Distance to edge (log)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot2

# CWM SLA
lm1 <- lm(CWM_SLA ~ log_dist_edge, data = plots_field_spectral)
smry1 = summary(lm1)
plot3 <- ggplot(data = plots_field_spectral, aes(x = dist_edge, y = CWM_SLA)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = plots_field_spectral, aes(x = dist_edge, y = CWM_SLA), formula = y ~ log(x), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x=300, y= 5,
           label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("CWM SLA") + xlab("Distance to edge (log)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot3

# CWM WD
lm1 <- lm(CWM_WD ~ log_dist_edge, data = plots_field_spectral)
smry1 = summary(lm1)
plot4 <- ggplot(data = plots_field_spectral, aes(x = dist_edge, y = CWM_WD)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = plots_field_spectral, aes(x = dist_edge, y = CWM_WD), formula = y ~ log(x), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x=300, y= .9,
           label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("CWM WD") + xlab("Distance to edge (log)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot4
# all plots 
library("gridExtra")
grid.arrange(plot1,  plot2, plot3, plot4 ,nrow=2) 

## export png ##
# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/plots_div_vs_dist_edge.png", width=150, height=130, units = 'mm', res = 300) 
grid.arrange(plot1,  plot2, plot3, plot4 ,nrow=2) 
dev.off()

####   ####   ####   ####   ####   ####   ####   #### 
#### spectral diversity VS. distance to edge from plots only  #### 
####   ####   ####   ####   ####   ####   ####   #### 
##### Linear model ##### 
summary(lm(plots_field_spectral$PCAS2.mean.PC.1 ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$PCAS2.mean.PC.1 ~ plots_field_spectral$log_dist_edge)

summary(lm(plots_field_spectral$PCAS2.mean.PC.6 ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$PCAS2.mean.PC.6 ~ plots_field_spectral$log_dist_edge)


summary(lm(plots_field_spectral$betadiv_PC1_S2 ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$betadiv_PC1_S2 ~ plots_field_spectral$log_dist_edge)

summary(lm(plots_field_spectral$Ric_S2 ~ plots_field_spectral$log_dist_edge))
plot(plots_field_spectral$Ric_S2 ~ plots_field_spectral$log_dist_edge)


#### nice plot ####
# alpha div (nb sp spectral sp in 3*3 cells)
library(ggplot2)
lm1 <- lm(Ric_S2 ~ log_dist_edge, data = plots_field_spectral)
smry1 = summary(lm1)
plot1 <- ggplot(data = plots_field_spectral, aes(x = dist_edge, y = Ric_S2)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = plots_field_spectral, aes(x = dist_edge, y = Ric_S2), formula = y ~ log(x), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x=300, y= 3,
           label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("Spectral richness") + xlab("Distance to edge (log)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot1

# beta div (PCoA axis 1)

library(ggplot2)
lm1 <- lm(betadiv_PC1_S2 ~ log_dist_edge, data = plots_field_spectral)
smry1 = summary(lm1)
plot2 <- ggplot(data = plots_field_spectral, aes(x = dist_edge, y = betadiv_PC1_S2)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = plots_field_spectral, aes(x = dist_edge, y = betadiv_PC1_S2), formula = y ~ log(x), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x=300, y= -.3,
           label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("Spectral β-diversity (PCoA 1)") + xlab("Distance to edge (log)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot2

# all plots 
library("gridExtra")
grid.arrange(plot1,  plot2,nrow = 1) 

## export png ##
# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/plots_spectral_vs_dist_edge.png", width=150, height=65, units = 'mm', res = 300) 
grid.arrange(plot1,  plot2,nrow = 1) 
dev.off()


##########################################################################################
###                                                                                    ###
####      Influence of fragmentation on spectral diversity for all cells in the PDL   ####
###                                                                                    ###
##########################################################################################


#### get fragmentation data for all pixels in the PDL ( see script "analyse_frag_betadiv_from_biodivmap_PDL.R") ####
# get fragmentation metrics and mean environmental variables for 50*50m cells 
centroids_cells_ok_landscape_env <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/centroids_cells_ok_landscape_env.rds")
# get fragmentation metrics and mean environmental variables for 30*30m cells with mean S2 bands and mean spectral PCA values
centroids_cells_30m_ok_landscape_env_S2 <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/centroids_cells_30m_ok_landscape_env_S2.rds")
centroids_cells_30m_ok_landscape_env_S2$log_dist_edge <- log(centroids_cells_30m_ok_landscape_env_S2$dist_edge)

#### get aplha and beta spectral diversity for all landscape ####
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL"

shannon_pdl_landscape <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/ALPHA/Shannon_3_MeanFilter_Fullres")
beta_pcoa <- stack("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/BETA/BetaDiversity_BCdiss_PCO_3_Fullres")
 #### aggregate 10*10m resulution raster (S2 original) to get a user define resolution fo compute fragmentation indices ####
# shannon_pdl_landscape.aggregate <- aggregate(shannon_pdl_landscape, fact=3)
# beta_pcoa.aggregate <- aggregate(beta_pcoa, fact=3)

# # saveRDS(shannon_pdl_landscape.aggregate, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/shannon_pdl_landscape_30x30.rds")
# # saveRDS(beta_pcoa.aggregate, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/beta_pcoa_30x30.rds")
shannon_pdl_landscape.aggregate <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/shannon_pdl_landscape_30x30.rds")

beta_pcoa.aggregate <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/beta_pcoa_30x30.rds")

#### get values for each centroid located in forest (use the centroid from the script analyse_frag_betadiv_from_biodivmap_PDL_on_server_new.R) ####
shannon_centroid <- extract( shannon_pdl_landscape.aggregate, centroids_cells_30m_ok_landscape_env_S2)
beta_pcoa_centroid1 <- extract( beta_pcoa.aggregate$PCoA.1, centroids_cells_30m_ok_landscape_env_S2)
beta_pcoa_centroid2 <- extract( beta_pcoa.aggregate$PCoA.2, centroids_cells_30m_ok_landscape_env_S2)
beta_pcoa_centroid3 <- extract( beta_pcoa.aggregate$PCoA.3, centroids_cells_30m_ok_landscape_env_S2)

centroids_cells_30m_ok_landscape_env_S2$spectral_shannon <- shannon_centroid

centroids_cells_30m_ok_landscape_env_S2$spectral_beta_PCoA1 <- beta_pcoa_centroid1
centroids_cells_30m_ok_landscape_env_S2$spectral_beta_PCoA2 <- beta_pcoa_centroid2
centroids_cells_30m_ok_landscape_env_S2$spectral_beta_PCoA3 <- beta_pcoa_centroid3

#### analysis ####

#### multivariate model selection on spectral PCA  ####

library(MuMIn)
library(parallel)
#### select variables ####
# choose set of explanatory variable

names(centroids_cells_30m_ok_landscape_env_S2)

expl_v_list_list <- list(c( "log_dist_edge", "prop.landscape_100_centroid",
                          "elevation","slope",                       
                        "curvature", "twi", ))

expl_v_list_list <- list(c( "log_dist_edge", "total.edge_500_centroid", "effective.mesh.size_500_centroid",
                            "elevation","slope",                       
                            "curvature", "twi", ))
"total.edge_500_centroid", "effective.mesh.size_500_centroid"

# set of response variables
resp_v_list <- c("mean.PC.1", "mean.PC.6","mean.PC.7" , "mean.PC.8", 
                 "spectral_shannon", "spectral_beta_PCoA1", "spectral_beta_PCoA2", "spectral_beta_PCoA3"     )

# total set of variable
data_for_mod_all <- data.frame(centroids_cells_30m_ok_landscape_env_S2)
#### model selection ####
list_list_model_tab <- list()
list_list_models_best <- list()
list_list_first_best_model <- list()
list_list_model_avg <- list()
list_list_model_full <- list()
for(l in 1:length(expl_v_list_list)){
  expl_v_list <- list(expl_v_list_list[[l]])
  
  list_model_tab <- list()
  list_models_best <- list()
  list_first_best_model <- list()
  list_model_avg <- list()
  list_model_full <- list()
  for(i in 1:length(resp_v_list)){
    data_for_mod <- plots_field_spectral
    print(paste("model selection for model", i, "on", length(resp_v_list), "response variables, for", l, "on", length(expl_v_list_list), "groups of predictors"))
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
  list_list_model_tab[[l]] <- list_model_tab
  list_list_models_best[[l]] <- list_models_best
  list_list_first_best_model[[l]] <- list_first_best_model
  list_list_model_avg[[l]] <- list_model_avg
  list_list_model_full[[l]] <- list_model_full
}

# list_model_tab
# list_first_best_model

##### export results of best models #####

list_coef_best <- list()
for(l in 1:length(list_list_first_best_model)){
  list_first_best_model <- list_list_first_best_model[[l]]
  expl_v_list <- list(expl_v_list_list[[l]])
  list_model_tab <-  list_list_model_tab[[l]]
  coef_best <- c()
  for(i in 1: length(list_first_best_model)){
    coef_tmp <- rep(NA, length( expl_v_list[[1]]))
    names(coef_tmp) <- expl_v_list[[1]]
    sm_tmp <- summary(list_first_best_model[[i]])
    if(nrow(sm_tmp$coefficients)>1){
      for(j in 2: nrow(sm_tmp$coefficients)){
        coef_tmp[names(coef_tmp) == rownames(sm_tmp$coefficients)[j]] <- 
          paste0(round(sm_tmp$coefficients[j,1], digits = 2), gtools::stars.pval(sm_tmp$coefficients[j,4]))
      }
    }
    # add R2
    R2 <- round(list_model_tab[[i]]$`R^2`[1], digits = 2)
    coef_tmp <- c(coef_tmp, R2 = R2)
    names(coef_tmp) <- c(expl_v_list[[1]],"R2")
    coef_best <- rbind(coef_best, coef_tmp)
    rownames(coef_best)[i] <- names(list_first_best_model)[i]
  }
  # order results
  coef_best <- data.frame(coef_best)
  list_coef_best <- coef_best
}

# order results
list_coef_best <- list()
for(l in 1:length(list_list_first_best_model)){
  list_coef_best[[l]] <- coef_best[rev(order(coef_best$R2)),]
}

# explore results 
list_coef_best

# export results 
# write.csv2(list_coef_best, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/results_models_spectral_frag_landscape.csv")

##### Linear model ##### 

# PCA from bands

summary(lm(centroids_cells_30m_ok_landscape_env_S2$mean.PC.1 ~ centroids_cells_30m_ok_landscape_env_S2$log_dist_edge))
plot(centroids_cells_30m_ok_landscape_env_S2$mean.PC.1 ~ centroids_cells_30m_ok_landscape_env_S2$log_dist_edge)

plot(centroids_cells_30m_ok_landscape_env_S2$mean.PC.1 ~ centroids_cells_30m_ok_landscape_env_S2$dist_edge)

# alpha div
summary(lm(centroids_cells_30m_ok_landscape_env_S2$spectral_shannon ~ centroids_cells_30m_ok_landscape_env_S2$log_dist_edge))

# beta div
summary(lm(centroids_cells_30m_ok_landscape_env_S2$spectral_beta_PCoA1 ~ centroids_cells_30m_ok_landscape_env_S2$log_dist_edge))

summary(lm(centroids_cells_30m_ok_landscape_env_S2$spectral_beta_PCoA2 ~ centroids_cells_30m_ok_landscape_env_S2$log_dist_edge))

summary(lm(centroids_cells_30m_ok_landscape_env_S2$spectral_beta_PCoA3 ~ centroids_cells_30m_ok_landscape_env_S2$log_dist_edge))

plot(centroids_cells_30m_ok_landscape_env_S2$mean.PC.1 ~ centroids_cells_30m_ok_landscape_env_S2$log_dist_edge)






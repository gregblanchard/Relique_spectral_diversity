 
library(SSDM)
library(raster)

##################################################################################################
#### load data ####
##################################################################################################

merge_all_KB <- readRDS( "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/merge_all_KB.rds")
merge_all_KB_ok <- merge_all_KB[!is.na(merge_all_KB$coordX),]
site_sp_mat_ok <- readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/site_sp_mat_ok.rds")
stack20 <- readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/stack20.rds")

##################################################################################################
#### create sp xy table for occurences ####
##################################################################################################
# occurences no absences
Occ <- data.frame()
for (i in 1:nrow(site_sp_mat_ok)){
  plot_tmp <- rownames(site_sp_mat_ok)[i]
  sp_tmp <- rep(colnames(site_sp_mat_ok), site_sp_mat_ok[i,])
  df_tmp <- cbind(SPECIES = sp_tmp,
                  LONGITUDE = rep(merge_all_KB_ok$coordX[merge_all_KB_ok$id_plot == plot_tmp],length(sp_tmp)),
                  LATITUDE = rep(merge_all_KB_ok$coordY[merge_all_KB_ok$id_plot == plot_tmp],length(sp_tmp)))
  # compil alla in a dataframe
  Occ <- rbind(Occ, df_tmp)
}
Occ$LONGITUDE <- as.numeric(Occ$LONGITUDE)
Occ$LATITUDE <- as.numeric(Occ$LATITUDE)

# occurences WITH absences
Occ_abs <- data.frame()
for (i in 1:nrow(site_sp_mat_ok)){
  plot_tmp <- rownames(site_sp_mat_ok)[i]
  sp_tmp <- rep(colnames(site_sp_mat_ok), site_sp_mat_ok[i,])
  abs_tmp <- colnames(site_sp_mat_ok)[site_sp_mat_ok[i,] == 0 ]
  sp_abs_tmp <- c(sp_tmp, abs_tmp)
  Pcol_tmp <- c(rep(1, length(sp_tmp)), rep(0, length(abs_tmp)))
  df_tmp <- cbind(SPECIES = sp_abs_tmp,
                  LONGITUDE = rep(merge_all_KB_ok$coordX[merge_all_KB_ok$id_plot == plot_tmp],length(sp_abs_tmp)),
                  LATITUDE = rep(merge_all_KB_ok$coordY[merge_all_KB_ok$id_plot == plot_tmp],length(sp_abs_tmp)),
                  PresAbs = Pcol_tmp)
  # compil alla in a dataframe
  Occ_abs <- rbind(Occ_abs, df_tmp)
}
Occ_abs$LONGITUDE <- as.numeric(Occ_abs$LONGITUDE)
Occ_abs$LATITUDE <- as.numeric(Occ_abs$LATITUDE)
Occ_abs$PresAbs <- as.numeric(Occ_abs$PresAbs)


# occurences WITH only presences and absences (not repeated presences for abundances)
site_sp_mat_ok_presabs <- site_sp_mat_ok
site_sp_mat_ok_presabs[site_sp_mat_ok_presabs>0] = 1

Occ_presabs <- data.frame()
for (i in 1:nrow(site_sp_mat_ok_presabs)){
  plot_tmp <- rownames(site_sp_mat_ok_presabs)[i]
  sp_tmp <- rep(colnames(site_sp_mat_ok_presabs), site_sp_mat_ok_presabs[i,])
  abs_tmp <- colnames(site_sp_mat_ok_presabs)[site_sp_mat_ok_presabs[i,] == 0 ]
  sp_abs_tmp <- c(sp_tmp, abs_tmp)
  Pcol_tmp <- c(rep(1, length(sp_tmp)), rep(0, length(abs_tmp)))
  df_tmp <- cbind(SPECIES = sp_abs_tmp,
                  LONGITUDE = rep(merge_all_KB_ok$coordX[merge_all_KB_ok$id_plot == plot_tmp],length(sp_abs_tmp)),
                  LATITUDE = rep(merge_all_KB_ok$coordY[merge_all_KB_ok$id_plot == plot_tmp],length(sp_abs_tmp)),
                  PresAbs = Pcol_tmp)
  # compil alla in a dataframe
  Occ_presabs <- rbind(Occ_presabs, df_tmp)
}
Occ_presabs$LONGITUDE <- as.numeric(Occ_presabs$LONGITUDE)
Occ_presabs$LATITUDE <- as.numeric(Occ_presabs$LATITUDE)
Occ_presabs$PresAbs <- as.numeric(Occ_presabs$PresAbs)
#### save objects ####

# saveRDS(stack20, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/predictors.rds")
# saveRDS(Occ, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/occurences.rds")
# saveRDS(Occ_abs, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/occurences_abs.rds")
# saveRDS(Occ_presabs, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/occurences_presabs.rds")

##################################################################################################
#### Individual species distribution models ####
##################################################################################################
Env <- stack20
SDM <- modelling('GLM', subset(Occ, Occ$SPECIES == 'Gymnostoma deplancheanum'), 
                 Env, Xcol = 'LONGITUDE', Ycol = 'LATITUDE', verbose = FALSE)



plot(SDM@projection, main = 'SDM\nfor Gymnostoma deplancheanum\nwith GLM algorithm')

##################################################################################################
#### Ensemble species distribution models ####
##################################################################################################
Occ <-  readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/occurences.rds")

ESDM <- ensemble_modelling(c('GLM','CTA', 'RF'), subset(Occ, Occ$SPECIES == 'Gymnostoma deplancheanum'), 
                           En = stack20_ok, rep = 10, ensemble.thresh = 0, Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                           Spcol = 'SPECIES', verbose = FALSE, cv.param = c(0.7, 5), 
                           cores = 3, uncertainty=FALSE)
plot(ESDM@projection, main = 'ESDM\nfor Gymnostoma deplancheanum\nwith CTA and MARS algorithms')
# saveRDS(ESDM, "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/ESDM_Gymnp.rds")

##################################################################################################
#### Stacked species distribution models ####
##################################################################################################

SSDM <- stack_modelling(c('GAM', 'SVM'), Occ, Env, rep = 1, ensemble.thresh = 0,
                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                        Spcol = 'SPECIES', method = "pSSDM", verbose = FALSE, cv.param = c(0.7, 1), 
                        cores = 3)

##################################################################################################
#### Stacked species distribution models ON SERVER ####
##################################################################################################

# ssh greg@niamoto.ird.nc

# scp -r  /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/predictors.rds   greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM
# scp -r  /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/occurences.rds   greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM
# scp -r  /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/occurences_abs.rds   greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM
# scp -r  /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/occurences_presabs.rds   greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM

library(SSDM)
library(raster)

stack20 <- readRDS("~/projects/reliques/kuebini/SSDM/predictors.rds")
Occ <- readRDS("~/projects/reliques/kuebini/SSDM/occurences.rds")
Occ_abs <- readRDS("~/projects/reliques/kuebini/SSDM/occurences_abs.rds")
Occ_presabs <- readRDS("~/projects/reliques/kuebini/SSDM/occurences_presabs.rds")

# sel env. variables
cor(as.matrix(stack20)[complete.cases(as.matrix(stack20)),])
stack20_ok <- stack20
# no gap fraction
stack20_ok <- raster::subset(stack20_ok, names(stack20_ok)[!names(stack20_ok) == "gap_fraction_utm"])
# no wetness index
stack20_ok <- raster::subset(stack20_ok, names(stack20_ok)[!names(stack20_ok) == "wetness_index"])
# no elevation
stack20_ok <- raster::subset(stack20_ok, names(stack20_ok)[!names(stack20_ok) == "KubeniASPECT5"])


# SSDM => set uncertainty=FALSE to avoid errors (see https://github.com/sylvainschmitt/SSDM/issues/89)
SSDM <- stack_modelling(algorithms = c('GLM', 'CTA', 'RF'), Occurrences = Occ_abs, Env = stack20_ok, rep = 10, ensemble.thresh = 0,
                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE', Pcol = "PresAbs",
                        Spcol = 'SPECIES', method = "pSSDM", verbose = FALSE, cv.param = c(0.7, 5), 
                        cores = 11, uncertainty=FALSE)

# saveRDS(SSDM, "~/projects/reliques/kuebini/SSDM/SSDM3.rds")
# scp -r  greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM/SSDM1.rds /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM 
# scp -r  greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM/SSDM2.rds /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM 
# scp -r  greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM/SSDM3.rds /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM 

# SSDM => set uncertainty=FALSE to avoid errors (see https://github.com/sylvainschmitt/SSDM/issues/89)
SSDM_presabs <- stack_modelling(algorithms = c('GLM', 'CTA', 'RF'), Occurrences = Occ_presabs, Env = stack20_ok, rep = 10, ensemble.thresh = 0,
                        Xcol = 'LONGITUDE', Ycol = 'LATITUDE', Pcol = "PresAbs",
                        Spcol = 'SPECIES', method = "pSSDM", verbose = T, cv.param = c(0.7, 5), 
                        cores = 11, uncertainty=FALSE)

# saveRDS(SSDM_presabs, "~/projects/reliques/kuebini/SSDM/SSDM2_presabs.rds")

# scp -r  greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM/SSDM1_presabs.rds /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM 
# scp -r  greg@niamoto.ird.nc:~/projects/reliques/kuebini/SSDM/SSDM1_presabs.rds /home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM 



##################################################################################################
#### Model accuracy assessment ####
##################################################################################################

knitr::kable(SSDM@evaluation)
knitr::kable(SSDM@variable.importance)
knitr::kable(SSDM_presabs@evaluation)
knitr::kable(SSDM_presabs@variable.importance)

##################################################################################################
#### Check model output : local ####
##################################################################################################
library(rgdal)

SSDM <- readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/SSDM1.rds")
SSDM <- readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/SSDM1_presabs.rds")
# last selection
SSDM <- readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/SSDM3.rds")
ESDM_Gymno <- readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/SSDM/ESDM_Gymnp.rds")

#### Model accuracy assessment ####
knitr::kable(SSDM@evaluation)
knitr::kable(SSDM@variable.importance)

#### Plot species richness ####

plot(SSDM@diversity.map, main = 'SSDM\nfor Kuebini diversity')

patch_selection <- readOGR("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/GIS/sample_points/patch_selection.shp")
patch_selection <- patch_selection[rev(order(patch_selection$area)),]
patch_selection <- spTransform(patch_selection, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(patch_selection, add = T)

####################################################
#### plot individual species habitat suitability ####
####################################################
library(raster)
#### sp for regeneration ####
plot(SSDM@esdms$`Alphitonia neocaledonica.Ensemble.SDM`@projection)
plot(SSDM@esdms$`Agathis ovata.Ensemble.SDM`@projection)
plot(SSDM@esdms$`Tristaniopsis macphersonii.Ensemble.SDM`@projection)
plot(SSDM@esdms$`Tristaniopsis calobuxus.Ensemble.SDM`@projection)
plot(SSDM@esdms$`Tristaniopsis yateensis.Ensemble.SDM`@projection)
plot(ESDM_Gymno@projection)

#### forest species for enrichment ####
plot(SSDM@esdms$`Calophyllum caledonicum.Ensemble.SDM`)

#### 5 sp clusters ####

res <- readRDS( "/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/cluster_sp_4grp_nmds1.rds")
sp_grps_5 <- res$cluster

# binary maps
dummy_empty_sdm <- ESDM_Gymno@binary
dummy_empty_sdm[] <- 0
SSDM_grp_bin <- list()
for(i in 1:max(sp_grps_5)){
  sp_tmp <- names(sp_grps_5[sp_grps_5 == i])
  esdm_bin_tmp <- list()
  for(j in 1:length(sp_tmp)){
    if(sp_tmp[j] == "Gymnostoma deplancheanum"){ esdm_bin_tmp[[j]] <- ESDM_Gymno@binary
    }else if(length(grep(sp_tmp[j], names(SSDM@esdms)))<1){ esdm_bin_tmp[[j]] <- dummy_empty_sdm
    }else{esdm_bin_tmp[[j]] <- SSDM@esdms[[grep(sp_tmp[j], names(SSDM@esdms))]]@binary}
  } 
  SSDM_grp_bin[[i]] <- Reduce("+", esdm_bin_tmp, accumulate = F)
}
fviz_dend(res, rect = TRUE, cex = 0.5)

plot(SSDM_grp_bin[[1]])
plot(SSDM_grp_bin[[2]])
plot(SSDM_grp_bin[[3]])
plot(SSDM_grp_bin[[4]])

# sum prob maps
SSDM_grp_sumprob <- list()
for(i in 1:max(sp_grps_5)){
  sp_tmp <- names(sp_grps_5[sp_grps_5 == i])
  esdm_sumprob_tmp <- list()
  for(j in 1:length(sp_tmp)){
    if(sp_tmp[j] == "Gymnostoma deplancheanum"){ esdm_sumprob_tmp[[j]] <- ESDM_Gymno@projection
    }else if(length(grep(sp_tmp[j], names(SSDM@esdms)))<1){ esdm_sumprob_tmp[[j]] <- dummy_empty_sdm
    }else{esdm_sumprob_tmp[[j]] <- SSDM@esdms[[grep(sp_tmp[j], names(SSDM@esdms))]]@projection}
  } 
  SSDM_grp_sumprob[[i]] <- Reduce("+", esdm_sumprob_tmp, accumulate = F)
  SSDM_grp_sumprob[[i]] <- SSDM_grp_sumprob[[i]]/max(values(SSDM_grp_sumprob[[i]]), na.rm = T)
}
fviz_dend(res, rect = TRUE, cex = 0.5)
sort(res$cluster)

plot(SSDM_grp_sumprob[[1]])
plot(SSDM_grp_sumprob[[2]])
plot(SSDM_grp_sumprob[[3]])
plot(SSDM_grp_sumprob[[4]])

#### check variable influence for individual sp ####

SSDM@esdms$`Calophyllum caledonicum.Ensemble.SDM`@variable.importance
SSDM@esdms$`Calophyllum caledonicum.Ensemble.SDM`@parameters
SSDM@esdms$`Calophyllum caledonicum.Ensemble.SDM`@sdms



#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
####  Analysis of the signature of forest fragmentation on spectral beta diversity (Sentinel 2)  ####
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################


library(spatial)
library(sp)
library(sf)
library(biodivMapR)
library(reticulate)
library(rgee)
library(googledrive)
library(spatial)
library(rgdal)
library(raster)
library(rgdal)
library(UScensus2010)
library(rgeos)
library(maptools)
library(RSAGA)
library(SDMTools)
library(mapview)
library(spex)
library(exactextractr)




#####################################################################################################
####  Spectral div from NC composite corrected with cloud masking algorythm from GEE ####
#####################################################################################################

path_dir <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all"

####  100sp ####  
# min sunlit 0.8
im_name <- "NC - year - 100sp"
beta_div <- stack(paste0(path_dir,"/RESULTS/S2_NC_int/SPCA/100sp_min80/BETA/BetaDiversity_BCdiss_PCO_10_Fullres"))
alpha_div <- raster(paste0(path_dir,"/RESULTS/S2_NC_int/SPCA/100sp_min80/ALPHA/Shannon_10_Fullres"))
funct_div <- stack(paste0(path_dir,"/RESULTS/S2_NC_int/SPCA/100sp_min80/FUNCTIONAL/FunctionalDiversity_Map_MeanFilter_Fullres"))
# get sampled points 100sp
path_dir_frag <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/results/test/"
centroids_cells_points_wgs <- readRDS(paste0(path_dir_frag,"centroids_cells_100_points_wgs_landscape_metrics_topovar.rds"))
centroids_cells_points_utm <- st_transform(centroids_cells_points_wgs, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

####  25sp #### 
# min sunlit 0.8
im_name <- "NC - year - 25sp"
beta_div <- stack(paste0(path_dir,"/RESULTS/S2_NC_int/SPCA/25sp_min80/BETA/BetaDiversity_BCdiss_PCO_5_Fullres"))
alpha_div <- raster(paste0(path_dir,"/RESULTS/S2_NC_int/SPCA/25sp_min80/ALPHA/Shannon_5_MeanFilter_Fullres"))
funct_div <- stack(paste0(path_dir,"/RESULTS/S2_NC_int/SPCA/25sp_min80/FUNCTIONAL/FunctionalDiversity_Map_Fullres"))
# get sampled points 25sp
path_dir_frag <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/results/test/"
centroids_cells_points_wgs <- readRDS(paste0(path_dir_frag,"centroids_cells_50_points_wgs_landscape_metrics_topovar.rds"))
centroids_cells_points_utm <- st_transform(centroids_cells_points_wgs, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# aggregate spectral div to real resolution (BiodivmapR did not!!!! maybe see the fullres option)
# alpha_div_ok <- aggregate(alpha_div, fact=5)
# writeRaster(alpha_div_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/ALPHA/Shannon_5_MeanFilter_50mres.tif", drivername = "GTiff", datatype="FLT4S", options="COMPRESS=LZW")
# beta_div_ok <- aggregate(beta_div, fact=5)
# writeRaster(beta_div_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/BETA/BetaDiversity_BCdiss_PCO_50mres.tif", drivername = "GTiff", datatype="FLT4S", options="COMPRESS=LZW")
# funct_div_ok <- aggregate(funct_div, fact=5)
# writeRaster(funct_div_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/FUNCTIONAL/FunctionalDiversity_Map_50mres.tif", drivername = "GTiff", datatype="FLT4S", options="COMPRESS=LZW")


#####################################################################################################
####  Spectral div from tile 58KFA corrected with Overland ####
#####################################################################################################

path_dir <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_sudNC_58KFA/58FA"

####  100sp ####  
# min sunlit 0.25
# beta_div <- stack(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/100sp/BETA/BetaDiversity_BCdiss_PCO_10_Fullres"))
# alpha_div <- raster(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/100sp/ALPHA/Shannon_10_Fullres"))
# funct_div <- stack(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/100sp/FUNCTIONAL/FunctionalDiversity_Map_Fullres"))

# min sunlit 0.8
im_name <- "58KFA- Overland - 100sp"
beta_div <- stack(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/100sp/BETA/BetaDiversity_BCdiss_PCO_10_80percent.envi"))
alpha_div <- raster(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/100sp/ALPHA/Shannon_10_Fullres_80percent.envi"))
funct_div <- stack(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/100sp/FUNCTIONAL/FunctionalDiversity_Map_80percent.envi"))
# get sampled points 100sp
path_dir_frag <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/results/test/"
centroids_cells_points_wgs <- readRDS(paste0(path_dir_frag,"centroids_cells_100_points_wgs_landscape_metrics_topovar.rds"))
centroids_cells_points_utm <- st_transform(centroids_cells_points_wgs, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

####  25sp #### 
# min sunlit 0.8
im_name <- "58KFA- Overland - 25sp"
beta_div <- stack(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/25sp/BETA/BetaDiversity_BCdiss_PCO_5_Fullres"))
alpha_div <- raster(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/25sp/ALPHA/Shannon_5_Fullres"))
funct_div <- stack(paste0(path_dir,"/RESULTS/image_dehazed/SPCA/25sp/FUNCTIONAL/FunctionalDiversity_Map_Fullres"))
# get sampled points 25sp
path_dir_frag <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/results/test/"
centroids_cells_points_wgs <- readRDS(paste0(path_dir_frag,"centroids_cells_50_points_wgs_landscape_metrics_topovar.rds"))
centroids_cells_points_utm <- st_transform(centroids_cells_points_wgs, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#####################################################################################################
####  Spectral div from tile 58KFA raw ####
#####################################################################################################

path_dir <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_sudNC_58KFA/58KFA_raw"

####  100sp ####  
# min sunlit 0.8
im_name <- "58KFA - Raw - 100sp"
beta_div <- stack(paste0(path_dir,"/RESULTS/S2_58KFA_int_raw/SPCA/100sp/BETA/BetaDiversity_BCdiss_PCO_10_Fullres"))
alpha_div <- raster(paste0(path_dir,"/RESULTS/S2_58KFA_int_raw/SPCA/100sp/ALPHA/Shannon_10_Fullres"))
funct_div <- stack(paste0(path_dir,"/RESULTS/S2_58KFA_int_raw/SPCA/100sp/FUNCTIONAL/FunctionalDiversity_Map_Fullres"))
# get sampled points 100sp
path_dir_frag <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/results/test/"
centroids_cells_points_wgs <- readRDS(paste0(path_dir_frag,"centroids_cells_100_points_wgs_landscape_metrics_topovar.rds"))
centroids_cells_points_utm <- st_transform(centroids_cells_points_wgs, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

####  25sp #### 
# min sunlit 0.8
im_name <- "58KFA - Raw - 25sp"
beta_div <- stack(paste0(path_dir,"/RESULTS/S2_58KFA_int_raw/SPCA/25sp/BETA/BetaDiversity_BCdiss_PCO_5_Fullres"))
alpha_div <- raster(paste0(path_dir,"/RESULTS/S2_58KFA_int_raw/SPCA/25sp/ALPHA/Shannon_5_Fullres"))
funct_div <- stack(paste0(path_dir,"/RESULTS/S2_58KFA_int_raw/SPCA/25sp/FUNCTIONAL/FunctionalDiversity_Map_Fullres"))
# get sampled points 25sp
path_dir_frag <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/results/test/"
centroids_cells_points_wgs <- readRDS(paste0(path_dir_frag,"centroids_cells_50_points_wgs_landscape_metrics_topovar.rds"))
centroids_cells_points_utm <- st_transform(centroids_cells_points_wgs, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

########################################################################################################################
########################################################################################################################
########################################################################################################################
#### filter sampled points ####
########################################################################################################################
########################################################################################################################
########################################################################################################################
#### get cloud mask for 58KFA corrected in overland (other images already have a cloudfilter) ####
clouds <- terra::rast("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/cloud_mask/T58KFA_cloud_final_mask.tif")

# #### save extent as polygon ####
# e <- as(extent(beta_div), 'SpatialPolygons')
# crs(e) <- "+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# e <- spTransform(e, CRS("+proj=longlat +datum=WGS84 +no_defs"))
# df = data.frame( ID=1:length(e))
# p <- SpatialPolygonsDataFrame(e, df) 
# writeOGR(p, layer = "58KFA_extent", path_dir, driver = "ESRI Shapefile") 


# only points in image extent 
e <- extent(beta_div)
if(as.character(crs(beta_div)) == "+proj=longlat +datum=WGS84 +no_defs") centroids_cells_points_ok <- st_crop(centroids_cells_points_wgs, e)
if(as.character(crs(beta_div)) == "+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs") centroids_cells_points_ok <- st_crop(centroids_cells_points_utm, e)

    
dim(centroids_cells_points_ok)

# extract values
beta_div_spl <- raster::extract(beta_div, centroids_cells_points_ok)
alpha_div_spl <- raster::extract(alpha_div, centroids_cells_points_ok)
funct_div_spl <- raster::extract(funct_div, centroids_cells_points_ok)
cloud_spl <- terra::extract(clouds, terra::vect(centroids_cells_points_ok))

# merge 
data_div_frag <- cbind(centroids_cells_points_ok, beta_div_spl, Shannon = alpha_div_spl, funct_div_spl, no_clouds = cloud_spl)

# no NA in betadiv
data_div_frag <- data_div_frag[complete.cases(beta_div_spl),]
dim(data_div_frag)

# no clouds for 58KFA from overland 
if(grepl("58FA", path_dir, fixed = TRUE)) data_div_frag <- data_div_frag[data_div_frag$no_clouds.cloudmask == 1,]
dim(data_div_frag)

##### arrange data ####
data_div_frag <- as.data.frame(data_div_frag)

#### get long lat data ####
library(tidyverse)
data_div_frag <- data_div_frag %>%
  mutate(long = unlist(map(data_div_frag$geometry,1)),
         lat = unlist(map(data_div_frag$geometry,2)))

#####################################################################################################
#####################################################################################################
#### model selection with glm for each variable (PCA axis, functional div, PCoA axis for betadiv etc...) ####
#####################################################################################################
#####################################################################################################

library(MuMIn)
library(parallel)
##### Choose variable to analyse ####

names_var <- "PCoA1" # for betadiv
names_var <- "PCoA.2" # for betadiv
names_var <- "PCoA.3" # for betadiv

names_var <- "Shannon" # for alphadiv
names_var <- "RICHNESS"
names_var <- "EVENNESS"
names_var <- "DIVERGENCE" # for funct div
# names_var <- "RICHNESS" # for funct div

##### function for simple model selection for landscape buffers of different sizes ####
msel_aic_uni = function(Y,X, empty = F){
  aic = c()
  for(i in 1:ncol(X)){
    lm = lm (Y~X[,i])
    aic = c(aic,AIC(lm))
    if(empty & i == ncol(X)) aic = c(aic, AIC(lm(Y~1)))
  }
  if(empty) tab = cbind(c(colnames(X), "intercept"), aic) else  tab = cbind(colnames(X), aic)
  return(list(pred = colnames(X)[order(aic)][1], tab = tab))
}  

##### get log-transformed distance to edge ####
data_div_frag$log_edge <- log(data_div_frag$dist_edge_centroid+1)

####  remove outliers  ####  
if(names_var == "Shannon")  probs=c(.01, 1)
if(names_var == "RICHNESS")  probs=c(0, 0.99)
if(names_var == "DIVERGENCE")  probs=c(0, 0.99)

Q <- quantile(data_div_frag[,names_var], probs=probs, na.rm = T)
data_div_frag_ok <- data_div_frag[data_div_frag[,names_var] > Q[1] & data_div_frag[,names_var] < Q[2] ,]

#### minimum distance from edge for points (10: cells are not a the edge) ####
####  no points < d m from edge #### 
# for 50m cells with 25sp
dmin = 35
dmax = 100

# for 100m cells with 100sp

dmin = 70
dmax = 300

dmax = 1000000

data_div_frag_ok <- data_div_frag_ok[data_div_frag_ok$dist_edge_centroid > dmin & data_div_frag_ok$dist_edge_centroid < dmax ,]
nrow(data_div_frag_ok)

#### get coordinates ####
data_div_frag_ok <- data_div_frag_ok[!is.na(data_div_frag_ok$long),]
coords_spl_pts <- data_div_frag_ok[,c("lat", "long")]

#### for loop : model selection for each pca axis ####
select_col <- grep(names_var, colnames(data_div_frag_ok))
expl_col <- c(1:ncol(data_div_frag_ok))[!grepl(names_var, colnames(data_div_frag_ok))]
if(length(names_var)>1 & !data_div_frag_ok == "PCoA") expl_col <- c(1:ncol(data_div_frag_ok))[!colnames(data_div_frag_ok) %in% names_var]
if(length(names_var)>1 & !names_var == "PCoA") select_col <- c(1:ncol(data_div_frag_ok))[colnames(data_div_frag_ok) %in% names_var]

#### model selection ####
list_model_tab <- list()
list_model_best <- list()
list_model_avg <- list()
for(i in 1:length(select_col)){
  print(paste("model selection for response variable", i, "on", length(select_col)))
  var_tmp <- colnames(data_div_frag_ok)[select_col[i]]
  data_tmp <- data_div_frag_ok[,c(select_col[i], expl_col)]
  
  # distance edge (log of standard)
  aic_sel_edge <- msel_aic_uni(data_tmp[,1], cbind(dist_edge = data_tmp$dist_edge_centroid, log_edge = data_tmp$log_edge))
  
  # habitat amount
  aic_sel_HA <- msel_aic_uni(data_tmp[,1], data_tmp[,c( paste0(rep("prop.landscape_", 4), c(250,500,1000,2000), "_centroid"))])
  aic_sel_HA

  # effective mesh size
  aic_sel_EMS <- msel_aic_uni(data_tmp[,1], data_tmp[,c( paste0(rep("effective.mesh.size_", 4), c(250,500,1000,2000), "_centroid"))])
  aic_sel_EMS
  
  colnames(data_tmp)
  land_var <- paste0(rep(c("total.edge_", "effective.mesh.size_", "prop.landscape_"), each = 4), c(250,500,1000,2000), "_centroid")
  # from selection
  land_var <- c(aic_sel_HA$pred, aic_sel_EMS$pred)
  land_var <- c(aic_sel_HA$pred)
  
  topo_var <- c("elevation","slope","aspect" ,"curvature","substrat_centroid")
  sel_var <- c(aic_sel_edge$pred, "area",land_var, topo_var)
  sel_var <- c("log_edge", land_var, topo_var) # need to log distance to edge !!! (see distribution)
  # no NA values
  data_for_mod <- data_tmp[, c(var_tmp, sel_var)]
  # no NA values
  complet_cases <- complete.cases(data_for_mod) 
  data_for_mod <- data_for_mod[complet_cases,]
  coords_spl_pts <- coords_spl_pts[complet_cases,]
  coordinates(coords_spl_pts) <- coords_spl_pts
  
  #### final number of points ####
  dim(data_for_mod)
  #### standardize variables ####
  scale_01 <- function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
 
  data_for_mod_toscale <- data_for_mod
  data_for_mod_toscale$substrat_centroid <- as.numeric(as.factor(data_for_mod_toscale$substrat_centroid))
  data_for_mod_scaled <- as.data.frame(apply(data_for_mod_toscale, 2, scale))

  #### model selection for each PCA axis ####
  # with gaussiian family distribution 
  # do not keep Variance Inflation Factors > 5
  
  glm1 <- glm( formula(paste0(var_tmp,  " ~ .")), data = data_for_mod_scaled, na.action = "na.fail")
  vif <- car::vif(glm1)
  data_for_mod_scaled_ok <- data_for_mod_scaled
  while(max(vif)>5){
    data_for_mod_scaled_ok <- data_for_mod_scaled_ok[, !colnames(data_for_mod_scaled_ok) %in% names(vif)[vif == max(vif)]]
    glm1 <- glm(formula(paste0(var_tmp,  " ~ .")), data = data_for_mod_scaled_ok, na.action = "na.fail")
    vif <- car::vif(glm1)
  }
  vif
  glm1 <- glm(formula(paste0(var_tmp,  " ~ aspect + I(aspect^2) + .")), data = data_for_mod_scaled_ok, na.action = "na.fail")
  
  rm(dd)
  cl <- makeCluster(3)
  clusterExport(cl, "data_for_mod_scaled_ok")
  dd <- pdredge(glm1,cluster=cl,rank = "AIC", beta = "partial.sd",
                extra = list(
                  "R^2", "*" = function(x) {
                    s <- summary(x)
                    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
                      F = s$fstatistic[[1]])
                  })
  )
  stopCluster(cl)
  dd
  list_model_tab[[i]] <- dd 
  # selection few best models based on delta AIC
  rm(models.list)
  models.list <- get.models(dd,subset =  delta < 2) 
  print(paste(i, length(models.list)))
  list_model_best[[i]] <- models.list
  # average models and coeficiants based on weigth
  if(length(models.list) > 1){   muminavg <- model.avg(models.list, revised.var = TRUE, fit = TRUE, beta = "partial.sd")
  list_model_avg[[i]] <- muminavg
  }else{list_model_avg[[i]] <- NA}
}

names(list_model_tab) =  names(list_model_best) =  names(list_model_avg) <- colnames(data_div_frag)[select_col]
list_model_tab
#### plot results ####
library(GGally)
plot_coef_list <- list()
for(i in 1:length(list_model_tab)){
  
  muminavg <- list_model_avg[[i]]
  model_best <- list_model_best[[i]]
  if(!is.na(muminavg)){
    weighted_confint <- confint(muminavg)
    weighted_coef <- coef(muminavg)
  }else{
    weighted_confint <- confint(model_best[[1]])
    weighted_coef <- coef(model_best[[1]])
  }
  
  df_coef <- data.frame(term = names(weighted_coef[2:length(weighted_coef)]),
                        estimate = weighted_coef[2:length(weighted_coef)],
                        conf.low = weighted_confint[2:length(weighted_coef),1],
                        conf.high = weighted_confint[2:length(weighted_coef),2])
  
  plot_coef_list[[i]] <- ggcoef(df_coef, sort = "ascending", vline_linetype =  "dotted", mapping = aes(x = estimate, y = term, colour = term),
                                errorbar_height = .1) + 
    theme_bw() + 
    ggtitle(paste(names(list_model_tab[i]), 'from', im_name)) +
    theme(legend.position="none",
          plot.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15)) +
    xlab('standardized coefficients') 
}

# change term names
  plot_coef_list[[1]]$data$term[plot_coef_list[[1]]$data$term == "I(aspect^2)"] <- "aspect^2"
  
  levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "substrat_centroid"] <- "UM substrat"
  levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "prop.landscape_250_centroid"] <- "Habitat amount (250m)"
  levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "prop.landscape_2000_centroid"] <- "Habitat amount (2000m)"
  
  plot_coef_list
  
    
best_model <- list_model_best[[1]][[1]]
summary(best_model)
# plot(best_model)
shapiro.test(best_model$residuals[seq(1,length(best_model$residuals), 5)])

################################################################
#### Moran's I test for residual spatial autocorrelation ####
################################################################
library(spdep)

nb <- knn2nb(knearneigh(coords_spl_pts, k=2)) 
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

################################################################
#### nice plots Betadiv   ####
################################################################

### plot with H.A.  ####

data_plot1 <- data.frame(Betadiv = data_for_mod$PCoA.1, HA= data_for_mod$prop.landscape_2000)

data_plot1 <- data.frame(Betadiv = data_for_mod_scaled_ok$PCoA.1, HA= data_for_mod$prop.landscape_2000)

lm1 <- lm(Betadiv ~ HA, data = data_plot1)
smry1 = summary(lm1)
smry1
plot <- ggplot(data = data_plot1, aes(x = HA, y = Betadiv)) + 
   # geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = HA, y = Betadiv), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  # annotate("text", x=.9, y= 1,
  #          label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
  #          color='red') +
  ylab("Betadiv") + xlab("H.A.") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

### plot with residuals from topo  ####

lmtopo <- lm(PCoA.1 ~ slope + aspect + substrat_centroid, data = data_for_mod)
smrytopo = summary(lmtopo)
smrytopo
data_plot1 <- data.frame(ResBeta = lmtopo$residuals, HA= data_for_mod$prop.landscape_2000)

lm1 <- lm(ResBeta ~ HA, data = data_plot1)
smry1 = summary(lm1)
smry1
plot <- ggplot(data = data_plot1, aes(x = HA, y = ResBeta)) + 
  # geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = HA, y = ResBeta), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  # annotate("text", x=.9, y= 1,
  #          label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
  #          color='red') +
  ylab("ResShannon") + xlab("H.A.") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

#### 2 variables ####
data_mod <- data.frame(Betadiv = data_for_mod$PCoA.1, HA= data_for_mod$prop.landscape_2000, log_edge= data_for_mod$log_edge)
lm1 <- lm(Betadiv ~ HA + log_edge, data = data_mod)
summary(lm1)

#### plot with distance to edge  ####

data_plot1 <- data.frame(Betadiv = data_for_mod$PCoA.1, log_edge= data_for_mod$log_edge)

data_plot1 <- data.frame(Betadiv = data_for_mod_scaled_ok$PCoA.1, log_edge= data_for_mod$log_edge)

lm1 <- lm(Betadiv ~ log_edge, data = data_plot1)
smry1 = summary(lm1)
plot <- ggplot(data = data_plot1, aes(x = log_edge, y = Betadiv)) + 
   geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = log_edge, y = Betadiv), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  # annotate("text", x=.9, y= 1,
  #          label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
  #          color='red') +
  ylab("Betadiv") + xlab("log_edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

### plot with residuals from topo  ####

lmtopo <- lm(PCoA.1 ~ slope + aspect + substrat_centroid, data = data_for_mod)
smrytopo = summary(lmtopo)
smrytopo
data_plot1 <- data.frame(ResBeta = lmtopo$residuals, log_edge= data_for_mod$log_edge)

lm1 <- lm(ResBeta ~ log_edge, data = data_plot1)
smry1 = summary(lm1)
smry1
plot <- ggplot(data = data_plot1, aes(x = log_edge, y = ResBeta)) + 
  # geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = log_edge, y = ResBeta), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  # annotate("text", x=.9, y= 1,
  #          label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
  #          color='red') +
  ylab("ResBeta") + xlab("log_edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot


################################################################
#### nice plots Shannon   ####
################################################################

### plot with HA  ####
data_mod1 <- data.frame(Shannon = data_for_mod$Shannon, HA= log(data_for_mod$prop.landscape_250_centroid))

data_mod1 <- data.frame(Shannon = data_for_mod_scaled_ok$Shannon, HA= data_for_mod$prop.landscape_250_centroid)

lm1 <- lm(Shannon ~ HA, data = data_mod1)
smry1 = summary(lm1)
smry1

#### plot model ####
data_plot1 =  data.frame(Shannon = data_for_mod$Shannon, HA= data_for_mod$prop.landscape_250_centroid)
#data_plot1 <- data.frame(FRic = data$FRic, log_edge= data$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = HA, y = Shannon)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = HA, y = Shannon), formula = y ~ log(x), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
   annotate("text", x=.5, y= 2.9,
            label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
            color='red') +
  ylab("Shannon") + xlab("H.A.") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

### plot with distance to edge  ####
data_mod1 <- data.frame(Shannon = data_for_mod$Shannon, dist_edge= data_for_mod$log_edge)

data_mod1 <- data.frame(Shannon = data_for_mod_scaled_ok$Shannon, dist_edge= data_for_mod$log_edge)

cor.test(data_mod1$Shannon,data_mod1$dist_edge, method = "kendall")

lm1 <- lm(Shannon ~ dist_edge, data = data_mod1)
smry1 = summary(lm1)
smry1
#### plot model ####
data_plot1 <- data.frame(Shannon = data_for_mod$Shannon, dist_edge= exp(data_for_mod$log_edge))
#data_plot1 <- data.frame(FRic = data$FRic, log_edge= data$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = Shannon)) + 
   geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = Shannon), formula = y ~ log(x), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x= 50, y= 2.9,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("Shannon") + xlab("Distance to edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

#### 2 variables ####
data_mod2 <- data.frame(Shannon = data_for_mod$Shannon, HA= data_for_mod$prop.landscape_250_centroid, log_edge= data_for_mod$log_edge)
lm1 <- lm(Shannon ~ HA + log_edge, data = data_mod2)
summary(lm1)

### plot with residuals from topo  ####

lmtopo <- lm(Shannon ~ curvature + aspect + substrat_centroid, data = data_for_mod)
smrytopo = summary(lmtopo)
smrytopo
data_plot1 <- data.frame(ResShannon = lmtopo$residuals, dist_edge= data_for_mod$log_edge, HA= data_for_mod$prop.landscape_250_centroid)

lm1 <- lm(ResShannon ~ dist_edge, data = data_plot1)
smry1 = summary(lm1)
smry1

lm2 <- lm(ResShannon ~ dist_edge + HA, data = data_plot1)
smry2 = summary(lm2)
smry2

#data_plot1 <- data.frame(FRic = data$FRic, log_edge= data$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = ResShannon)) + 
   geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = ResShannon), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
   annotate("text", x=4, y= 0.5,
            label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
            color='red') +
  ylab("ResShannon") + xlab("Distance to edge (log)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

### plot with aspect ####
data_mod1 <- data.frame(Shannon = data_for_mod$Shannon, aspect= data_for_mod$aspect)

lm1 <- lm(Shannon ~ aspect + I(aspect^2), data = data_mod1)
smry1 = summary(lm1)
smry1
#### plot model ####
data_plot1 <- data.frame(Shannon = data_for_mod$Shannon, aspect= data_for_mod$aspect)
#data_plot1 <- data.frame(FRic = data$FRic, log_edge= data$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = aspect, y = Shannon)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = aspect, y = Shannon), formula = y ~ x + I(x^2), method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x= 50, y= 2.9,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[11])),
           color='red') +
  ylab("Shannon") + xlab("Aspect (degree)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot


### plot with Substrat ####
data_mod1 <- data.frame(Shannon = data_for_mod$Shannon, substrat= data_for_mod$substrat_centroid)
data_mod1$substrat <- as.factor((ifelse(data_mod1$substrat == "non_UM","VS", "UM")))
lm1 <- lm(Shannon ~ substrat, data = data_mod1)
smry1 = summary(lm1)
smry1
#### plot model ####
data_plot1 <- data_mod1
#data_plot1 <- data.frame(FRic = data$FRic, log_edge= data$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = substrat, y = Shannon)) + 
  geom_point(color='navy', size = 1, alpha =.3) +

    # annotate("text", x= 50, y= 2.9,
    #        label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[11])),
    #        color='red') +
  ylab("Shannon") + xlab("Aspect (degree)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot


################################################################
#### nice plots Funct div   ####
################################################################

### plot with distance to edge  ####

data_plot1 <- data.frame(Functdiv = data_for_mod$RICHNESS, dist_edge = data_for_mod$log_edge)

data_plot1 <- data.frame(Functdiv = data_for_mod_scaled_ok$RICHNESS, dist_edge = data_for_mod$log_edge)

lm1 <- lm(Functdiv ~ dist_edge, data = data_plot1)
smry1 = summary(lm1)
plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = Functdiv)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = Functdiv), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x=.9, y= 1,
           label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("Functdiv") + xlab("Distance to edge (log)") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

### plot with residuals from topo  ####
data_plot1 <- data.frame(Functdiv = data_for_mod$RICHNESS, dist_edge = data_for_mod$log_edge)

lmtopo <- lm(RICHNESS ~ slope + aspect + substrat_centroid , data = data_for_mod)
smrytopo = summary(lmtopo)
smrytopo
data_plot1 <- data.frame(ResFunct = lmtopo$residuals, dist_edge= data_for_mod$log_edge)

lm1 <- lm(ResFunct ~ dist_edge, data = data_plot1)
smry1 = summary(lm1)
smry1
plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = ResFunct)) + 
  geom_point(color='navy', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = ResFunct), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
  annotate("text", x=.9, y= 1,
           label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
           color='red') +
  ylab("ResShannon") + xlab("Distance to edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

############################################################################################################
############################################################################################################
############################################################################################################
#### analyse per altitudinal bands ####
############################################################################################################
############################################################################################################
############################################################################################################
library(tidyverse)
library(gapminder)
library(ggplot2)

names(data_div_frag)
spdf <- SpatialPointsDataFrame(coords = data_div_frag[,c("long","lat")], data = data_div_frag,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# add altitudinal bands to data.frame of sampluing values
dem_class <- raster("/home/thesardfou/Documents/GIS/MNT/100mbands_100m.tif")

elev_band <- raster::extract(dem_class, spdf)
data_div_frag$elev_band <- elev_band
data_div_frag_to_elev <- data_div_frag[!is.na(data_div_frag$elev_band),]

# number of point per elev
nb_pts_elev <- c()
for(i in 1:max(data_div_frag_to_elev$elev_band)) nb_pts_elev <- c(nb_pts_elev,sum(data_div_frag_to_elev$elev_band == i)) 

data_div_frag_to_elev %>% 
  ggplot(aes(x=elev_band,y=Shannon, fill=substrat_centroid)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) 


data_div_frag_to_elev %>% 
  filter(substrat_centroid %in% unique(data_div_frag$substrat_centroid)) %>%
  ggplot(aes(x=factor(elev_band),y=Shannon, fill=factor(substrat_centroid))) +
  geom_point(position=position_jitterdodge(),alpha=0.03) +
  geom_boxplot() + 
  labs(fill = "substrat_centroid") + 
  theme_bw(base_size = 16)

#### altitudinal bands ####

# Generate a reclass raster (begin, end, value ...)
# dem <- raster("/home/thesardfou/Documents/GIS/MNT/mnt501.tif")
# rclmat = matrix( c(seq(0,1600, 100), seq(100,1700, 100), seq(1,17)), ncol=3, byrow=F)
reclassified_dem = reclassify(dem, rclmat)
# writeRaster(reclassified_dem, "/home/thesardfou/Documents/GIS/MNT/100mbands_100m.tif", drivername = "GTiff", datatype="INT2U", options="COMPRESS=LZW")
dem_class <- raster("/home/thesardfou/Documents/GIS/MNT/100mbands_100m.tif")

# aggregate spectral div to real resolution (BiodivmapR did not!!!!)
# alpha_div_ok <- aggregate(alpha_div, fact=5)
# writeRaster(alpha_div_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/ALPHA/Shannon_5_MeanFilter_50mres.tif", drivername = "GTiff", datatype="FLT4S", options="COMPRESS=LZW")
alpha_div_ok <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/ALPHA/Shannon_5_MeanFilter_50mres.tif")
beta_div_ok <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/BETA/BetaDiversity_BCdiss_PCO_50mres.tif")
funct_div_ok <- stack("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/FUNCTIONAL/FunctionalDiversity_Map_50mres.tif")

# get elevation raster
# resample DEM to match with response layers
# e <- extent(alpha_div_ok)
# dem <- crop(dem, e)
# dem_ok <- raster::resample(dem, alpha_div_ok, method="bilinear")
# writeRaster(dem_ok, "/home/thesardfou/Documents/GIS/MNT/dem_ok_50m.tif", drivername = "GTiff", datatype="INT2U", options="COMPRESS=LZW", overwrite=T)
dem_ok <- raster("/home/thesardfou/Documents/GIS/MNT/dem_ok_50m.tif")

# get elevation classes raster
# resample classified DEM to match with response layers
# e <- extent(alpha_div_ok)
# dem_class <- crop(dem_class, e)
# dem_class_ok <- raster::resample(dem_class, alpha_div_ok, method="bilinear")
# writeRaster(dem_class_ok, "/home/thesardfou/Documents/GIS/MNT/100mbands_50m.tif", drivername = "GTiff", datatype="INT2U", options="COMPRESS=LZW", overwrite=T)
dem_class_ok <- raster("/home/thesardfou/Documents/GIS/MNT/100mbands_50m.tif")

# add substrat raster
# UM_raster <- raster("/home/thesardfou/Documents/GIS/divers/NC_GT_UM_areas.tif")
# e <- extent(alpha_div_ok)
# UM_raster <- crop(UM_raster, e)
# UM_raster_ok <- raster::resample(UM_raster, alpha_div_ok, method="bilinear")
# writeRaster(UM_raster_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/UM_raster_ok_50m.tif", drivername = "GTiff", datatype="INT2U", options="COMPRESS=LZW", overwrite=T)
UM_raster_ok <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/UM_raster_ok_50m.tif")

# add edge distance raster
# edge_dist_raster_utm <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/distance_to_forest_edge_raster_utm.tif")
# edge_dist_raster <- projectRaster(edge_dist_raster_utm, crs = "+proj=longlat +datum=WGS84 +no_defs")
# e <- extent(alpha_div_ok)
# edge_dist_raster <- crop(edge_dist_raster, e)
# edge_dist_raster_ok <- raster::resample(edge_dist_raster, alpha_div_ok, method="bilinear")
# writeRaster(edge_dist_raster_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/edge_dist_raster_ok.tif", drivername = "GTiff", datatype="INT2U", options="COMPRESS=LZW", overwrite=T)
edge_dist_raster_ok <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/edge_dist_raster_ok.tif")

# add SSDM raster
# ssdm_raster <- raster("/home/thesardfou/Documents/GIS/NC_species_diversity_ssdm/a.tif")
# e <- extent(alpha_div_ok)
# ssdm_raster <- crop(ssdm_raster, e)
# ssdm_raster_ok <- raster::resample(ssdm_raster, alpha_div_ok, method="bilinear")
# writeRaster(ssdm_raster_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/ssdm_raster_ok.tif",
#             drivername = "GTiff", datatype="INT2U", options="COMPRESS=LZW", overwrite=T)
 ssdm_raster_ok <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/ssdm_raster_ok.tif")

 # add slope raster
 # slope <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/slope_GT.tif")
 # slope <- projectRaster(slope, crs = "+proj=longlat +datum=WGS84 +no_defs")
 # e <- extent(alpha_div_ok)
 # slope <- crop(slope, e)
 # slope_ok <- raster::resample(slope, alpha_div_ok, method="bilinear")
 # writeRaster(slope_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/slope_ok.tif",
 #             drivername = "GTiff", datatype="INT2U", options="COMPRESS=LZW", overwrite=T)
 slope_ok <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/slope_ok.tif")
 
 # add aspect raster
 # aspect <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/aspect_GT.tif")
 # aspect <- projectRaster(aspect, crs = "+proj=longlat +datum=WGS84 +no_defs")
 # e <- extent(alpha_div_ok)
 # aspect <- crop(aspect, e)
 # aspect_ok <- raster::resample(aspect, alpha_div_ok, method="bilinear")
 # writeRaster(aspect_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/aspect_ok.tif",
 #             drivername = "GTiff", datatype="INT2U", options="COMPRESS=LZW", overwrite=T)
 aspect_ok <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/aspect_ok.tif")
 
 # add curvature raster
 # curv <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/curvature_GT.tif")
 # curv <- projectRaster(curv, crs = "+proj=longlat +datum=WGS84 +no_defs")
 # e <- extent(alpha_div_ok)
 # curv <- crop(curv, e)
 # curv_ok <- raster::resample(curv, alpha_div_ok, method="bilinear")
 # writeRaster(curv_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/curv_ok.tif",
 # drivername = "GTiff", datatype="FLT4S", options="COMPRESS=LZW", overwrite=T)
 curv_ok <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_NC_all/RESULTS/S2_NC_int/SPCA/25sp_min80/curv_ok.tif")
 
 #### create a table with all data  ####
 div_um_edge_stack <- stack(alpha_div_ok,
                            beta_div_ok,
                            funct_div_ok,
                            UM_raster_ok,
                            edge_dist_raster_ok,
                            dem_ok, 
                            ssdm_raster_ok,
                            slope_ok,
                            aspect_ok,
                            curv_ok)

 
 vals <- getValues(div_um_edge_stack)
 div_um_edge_stack_df <- vals[!is.na(vals[,1]),]
 div_um_edge_stack_df <- data.frame(div_um_edge_stack_df)
 colnames(div_um_edge_stack_df) <- c("alpha", "beta",
                                     "funct_rich", "funct_eveness", "funct_div",
                                     "substrat", "dist_edge", "elev", "ssdm",
                                     "slope", "aspect", "curvature")
 div_um_edge_stack_df$substrat[is.na(div_um_edge_stack_df$substrat)] <- "VS"
 div_um_edge_stack_df$substrat[div_um_edge_stack_df$substrat == "1"] <- "UM"

 div_um_edge_stack_df[1:100,]
 
 div_um_edge_stack_df_ok <- div_um_edge_stack_df[!is.na(div_um_edge_stack_df$dist_edge),]
 div_um_edge_stack_df_ok <- div_um_edge_stack_df_ok[!is.na(div_um_edge_stack_df_ok$elev),]
 div_um_edge_stack_df_ok <- div_um_edge_stack_df_ok[!is.na(div_um_edge_stack_df_ok$ssdm),]
 
 # saveRDS(div_um_edge_stack_df_ok, "div_um_edge_stack_df_ok")
 
 plot(div_um_edge_stack_df_ok$alpha~div_um_edge_stack_df_ok$ssdm)
 
 test_df <- div_um_edge_stack_df_ok[sample(1:nrow(div_um_edge_stack_df_ok), 10000),]
 
 plot(test_df$alpha~test_df$ssdm)
 plot(test_df$alpha~test_df$elev)
 
 summary(lm(test_df$alpha~test_df$ssdm))
 summary(lm(test_df$alpha~test_df$elev))
 
 plot(test_df$alpha~test_df$dist_edge)
 
 summary(lm(test_df$alpha~test_df$dist_edge))
 
 boxplot(test_df$alpha~test_df$substrat)
 
 summary(lm(alpha ~ ssdm + elev + substrat + dist_edge, test_df))
 
 #### create a table with div data and substrat for each elevation band ####
 
div_um_edge_stack <- stack(alpha_div_ok,beta_div_ok, funct_div_ok, UM_raster_ok,edge_dist_raster_ok, ssdm_raster_ok)
raster(funct_div_ok, band=1)
div_um_stack_elev <- list()
for(i in 0:maxValue(dem_class_ok)){
  bnd_tmp <- dem_class_ok == i
  bnd_tmp[!bnd_tmp == 1] <- NA
  div_um_stack_bnd_tmp <- mask(div_um_stack, bnd_tmp)
  div_um_stack_elev[[i+1]] <- div_um_stack_bnd_tmp
  rm(div_um_stack_bnd_tmp)
}

# saveRDS(div_um_stack_elev, "/home/thesardfou/Documents/projets/Reliques/R/Reliques/synthese_spectral_div/div_um_stack_elev.rds")
div_um_stack_elev <- readRDS("/home/thesardfou/Documents/projets/Reliques/R/Reliques/synthese_spectral_div/div_um_stack_elev.rds")


div_um_elev_df <- c()
for(i in 1:length(div_um_stack_elev)){
  vals <- getValues(div_um_stack_elev[[i]])
  vals <- vals[!is.na(vals[,1]),]
  div_um_elev_df <- rbind(div_um_elev_df, cbind(vals, rep(i, nrow(vals))))
}

div_um_elev_df <- data.frame(div_um_elev_df)
colnames(div_um_elev_df) <- c("alpha", "beta", "funct_rich", "funct_eveness", "funct_div", "Substrat", "elev_band", "ssdm")
div_um_elev_df$Substrat[is.na(div_um_elev_df$Substrat)] <- "VS"
div_um_elev_df$Substrat[div_um_elev_df$Substrat == 1] <- "UM"
div_um_elev_df$Substrat <- as.factor(div_um_elev_df$Substrat)

div_um_elev_df %>% 
  filter(Substrat %in% unique(div_um_elev_df$Substrat)) %>%
  ggplot(aes(x=factor(elev_band),y=Shannon, fill=factor(Substrat))) +
  geom_point(position=position_jitterdodge(),alpha=0.03) +
  geom_boxplot() + 
  labs(fill = "Substrat") + 
  theme_bw(base_size = 16) +
  xlab("elevation band")

lm(div_um_elev_df)

############################################################################################################
############################################################################################################
############################################################################################################
#### BRT #### ===> from glm script, add selection of subset (distance from edge) and delete outliers 
############################################################################################################
############################################################################################################
############################################################################################################
library(gbm)
library(dismo)

##### Choose variable to analyse ####

names_var <- "PCoA.1" # for betadiv
names_var <- "PCoA.2" # for betadiv
names_var <- "PCoA.3" # for betadiv

names_var <- "Shannon" # for alphadiv
names_var <- "RICHNESS" # for funct div
names_var <-"EVENNESS" # for funct div
names_var <- "DIVERGENCE" # for funct div

##### function for simple model selection for landscape buffers of different sizes ####
msel_aic_uni = function(Y,X, empty = F){
  aic = c()
  for(i in 1:ncol(X)){
    lm = lm (Y~X[,i])
    aic = c(aic,AIC(lm))
    if(empty & i == ncol(X)) aic = c(aic, AIC(lm(Y~1)))
  }
  if(empty) tab = cbind(c(colnames(X), "intercept"), aic) else  tab = cbind(colnames(X), aic)
  return(list(pred = colnames(X)[order(aic)][1], tab = tab))
}  

##### get log-transformed distance to edge ####
data_div_frag$log_edge <- log(data_div_frag$dist_edge_centroid+1)

#### for loop : variable selection for model ####
rep_col <- grep(names_var, colnames(data_div_frag))
expl_col <- c(1:ncol(data_div_frag))[!grepl(names_var, colnames(data_div_frag))]
                                     
# minimum distance from edge for points (10: cells are not a the edge) #
#  no points < d m from edge #
d = 0
data_div_frag_ok <- data_div_frag[data_div_frag$dist_edge_centroid > d,]
nrow(data_div_frag_ok)

var_tmp <- colnames(data_div_frag_ok)[rep_col[i]]
data_tmp <- data_div_frag_ok[,c(rep_col[i], expl_col)]
# distance edge (log of standard)
aic_sel_edge <- msel_aic_uni(data_tmp[,1], cbind(dist_edge = data_tmp$dist_edge_centroid, log_edge = data_tmp$log_edge))

# habitat amount
aic_sel_HA <- msel_aic_uni(data_tmp[,1], data_tmp[,c( paste0(rep("prop.landscape_", 4), c(250,500,1000,2000), "_centroid"))])
aic_sel_HA

# effective mesh size
aic_sel_EMS <- msel_aic_uni(data_tmp[,1], data_tmp[,c( paste0(rep("effective.mesh.size_", 4), c(250,500,1000,2000), "_centroid"))])
aic_sel_EMS

colnames(data_tmp)
land_var <- paste0(rep(c("total.edge_", "effective.mesh.size_", "prop.landscape_"), each = 4), c(250,500,1000,2000), "_centroid")
# from selection
land_var <- c(aic_sel_HA$pred, aic_sel_EMS$pred)
land_var <- c(aic_sel_HA$pred)

topo_var <- c("elevation","slope","aspect" ,"curvature","substrat_centroid")
sel_var <- c("log_edge", land_var, topo_var) # need to log distance to edge !!! (see distribution)
# no NA values
data_for_mod <- data_tmp[, c(var_tmp, sel_var)]
# no NA values
data_for_mod <- data_for_mod[complete.cases(data_for_mod),]

#### final number of points ####
dim(data_for_mod)

#### change variable names ####
colnames(data_for_mod)
colnames(data_for_mod) <- c(colnames(data_for_mod)[1], "edge", "HA","elev","slop","asp","curv","subs")

data_for_mod_brt <- data_for_mod
data_for_mod_brt$subs <- as.factor(data_for_mod_brt$subs)

#### fit model ####
gbm_mod1 <- gbm.step(data= data_for_mod_brt, gbm.x = 2:ncol(data_for_mod_brt), gbm.y = 1,
                     family = "gaussian", tree.complexity = 5,
                     learning.rate = 0.01, bag.fraction = 0.5)

summary(gbm_mod1)

gbm.plot(gbm_mod1)
gbm.plot.fits(gbm_mod1)

gbm_mod1_simp <- gbm.simplify(gbm_mod1, n.drops = 5)

find.int <- gbm.interactions(gbm_mod1)
find.int$rank.list

gbm.perspec(gbm_mod1, 5, 4,z.range=c(4,7))
gbm.perspec(gbm_mod1, 7, 1,z.range=c(6,8))
gbm.perspec(gbm_mod1, 7, 1,z.range=c(6,8))

################################################
#### deviance explained ####
################################################
# D2 = 1 – (residual deviance/total deviance) (Nieto and Mélin, 2017, https://doi.org/10.1016/j.pocean.2016.11.009)
# D2 = (total deviance - cross validated residual deviance)/total deviance (Leathwick et al., 2006, doi:10.3354/meps321267)

#### results NC raw 100sp Beta PCoa1 : ####
# mean total deviance = 0.096 
# mean residual deviance = 0.077 
# estimated cv deviance = 0.082 ; se = 0 
expl_dev_Beta1 <- (0.096 - 0.077) /  0.096
expl_dev_Beta1

#### results NC raw 100sp Beta PCoa2 : ####
# mean total deviance = 0.062 
# mean residual deviance = 0.056 
# estimated cv deviance = 0.059 ; se = 0 
expl_dev_Beta2 <- (0.062 - 0.056) /  0.062
expl_dev_Beta2

#### results 58KFA raw 100sp Beta PCoa1 : ####
# mean total deviance = 0.066 
# mean residual deviance = 0.044 
# estimated cv deviance = 0.052 ; se = 0.001 

expl_dev_Beta1 <- (0.066 - 0.052) /  0.066

#### results 58KFA raw 100sp Shannon : ####
# mean total deviance = 0.152 
# mean residual deviance = 0.097 
# estimated cv deviance = 0.109 ; se = 0.004 

expl_dev_shannon <- (0.152 - 0.109) /  0.152

#### results 58KFA raw 100sp Richness : ####
# mean total deviance = 0.097 
# mean residual deviance = 0.065 
# estimated cv deviance = 0.074 ; se = 0.004 

expl_dev_richness <- (0.097 - 0.074) /  0.097


############################################################################################################
############################################################################################################
#### GAM ####
############################################################################################################
############################################################################################################
# see e.g. https://m-clark.github.io/generalized-additive-models/application.html
library(mgcv)
library(visreg)
data_for_mod_gam <- data_for_mod_brt
var_rep <- names(data_for_mod_gam)[1]
names(data_for_mod_gam)[1] <- "DIV"
names(data_for_mod_gam)
# need to transform 
#### one predictor ####
# linear model 
mod_lm <- gam(DIV ~ edge, data=data_for_mod_gam)
summary(mod_lm)

# GAM 
mod_gam1 <- gam(DIV ~ s(edge, bs="cr"), data=data_for_mod_gam)
summary(mod_gam1)
plot(mod_gam1)
visreg(mod_gam1)

# model comparison
AIC(mod_lm)
summary(mod_lm)$sp.criterion
summary(mod_lm)$r.sq 

AIC(mod_gam1)
summary(mod_gam1)$sp.criterion
summary(mod_gam1)$r.sq 

#### multiple predictors ####
# linear model 
mod_lm2 <- gam(DIV ~ edge + HA + subs + curv + elev + slop, data=data_for_mod_gam)
summary(mod_lm2)

# GAM 
mod_gam2 <- gam(DIV ~ s(edge) + s(HA) + s(curv) + s(elev) + s(slop) + subs, data=data_for_mod_gam)
summary(mod_gam2)
#  visreg(mod_gam2)

testdata = data.frame(edge = seq(min(data_for_mod_gam$edge), max(data_for_mod_gam$edge), length = 100),
                      HA = mean(mod_gam2$model$HA),
                      curv = mean(mod_gam2$model$curv),
                      elev = mean(mod_gam2$model$elev),
                      slop = mean(mod_gam2$model$slop),
                      subs = "UM")

fits = predict(mod_gam2, newdata=testdata, type='response', se=T)
predicts = data.frame(testdata, fits) %>% 
  mutate(lower = fit - 1.96*se.fit,
         upper = fit + 1.96*se.fit)
plot_mod_gam2_response = ggplot(aes(x=edge,y=fit), data=predicts) +
  geom_ribbon(aes(ymin = lower, ymax=upper), fill='gray90') +
  geom_line(color='#00aaff') +
  theme_bw()

plot_mod_gam2_response

#### test othe model ####
mod_gam3 <- gam(DIV ~ te(edge,HA) + subs, data=data_for_mod_gam)
summary(mod_gam3)
vis.gam(mod_gam3, type='response', plot.type='persp',
              phi=30, theta=40, n.grid=500, border=NA)

visreg2d(mod_gam3, xvar='edge', yvar='HA', scale='response')
#### nice plot draft ####

map(vars, function(x){
  p <- plotGAM(mod_gam1, smooth.cov = x, groupCovs = "Private") +
    geom_point(data = train.college, aes_string(y = "Outstate", x = x, color= "Private"), alpha = 0.2) +
    geom_rug(data = train.college, aes_string(y = "Outstate", x = x, color= "Private"  ), alpha = 0.2) +
    scale_color_manual("Private", values = c("#868686FF", "#0073C2FF")) +
    theme(legend.position="none")
  g <- ggplotGrob(p)
}) %>%
  {grid.arrange(grobs = (.), ncol = 3, nrow = 2)}

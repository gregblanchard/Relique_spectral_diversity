
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
####################################################################################
#### load files and data ####
####################################################################################
#### spatial data forest ####
# forest 
forest <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.shp")
zone <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/zone_plaine_des_lacs/zone.shp")
forest_map_raster <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT.tif")

# forest in the zone
CP <- as(extent(zone), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(forest))
forest_zone <-  gIntersection(forest, CP, byid=TRUE)
forest_zone$id_patch <- 1:length(forest_zone)
forest_zone_utm <- spTransform(forest_zone, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#####################################################################################################
####  Spectral div from tile 58KFA 2021 raw ####
#####################################################################################################

path_dir <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL"

####  25sp: 5*5 S2 cells ####  
# min sunlit 0.8
beta_div <- stack(paste0(path_dir,"/SPCA/BETA/BetaDiversity_BCdiss_PCO_5"))
alpha_div <- raster(paste0(path_dir,"/SPCA/ALPHA/Shannon_5_MeanFilter_Fullres"))
funct_div <- stack(paste0(path_dir,"/SPCA/FUNCTIONAL/FunctionalDiversity_Map_MeanFilter_Fullres"))

#### aggregate 10*10m resulution raster (S2 original) to get a user define resolution fo compute fragmentation indices ####
# alpha_div.aggregate <- aggregate(alpha_div, fact=3)
# res(alpha_div.aggregate)
# saveRDS(alpha_div.aggregate, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/rasterPDL_30x30.rds")
alpha_div.aggregate <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/rasterPDL_30x30.rds")

# 
# 
# # get sampled points 100sp
# path_dir_frag <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/results/test/"
# centroids_cells_points_wgs <- readRDS(paste0(path_dir_frag,"centroids_cells_100_points_wgs_landscape_metrics_topovar.rds"))
# centroids_cells_points_utm <- st_transform(centroids_cells_points_wgs, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#####################################################################################################
####  Get environemental variables and landscape metrics for all cell centroids ####
#####################################################################################################
raster_to_use <- beta_div$PCoA.1 

raster_to_use <- alpha_div.aggregate
# get cell centroids
xy_cells_centroids_utm <- SpatialPoints(xyFromCell(raster_to_use, 1:ncell(raster_to_use)), 
                          proj4string=CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
xy_cells_centroids_utm <- SpatialPointsDataFrame(xy_cells_centroids_utm, data = data.frame(values(raster_to_use)))
# keep only values defined in the raster (i.e. in forest)
xy_cells_centroids_utm <- xy_cells_centroids_utm[!is.na(xy_cells_centroids_utm$values.raster_to_use.),]

#####################################################################################################
####  Patch metrics from cell centroids ####
#####################################################################################################
# patch area 
forest_zone_utm$area <- st_area(st_as_sf(forest_zone_utm))
# ID patch
patchs_over_cells_utm <- sp::over(xy_cells_centroids_utm,forest_zone_utm) 
xy_cells_centroids_utm$id_patch <- patchs_over_cells_utm$id_patch
# patch area
xy_cells_centroids_utm$area <- patchs_over_cells_utm$area

# keep only centroids within forest 
centroids_cells_ok <- xy_cells_centroids_utm[!is.na(xy_cells_centroids_utm$id_patch),]
#### distance to edge #### 
n_pts <- nrow(centroids_cells_ok)
library(doSNOW)

cl <- makeCluster(3)
registerDoSNOW(cl)
clusterEvalQ(cl, library(sp))
clusterEvalQ(cl, library(rgeos))
pb <- txtProgressBar(max = n_pts, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
system.time(
  dist_edge_centroid <- foreach(i = 1:n_pts, .combine = c, 
                                .options.snow = opts) %dopar%
    {
      pts_tmp <- centroids_cells_ok[i,]
      if(!is.na(pts_tmp$id_patch)){
        patch_tmp <- forest_zone_utm[forest_zone_utm$id_patch == pts_tmp$id_patch,]
        forest_edge <- as(patch_tmp, "SpatialLines")
        dist_tmp <- gDistance(pts_tmp, forest_edge, byid=TRUE)
      }else{
        dist_tmp <- NA
      }
    }
)
close(pb)
stopCluster(cl) 

centroids_cells_ok$dist_edge <- dist_edge_centroid

dir_to_save <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/"
# saveRDS(centroids_cells_ok, paste0(dir_to_save,"centroids_cells_50m_PDL_utm.rds"))
# saveRDS(centroids_cells_ok, paste0(dir_to_save,"centroids_cells_30m_PDL_utm.rds"))

centroids_cells_ok <- readRDS(paste0(dir_to_save,"centroids_cells_50m_PDL_utm.rds"))
centroids_cells_ok <- readRDS(paste0(dir_to_save,"centroids_cells_30m_PDL_utm.rds"))

#####################################################################################################
##### landscape metrics for different landscape buffers ##### 
#####################################################################################################
# buffer values
buffers = c(50, 100, 250,500,1000)

# buffers = c(250,500)
n_pts <- nrow(centroids_cells_ok)
library(doSNOW)

cl <- makeCluster(3)
registerDoSNOW(cl)
pb <- txtProgressBar(max = n_pts, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
system.time(
  local_stats_centroid <- foreach(i = 1:n_pts, .packages= c("foreach","raster","sp","rgeos", "SDMTools"), .combine = rbind, 
                                  .options.snow = opts) %dopar%
    {
      foreach(j = length(buffers):1, .combine = cbind) %do% {
        points_spl_buffer_tmp <- gBuffer(centroids_cells_ok[i,], width=buffers[j], byid=TRUE)
        points_spl_buffer_tmp <- spTransform(points_spl_buffer_tmp, CRS("+proj=longlat +datum=WGS84 +no_defs"))
        if(j == length(buffers)){
          local_habitat =  crop(forest_map_raster, points_spl_buffer_tmp)
        }else{
          local_habitat =  crop(local_habitat, points_spl_buffer_tmp)
        }
        local_habitat =  mask(local_habitat, points_spl_buffer_tmp)
        local_stats_tmp = ClassStat(local_habitat)
        colnames(local_stats_tmp) = paste0( colnames(local_stats_tmp), "_", buffers[j], "_centroid")
        local_stats_tmp = local_stats_tmp[nrow(local_stats_tmp),]
      }
    }
)
close(pb)
stopCluster(cl) 

centroids_cells_ok_landscape <- cbind(centroids_cells_ok, local_stats_centroid)
#  saveRDS(centroids_cells_ok_landscape, paste0(dir_to_save,"centroids_cells_ok_landscape.rds"))
#  saveRDS(centroids_cells_ok_landscape, paste0(dir_to_save,"centroids_cells_30m_ok_landscape.rds"))

centroids_cells_ok_landscape <- readRDS(paste0(dir_to_save,"centroids_cells_30m_ok_landscape.rds"))

#####################################################################################################
#### Extract environmental variables for each cell ####
#####################################################################################################
#### raster with dimension = windowsize ####
S2_product <- alpha_div.aggregate

#### polygonize raster ####
poly <- polygonize(S2_product, na.rm = TRUE)
poly$id_cell_poly <- 1:nrow(poly)
# keep cells with defined controid

poly_ok <- st_intersection(st_as_sf(centroids_cells_ok_landscape),st_as_sf(poly))
poly_ok <- poly[poly$id_cell_poly %in% poly_ok$id_cell_poly,]

#### load environmental variables ####
# get values from gdal script (/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/gis_script/var_topo)

dem <- raster("/home/thesardfou/Documents/GIS/MNT/mnt10_GT.tif")
slope <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/slope_GT.tif")
aspect <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/aspect_GT.tif")
curv <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/curvature_GT.tif")
twi <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/twi_GT.tif")
prec <- raster("/home/thesardfou/Documents/GIS/precipitations/precipitations_utm58s.tif")

cells_alt <-  exact_extract(dem, poly_ok, "mean")
cells_slope <-  exact_extract(slope, poly_ok, "mean")
cells_aspect <-  exact_extract(aspect, poly_ok, "mean")
cells_curv <-  exact_extract(curv, poly_ok, "mean")
cells_twi <-  exact_extract(twi, poly_ok, "mean")
cells_prec <-  exact_extract(prec, poly_ok, "mean")

centroids_cells_ok_landscape$elevation <- cells_alt
centroids_cells_ok_landscape$slope <- cells_slope
centroids_cells_ok_landscape$aspect <- cells_aspect
centroids_cells_ok_landscape$curvature <- cells_curv
centroids_cells_ok_landscape$twi <- cells_twi
centroids_cells_ok_landscape$prec <- cells_prec

# saveRDS(centroids_cells_ok_landscape, paste0(dir_to_save,"centroids_cells_ok_landscape_env.rds"))
# saveRDS(centroids_cells_ok_landscape, paste0(dir_to_save,"centroids_cells_30m_ok_landscape_env.rds"))

centroids_cells_30m_ok_landscape_env <- readRDS(paste0(dir_to_save,"centroids_cells_30m_ok_landscape_env.rds"))

#####################################################################################################
#### compute the mean raw bands and mean PCs values for each 30*30 cell ####
#####################################################################################################
#### load raw bands and PCA values from S2 images ####
#### get raw spectral bands from S2 image  ####
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021"
S2_raw <- stack(paste0(file_path,"/T58KFA_raw_int_forest_zone_PDL.tif"))
# rename bands 
names(S2_raw) <- c("B2",     "B3",     "B4",     "B5",     "B6",     "B7",     "B8",    "B8A",    "B11",    "B12")
# sload PCA from selected S2 image
file_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL"
PCA_map <- stack(paste0(file_path,"/SPCA/PCA/OutputPCA_8_PCs"))

#### get polygons to extract
alpha_div.aggregate <- readRDS("/home/greg/projects/reliques/signature_spectrale_fragmentation/maps/biodivmapr/58FA_2021/rasterPDL_30x30.rds")

#### extract values for polygons ####
cells_mean_S2_spectral_raw <- exact_extract(S2_raw, poly_ok, "mean")
cells_mean_PCA_spectral_raw <- exact_extract(PCA_map, poly_ok, "mean")

####  add this to spatial points centroid of cells #### 
centroids_cells_30m_ok_landscape_env_S2 <- cbind(centroids_cells_30m_ok_landscape_env, cells_mean_S2_spectral_raw, cells_mean_PCA_spectral_raw)
# saveRDS(centroids_cells_30m_ok_landscape_env_S2, paste0(dir_to_save,"centroids_cells_30m_ok_landscape_env_S2.rds"))

####  add this to spatial polygons of cells #### 
poly_cells_30m_ok_landscape_env_S2 <- cbind(poly_ok, centroids_cells_30m_ok_landscape_env@data, cells_mean_S2_spectral_raw, cells_mean_PCA_spectral_raw)
# saveRDS(poly_cells_30m_ok_landscape_env_S2, paste0(dir_to_save,"poly_cells_30m_ok_landscape_env_S2.rds"))

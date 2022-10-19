

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
####  make / choose dir  ####
#####################################################################################################


path_dir <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/results/"
name_dir <- "test"
dir.create(paste0(path_dir,name_dir))


dir_to_save <- paste0(path_dir,name_dir,"/")

#####################################################################################################
####  load forest map  ####
#####################################################################################################

# forest polygones
# forest_map <- readOGR(
#   "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.shp")
# forest_map@data$id_tmp <- 1:dim(forest_map@data)[1]
# forest_map_utm <- spTransform(forest_map, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# all_patch_areas <- gArea(forest_map_utm,  byid=TRUE)
# forest_map@data$area <- forest_map_utm@data$area <- all_patch_areas
# saveRDS(forest_map, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.rds")

forest_map <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_NC.rds")
forest_map_utm <- spTransform(forest_map, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# plot(forest_map)

forest_map_utm$id_tmp <- 1:nrow(forest_map_utm)

forest_map_raster <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT.tif")
forest_map_raster_withna <- raster(
        "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT_withna.tif")
#####################################################################################################
####  Sample random points in forest  ####
#####################################################################################################

#### 20 000 pts from random raster cells ####
n <- 20000
## alternative (longer)
# points_spl <- spsample(forest_map_utm, n, type = "stratified")
 ## from raster cells centroids
 r <- forest_map_raster_withna
 # r is a raster, n is the number of points requested
 v <- raster::getValues(r)
 
 # v.ok <- which( v == 1 & !is.na(v)) # not na and only in forest
 # need to divid the number of cells by 2 becaus resulting vector is too large
 sub_v1 <- v[1:(ncell(r) / 2)]
 sub_v1.ok <- which( sub_v1 == 1 & !is.na(sub_v1)) # not na and only in forest

 sub_v2 <- v[((ncell(r) / 2)+1):ncell(r)]
 sub_v2.ok <- which( sub_v2 == 1 & !is.na(sub_v2)) + (ncell(r) / 2) # not na and only in forest
 
 v.ok <- c(sub_v1.ok,sub_v2.ok)
 # sample n cells 
 x <- sample(v.ok, n)
 samp_pts <- raster::xyFromCell(r, x, spatial=TRUE)
 plot(samp_pts)
 
points_spl <- SpatialPointsDataFrame(samp_pts, data.frame(ID=1:length(samp_pts)))
 points_spl_utm <- spTransform(points_spl, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
 
 # saveRDS(points_spl_utm, paste0(dir_to_save,"points_spl_utm_from_raster.rds"))

 points_spl_utm <- readRDS(paste0(dir_to_save,"points_spl_utm_from_raster.rds"))

 
 #####################################################################################################
 ####  Patch metrics from random points in forest  ####
 #####################################################################################################
 

 # ID patch
 patchs_over_points_utm <- over(points_spl_utm,forest_map_utm) 
 points_spl_utm$id_patch <- patchs_over_points_utm$id_tmp
 # patch area
 points_spl_utm$area <- patchs_over_points_utm$area
 
 # keep only points in forest polygons
 points_spl_utm <- points_spl_utm[!is.na(points_spl_utm$id_patch),]
 
 #### distance to edge #### 
 n_pts <- nrow(points_spl_utm)
 library(doSNOW)
 
 cl <- makeCluster(3)
 registerDoSNOW(cl)
 clusterEvalQ(cl, library(sp))
 clusterEvalQ(cl, library(rgeos))
 pb <- txtProgressBar(max = n_pts, style = 3)
 progress <- function(n) setTxtProgressBar(pb, n)
 opts <- list(progress = progress)
 system.time(
         dist_edge <- foreach(i = 1:n_pts, .combine = c, 
                              .options.snow = opts) %dopar%
                 {
                         pts_tmp <- points_spl_utm[i,]
                         patch_tmp <- forest_map_utm[forest_map_utm$id_tmp == pts_tmp$id_patch,]
                         forest_edge <- as(patch_tmp, "SpatialLines") 
                         dist_tmp <- gDistance(pts_tmp, forest_edge, byid=TRUE)
                 }
 )
 close(pb)
 stopCluster(cl) 
 
 points_spl_utm$dist_edge <- dist_edge
 
 # saveRDS(points_spl_utm, paste0(dir_to_save,"points_spl_utm_patch_metrics.rds"))
 
 points_spl_utm <- readRDS(paste0(dir_to_save,"points_spl_utm_patch_metrics.rds"))
 
 #####################################################################################################
 #####  landscape metrics for diferent landscape buffers ##### 
 #####################################################################################################
 # buffer values 
 buffers = c(250,500,1000,2000)
 # buffers = c(250,500)
 n_pts <- nrow(points_spl_utm)
 library(doSNOW)
 
 
cl <- makeCluster(3)
registerDoSNOW(cl)
pb <- txtProgressBar(max = n_pts, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
system.time(
    local_stats <- foreach(i = 1:n_pts, .packages= c("foreach","raster","sp","rgeos", "SDMTools"), .combine = rbind, 
                                    .options.snow = opts) %dopar%
                         {
                                 foreach(j = length(buffers):1, .combine = cbind) %do% {
                                         points_spl_buffer_tmp <- gBuffer(points_spl_utm[i,], width=buffers[j], byid=TRUE)
                                         points_spl_buffer_tmp <- spTransform(points_spl_buffer_tmp, CRS("+proj=longlat +datum=WGS84 +no_defs"))
                                        if(j == length(buffers)){
                                                local_habitat =  crop(forest_map_raster, points_spl_buffer_tmp)
                                        }else{
                                                local_habitat =  crop(local_habitat, points_spl_buffer_tmp)
                                                }
                                         local_habitat =  mask(local_habitat, points_spl_buffer_tmp)
                                         local_stats_tmp = ClassStat(local_habitat)
                                         colnames(local_stats_tmp) = paste0( colnames(local_stats_tmp), "_", buffers[j])
                                         local_stats_tmp = local_stats_tmp[nrow(local_stats_tmp),]
                                 }
                       
                         }
                 
                 
)
close(pb)
stopCluster(cl) 

#  saveRDS(local_stats, paste0(dir_to_save,"landscape_metrics.rds"))
 
local_stats <- readRDS(paste0(dir_to_save,"landscape_metrics.rds"))

 points_spl_utm <- cbind(points_spl_utm, local_stats)
 
 #####################################################################################################
 #### extract environmental variables ####
 #####################################################################################################
 # get values from gdal script (/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/gis_script/var_topo)
 
 dem <- raster("/home/thesardfou/Documents/GIS/MNT/mnt10_GT.tif")
 slope <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/slope_GT.tif")
 aspect <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/aspect_GT.tif")
 curv <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/curvature_GT.tif")
 twi <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/twi_GT.tif")
 
 pts_alt <- extract(dem, points_spl_utm)
 pts_slope <- extract(slope, points_spl_utm)
 pts_aspect <- extract(aspect, points_spl_utm)
 pts_cur <- extract(curv, points_spl_utm)
 pts_twi <- extract(twi, points_spl_utm)
 
 points_spl_utm$elevation <- pts_alt
 points_spl_utm$slope <- pts_slope
 points_spl_utm$aspect <- pts_aspect
 points_spl_utm$curvature <- pts_cur
 
 # distance to river network ==> récuérer le reseau hydrographique sur le nas! 
 # channel <- readOGR(".shp")
 # all_dist_to_channel <- gDistance(points_spl, channel, byid=TRUE)
 # dist_to_channel <- apply(points_spl, 2,min)

 ##### UM and non-UM areas ####
 
 UM_GT <- readOGR("/home/thesardfou/Documents/GIS/divers/NC_GT_UM_areas.shp")
 UM_GT_UTM <-  spTransform(UM_GT, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs "))
 
 pts_in_UM <- over( points_spl_utm, UM_GT_UTM)
 substrat <- ifelse(is.na(pts_in_UM$gid), "non_UM", "UM")
 points_spl_utm$substrat <- substrat
 
 # saveRDS(points_spl_utm, paste0(dir_to_save,"points_spl_utm_landscape_metrics_topovar.rds"))
 
 #####################################################################################################
 #####################################################################################################
 #####################################################################################################
 ####  get a grid with resolution * 10 (10m -> 100m) to get spectral communities (100 species)  ####
 ####  get a grid with resolution * 5 (10m -> 50m) to get spectral communities (25 species)  ####
 #####################################################################################################
 #####################################################################################################
 #####################################################################################################
 # load image from GEE to get raster dimensions to consider
 im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons"
 im_name <- "pca_S2_NC_SS_int_new"
 im_ex <- stack( paste0(im_path, "/", im_name, ".tif") )
 
 # get sampled points
 points_spl_utm <- readRDS(paste0(dir_to_save,"points_spl_utm_landscape_metrics_topovar.rds"))
 points_spl_wgs <- spTransform(points_spl_utm, CRS("+proj=longlat +datum=WGS84 +no_defs"))
 
 # raster with dimension = windowsize
 windowsize <- 10
 windowsize <- 5
 res_S2 <- res(im_ex)
 r <- raster(ext = extent(im_ex), res=c(res_S2*windowsize,res_S2*windowsize))
 # set all celles to NA
 r[] <- NA
 # set sampled cells to 1 (in order to polygonize only sampled cells)
 cells_over_spl_pts <- cellFromXY(r,points_spl_wgs)
 r[cells_over_spl_pts] <- 1

 # #### polygonize all raster ####
 # polygon_grid_50m <- polygonize(r, na.rm = TRUE)
 # ID_cell <- rownames(polygon_grid_50m)
 # polygon_grid_50m <- as(polygon_grid_50m, "Spatial")
 # polygon_grid_50m$ID_cell <- as.numeric(ID_cell)
 
 # #polygon_grid_50m <- as(polygon_grid_50m, "Spatial")
 # polygon_grid_50m$ID_cell <- as.numeric(ID_cell)
 # 
 # # merge with points
 # points_spl_wgs$ID_cell <- cells_over_spl_pts
 # polygon_grid_50m <- st_join(st_as_sf(polygon_grid_50m), st_as_sf(points_spl_wgs, by = "ID_cell"))
 # polygon_grid_50m <- as(polygon_grid_50m, "Spatial")
 # polygon_grid_50m$ID_cell <- polygon_grid_50m$ID_cell.x
 # # saveRDS(polygon_grid_50m, paste0(dir_to_save,"cells_50_wgs.rds"))
 
 # #### split raster in several parts to perform poligonize ####
  library(SpaDES)
 sub_r <- splitRaster(r,3,3)
 
 n_tiles <- length(sub_r)
 library(doSNOW)
 cl <- makeCluster(3)
 registerDoSNOW(cl)
 clusterEvalQ(cl, library(spex))
 clusterEvalQ(cl, library(raster))
 pb <- txtProgressBar(max = n_tiles, style = 3)
 progress <- function(n) setTxtProgressBar(pb, n)
 opts <- list(progress = progress)
 system.time(
    polygon_grid_list <- foreach(i = 1:n_tiles, .combine = 'c', 
                                  .options.snow = opts) %dopar%
       {
          r_tmp <- sub_r[[i]]
          poly <- polygonize(r_tmp, na.rm = TRUE)
       }
 )
 close(pb)
 stopCluster(cl) 
 
 polygon_grid_list_ok <- list()
 for(i in 1:(length(polygon_grid_list)/2)) polygon_grid_list_ok[[i]] <- polygon_grid_list[[ seq(2,18,2)[i] ]]
 polygon_grid_list_ok[sapply(polygon_grid_list_ok, function(x) length(x) == 0)] <- NULL

 for(i in 1:length(polygon_grid_list_ok)) polygon_grid_list_ok[[i]] <- st_as_sf(polygon_grid_list_ok[[i]])
 polygon_grid <- mapedit:::combine_list_of_sf(polygon_grid_list_ok)
 
 polygon_grid_50m <- polygon_grid
 
  # saveRDS(polygon_grid_50m, paste0(dir_to_save,"cells_50_wgs.rds"))
 

 cells_over_points_wgs <- readRDS( paste0(dir_to_save,"cells_100_over_points_wgs.rds"))
 cells_over_points_wgs <- readRDS( paste0(dir_to_save,"cells_50_wgs.rds"))
 # # this is no more needed for new method to generate larges cells
 # cells_over_points_wgs_unique <- cells_over_points_wgs[!duplicated(cells_over_points_wgs$ID_cell),]
 # saveRDS(cells_over_points_wgs_unique, paste0(dir_to_save,"cells_over_points_wgs_unique.rds"))

 centroids_cells_points_wgs <- st_centroid(st_as_sf(cells_over_points_wgs))
 centroids_cells_points_wgs$ID_cell <- exact_extract(r, centroids_cells_points_wgs)
 
 # saveRDS(centroids_cells_points_wgs, paste0(dir_to_save,"centroids_cells_100_points_wgs.rds"))
 centroids_cells_points_wgs <- readRDS( paste0(dir_to_save,"centroids_cells_100_points_wgs.rds"))

  #####################################################################################################
 #### get 100*100 grid centroids ####
 #####################################################################################################
 centroids_cells_points_utm <- st_transform(centroids_cells_points_wgs, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs "))
 centroids_cells_points_wgs_sp <- as(centroids_cells_points_wgs, "Spatial")
 centroids_cells_points_utm_sp <- as(centroids_cells_points_utm, "Spatial")
  #####################################################################################################
 ####  FOR 100*100 grid: Patch metrics from random points in forest  ####
 #####################################################################################################

 # ID patch
 patchs_over_points_utm <- sp::over(centroids_cells_points_utm_sp,forest_map_utm) 
 centroids_cells_points_utm_sp$id_patch <- patchs_over_points_utm$id_tmp
 # patch area
 centroids_cells_points_utm_sp$area <- patchs_over_points_utm$area
 
 # keep only centroids within forest 
 centroids_cells_points_utm_sp_ok <- centroids_cells_points_utm_sp[!is.na(centroids_cells_points_utm_sp$id_patch),]
 centroids_cells_points_wgs <- centroids_cells_points_wgs[!is.na(centroids_cells_points_utm_sp$id_patch),]
 centroids_cells_points_sp_utm <- centroids_cells_points_sp_utm_ok
 #### distance to edge #### 
 n_pts <- nrow(centroids_cells_points_wgs)
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
       pts_tmp <- centroids_cells_points_sp_utm[i,]
       patch_tmp <- forest_map_utm[forest_map_utm$id_tmp == pts_tmp$id_patch,]
       forest_edge <- as(patch_tmp, "SpatialLines")
       dist_tmp <- gDistance(pts_tmp, forest_edge, byid=TRUE)
     }
 )
 close(pb)
 stopCluster(cl) 
 
 centroids_cells_points_wgs$dist_edge_centroid <- dist_edge_centroid
 
 # saveRDS(centroids_cells_points_wgs, paste0(dir_to_save,"centroids_cells_100_points_wgs.rds"))
 # saveRDS(centroids_cells_points_wgs, paste0(dir_to_save,"centroids_cells_50_points_wgs.rds"))
 
 centroids_cells_points_wgs <- readRDS(paste0(dir_to_save,"centroids_cells_100_points_wgs.rds"))
 centroids_cells_points_wgs <- readRDS(paste0(dir_to_save,"centroids_cells_50_points_wgs.rds"))
 
 #####################################################################################################
 ##### FOR 100*100 grid: landscape metrics for diferent landscape buffers ##### 
 #####################################################################################################
 # buffer values 
 buffers = c(250,500,1000,2000)
 # buffers = c(250,500)
 n_pts <- nrow(centroids_cells_points_wgs)
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
         points_spl_buffer_tmp <- gBuffer(centroids_cells_points_sp_utm[i,], width=buffers[j], byid=TRUE)
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
 
 #  saveRDS(local_stats_centroid, paste0(dir_to_save,"landscape_metrics_centroid.rds"))
 
 local_stats_centroid <- readRDS(paste0(dir_to_save,"landscape_metrics_centroid.rds"))
 
 centroids_cells_points_wgs <- cbind(centroids_cells_points_wgs, local_stats_centroid)
 
 #####################################################################################################
 #### FOR 100*100 grid: extract environmental variables ####
 #####################################################################################################
 # get values from gdal script (/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/gis_script/var_topo)
 
 dem <- raster("/home/thesardfou/Documents/GIS/MNT/mnt10_GT.tif")
 slope <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/slope_GT.tif")
 aspect <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/aspect_GT.tif")
 curv <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/curvature_GT.tif")
 twi <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/twi_GT.tif")
 prec <- raster("/home/thesardfou/Documents/GIS/precipitations/precipitations_utm58s.tif")
 
 ## !!!!! here get the mean of the 100*100 cell! 
 cells_over_points_wgs_sf <- st_as_sf(readRDS( paste0(dir_to_save,"cells_50_wgs.rds")))
 
 cells100_alt <-  exact_extract(dem, cells_over_points_wgs_sf, "mean")
 cells100_slope <-  exact_extract(slope, cells_over_points_wgs_sf, "mean")
 cells100_aspect <-  exact_extract(aspect, cells_over_points_wgs_sf, "mean")
 cells100_curv <-  exact_extract(curv, cells_over_points_wgs_sf, "mean")
 cells100_twi <-  exact_extract(twi, cells_over_points_wgs_sf, "mean")
 cells100_prec <-  exact_extract(prec, cells_over_points_wgs_sf, "mean")
 
 centroids_cells_points_wgs$elevation <- cells100_alt
 centroids_cells_points_wgs$slope <- cells100_slope
 centroids_cells_points_wgs$aspect <- cells100_aspect
 centroids_cells_points_wgs$curvature <- cells100_curv
 centroids_cells_points_wgs$twi <- cells100_twi
 centroids_cells_points_wgs$prec <- cells100_prec
 
 # distance to river network ==> récuérer le reseau hydrographique sur le nas! 
 # channel <- readOGR(".shp")
 # all_dist_to_channel <- gDistance(points_spl, channel, byid=TRUE)
 # dist_to_channel <- apply(points_spl, 2,min)
 
 ##### UM and non-UM areas ####
 
 UM_GT <- readOGR("/home/thesardfou/Documents/GIS/divers/NC_GT_UM_areas.shp")
 centroids_cells_points_wgs_sp <- as(centroids_cells_points_wgs, "Spatial")
 
 cells100_in_UM <- over( centroids_cells_points_wgs_sp, UM_GT)
 
 cells100_substrat <- ifelse(is.na(cells100_in_UM$gid), "non_UM", "UM")
 centroids_cells_points_wgs$substrat_centroid <- cells100_substrat
 
 # saveRDS(centroids_cells_points_wgs, paste0(dir_to_save,"centroids_cells_100_points_wgs_landscape_metrics_topovar.rds"))
 # saveRDS(centroids_cells_points_wgs, paste0(dir_to_save,"centroids_cells_50_points_wgs_landscape_metrics_topovar.rds"))
 
  #####################################################################################################
 #####################################################################################################
 ####  GET SAMPLED POINTS WITH ALL LANDSACAPE & TOPO METRICS FOR FRAGMENTATION ANALYSES  ####
 #####################################################################################################
 #####################################################################################################
 
 points_spl_utm_ok <- readRDS(paste0(dir_to_save,"points_spl_utm_landscape_metrics_topovar.rds"))
 
 points_spl_utm_ok$log_edge <- log(points_spl_utm_ok$dist_edge+1)
 #####################################################################################################
 #####################################################################################################
 ####  Spectral images from Sentinel 2  ####
 #####################################################################################################
 #####################################################################################################
 
 #####################################################################################################
 ####  load PCA from GEE and apply forest mask  ####
 #####################################################################################################
 
 #### select PCA image of dry season ####
 
 im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons"
 im_name <- "pca_S2_NC_SS_int_new"
 
 pca_S2 <- stack( paste0(im_path, "/", im_name, ".tif") )
 
 #### mask raster with forest #### 
 # get forest raster
 forest_map_raster_withna <- raster(
         "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT_withna.tif")
 
 pca_S2_forest <- mask(pca_S2, forest_map_raster_withna)
 
 # writeRaster(pca_S2_forest,  paste0(dir_to_save, "pca_S2_forest.tif"),
 #             datatype="INT1U", options="COMPRESS=LZW")
 
 
 #####################################################################################################
 ####  load PCA from GEE with forest mask (if already done)  ####
 #####################################################################################################
 im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons"
 im_name <- "pca_S2_NC_SH_int_forest"
 
 pca_S2_NC_SH_forest <- stack( paste0(im_path, "/", im_name, ".tif") )
 
 im_name <- "pca_S2_NC_SS_int_forest"
 pca_S2_NC_SS_forest <- stack( paste0(im_path, "/", im_name, ".tif") )
 
 
 #####################################################################################################
 ####  Variance explained by PCA (get this from GEE, change file if needed)  ####
 #####################################################################################################
 var_explained_SH <- read.delim("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons/variance_explained_SH.txt",
                                header = F)
 var_explained_SH <- as.numeric(substring(var_explained_SH$V1, 3,))

 #####################################################################################################
 ####  export PCA values for sampled points  ####
 #####################################################################################################
 
 pca_SH_pts <-  terra::extract(pca_S2_NC_SH_forest, st_as_sf(points_spl_utm_ok))
 pca_SS_pts <- terra::extract(pca_S2_NC_SS_forest, st_as_sf(points_spl_utm_ok))
 
  # saveRDS(pca_SH_pts, paste0(dir_to_save,"pca_SH_pts.rds"))
  # saveRDS(pca_SS_pts, paste0(dir_to_save,"pca_SS_pts.rds"))

 pca_SH_pts <- readRDS(paste0(dir_to_save,"pca_SH_pts.rds"))
 pca_SS_pts <- readRDS(paste0(dir_to_save,"pca_SH_pts.rds"))

 points_spl_pca_SH <- cbind(data.frame(pca_SH_pts), points_spl_utm_ok)
 points_spl_pca_SS <- cbind(data.frame(pca_SS_pts), points_spl_utm_ok)
 
 #####################################################################################################
 #### any differences between dry and humid seasons for each PCA axis ####
 #####################################################################################################
 
 for(i in 1:10){
   points_spl_pca_SH_tmp <- points_spl_pca_SH[,i]
   points_spl_pca_SS_tmp <- points_spl_pca_SS[,i]
   points_spl_pca_SH_tmp <- points_spl_pca_SH_tmp[!is.na(points_spl_pca_SH_tmp)]
   points_spl_pca_SS_tmp <- points_spl_pca_SS_tmp[!is.na(points_spl_pca_SS_tmp)]
   plot(points_spl_pca_SH[,i] ~ points_spl_pca_SS[,i], main = paste("PC",i))
 }
 
 #####################################################################################################
 #### model selection with glm for each PCA axis ####
 #####################################################################################################
 library(MuMIn)
 library(parallel)
 
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
 
 ##### Choose dry or humid season! ####
 # SH
 data_pca <- as.data.frame(points_spl_pca_SH)
 names_pca <- "pca_S2_NC_SH_int_forest"
 
 # SS
 data_pca <- as.data.frame(points_spl_pca_SS)
 names_pca <- "pca_S2_NC_SS_int_forest"
 
 # arrange data 
 dim(data_pca)
 data_pca <- data_pca[complete.cases(data_pca[,1:10]),]
 dim(data_pca)
 colnames(data_pca)
 
 #### for loop : model selection for each pca axis ####
 pca_axis_col <- grep(names_pca, colnames(data_pca))
 expl_col <- c(1:ncol(data_pca))[!grepl(names_pca, colnames(data_pca))]
 
 #### minimum distance from edge for points (10: cells are not a the edge) ####
 d = 0
 
 list_model_tab <- list()
 list_model_best <- list()
 list_model_avg <- list()
 for(i in 1:length(pca_axis_col)){
   print(paste("model selection for PCA axis", i, "on", length(pca_axis_col)))
   pca_axis_tmp <- colnames(data_pca)[pca_axis_col[i]]
   data_tmp <- data_pca[,c(pca_axis_col[i], expl_col)]
   ####  no points < d m from edge #### 
   data_tmp <- data_tmp[data_tmp$dist_edge > d,]
   # distance edge (log of standard)
   aic_sel_edge <- msel_aic_uni(data_tmp[,1], cbind(dist_edge = data_tmp$dist_edge, log_edge = data_tmp$log_edge))
   
   # total habitat
   aic_sel_HA <- msel_aic_uni(data_tmp[,1], data_tmp[,c( paste0(rep("prop.landscape_", 4), c(250,500,1000,2000)))])
   aic_sel_HA
   # habitat amount
   aic_sel_totalarea <- msel_aic_uni(data_tmp[,1], data_tmp[,c( paste0(rep("total.area_", 4), c(250,500,1000,2000)))])
   aic_sel_totalarea
   
   # effective mesh size
   aic_sel_EMS <- msel_aic_uni(data_tmp[,1], data_tmp[,c( paste0(rep("effective.mesh.size_", 4), c(250,500,1000,2000)))])
   aic_sel_EMS
   
   colnames(data_tmp)
   land_var <- paste0(rep(c("total.edge_", "total.area_", "effective.mesh.size_", "prop.landscape_"), each = 4), c(250,500,1000,2000))
   # from selection
   land_var <- c(aic_sel_HA$pred, aic_sel_EMS$pred)
   
   topo_var <- c("elevation","slope","aspect" ,"curvature","substrat")
   sel_var <- c(aic_sel_edge$pred, "area",land_var, topo_var)
   # no NA values
   data_for_mod <- data_tmp[, c(pca_axis_tmp, sel_var)]
   # no NA values
   data_for_mod <- data_for_mod[complete.cases(data_for_mod),]
   
   #### final number of points ####
   dim(data_for_mod)
   #### standardize variables ####
   scale_01 <- function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
   
   data_for_mod_scaled <- as.data.frame(apply(data_for_mod, 2, function(x) if(!any(is.numeric(x))){
     scale(as.numeric(as.factor(x)))
   }else{
     scale(is.numeric(x))
   }))
   
   #### model selection for each PCA axis ####
   # with gaussiian family distribution 
   # do not keep Variance Inflation Factors > 5
   
   glm1 <- glm( formula(paste0(pca_axis_tmp,  " ~ .")), data = data_for_mod_scaled, na.action = "na.fail")
   vif <- car::vif(glm1)
   data_for_mod_scaled_ok <- data_for_mod_scaled
   while(max(vif)>5){
     data_for_mod_scaled_ok <- data_for_mod_scaled_ok[, !colnames(data_for_mod_scaled_ok) %in% names(vif)[vif == max(vif)]]
     glm1 <- glm(formula(paste0(pca_axis_tmp,  " ~ .")), data = data_for_mod_scaled_ok, na.action = "na.fail")
     vif <- car::vif(glm1)
   }
   vif
   glm1 <- glm(formula(paste0(pca_axis_tmp,  " ~ .")), data = data_for_mod_scaled_ok, na.action = "na.fail")
   
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
 
 names(list_model_tab) =  names(list_model_best) =  names(list_model_avg) <- colnames(data_pca)[pca_axis_col]
 
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
     ggtitle(names(list_model_tab[i])) +
     theme(legend.position="none",
           plot.title = element_text(size = 10, face = "bold"),
           axis.text.x = element_text(size=15),
           axis.text.y = element_text(size=15),
           axis.title.x = element_text(size=15),
           axis.title.y = element_text(size=15)) +
     xlab('standardized coefficients') 
 }
 
 plot_coef_list
 
 
 ##############################################################################
 #### work on residuals from regression on topographic variables  ####
 ##############################################################################
 
 data_pca_nona <- data_pca[!is.na(data_pca$elevation),]
 data_pca_nona <- data_pca_nona[!is.na(data_pca_nona$substrat),]
 data_pca_nona <- data_pca_nona[!is.na(data_pca_nona$slope),]
 data_pca_nona <- data_pca_nona[!is.na(data_pca_nona$aspect),]
 data_pca_nona <- data_pca_nona[!is.na(data_pca_nona$curvature),]
 
 lm <- lm(as.formula(paste0(names_pca,".1  ~ elevation + substrat + slope + aspect")), data = data_pca_nona)
 lm <- lm(as.formula(paste0(names_pca,".2  ~ elevation + substrat + slope + aspect + curvature")), data = data_pca_nona)
 
 stand_res <- scale(lm$residuals)
 
 plot(stand_res~data_pca_nona$log_edge)
 
 lm_res <- lm(stand_res ~ data_pca_nona$log_edge)
 summary(lm_res)
 abline(lm_res, col = 'blue')
 smry1 = summary(lm_res)
 
 plot(stand_res~data_pca_nona$prop.landscape_250)
 lm_res <- lm(stand_res ~ data_pca_nona$prop.landscape_250)
 summary(lm_res)
 abline(lm_res, col = 'blue')
 
 lm_res <- lm(stand_res ~ data_pca_nona$prop.landscape_250 + data_pca_nona$log_edge)
 summary(lm_res)
 
 #### nice plot ####
 
 data_plot1 <- data.frame(stand_res, log_edge= data_pca_nona$log_edge)
 
 plot <- ggplot(data = data_plot1, aes(x = log_edge, y = stand_res)) + 
   geom_point(color='navy', size = 1, alpha =.3) +
   geom_smooth(data = data_plot1, aes(x = log_edge, y = stand_res), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
   
   annotate("text", x=6.5, y= -5,
            label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
            color='red') +
   ylab("Standardized residuals\nfrom PCA2 ~ topography") + xlab("Log distance edge") +
   theme(
     axis.text.x = element_text(size=12),
     axis.text.y = element_text(size=12)
   ) +
   theme_bw()
 
 plot
 
 ##############################################################################
 #### Conclusions  ####
 ##############################################################################
 
 # signature of fragmentation is higher on the dry season 
 
 # PCA axes to select : 
 # 2,3,5
 
 #####################################################################################################
 #####################################################################################################
 ####  Extract PCA values for spectral communities (from PCA axes (100*100 cells) ####
 #####################################################################################################
 #####################################################################################################
 
 #### get polygon grid 100*100m ####
 cells_over_points_wgs <- readRDS( paste0(dir_to_save,"cells_over_points_wgs_unique.rds"))
 #### get polygon grid 50*50m ####
 cells_over_points_wgs <- readRDS( paste0(dir_to_save,"cells_50_wgs.rds"))
 cells_over_points_wgs$ID_cell <- 1:length(cells_over_points_wgs)

sqrt(area(cells_over_points_wgs[1,]))

#### get cells in forest  #### 
cells_over_points_wgs_utm_tmp <- st_centroid(st_transform(st_as_sf(cells_over_points_wgs), "+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs"))
cells_over_points_wgs_utm_tmp <- as(cells_over_points_wgs_utm_tmp, "Spatial")
cells_over_points_wgs_utm_tmp$ID_cell <- 1:nrow(cells_over_points_wgs_utm_tmp)

cells_ok <- sp::over(cells_over_points_wgs_utm_tmp,forest_map_utm)
cells_ok$id_tmp
cells_over_points_wgs_inforest <- cells_over_points_wgs[!is.na(cells_ok$id_tmp),]
# saveRDS(cells_over_points_wgs_inforest, paste0(dir_to_save,"cells50_inforest.rds"))
cells_over_points_wgs <- readRDS( paste0(dir_to_save,"cells50_inforest.rds"))

cells_over_points_wgs_sf <- st_as_sf(cells_over_points_wgs) 

#### get PCA value from cells "community" of PCA raster for polygons over sampling points on the 100*100m grid ####
 # Select PCA axes to keep (selection from previous analyses PCA~ fragmentation)
 sel_axes <- c(2,3,5)
 pca_S2_NC_SH_forest_select <- subset(pca_S2_NC_SH_forest, sel_axes)
 pca_S2_NC_SS_forest_select <- subset(pca_S2_NC_SS_forest, sel_axes)
 
 # Best solution :exactextractr package with parallel computing : ~ 3h for > 11 000 for 100*100 cells on a 10 PCA layer stacks 
 n_cells <- nrow(cells_over_points_wgs)
 # n_cells <- 10
 library(doSNOW)
 
 cl <- makeCluster(3)
 registerDoSNOW(cl)
 clusterEvalQ(cl, library(raster))
 clusterEvalQ(cl, library(sf))
 clusterEvalQ(cl, library(exactextractr))
 pb <- txtProgressBar(max = n_cells, style = 3)
 progress <- function(n) setTxtProgressBar(pb, n)
 opts <- list(progress = progress)
 system.time(
   pca_spectral_com_SH <- foreach(i = 1:n_cells, .combine = c, 
                                  .options.snow = opts) %dopar%
     {
       pca_spectral_com_SH_tmp <- exact_extract(pca_S2_NC_SH_forest_select, cells_over_points_wgs_sf[i,])
       pca_spectral_com_SH_tmp[[1]]$ID_cell <- cells_over_points_wgs_sf[i,]$ID_cell
       return(pca_spectral_com_SH_tmp)
     }
 )
 close(pb)
 stopCluster(cl)
 # saveRDS(pca_spectral_com_SH, paste0(dir_to_save,"pca_235_spectral_com_SH.rds"))
 # saveRDS(pca_spectral_com_SH, paste0(dir_to_save,"pca_235_spectral_com_SH_cell50m.rds"))
 
 cl <- makeCluster(3)
 registerDoSNOW(cl)
 clusterEvalQ(cl, library(raster))
 clusterEvalQ(cl, library(sf))
 clusterEvalQ(cl, library(exactextractr))
 pb <- txtProgressBar(max = n_cells, style = 3)
 progress <- function(n) setTxtProgressBar(pb, n)
 opts <- list(progress = progress)
 system.time(
   pca_spectral_com_SS <- foreach(i = 1:n_cells, .combine = c, 
                                  .options.snow = opts) %dopar%
     {
       pca_spectral_com_SS_tmp <- exact_extract(pca_S2_NC_SS_forest_select, cells_over_points_wgs_sf[i,])
       pca_spectral_com_SS_tmp[[1]]$ID_cell <- cells_over_points_wgs_sf[i,]$ID_cell
       return(pca_spectral_com_SS_tmp)
     }
 )
 close(pb)
 stopCluster(cl)
 # saveRDS(pca_spectral_com_SS, paste0(dir_to_save,"pca_235_spectral_com_SS.rds"))
 # saveRDS(pca_spectral_com_SS, paste0(dir_to_save,"pca_235_spectral_com_SS_cell50m.rds"))
 
 #####################################################################################################
 #####################################################################################################
 ####  use spectral communities to compute FD indices from PCA axes  ####
 #####################################################################################################
 #####################################################################################################
 
 #### choose PCA raster nere ####
 pca_spectral_com_SS <- readRDS(paste0(dir_to_save,"pca_235_spectral_com_SS.rds"))
 pca_spectral_com_SS <- readRDS(paste0(dir_to_save,"pca_235_spectral_com_SS_cell50m.rds"))
 
 ###  SCALE PCA VALUES HERE! ###
 pca_spectral_com_SS_df <- do.call(rbind, pca_spectral_com_SS )
 means_pca <- apply(pca_spectral_com_SS_df[1:(ncol(pca_spectral_com_SS_df)-2)],2, mean, na.rm = T)
 sd_pca <- apply(pca_spectral_com_SS_df[1:(ncol(pca_spectral_com_SS_df)-2)],2, sd, na.rm = T)
 for(i in 1:length(pca_spectral_com_SS)) {
   tmp <- pca_spectral_com_SS[[i]][,1:(ncol(pca_spectral_com_SS[[i]])-2)] 
   for(j in 1:ncol(tmp))  tmp[,j] <- (tmp[,j] - means_pca[j]) / sd_pca[j]
   pca_spectral_com_SS[[i]][,1:(ncol(pca_spectral_com_SS[[i]])-2)] <- tmp
 }
 
 #### KEEP ONLY CELLS WITH coverage_fraction == 1 !!!! ### 
 
 for (i in 1:length(pca_spectral_com_SS)) pca_spectral_com_SS[[i]] <-  pca_spectral_com_SS[[i]][pca_spectral_com_SS[[i]]$coverage_fraction == 1,]
 
 #### Compute "functional" diversity indices  ####
 pca_spectral_com_SS_ok <- pca_spectral_com_SS
 
 library(FD)
 library(doSNOW)
 
 n_spectral_coms <- length(pca_spectral_com_SS_ok)
 cl <- makeCluster(3)
 registerDoSNOW(cl)
 clusterEvalQ(cl, library(FD))
 pb <- txtProgressBar(max = n_spectral_coms, style = 3)
 progress <- function(n) setTxtProgressBar(pb, n)
 opts <- list(progress = progress)
 system.time(
   FD_spectral_SS <- foreach(i = 1:n_spectral_coms,
                             .options.snow = opts) %dopar%
     {FD_tmp <- tryCatch(
       {
         com_tmp <- pca_spectral_com_SS_ok[[i]]
         # keep only the 100 cells with a coverage fraction of ~1 
         com_tmp <- com_tmp[com_tmp$coverage_fraction == 1,]
         if(!nrow(com_tmp) == 100) print(paste0("not 100 cells for cell ", i, "!!!"))
         ID_point <-  com_tmp$ID_point[1]
         ID_cell <-  com_tmp$ID_cell[1]
         # keep only pca colunms 
         com_tmp <- com_tmp[,1:(ncol(com_tmp)-2)]
         # no NA (but if thare are NA the spectral community is not complete...) 
         com_tmp <- com_tmp[complete.cases(com_tmp),]
         FD_tmp <- dbFD(com_tmp, stand.x = F)
       }, error = function(e) {NA})
     if(is.na(FD_tmp)) names(FD_tmp) <- "nbsp"
     FD_tmp$ID_point <-  ID_point
     FD_tmp$ID_cell <-  ID_cell
     names(FD_tmp$nbsp) <- "Community1"
     return(FD_tmp)
     }
 )
 close(pb)
 stopCluster(cl)
 
 # saveRDS(FD_spectral_SS, paste0(dir_to_save,"pca_235_FD_spectral_com_SS.rds"))
 # saveRDS(FD_spectral_SS, paste0(dir_to_save,"pca_235_FD_spectral_com_50m_SS.rds"))
 
 
 #####################################################################################################
 #####################################################################################################
 ####  Analyze"functional" diversity indices  ####
 #####################################################################################################
 #####################################################################################################
 
 FD_spectral_SS <- readRDS( paste0(dir_to_save,"pca_235_FD_spectral_com_SS.rds"))
 FD_spectral_SS <- readRDS( paste0(dir_to_save,"pca_235_FD_spectral_com_50m_SS.rds"))
 
 FD_spectral <- FD_spectral_SS
 # transform to data.frame
 FD_spectral <- do.call(dplyr::bind_rows,lapply(FD_spectral, function(x) unlist(lapply(x, as.numeric))))
 FD_spectral <- as.data.frame(FD_spectral)
#### get cells in forest ####
 cells_over_points_wgs <- readRDS( paste0(dir_to_save,"cells50_inforest.rds"))
 FD_spectral <- FD_spectral[FD_spectral$ID_cell %in% cells_over_points_wgs$ID_cell,]
 
 #### merge with data on fragmentation ####
 # 100*100m cells 
 centroids_cells_100_points_wgs <- readRDS(paste0(dir_to_save,"centroids_cells_100_points_wgs_landscape_metrics_topovar.rds"))
 # 50*50m cells 
 centroids_cells_100_points_wgs <- readRDS(paste0(dir_to_save,"centroids_cells_50_points_wgs_landscape_metrics_topovar.rds"))
 centroids_cells_100_points_wgs_df <- as.data.frame(centroids_cells_100_points_wgs)
 
 # FD_spectral_frag <- merge(FD_spectral_ok,centroids_cells_100_points_wgs_df, by = "ID_cell")
 FD_spectral_frag <- cbind(FD_spectral,centroids_cells_100_points_wgs_df)
 
 # keep only spectral communities with 70 sp.
 #  FD_spectral_ok <- FD_spectral[FD_spectral$nbsp >= 70,] # for 100sp
 #  FD_spectral_ok <- FD_spectral[FD_spectral$nbsp >= 20,] # for 50sp
 
 dim(FD_spectral_ok)
 
 ####################################################
 #### model selection ####
 ####################################################
 library(MuMIn)
 ##### simple model selection for landscape buffers ####
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
 
 ## choose data to analyse 
 # 
 # more than 80 calls/100 with values
 data <- FD_spectral_frag
 data <- data[, !colnames(data) %in% c("geometry.x",  "geometry.y")]
 data <- data.frame(data)
 # no points < d m from edge
 d =50
 d =25
 
 data <- data[data$dist_edge_centroid > d,]
 dim(data)
 colnames(data)
 
 # add log edge
 data$log_edge_centroid <- log(data$dist_edge_centroid+1)

 ##### preliminary tests #####
 data_test <- data[complete.cases(data$FDis), ]
 data_test <- data_test[complete.cases(data_test$RaoQ), ]
 
 plot(data_test$FDis~data_test$RaoQ)
 cor(data_test$FDis,data_test$RaoQ)
 
 data_scaled <- as.data.frame(data)
 for(i in colnames(data_scaled)[!colnames(data_scaled) %in% c("nbsp", "sing.sp", "ID_cell", "ID_point", "id_patch", "geometry")]){
    if(is.numeric( data_scaled[,i])){data_scaled[,i] <- scale(data_scaled[,i])
    }else{
       data_scaled[,i] <- scale(as.numeric(as.factor(data_scaled[,i])))
    }
 }
 
 plot(data_scaled$FDiv~data_scaled$log_edge_centroid)
 summary(lm(data_scaled$FDiv~data_scaled$log_edge_centroid))
 
 plot(data_scaled$FDis~data_scaled$log_edge_centroid)
 summary(lm(data_scaled$FDis~data_scaled$log_edge_centroid))
 
 plot(data_scaled$FRic~data_scaled$log_edge_centroid)
 summary(lm1 <- lm(data_scaled$FRic~data_scaled$log_edge_centroid))
 abline(lm1, col = "blue")
 
 plot(data_scaled$RaoQ~data_scaled$log_edge_centroid)
 summary(lm1 <- lm(data_scaled$RaoQ~data_scaled$log_edge_centroid))
 abline(lm1, col = "blue")
 
 plot(data_scaled$CWM1~data_scaled$log_edge_centroid)
 summary(lm1 <- lm(data_scaled$CWM1~data_scaled$log_edge_centroid))
 abline(lm1, col = "blue")
 
 plot(data_scaled$CWM2~data_scaled$log_edge_centroid)
 summary(lm1 <- lm(data_scaled$CWM2~data_scaled$log_edge_centroid))
 abline(lm1, col = "blue")
 
 plot(data_scaled$CWM3~data_scaled$log_edge_centroid)
 summary(lm1 <- lm(data_scaled$CWM3~data_scaled$log_edge_centroid))
 abline(lm1, col = "blue")
 
 ##### model selection #####
 # total habitat
 aic_sel_HA <- msel_aic_uni(data$FDiv, data[,c( paste0(rep("prop.landscape_", 4), c(250,500,1000,2000), rep("_centroid", 4)))])
 aic_sel_HA
 
 # effective mesh size
 aic_sel_EMS <- msel_aic_uni(data$FDiv, data[,c( paste0(rep("effective.mesh.size_", 4), c(250,500,1000,2000), rep("_centroid", 4)))])
 aic_sel_EMS
 
 colnames(data)
 land_var <- paste0(rep(c("total.edge_", "total.area_", "effective.mesh.size_", "prop.landscape_"), each = 4), c(250,500,1000,2000), rep("_centroid", 4))
 # from selection
 land_var <- c(aic_sel_HA$pred, aic_sel_EMS$pred)
 
 topo_var <- c("elevation","slope","aspect" ,"curvature","substrat_centroid")
 sel_var <- c("dist_edge_centroid","log_edge_centroid",land_var, topo_var)
 sel_var1 <- c("dist_edge_centroid",land_var, topo_var)
 sel_var2 <- c("log_edge_centroid",land_var, topo_var)
 
 ##### choose model here #####
 ## FRic ##
 data_for_mod <- data[, c("FRic", sel_var1)]
 data_for_mod <- data[, c("FRic", sel_var2)]
 
 # standardize variables 
 data_for_mod_scaled <- data_for_mod
 for(i in colnames(data_for_mod_scaled)[!colnames(data_for_mod_scaled) %in% c("nbsp", "sing.sp", "ID_cell", "ID_point", "id_patch")]){
    if(is.numeric( data_for_mod_scaled[,i])){data_for_mod_scaled[,i] <- scale(data_for_mod_scaled[,i])
    }else{
       data_for_mod_scaled[,i] <- scale(as.numeric(as.factor(data_for_mod_scaled[,i])))
    }
 }
 # no na
 data_for_mod_scaled <- data_for_mod_scaled[complete.cases(data_for_mod_scaled),]
 dim(data_for_mod_scaled)
 # total model
 glm1 <- glm(FRic~ . , data = data_for_mod_scaled, na.action = "na.fail")
 # Variance Inflation Factors
 car::vif(glm1)
 VIF <- car::vif(glm1)
 data_for_mod_scaled_new <- data_for_mod_scaled
 while(any(VIF>5)){
    omit_var <- names(sort(VIF, decreasing = T))[1]
    data_for_mod_scaled_new <- data_for_mod_scaled_new[! colnames(data_for_mod_scaled_new) %in% omit_var]
    glm1 <- glm(FRic~ . , data = data_for_mod_scaled_new, na.action = "na.fail")
    VIF <- car::vif(glm1)
 }
 
 ## FDis ##
 data_for_mod <- data[, c("FDis", sel_var2)]
 # standardize variables 
 data_for_mod_scaled <- data_for_mod
 for(i in colnames(data_for_mod_scaled)[!colnames(data_for_mod_scaled) %in% c("nbsp", "sing.sp", "ID_cell", "ID_point", "id_patch")]){
    if(is.numeric( data_for_mod_scaled[,i])){data_for_mod_scaled[,i] <- scale(data_for_mod_scaled[,i])
    }else{
       data_for_mod_scaled[,i] <- scale(as.numeric(as.factor(data_for_mod_scaled[,i])))
    }
 }
 # no na
 data_for_mod_scaled <- data_for_mod_scaled[complete.cases(data_for_mod_scaled),]
 dim(data_for_mod_scaled)
 # total model
 glm1 <- glm(FDis~ . , data = data_for_mod_scaled, na.action = "na.fail")
 # Variance Inflation Factors
 car::vif(glm1)
 VIF <- car::vif(glm1)
 data_for_mod_scaled_new <- data_for_mod_scaled
 while(any(VIF>5)){
    omit_var <- names(sort(VIF, decreasing = T))[1]
    data_for_mod_scaled_new <- data_for_mod_scaled_new[! colnames(data_for_mod_scaled_new) %in% omit_var]
    glm1 <- glm(FDis~ . , data = data_for_mod_scaled_new, na.action = "na.fail")
    VIF <- car::vif(glm1)
 }
 
 ## RaoQ ##
 data_for_mod <- data[, c("RaoQ", sel_var2)]
 # standardize variables 
 data_for_mod_scaled <- data_for_mod
 for(i in colnames(data_for_mod_scaled)[!colnames(data_for_mod_scaled) %in% c("nbsp", "sing.sp", "ID_cell", "ID_point", "id_patch")]){
    if(is.numeric( data_for_mod_scaled[,i])){data_for_mod_scaled[,i] <- scale(data_for_mod_scaled[,i])
    }else{
       data_for_mod_scaled[,i] <- scale(as.numeric(as.factor(data_for_mod_scaled[,i])))
    }
 }
 # no na
 data_for_mod_scaled <- data_for_mod_scaled[complete.cases(data_for_mod_scaled),]
 dim(data_for_mod_scaled)
 # total model
 glm1 <- glm(RaoQ~ . , data = data_for_mod_scaled, na.action = "na.fail")
 # Variance Inflation Factors
 car::vif(glm1)
 VIF <- car::vif(glm1)
 data_for_mod_scaled_new <- data_for_mod_scaled
 while(any(VIF>5)){
    omit_var <- names(sort(VIF, decreasing = T))[1]
    data_for_mod_scaled_new <- data_for_mod_scaled_new[! colnames(data_for_mod_scaled_new) %in% omit_var]
    glm1 <- glm(RaoQ~ . , data = data_for_mod_scaled_new, na.action = "na.fail")
    VIF <- car::vif(glm1)
 }
 
 ## FDiv ##
 data_for_mod <- data[, c("FDiv", sel_var2)]
 # standardize variables 
 data_for_mod_scaled <- data_for_mod
 for(i in colnames(data_for_mod_scaled)[!colnames(data_for_mod_scaled) %in% c("nbsp", "sing.sp", "ID_cell", "ID_point", "id_patch")]){
    if(is.numeric( data_for_mod_scaled[,i])){data_for_mod_scaled[,i] <- scale(data_for_mod_scaled[,i])
    }else{
       data_for_mod_scaled[,i] <- scale(as.numeric(as.factor(data_for_mod_scaled[,i])))
    }
 }
 # no na
 data_for_mod_scaled <- data_for_mod_scaled[complete.cases(data_for_mod_scaled),]
 dim(data_for_mod_scaled)
 # total model
 glm1 <- glm(FDiv~ . , data = data_for_mod_scaled, na.action = "na.fail")
 # Variance Inflation Factors
 car::vif(glm1)
 VIF <- car::vif(glm1)
 data_for_mod_scaled_new <- data_for_mod_scaled
 while(any(VIF>5)){
    omit_var <- names(sort(VIF, decreasing = T))[1]
    data_for_mod_scaled_new <- data_for_mod_scaled_new[! colnames(data_for_mod_scaled_new) %in% omit_var]
    glm1 <- glm(FDiv~ . , data = data_for_mod_scaled_new, na.action = "na.fail")
    VIF <- car::vif(glm1)
 }
 
 ##### Model selection #####
 library(MuMIn)
 library(parallel)
 cl <- makeCluster(3)
 clusterExport(cl, "data_for_mod_scaled_new")
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
 # selection few best models based on delta AIC
 models.list <- get.models(dd,subset =  delta < 2) 
 
 # average models and coeficiants based on weigth
 if(length(models.list)>1){
    muminavg <- model.avg(models.list, revised.var = TRUE, fit = TRUE, beta = "partial.sd")
 } else {muminavg <- models.list[[1]]}
 summary(muminavg)
 
 # transform AIC-weigthed standardized coeficients
 weighted_coef <- coef(muminavg)
 sign_coef <- weighted_coef
 # transform to relative influence 
 relat_inf <- abs(weighted_coef)/sum(abs(weighted_coef))
 relat_inf <- relat_inf[2:length(relat_inf)]
 
 # get conf. intervals for conditional average (selected models) 
 weighted_confint <- confint(muminavg)
 
 #### plot ####
 library(GGally)
 
 df_coef <- data.frame(term = names(weighted_coef[2:length(weighted_coef)]),
                       estimate = weighted_coef[2:length(weighted_coef)],
                       conf.low = weighted_confint[2:length(weighted_coef),1],
                       conf.high = weighted_confint[2:length(weighted_coef),2])
 
 # change names 
 # df_coef$term <- c("aspect","effective mesh size", "elevation","distance edge","habitat amount", "slope",                   
 #                    "substrat","curvature")
 
 ggcoef(df_coef, sort = "ascending", vline_linetype =  "dotted", mapping = aes(x = estimate, y = term, colour = term),
        errorbar_height = .1) + 
    theme_bw() + 
    theme(legend.position="none",
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15)) +
    xlab('standardized coefficients')
 
 #### nice plot edge  ####
 lm1 <- lm(FDis ~ log_edge_centroid, data = data_for_mod_scaled_new)

 lm_topo <- lm(FDis ~ slope + aspect + elevation + substrat_centroid + curvature, data = data_for_mod_scaled_new)
 smry_topo = summary(lm_topo)
 smry_topo
 lm_topo_edge <- lm(FDis ~ slope + aspect + elevation + substrat_centroid + curvature + log_edge_centroid, data = data_for_mod_scaled_new)
 smry_topo_edge = summary(lm_topo_edge)
 smry_topo_edge
 
 data_res <- data.frame(FDis_res = lm_topo$residuals, log_edge =  data_for_mod_scaled_new$log_edge_centroid)
 lm_res <- lm(FDis_res ~ log_edge, data = data_res)
 smry_res = summary(lm_res)
 smry_res
 data_plot1 <- data_res
 plot <- ggplot(data = data_plot1, aes(x = log_edge, y = FDis_res)) + 
    geom_point(color='navy', size = 1, alpha =.3) +
    geom_smooth(data = data_plot1, aes(x = log_edge, y = FDis_res), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
    annotate("text", x=2, y= 10,
             label= paste0("R²=", round(smry_res$r.squared, digits = 3), gtools::stars.pval(smry_res$coefficients[8])),
             color='red') +
    ylab("FDis_res") + xlab("Log distance edge") +
    theme(
       axis.text.x = element_text(size=12),
       axis.text.y = element_text(size=12)
    ) +
    theme_bw()
 
 plot
 #####################################################################################################
 #####################################################################################################
 ####  load k-means image from GEE ==> get previously saved masked image just below if already done!!!!  ####
 #####################################################################################################
 ##################################################################################################### 
 #### select clustered image of dry season ####
 
 im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons"
 im_name <- "kmeans_S2_NC_SH_50grp_pca2to5"
 im_name <- "kmeans_S2_NC_SH_50grp_pca235"
 
 # kmeans groups
 kmeans_S2 <- raster( paste0(im_path, "/", im_name, ".tif") )
 
 #### mask raster with forest #### 
 # get forest raster
 forest_map_raster_withna <- raster(
         "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT_withna.tif")
 
 kmeans_S2_forest <- mask(kmeans_S2, forest_map_raster_withna)
 
 # writeRaster(kmeans_S2_forest,  paste0(dir_to_save,"kmeans_S2_SH_50grp_pca235_forest.tif"),
 #             datatype="INT1U", options="COMPRESS=LZW", overwrite=TRUE)
 
 im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons"
 im_name <- "kmeans_S2_NC_SS_50grp_pca2to5"
 im_name <- "kmeans_S2_NC_SS_50grp_pca235"
 
 # kmeans groups
 kmeans_S2 <- raster( paste0(im_path, "/", im_name, ".tif") )
 
 #### mask raster with forest #### 
 # get forest raster
 forest_map_raster_withna <- raster(
    "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT_withna.tif")
 
 kmeans_S2_forest <- mask(kmeans_S2, forest_map_raster_withna)
 
 # writeRaster(kmeans_S2_forest,  paste0(dir_to_save,"kmeans_S2_SS_50grp_pca235_fores_forest.tif"),
 #             datatype="INT1U", options="COMPRESS=LZW", overwrite=TRUE)
 
 
 #####################################################################################################
 ####  load K-means Spectral Species from GEE with forest mask  ####
 #####################################################################################################
 im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/2_seasons"
 
 im_name <- "kmeans_S2_NC_SH_50grp_pca2to5_forest"
 
 kmeans_S2_NC_SH_50grp_forest <- raster( paste0(im_path, "/", im_name, ".tif") )
 
 im_name <- "kmeans_S2_NC_SS_50grp_pca2to5_forest"
 
 kmeans_S2_NC_SS_50grp_forest <- raster( paste0(im_path, "/", im_name, ".tif") )
 
 ##########################################################################################
 #### get cells over points in the grid with resolution * 10 (10m -> 100m), generated in PCA analyses  ####
 ##########################################################################################
 cells_over_points_wgs <- readRDS( paste0(dir_to_save,"cells_over_points_wgs_unique.rds"))
 
 centroids_cells_100_points_wgs <- readRDS( paste0(dir_to_save,"centroids_cells_100_points_wgs_landscape_metrics_topovar.rds"))
 
 ##############################################################################
 #### extract spectral species communities  ####
 ##############################################################################
 
 n_cells <- nrow(cells_over_points_wgs)
 # n_cells <- 10
 library(doSNOW)
 
 cl <- makeCluster(3)
 registerDoSNOW(cl)
 clusterEvalQ(cl, library(raster))
 clusterEvalQ(cl, library(sf))
 clusterEvalQ(cl, library(exactextractr))
 pb <- txtProgressBar(max = n_cells, style = 3)
 progress <- function(n) setTxtProgressBar(pb, n)
 opts <- list(progress = progress)
 system.time(
   species_spectral_com_SH <- foreach(i = 1:n_cells, 
                                      .options.snow = opts) %dopar%
     {
       species_spectral_com_SH_tmp <- exact_extract(kmeans_S2_NC_SH_50grp_forest, cells_over_points_wgs[i,])
       species_spectral_com_SH_tmp <- species_spectral_com_SH_tmp[[1]]
       species_spectral_com_SH_tmp <- data.frame(spectral_species = species_spectral_com_SH_tmp$value[species_spectral_com_SH_tmp$coverage_fraction == 1])
       species_spectral_com_SH_tmp$ID_cell <- cells_over_points_wgs[i,]$ID_cell
       species_spectral_com_SH_tmp$ID_point <- cells_over_points_wgs[i,]$ID
       return(species_spectral_com_SH_tmp)
     }
 )
 close(pb)
 stopCluster(cl)
 # saveRDS(species_spectral_com_SH, paste0(dir_to_save,"species_spectral_pca2to5_com_SH.rds"))
 
 # saveRDS(species_spectral_com_SS, paste0(dir_to_save,"species_spectral_pca235_com_SS.rds"))
 
 cl <- makeCluster(3)
 registerDoSNOW(cl)
 clusterEvalQ(cl, library(raster))
 clusterEvalQ(cl, library(sf))
 clusterEvalQ(cl, library(exactextractr))
 pb <- txtProgressBar(max = n_cells, style = 3)
 progress <- function(n) setTxtProgressBar(pb, n)
 opts <- list(progress = progress)
 system.time(
   species_spectral_com_SS <- foreach(i = 1:n_cells, 
                                      .options.snow = opts) %dopar%
     {
       species_spectral_com_SS_tmp <- exact_extract(kmeans_S2_NC_SS_50grp_forest, cells_over_points_wgs[i,])
       species_spectral_com_SS_tmp <- species_spectral_com_SS_tmp[[1]]
       species_spectral_com_SS_tmp <- data.frame(spectral_species = species_spectral_com_SS_tmp$value[species_spectral_com_SS_tmp$coverage_fraction == 1])
       species_spectral_com_SS_tmp$ID_cell <- cells_over_points_wgs[i,]$ID_cell
       species_spectral_com_SS_tmp$ID_point <- cells_over_points_wgs[i,]$ID
       return(species_spectral_com_SS_tmp)
     }
 )
 close(pb)
 stopCluster(cl)
 # saveRDS(species_spectral_com_SS, paste0(dir_to_save,"species_spectral_pca2to5_com_SS.rds"))
 # saveRDS(species_spectral_com_SS, paste0(dir_to_save,"species_spectral_pca235_com_SS.rds"))
 

 ##############################################################################
 ##############################################################################
 #### alpha & beta diversity  spectral community matrix ####
 ##############################################################################
 ##############################################################################
 

 #### get spectral communities ####
 # species_spectral_com_SH <- readRDS( "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/data_frag/species_spectral_pca2to5_com_SH.rds")
 
 species_spectral_com_SS <- readRDS(paste0(dir_to_save,"species_spectral_pca2to5_com_SS.rds")); season = "SS"
 species_spectral_com_SS <- readRDS(paste0(dir_to_save,"species_spectral_pca235_com_SS.rds")); season = "SS"
 
 #### get a spectral community matrix ####
 species_spectral_com_df <- do.call(rbind.data.frame, species_spectral_com_SS)
 
 com_mat_0 <- t(as.matrix(table(species_spectral_com_df[,2:1])))
 com_mat <- as.matrix.data.frame(com_mat_0)
 rownames(com_mat) <- rownames(com_mat_0)
 colnames(com_mat)<- colnames(com_mat_0)
 
 ##############################################################################
 ##############################################################################
 #### alpha diversity  ####
 ##############################################################################
 ##############################################################################
 
 #### species richness, Shannon, Simpson (Hill numbers) ####
 library(iNEXT)
 
 # keep only spectral communities with > 20 species
 # com_mat <- com_mat[rowSums(com_mat) > 20,]
 # com_mat <- com_mat[,colSums(com_mat) > 20]

 # set a series of sample sizes (m) for R/E computation
 m <- c(50, 100)
 
 spectral_alphadiv <- iNEXT(com_mat, q=c(0,1,2), datatype="abundance", size=m)
 
 # saveRDS(spectral_alphadiv, paste0(dir_to_save,"spectral_alphadiv", season, ".rds"))
 
 
 df_alpha_div <- c()
 for (i in 1:length(spectral_alphadiv$iNextEst)){
   tmp <- spectral_alphadiv$iNextEst[[i]]
   tmp50Q0 <- tmp[tmp$m == 50 & tmp$order == 0,]
   tmp50Q1 <- tmp[tmp$m == 50 & tmp$order == 1,]
   tmp50Q2 <- tmp[tmp$m == 50 & tmp$order == 2,]
   df_alpha_div <- rbind(df_alpha_div,
                            c(tmp50Q0$qD[1], tmp50Q1$qD[1], tmp50Q2$qD[1], names(spectral_alphadiv$iNextEst)[[i]]))
 }
 ##############################################################################
 ##############################################################################
 ####  beta diversity ####
 ##############################################################################
 ##############################################################################
 library(vegan)
 
 com_mat_for_vegan <- t(com_mat)
 
 BC_dissimilarity <- vegdist(com_mat_for_vegan)
 dim(BC_dissimilarity)
 
 # saveRDS(BC_dissimilarity, paste0(dir_to_save,"BC_dissimilarity_pca2to5_com", season, ".rds"))
 #  BC_dissimilarity <- readRDS(paste0(dir_to_save,"BC_dissimilarity_pca2to5_com", season, ".rds"))
 
 # saveRDS(BC_dissimilarity, paste0(dir_to_save,"BC_dissimilarity_pca235_com", season, ".rds"))
 #  BC_dissimilarity <- readRDS(paste0(dir_to_save,"BC_dissimilarity_pca235_com", season, ".rds"))
 
 #### get mean betadiv per spectral com ####
 BC_dissimilarity <- as.matrix(BC_dissimilarity)
 diag(BC_dissimilarity) <- 0
 
 mean_BC_com <- rowMeans(BC_dissimilarity)
 mean_BC_com_df <- data.frame(ID_cell = names(mean_BC_com), mean_BC = mean_BC_com)

  ##############################################################################
 ####  analyses beta diversity ~ fragmentation ####
 ##############################################################################
 # use centroids of 100*100 cells
 beta_div_SH_frag <- merge(mean_BC_com_df, centroids_cells_100_points_wgs, by = "ID_cell")
 
 data <- beta_div_SH_frag
 #### keep only points > distmin from edge ####
 #### IMPORTANT: as cells are 100*100m (100m wide) they might include open areas (not forest)
 #### set minimum to 50m to get only cells with mostly forest areas 
 distmin <- 100
 distmin <- 60
 
 data <- data[data$dist_edge_centroid>distmin,]
 dim(data)
 data$log_edge_centroid <- log(data$dist_edge_centroid)
 lm1 <- lm(data$mean_BC~data$log_edge_centroid)
 summary(lm1)
 plot(data$mean_BC~data$log_edge_centroid)
 abline(lm1, col = "blue")
 
 ####################################################
 #### model selection ####
 ####################################################
 
 library(MuMIn)
 
 ##### simple model selection for landscape buffers ####
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
 # total habitat
 aic_sel_HA <- msel_aic_uni(data$mean_BC, data[,c( paste0(rep("prop.landscape_", 4), c(250,500,1000,2000), rep("_centroid", 4)))])
 aic_sel_HA
 
 # effective mesh size
 aic_sel_EMS <- msel_aic_uni(data$mean_BC, data[,c( paste0(rep("effective.mesh.size_", 4), c(250,500,1000,2000), rep("_centroid", 4)))])
 aic_sel_EMS
 
 colnames(data)
 land_var <- paste0(rep(c("total.edge_", "total.area_", "effective.mesh.size_", "prop.landscape_"), each = 4), c(250,500,1000,2000), rep("_centroid", 4))
 # from selection
 land_var <- c(aic_sel_HA$pred, aic_sel_EMS$pred)
 
 topo_var <- c("elevation","slope","aspect" ,"curvature","substrat_centroid")
 sel_var <- c("dist_edge_centroid","log_edge_centroid",land_var, topo_var)
 
 sel_var1 <- c("dist_edge_centroid",land_var, topo_var)
 
 sel_var2 <- c("log_edge_centroid",land_var, topo_var)
 
 
 ##### choose model here #####
 ## R ##
 data_for_mod <- data[, c("mean_BC", sel_var1)]
 
 data_for_mod <- data[, c("mean_BC", sel_var2)]
 
 # standardize variables 
 data_for_mod_scaled <- data_for_mod
 for(i in colnames(data_for_mod_scaled)[!colnames(data_for_mod_scaled) %in% c("nbsp", "sing.sp", "ID_cell", "ID_point", "id_patch")]){
    if(is.numeric( data_for_mod_scaled[,i])){data_for_mod_scaled[,i] <- scale(data_for_mod_scaled[,i])
    }else{
       data_for_mod_scaled[,i] <- scale(as.numeric(as.factor(data_for_mod_scaled[,i])))
    }
 }
 # no na
 data_for_mod_scaled <- data_for_mod_scaled[complete.cases(data_for_mod_scaled),]
 dim(data_for_mod_scaled)
 # total model
 glm1 <- glm(mean_BC ~ . , data = data_for_mod_scaled, na.action = "na.fail")
 # Variance Inflation Factors
 car::vif(glm1)
 VIF <- car::vif(glm1)
 data_for_mod_scaled_new <- data_for_mod_scaled
 while(any(VIF>5)){
    omit_var <- names(sort(VIF, decreasing = T))[1]
    data_for_mod_scaled_new <- data_for_mod_scaled_new[! colnames(data_for_mod_scaled_new) %in% omit_var]
    glm1 <- glm(mean_BC~ . , data = data_for_mod_scaled_new, na.action = "na.fail")
    VIF <- car::vif(glm1)
 }
 
 ##### Model selection #####
 library(MuMIn)
 library(parallel)
 cl <- makeCluster(3)
 clusterExport(cl, "data_for_mod_scaled_new")
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
 
 # selection few best models based on delta AIC
 models.list <- get.models(dd,subset =  delta < 2) 
 
 # average models and coeficiants based on weigth
 if(length(models.list)>1){
    muminavg <- model.avg(models.list, revised.var = TRUE, fit = TRUE, beta = "partial.sd")
 } else {muminavg <- models.list[[1]]}
 summary(muminavg)
 
 # transform AIC-weigthed standardized coeficients
 weighted_coef <- coef(muminavg)
 sign_coef <- weighted_coef
 # transform to relative influence 
 relat_inf <- abs(weighted_coef)/sum(abs(weighted_coef))
 relat_inf <- relat_inf[2:length(relat_inf)]
 
 # get conf. intervals for conditional average (selected models) 
 weighted_confint <- confint(muminavg)
 
 #### plot ####
 library(GGally)
 
 df_coef <- data.frame(term = names(weighted_coef[2:length(weighted_coef)]),
                       estimate = weighted_coef[2:length(weighted_coef)],
                       conf.low = weighted_confint[2:length(weighted_coef),1],
                       conf.high = weighted_confint[2:length(weighted_coef),2])
 
 # change names 
 # df_coef$term <- c("aspect","effective mesh size", "elevation","distance edge","habitat amount", "slope",                   
 #                    "substrat","curvature")
 
 ggcoef(df_coef, sort = "ascending", vline_linetype =  "dotted", mapping = aes(x = estimate, y = term, colour = term),
        errorbar_height = .1) + 
    theme_bw() + 
    theme(legend.position="none",
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15)) +
    xlab('standardized coefficients')
 
 #### nice plot   ####
 lm1 <- lm(mean_BC ~ prop.landscape_2000_centroid, data = data)
 smry1 = summary(lm1)
 data_plot1 <- data.frame(BC = data$mean_BC, HA= data$prop.landscape_2000_centroid)
 #data_plot1 <- data.frame(FRic = data$FRic, log_edge= data$dist_edge)
 plot <- ggplot(data = data_plot1, aes(x = HA, y = BC)) + 
    geom_point(color='navy', size = 1, alpha =.3) +
    geom_smooth(data = data_plot1, aes(x = HA, y = BC), formula = y ~ x, method = "lm", color='red', se=T, fill = "red" , alpha = 0.12) +
    annotate("text", x=.9, y= 1,
             label= paste0("R²=", round(smry1$r.squared, digits = 2), gtools::stars.pval(smry1$coefficients[8])),
             color='red') +
    ylab("BC") + xlab("Habitat amount (2000m)") +
    theme(
       axis.text.x = element_text(size=12),
       axis.text.y = element_text(size=12)
    ) +
    theme_bw()
 
 plot
 
 
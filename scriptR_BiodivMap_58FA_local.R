
####### Connect to server ####### 

# ssh -L 3499:localhost:3499 greg@niamoto


####### Packages ####### 
library("raster")
library("devtools")
# devtools::install_github('jbferet/biodivMapR')
library("biodivMapR")
library("mapview")


## rasterOptions(tmpdir="D:/A_traiter/R_temp")
system.file(package = "biodivMapR")
############################################# 
############### 1/ Get S2 images ############### 
############################################# 

# note that there is a function in BiodivMapR to merge multiple images : build_image_from_list 

# S2 images
im_path <- "/home/greg/projects/reliques/signature_spectrale_fragmentation/maps/biodivmapr/UM"
im_path <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021"

setwd(im_path)
list.files()

# this is an image with topographic correction in GEE
im_name <- "T58KFA_int_forest_zone_PDL"
# this is an image with raw data
im_name <- "T58KFA_raw_int_forest_zone_PDL"

sufix <- ".tif"
im <- stack(paste0(im_path, "/", im_name, sufix))
names(im) <- c("B2",     "B3",     "B4",     "B5",     "B6",     "B7",     "B8",    "B8A",    "B11",    "B12")

# #### generate mask from the first layer (using gdal rasterize with original S2 image) ####
# mask <- raster("/home/greg/projects/reliques/signature_spectrale_fragmentation/maps/biodivmapr/mask/image_dehazed_ok.tif")
# get the NIR band (B8) to make the mask, as it seems to never have zero value on forest.
mask <- im[["B8"]]
# mask[!is.na(mask[])] <- 1 
# mask[is.na(mask[])] <- 0 
# mask[]<- 1
mask[!mask[] == 0] <- 1

### crop image to exclude large masked areas (as this seems to generate zeros in the PCA) ####
e <- extent(im)
e <- as.vector(e)
# e[3] <- 7513000
e <- as(extent(e), 'SpatialPolygons')
crs(e) <- "+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

crop_im =  crop(im, e)
crop_mask =  crop(mask, e)

im <- crop_im
mask <- crop_mask
#### select S2 spectral bands to match with HDR template (exclude band with 60m resolution) ####
sel_bands <- c("B2",     "B3",     "B4",     "B5",     "B6",     "B7",     "B8",    "B8A",    "B11",    "B12")
im_ok <- im[[sel_bands]]
im <- im_ok

# convert image to ENVI format (keep INT2U format INT2S (-32,767	32,767) or INT2U (0 - 65534))
if(!file.exists(paste0(im_path, "/", im_name, ".envi"))){
  writeRaster(im, paste0(im_path, "/", im_name), format="ENVI", datatype = "INT2S", options="INTERLEAVE=BSQ", overwrite=TRUE) 
}

# convert mask to ENVI format (keep INT2U format (0 - 65534))
if(!file.exists(paste0(im_path, "/", im_name, "_mask.envi"))){
  writeRaster(mask, paste0(im_path, "/", im_name, "_mask"), format="ENVI", datatype="INT2S", options="COMPRESS=LZW", overwrite=TRUE)
}
# get .envi file mask
mask <- raster(paste0(im_path, "/", im_name, "_mask.envi"))

# complete the HDR file from package templates (own function)
# source("/home/greg/projects/reliques/signature_spectrale_fragmentation/R_scripts/write_hdr_modif.R")
# HDRpath <- paste0(im_path, "/", im_name, ".hdr")
# Sensor <- 'SENTINEL_2A'
# comp_ENVI_header_from_template(HDRpath,Sensor)

# select image with hdr file
sufix <- ".envi"

##########################################
##########################################
############## PCA ############## 
##########################################
##########################################

###### Computing options ##### 
# where to store results
Output_Dir = paste0(im_path,"/biodivmapr_local/RESULTS")

Input_Image_File <-  paste0(im_path, "/", im_name, sufix)
file.exists(Input_Image_File)

Input_Mask_Path <- paste0(im_path, "/", im_name, "_mask.envi")
file.exists(Input_Mask_Path)

im <- stack(Input_Image_File)
mask <- raster(Input_Mask_Path)

Output_Dir = paste0(im_path,"/biodivmapr_local/RESULTS")

############## Radiometric filterings before PCA ############## 
# NDVI_Thresh allows filtering to eliminate non-vegetated pixels. Nothing fancy so you may need to deal with mixed pixels…
NDVI_Thresh <- 0.25
# Blue_Thresh allows filtering of clouds, based on the hypothesis that atmospheric scattering will lead to higher reflectance in the blue domain. 
# Blue_Thresh defines the maximum Blue reflectance to be kept.
Blue_Thresh <- NULL
# NIR_Thresh allows filtering of shadows and pixels with very low signal. NIR_Thresh defines the minimum NIR value to be kept.
NIR_Thresh <- 1500

# print("PERFORM RADIOMETRIC FILTERING")
# Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File, Mask_Path = Input_Mask_File,
#                                                  Output_Dir = Output_Dir, TypePCA = TypePCA,
#                                                  NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
#                                                  NIR_Thresh = NIR_Thresh)

############## compute PCA ############## 

parallel::detectCores()

nbCPU         = 4
MaxRAM        = 0.5
TypePCA = 'SPCA'

print("PERFORM PCA ON RASTER")
PCA_Output    <- perform_PCA(Input_Image_File, Input_Mask_Path, Output_Dir,
                             TypePCA = TypePCA, FilterPCA = T, nbCPU = nbCPU,MaxRAM = MaxRAM, Continuum_Removal = T)
#### you may have to add the sufix ".envi" manually ####


# correct pca values ?
# 
# pca <- stack(paste0(im_path, "/RESULTS/PCA_", im_name, "/", im_name, "/SPCA/PCA/OutputPCA_8_PCs"))
# 
# pca_correct <- pca
# system.time({for(i in 1:nlayers(pca_correct)){pca_correct[[i]][pca_correct[[i]] == 0] <- NA}})
# writeRaster(pca_correct,  paste0(im_path, "/RESULTS/PCA_", im_name, "/", im_name, "/SPCA/PCA/OutputPCA_8_PCs"),
#             format="ENVI", datatype="INT2U", overwrite=TRUE)
# 
# pca[pca == 0] <- NA
# writeRaster(pca,  paste0(im_path, "/RESULTS/PCA_", im_name, "/", im_name, "/SPCA/PCA/OutputPCA_8_PCs"),
#           format="ENVI", datatype="INT2U", overwrite=TRUE)

##########################################
############## Spectral species ############## 
##########################################
# path_pca_info <- paste0(im_path, "/RESULTS/PCA_", im_name, "/", im_name,"/SPCA/PCA/PCA_Info.RData")

path_pca_info <- paste0(im_path, "/biodivmapr_local/RESULTS/", im_name,"/SPCA/PCA/PCA_Info.RData")
file.exists(path_pca_info)
load(paste0(path_pca_info))

path_pca_info
# Path for the PCA raster
PCA_Files
# change PCA_Files
PCA_Files <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/OutputPCA_8_PCs"
# number of pixels used for each partition used for k-means clustering
Pix_Per_Partition 
# number of partitions used for k-means clustering
nb_partitions 
# Path for the updated mask
Input_Mask_File <- MaskPath
# change PCA_Files
Input_Mask_File <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/ShadeMask_Update_PCA"
# parameters of the PCA model
PCA_model   
# variance explained
summary(PCA_model)
# definition of spectral bands to be excluded from the analysis
SpectralFilter
# nuber of spectral species (clusters)
nbclusters    = 20

nbCPU         = 4
MaxRAM        = 0.5
TypePCA = 'SPCA'

Input_Image_File <-  paste0(im_path, "/", im_name, sufix)
Output_Dir = paste0(im_path,"/biodivmapr_local/RESULTS")

print("Select PCA components for diversity estimations")
# # CHOICE:  open Vim text editor (edit and quit using :q)
# select_PCA_components(Input_Image_File, Output_Dir, PCA_Files, File_Open = TRUE)
# OR: create and export table 
sel_PCs <- data.frame(c(1,2,6,7,8), check.names = F)
sel_PCs_path <- paste0(c(unlist(strsplit(PCA_Files, split="/"))[1:(length(unlist(strsplit(PCA_Files, split="/")))-1)],
                         "Selected_Components.txt"), collapse = "/")
write.table(sel_PCs, file = sel_PCs_path, row.names = F, col.names = F)

# check selected component .txt file
read.table(select_PCA_components(Input_Image_File, Output_Dir, PCA_Files, File_Open = F))

print("MAP SPECTRAL SPECIES")
map_spectral_species(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, PCA_Files = PCA_Files,
                     PCA_model = PCA_model, SpectralFilter = SpectralFilter, Input_Mask_File = Input_Mask_File,
                     Pix_Per_Partition = Pix_Per_Partition, nb_partitions = nb_partitions, nbCPU=nbCPU, MaxRAM=MaxRAM, 
                     nbclusters = nbclusters, TypePCA = TypePCA, Continuum_Removal = TRUE)

####################################################################################
####################################################################################
############## crop images if some zero are generated ############## 
####################################################################################
####################################################################################



####################################################################################
####################################################################################
############## Spectral communities ############## 
####################################################################################
####################################################################################

#### set a window size in pixel unit ####
window_size = 5 # 25 spectral species, 50*50m 
# window_size = 10 # 100 spectral species, 100*100m

################################################################################
##              MAP FUNCTIONAL DIVERSITY METRICS FRic, FEve, FDiv             ##
##          (Villeger et al, 2008 https://doi.org/10.1890/07-1206.1)          ##
################################################################################
## read selected features from dimensionality reduction
file.exists(paste0(im_path, "/biodivmapr_local/RESULTS/", im_name, "/SPCA/PCA/Selected_Components.txt"))
# Define PCs to be selected. Set to FALSE if you want to use the "Selected_Components.txt" file
Sel_PC <- paste0(im_path, "/biodivmapr_local/RESULTS/", im_name, "/SPCA/PCA/Selected_Components.txt")

Selected_Features <- read.table(Sel_PC)[[1]]
## path for selected components
map_functional_div(Original_Image_File = Input_Image_File, Functional_File = PCA_Files,
                   Selected_Features = Selected_Features, Output_Dir = Output_Dir,
                   window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM,TypePCA = TypePCA,
                   MinSun = 0.8)

# pca <- stack(paste0(im_path, "/RESULTS/PCA_", im_name, "/", im_name, "/SPCA/PCA/OutputPCA_8_PCs"))
# specsp <- raster(paste0(im_path, "/RESULTS/PCA_", im_name, "/", im_name, "/SPCA/SpectralSpecies/SpectralSpecies"))
# S2 <- stack(paste0(im_path,"/", im_name, ".envi"))
# ncell(pca)
# ncell(S2)
# ncell(specsp)
# spl_cell <- sample(1:ncell(pca), 1000)
# spl_pca <- pca[spl_cell]
# spl_S2 <- S2[spl_cell]
# spl_specsp <- specsp[spl_cell]
# cbind(spl_pca[1:1000,1],
# spl_S2[1:1000,1],
# spl_specsp[1:1000])
# spl_pca[complete.cases(spl_pca),]
# spl_S2[complete.cases(spl_S2),]
# spl_pca1 <- sample(pca,1000)
# spl_pca1[complete.cases(spl_pca1),]

print("MAP ALPHA DIVERSITY")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')
map_alpha_div(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, TypePCA = TypePCA,
              window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha, nbclusters = nbclusters,
              MinSun = 0.8)

print("MAP BETA DIVERSITY")
map_beta_div(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, TypePCA = TypePCA,
             window_size = window_size, nb_partitions=nb_partitions, nbCPU = nbCPU, MaxRAM = MaxRAM,
             nbclusters = nbclusters,
             MinSun = 0.8)


####################################################################################
####################################################################################
############## Validation ############## 
####################################################################################
####################################################################################


#### load data from biodivmapR pipeline ####
# infos from PCA and Kmeans imported from server 
load("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/PCA_Info.RData")
load("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/Kmeans_Info.RData")
#### load data from biodivmapR pipeline ####
# infos from PCA and Kmeans generated on local 
load("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/PCA_Info.RData")
load("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/Kmeans_Info.RData")


# selected PCs imported from server 
Sel_PC <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/Selected_Components.txt"
Selected_Features <- read.table(Sel_PC)[[1]]
# selected PCs from local
Sel_PC <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/Selected_Components.txt"
Selected_Features <- read.table(Sel_PC)[[1]]
# change selected PC
 # Selected_Features <- c(1,6)
####  IF THERE ARE MULTIPLE SHAPEFILES :  location of the directory where shapefiles used for validation are saved ####
VectorDir <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data"
# list vector data
Path_Vector <- list_shp(VectorDir)
Name_Vector <- tools::file_path_sans_ext(basename(Path_Vector))
# location of the spectral species raster needed for validation
Path_SpectralSpecies <- Kmeans_info$SpectralSpecies

#### modify PCA_Files ####
# imported from server  (server path => local path)
PCA_Files <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/OutputPCA_8_PCs"
Path_SpectralSpecies <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr/RESULTS/T58KFA_raw_int_forest_zone_PDL/T58KFA_raw_int_forest_zone_PDL/SPCA/SpectralSpecies/SpectralSpecies"
# generated on local 
PCA_Files <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/PCA/OutputPCA_8_PCs"
Path_SpectralSpecies <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/biodivmapr_local/RESULTS/T58KFA_raw_int_forest_zone_PDL/SPCA/SpectralSpecies/SpectralSpecies"

#### filter field plots that have values in PCA raster ####
# get plot spatial data
library(rgdal)
plot_data <- readOGR(Path_Vector, Name_Vector)
# get PCA raster
PCA_map <- stack(PCA_Files)
# extract values
library(exactextractr)
test_PCA <- exact_extract(PCA_map, plot_data)
# keep fields plots with at more than n pixel values from PCA raster based on n selected PCs.
nb_PCs <- length(Selected_Features)
keep_plots <- unlist(lapply(test_PCA, function(x) sum(!is.na(x[,1])) )) > nb_PCs
plot_data_ok <- plot_data[keep_plots,]

#### keep only plots with aproximatly the same size (314m2 to 400m2) ####
plot_data_ok_same <- plot_data_ok[plot_data_ok$plot %in% c("circle_R10","circle_R11.28","plot_20x20"),]

#### create polygons that perfecly match with raster dimensions #### 
#### (here 9 cells from plot centroid, but see "adjacent function for other possibilities) ####
xy_plots_centers = rgeos::gCentroid(plot_data_ok_same,byid=TRUE)
r <- raster(PCA_Files)
cells <- cellFromXY(r, xy_plots_centers)
adj <- adjacent(r, cells, 8, include=TRUE)
r[] <- 0
r[adj[,2]] <- 1
plot_data_ok_rast <- sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(r), 
                                         as_points = FALSE, merge = TRUE)) 
#### export in a new file ####
dir.create("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_ok")
# export
library(sf)
plot_data_ok <- st_as_sf(plot_data_ok_rast)
st_write(plot_data_ok_rast, dsn = '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_ok',
          layer = 'geom_plots_PDL_div_indices_utm_ok', driver = "ESRI Shapefile", append=FALSE)

#### get the updated shapefile ####
VectorDir <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_ok"
# list vector data
Path_Vector <- list_shp(VectorDir)
Name_Vector <- tools::file_path_sans_ext(basename(Path_Vector))

#### get the test shapefile ####
# VectorDir <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_test"
# # list vector data
# Path_Vector <- list_shp(VectorDir)
# Name_Vector <- tools::file_path_sans_ext(basename(Path_Vector))


#### Compute diversity indicators ####
# get diversity indicators corresponding to shapefiles (no partitioning of spectral diversity based on field plots so far...)
Biodiv_Indicators <- diversity_from_plots(Raster_SpectralSpecies = Path_SpectralSpecies, Plots = Path_Vector,
                                          Raster_Functional = PCA_Files, Selected_Features = Selected_Features, nbclusters = nbclusters)

Shannon_RS <- c(Biodiv_Indicators$Shannon)[[1]]
FRic <- c(Biodiv_Indicators$FunctionalDiversity$FRic)
FEve <- c(Biodiv_Indicators$FunctionalDiversity$FEve)
FDiv <- c(Biodiv_Indicators$FunctionalDiversity$FDiv)
# if no name for plots
Biodiv_Indicators$Name_Plot = seq(1,length(Biodiv_Indicators$Shannon[[1]]),by = 1)

#### The diversity indices corresponding to the plots can then be written in CSV files. ####

Path_Results <- file.path(Output_Dir, im_name, TypePCA, 'VALIDATION')
dir.create(Path_Results, showWarnings = FALSE,recursive = TRUE)

# save RDS object wih all indices
# saveRDS(Biodiv_Indicators, paste0(Path_Results, "Biodiv_Indicators.rds"))

# write a table for Shannon index
write.table(Shannon_RS, file = file.path(Path_Results,"ShannonIndex.csv"),
            sep="\t", dec=".", na=" ", row.names = Biodiv_Indicators$Name_Plot, col.names= F,quote=FALSE)

# write a table for all spectral diversity indices corresponding to alpha diversity
Results <- data.frame(Name_Vector, Biodiv_Indicators$Richness, Biodiv_Indicators$Fisher,
                      Biodiv_Indicators$Shannon, Biodiv_Indicators$Simpson,
                      Biodiv_Indicators$FunctionalDiversity$FRic,
                      Biodiv_Indicators$FunctionalDiversity$FEve,
                      Biodiv_Indicators$FunctionalDiversity$FDiv)

names(Results)  = c("ID_Plot", "Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
write.table(Results, file = file.path(Path_Results,"AlphaDiversity.csv"),
            sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)

# write a table for Bray Curtis dissimilarity
BC_mean <- Biodiv_Indicators$BCdiss
colnames(BC_mean) <- rownames(BC_mean) <- Biodiv_Indicators$Name_Plot
write.table(BC_mean, file = file.path(Path_Results,"BrayCurtis.csv"),
            sep="\t", dec=".", na=" ", row.names = F, col.names= T,quote=FALSE)


#############################################################################################################################################################







#### DRAFTS ####








# # convert the image using Convert.Raster2BIL if not in the proper format
Input_Image_File  <- raster2BIL(Raster_Path = paste0(im_path, "/", im_name, ".envi"),
                                Sensor = "unknown",
                                Convert_Integer = FALSE,
                                Multiplying_Factor= 10000,
                                #Multiplying_Factor_Last =1,
                                Output_Dir = FALSE)

im_bil=stack(paste0(im_path, "/biodivMapR_Convert_BIL/", im_name))

im_envi=stack(paste0(im_path, "/", im_name, ".envi"))
plot(im_envi)



Input_Image_File=(paste0(im_path, "/biodivMapR_Convert_BIL/", im_name))
Input_Image_File
#Input_Mask_File   ="D:/Mes Donnees/Google_Cloud/mask_NC_sud"
Output_Dir        = paste0(im_path, "/biodivmapr_local/RESULTS/", im_name)
Input_Mask_File   = 
  #
  #red=raster('J:/Sentinel2/Gilles/biodivMapR_Convert_BIL/Saotome_S2', band=3)
  #nir=raster('J:/Sentinel2/Gilles/biodivMapR_Convert_BIL/Saotome_S2', band=7)
  #nir=calc(nir,function(x)ifelse(x>0,x,NA))
  #range(values(nir), na.rm = T)
  #ndvi=(nir-red)/(nir+red)
  #toto=raster('J:/Sentinel2/Gilles/biodivMapR_Convert_BIL/Saotome_S2', band=7)
  #mapview(toto)
  
  
  window_size = 10
TypePCA = 'SPCA'
FilterPCA = TRUE

########Computing options
nbCPU         = 3
MaxRAM        = 0.5
nbclusters    = 80

# 1- MASK : Filter data in order to discard non vegetated / shaded / cloudy pixels
NDVI_Thresh = 0.53
Blue_Thresh = 500
NIR_Thresh  = 1500
print("PERFORstM RADIOMETRIC FILTERING")
Input_Mask_File = perform_radiometric_filtering(Input_Image_File,Input_Mask_File,Output_Dir,
                                                NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                NIR_Thresh = NIR_Thresh)


############## PCA
print("PERFORM PCA ON RASTER")
PCA_Output        = perform_PCA(Input_Image_File, Input_Mask_File, Output_Dir,
                                TypePCA = TypePCA, FilterPCA = FilterPCA, nbCPU = nbCPU,MaxRAM = MaxRAM, Continuum_Removal = TRUE)


# Path for the PCA raster
PCA_Files         = PCA_Output$PCA_Files
# number of pixels used for each partition used for k-means clustering
Pix_Per_Partition = PCA_Output$Pix_Per_Partition
# number of partitions used for k-means clustering
nb_partitions     = PCA_Output$nb_partitions

# Path for the updated mask
Input_Mask_File   = PCA_Output$MaskPath
# parameters of the PCA model
PCA_model         = PCA_Output$PCA_model
# definition of spectral bands to be excluded from the analysis
SpectralFilter    = PCA_Output$SpectralFilter


print("Select PCA components for diversity estimations")
select_PCA_components(Input_Image_File, Output_Dir, PCA_Files, File_Open = TRUE)

print("MAP SPECTRAL SPECIES")
map_spectral_species(Input_Image_File, Output_Dir, PCA_Files, PCA_model, SpectralFilter, Input_Mask_File,
                     Pix_Per_Partition, nb_partitions, nbCPU=nbCPU, MaxRAM=MaxRAM, 
                     nbclusters = nbclusters, TypePCA = TypePCA, Continuum_Removal = TRUE)
###################################
####"#alpha an beta diversity maps

print("MAP ALPHA DIVERSITY")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha   = c('Shannon')
map_alpha_div(Input_Image_File, Output_Dir, window_size, nbCPU=nbCPU,
              MaxRAM=MaxRAM, Index_Alpha = Index_Alpha, nbclusters = nbclusters)

print("MAP BETA DIVERSITY")
map_beta_div(Input_Image_File, Output_Dir, window_size, nb_partitions=nb_partitions,
             nbCPU=nbCPU, MaxRAM=MaxRAM, nbclusters = nbclusters)


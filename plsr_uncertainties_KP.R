######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
####                                    PLSR for uncertainties                                    ####
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

library(pls)
library(spatial)
library(sp)
library(sf)
library(rgdal)
library(raster)
library(rgeos)
library(RSAGA)
#### load data ####
# get fragmentation metrics and mean environmental variables for 30*30m cells with mean S2 bands and mean spectral PCA values
poly_cells_30m_ok_landscape_env_S2 <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58KDB_2020/fragmentation_data/poly_cells_30m_ok_landscape_env_S2.rds")
poly_cells_30m_ok_landscape_env_S2$log_dist_edge <- log(poly_cells_30m_ok_landscape_env_S2$dist_edge)
landscape_data_utm <- poly_cells_30m_ok_landscape_env_S2[,c("mean.B2","mean.B3","mean.B4",                              
                                                            "mean.B5","mean.B6","mean.B7",                              
                                                            "mean.B8","mean.B8A","mean.B11","mean.B12")]

landscape_data <- data.frame(landscape_data_utm)[,c("mean.B2","mean.B3","mean.B4",                              
                                                    "mean.B5","mean.B6","mean.B7",                              
                                                    "mean.B8","mean.B8A","mean.B11","mean.B12")]
colnames(landscape_data) <- paste0("S2.", colnames(landscape_data))

#### diversity indices form biological data ####
plots_field <- readRDS('/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58KDB_2020/plot_data_matrices/plots_KP_div_indices.rds')

#### plot true area around true plots ####
plots_field_ok <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58KDB_2020/plot_data/geom_plots_KP_div_indices_utm.shp")
plots_field_ok <- plots_field[plots_field$locality %in% plot_data_ok$locality,]
cbind(plots_field_ok$locality, plot_data_ok$locality)
#### create polygons that perfecly match with raster dimensions #### 
# extract S2 values for true plots
S2_image <- stack("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58KDB_2020/T58KDB_raw_int_forest_toDrive_ok.tif")

library(exactextractr)
test_S2 <- exact_extract(S2_image, plot_data_ok)

# keep fields plots intersecting more than n pixel values from S2 image in forest (not 0).
n = 5
keep_plots <- unlist(lapply(test_S2, function(x) sum(!x[,1] == 0) )) >= n
plot_data_ok <- plot_data_ok[keep_plots,]
#### (here 9 cells from plot centroid, but see "adjacent function for other possibilities) ####
xy_plots_centers = rgeos::gCentroid(plot_data_ok,byid=TRUE)
r <- raster(S2_image)
r[] <- 0
cells <- cellFromXY(r, xy_plots_centers)
plot_data_ok_rast_list <- list()
for(i in 1:length(cells)){
  adj <- adjacent(r, cells[i], 8, include=TRUE)
  r_tmp <- r
  r_tmp[adj[,2]] <- 1
  plot_data_ok_rast_tmp <- sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(r_tmp), 
                                                   as_points = FALSE, merge = TRUE)) 
  plot_data_ok_rast_tmp$locality <- plot_data_ok$locality[i]
  plot_data_ok_rast_list[[i]] <-  plot_data_ok_rast_tmp
}
plot_data_ok_rast <- bind(unlist(lapply(plot_data_ok_rast_list, function(x) x[1,])))
plot_data_ok_rast

#### export in a new file ####
dir.create("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58KDB_2020/plot_data_ok")
# export
library(sf)
plot_data_ok_rast <- st_as_sf(plot_data_ok_rast)
# st_write(plot_data_ok_rast, dsn = '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58KDB_2020/plot_data_ok',
#           layer = 'geom_plots_KP_div_indices_utm_ok', driver = "ESRI Shapefile", append=FALSE)
#### extract S2 values for plots for area that match S2 resolution ####
names(S2_image) <- c("B2",     "B3",     "B4",     "B5",     "B6",     "B7",     "B8",    "B8A",    "B11",    "B12")
plot_CWM_raw_spectral <- exact_extract(S2_image, plot_data_ok_rast, fun = "mean")
colnames(plot_CWM_raw_spectral) <- paste0("S2.",colnames(plot_CWM_raw_spectral))

#### compil with plot field data ####
cbind(plots_field_ok$locality, plot_data_ok_rast$locality)
plots_field_spectral <- cbind(plots_field, plot_CWM_raw_spectral)
# saveRDS(plots_field_spectral, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58KDB_2020/plot_data_ok/plots_field_spectral.rds")

######################################################################################################
####                                           test PLSR                                          ####
######################################################################################################
library(pls)
plots_field_spectral <- readRDS( "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58KDB_2020/plot_data_ok/plots_field_spectral.rds")

names(plots_field_spectral)

#### remove outliers taht include too much open areas in 30*30 cells ####
outliers <- c("FBA_070_KP", "FBA_071_KP", "FBA_072_KP", "FBA_103_KP")
plots_field_spectral <- plots_field_spectral[!plots_field_spectral$locality %in% outliers,]

# K-fold CV
model <- plsr(PCoA_PC1_taxo ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
pls::R2(model, estimate = "all") 

model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")


summary(model)
plot(RMSEP(model), legendpos = "topright")
pls::R2(model, estimate = "all") 

model <- plsr(CWM_trans_synth_CA_Axis1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")


summary(model)
plot(RMSEP(model), legendpos = "topright")
cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
plot(model, ncomp = best.dims, asp = 1, line = TRUE)
pls::R2(model, estimate = "all") 

pred <- drop(predict(model, ncomp = best.dims))

plots_field_spectral[,c("CWM_SLA", "locality")] [order(plots_field_spectral$CWM_SLA),]

a <- data.frame(cbind(pred, plots_field_spectral$locality))
a$pred <- as.numeric(a$pred)
a[order(a$pred),]

######################################################################################################
####                                         PLSR for SLA                                         ####
######################################################################################################

#### select best model with all data with a cross validation procedure ####

model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="none")
# leave one out 
model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="LOO")
# K-fold CV
model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
plot(model, ncomp = best.dims, asp = 1, line = TRUE)
pls::R2(model, estimate = "all") 

#### get variable importance as the standardized coefficent of predictors ####
library(plsVarSel)
# The regression coefficient represents a predictor variables’ importance in predicting the response
rc <- rev(sort(RC(model, best.dims)))
barplot(rc)
rc_SLA <-  rc/sum(rc)*100
# using caret?
var_imp <- varImp(model)
var_imp <- var_imp - min(var_imp)
var_imp <- var_imp / max(var_imp)*100
var_imp <- var_imp[rev(order(var_imp$Overall)),, drop=FALSE]
var_imp
plot(var_imp)
varimp_SLA <- var_imp
#### make spatial predictions at the landscape scale ####
landscape_pred <- drop(predict(model, landscape_data, ncomp = best.dims))
#### get spatial object from predictions ####
library(dplyr)
landscape_pred_utm <- landscape_data_utm
landscape_pred_CWM_SLA_utm <- cbind(landscape_pred_utm, pred_CWM_SLA = landscape_pred)

#### make models for generating uncertainties ####

# calibration and validation datasets to compute uncertainties
cal_prop <- 0.75

# for loop for iterations
iter = 100
valid_smpl <- c()
valid_data <- c()
preds <- c()
rmsep <- c()
models <- list()
for (i in 1:iter){
  calib_smpl <- sample(nrow(plots_field_spectral),nrow(plots_field_spectral) * cal_prop)
  hist(calib_smpl)
  valid_smpl_tmp <- (1:nrow(plots_field_spectral))[!1:nrow(plots_field_spectral) %in% calib_smpl]
  calib_data_tmp <- plots_field_spectral[calib_smpl,]
  valid_data_tmp <- plots_field_spectral[valid_smpl_tmp,]
  model_tmp <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                  S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                  S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
  
  
  preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
  rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
  hist(valid_data_tmp$CWM_SLA)
  plot(valid_data_tmp$CWM_SLA, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
  abline(0, 1, col="red")
  valid_smpl <- cbind(valid_smpl, valid_smpl_tmp)
  valid_data <- cbind(valid_data, valid_data_tmp$CWM_SLA)
  preds <- cbind(preds, preds_tmp)
  rmsep <- c(rmsep, data.frame(rmsep_tmp$val))
  models[[i]] <- model_tmp
}

# get all predictions in one table
preds_df <- data.frame(sample_data = unlist(data.frame(valid_smpl)),
                       observed = unlist(data.frame(valid_data)), 
                       predictions = unlist(data.frame(preds)))
plot(preds_df$predictions ~ preds_df$observed,main="Test Dataset")

# get statistics from all trhe predictions from iterations
models_stats_SLA <- caret::postResample(preds_df$predictions, preds_df$observed)
models_stats_SLA <- c(models_stats_SLA ,models_stats_SLA[1]/mean(plots_field_spectral$CWM_SLA)*100)
names(models_stats_SLA)[4] <- "%RMSE"
# get mean predictions and standard deviation in one table
mean_sd_preds_df <- data.frame(aggregate(predictions~sample_data, preds_df, function(x) c(mean = mean(x), sd = sd(x))))
mean_sd_preds_df <- cbind(mean_sd_preds_df[-ncol(mean_sd_preds_df)], mean_sd_preds_df[[ncol(mean_sd_preds_df)]])
mean_sd_preds_df$observed <- plots_field_spectral$CWM_SLA
mean_sd_preds_df$plow <- mean_sd_preds_df$mean - mean_sd_preds_df$sd
mean_sd_preds_df$pup <- mean_sd_preds_df$mean + mean_sd_preds_df$sd

# plot
plot_CV_SLA <- ggplot(data=mean_sd_preds_df)+aes(x=observed) +
    geom_point(aes(y=mean),color='blue',shape=19) +
    geom_errorbar(aes(ymin=plow,ymax=pup)) +
      geom_abline(slope = 1, intercept = 0, colour = "red") +
      xlab("Observed") +
      ylab("Predicted") +
      theme_bw()

#### make predictions from each models from iterations ####
pb = utils::txtProgressBar(min = 0, max = length(models), initial = 0) 
landsacpe_pred <- c()
for(i in 1:length(models)){
  utils::setTxtProgressBar(pb,i)
  model_tmp <- models[[i]]
  landsacpe_pred_tmp <- drop(predict(model_tmp, landscape_data, ncomp = best.dims))
  landsacpe_pred <- cbind(landsacpe_pred,landsacpe_pred_tmp)
}
close(pb)

landscape_pred_sd <- apply(landsacpe_pred, 1, sd)

#### get spatial object from predictions uncertainties ####
library(dplyr)
landscape_pred_CWM_SLA_utm <- cbind(landscape_pred_CWM_SLA_utm, pred_CWM_SLA_sd = landscape_pred_sd)
#### export spatial predictions with uncertainties ####
st_write(landscape_pred_CWM_SLA_utm, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_CWM_SLA_utm.shp")

######################################################################################################
####                                         PLSR for WD                                         ####
######################################################################################################

#### select best model with all data with a cross validation procedure ####

model <- plsr(CWM_WD ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="none")
# leave one out 
model <- plsr(CWM_WD ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="LOO")
# K-fold CV
model <- plsr(CWM_WD ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
plot(model, ncomp = best.dims, asp = 1, line = TRUE)
pls::R2(model, estimate = "all") 

#### get variable importance as the standardized coefficent of predictors ####
library(plsVarSel)
# The regression coefficient represents a predictor variables’ importance in predicting the response
rc <- rev(sort(RC(model, best.dims)))
barplot(rc)
rc_WD <- rc/sum(rc)*100
# using caret?
library(caret)
var_imp <- varImp(model)
var_imp <- var_imp - min(var_imp)
var_imp <- var_imp / max(var_imp)*100
var_imp <- var_imp[rev(order(var_imp$Overall)),, drop=FALSE]
var_imp
plot(var_imp)
varimp_WD <- var_imp

#### make spatial predictions at the landscape scale ####
landsacpe_pred <- drop(predict(model, landscape_data, ncomp = best.dims))
#### get spatial object from predictions ####
library(dplyr)
landscape_pred_utm <- landscape_data_utm
landscape_pred_CWM_WD_utm <- cbind(landscape_pred_utm, pred_CWM_WD = landsacpe_pred)

#### make models for generating uncertainties ####

# calibration and validation datasets to compute uncertainties
cal_prop <- 0.75

# for loop for iterations
iter = 100
valid_smpl <- c()
valid_data <- c()
preds <- c()
rmsep <- c()
models <- list()
for (i in 1:iter){
  calib_smpl <- sample(nrow(plots_field_spectral),nrow(plots_field_spectral) * cal_prop)
  hist(calib_smpl)
  valid_smpl_tmp <- (1:nrow(plots_field_spectral))[!1:nrow(plots_field_spectral) %in% calib_smpl]
  calib_data_tmp <- plots_field_spectral[calib_smpl,]
  valid_data_tmp <- plots_field_spectral[valid_smpl_tmp,]
  model_tmp <- plsr(CWM_WD ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                      S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                      S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                    data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
  
  
  preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
  rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
  hist(valid_data_tmp$CWM_WD)
  plot(valid_data_tmp$CWM_WD, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
  abline(0, 1, col="red")
  valid_smpl <- cbind(valid_smpl, valid_smpl_tmp)
  valid_data <- cbind(valid_data, valid_data_tmp$CWM_WD)
  preds <- cbind(preds, preds_tmp)
  rmsep <- c(rmsep, data.frame(rmsep_tmp$val))
  models[[i]] <- model_tmp
}

# get all predictions in one table
preds_df <- data.frame(sample_data = unlist(data.frame(valid_smpl)),
                       observed = unlist(data.frame(valid_data)), 
                       predictions = unlist(data.frame(preds)))
plot(preds_df$predictions ~ preds_df$observed,main="Test Dataset")

# get statistics from all trhe predictions from iterations
models_stats_WD <- caret::postResample(preds_df$predictions, preds_df$observed)
models_stats_WD <- c(models_stats_WD ,models_stats_WD[1]/mean(plots_field_spectral$CWM_WD)*100)
names(models_stats_WD)[4] <- "%RMSE"
# get mean predictions and standard deviation in one table
mean_sd_preds_df <- data.frame(aggregate(predictions~sample_data, preds_df, function(x) c(mean = mean(x), sd = sd(x))))
mean_sd_preds_df <- cbind(mean_sd_preds_df[-ncol(mean_sd_preds_df)], mean_sd_preds_df[[ncol(mean_sd_preds_df)]])
mean_sd_preds_df$observed <- plots_field_spectral$CWM_WD
mean_sd_preds_df$plow <- mean_sd_preds_df$mean - mean_sd_preds_df$sd
mean_sd_preds_df$pup <- mean_sd_preds_df$mean + mean_sd_preds_df$sd

# plot
plot_CV_WD <- ggplot(data=mean_sd_preds_df)+aes(x=observed) +
  geom_point(aes(y=mean),color='blue',shape=19) +
  geom_errorbar(aes(ymin=plow,ymax=pup)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  xlab("Observed") +
  ylab("Predicted") +
  theme_bw()

#### make predictions from each models from iterations ####
pb = utils::txtProgressBar(min = 0, max = length(models), initial = 0) 
landsacpe_pred <- c()
for(i in 1:length(models)){
  utils::setTxtProgressBar(pb,i)
  model_tmp <- models[[i]]
  landsacpe_pred_tmp <- drop(predict(model_tmp, landscape_data, ncomp = best.dims))
  landsacpe_pred <- cbind(landsacpe_pred,landsacpe_pred_tmp)
}
close(pb)

landscape_pred_sd <- apply(landsacpe_pred, 1, sd)

#### get spatial object from predictions uncertainties ####
library(dplyr)
landscape_pred_CWM_WD_utm <- cbind(landscape_pred_CWM_WD_utm, pred_CWM_WD_sd = landscape_pred_sd)
#### export spatial predictions with uncertainties ####
st_write(landscape_pred_CWM_WD_utm, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_CWM_WD_utm.shp")



######################################################################################################
####                                         PLSR for PCoA 1 taxo                                 ####
######################################################################################################
#### see script "updated_pcoa_taxo.R"
#### select best model with all data with a cross validation procedure ####

model <- plsr(PCoA_PC1_taxo ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="none")
# leave one out 
model <- plsr(PCoA_PC1_taxo ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="LOO")
# K-fold CV
model <- plsr(PCoA_PC1_taxo ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
plot(model, ncomp = best.dims, asp = 1, line = TRUE)
pls::R2(model, estimate = "all") 

#### get variable importance as the standardized coefficent of predictors ####
library(plsVarSel)
# The regression coefficient represents a predictor variables’ importance in predicting the response
rc <- rev(sort(RC(model, best.dims)))
barplot(rc)
rc_PCoA_1 <- rc/sum(rc)*100
# using caret?
var_imp <- varImp(model)
var_imp <- var_imp - min(var_imp)
var_imp <- var_imp / max(var_imp)*100
var_imp <- var_imp[rev(order(var_imp$Overall)),, drop=FALSE]
var_imp
plot(var_imp)
varimp_PCoA_1 <- var_imp

#### make spatial predictions at the landscape scale ####
landsacpe_pred <- drop(predict(model, landscape_data, ncomp = best.dims))
#### get spatial object from predictions ####
library(dplyr)
landscape_pred_utm <- landscape_data_utm
landscape_pred_PCoA1_utm <- cbind(landscape_pred_utm, pred_PCoA1 = landsacpe_pred)

#### make models for generating uncertainties ####

# calibration and validation datasets to compute uncertainties
cal_prop <- 0.75

# for loop for iterations
iter = 100
valid_smpl <- c()
valid_data <- c()
preds <- c()
rmsep <- c()
models <- list()
for (i in 1:iter){
  calib_smpl <- sample(nrow(plots_field_spectral),nrow(plots_field_spectral) * cal_prop)
  hist(calib_smpl)
  valid_smpl_tmp <- (1:nrow(plots_field_spectral))[!1:nrow(plots_field_spectral) %in% calib_smpl]
  calib_data_tmp <- plots_field_spectral[calib_smpl,]
  valid_data_tmp <- plots_field_spectral[valid_smpl_tmp,]
  model_tmp <- plsr(PCoA_PC1_taxo ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                      S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                      S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                    data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
  
  
  preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
  rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
  hist(valid_data_tmp$PCoA_PC1_taxo)
  plot(valid_data_tmp$PCoA_PC1_taxo, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
  abline(0, 1, col="red")
  valid_smpl <- cbind(valid_smpl, valid_smpl_tmp)
  valid_data <- cbind(valid_data, valid_data_tmp$PCoA_PC1_taxo)
  preds <- cbind(preds, preds_tmp)
  rmsep <- c(rmsep, data.frame(rmsep_tmp$val))
  models[[i]] <- model_tmp
}

# get all predictions in one table
preds_df <- data.frame(sample_data = unlist(data.frame(valid_smpl)),
                       observed = unlist(data.frame(valid_data)), 
                       predictions = unlist(data.frame(preds)))
plot(preds_df$predictions ~ preds_df$observed,main="Test Dataset")

# get statistics from all trhe predictions from iterations
models_stats_PCoA1 <- caret::postResample(preds_df$predictions, preds_df$observed)
models_stats_PCoA1 <- c(models_stats_PCoA1 ,models_stats_PCoA1[1]/mean(plots_field_spectral$PCoA_PC1_taxo)*100)
names(models_stats_PCoA1)[4] <- "%RMSE"
# get mean predictions and standard deviation in one table
mean_sd_preds_df <- data.frame(aggregate(predictions~sample_data, preds_df, function(x) c(mean = mean(x), sd = sd(x))))
mean_sd_preds_df <- cbind(mean_sd_preds_df[-ncol(mean_sd_preds_df)], mean_sd_preds_df[[ncol(mean_sd_preds_df)]])
mean_sd_preds_df$observed <- plots_field_spectral$PCoA_PC1_taxo
mean_sd_preds_df$plow <- mean_sd_preds_df$mean - mean_sd_preds_df$sd
mean_sd_preds_df$pup <- mean_sd_preds_df$mean + mean_sd_preds_df$sd

# plot
plot_CV_PCoA1 <- ggplot(data=mean_sd_preds_df)+aes(x=observed) +
  geom_point(aes(y=mean),color='blue',shape=19) +
  geom_errorbar(aes(ymin=plow,ymax=pup)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  xlab("Observed") +
  ylab("Predicted") +
  theme_bw()
plot_CV_PCoA1
#### make predictions from each models from iterations ####
pb = utils::txtProgressBar(min = 0, max = length(models), initial = 0) 
landsacpe_pred <- c()
for(i in 1:length(models)){
  utils::setTxtProgressBar(pb,i)
  model_tmp <- models[[i]]
  landsacpe_pred_tmp <- drop(predict(model_tmp, landscape_data, ncomp = best.dims))
  landsacpe_pred <- cbind(landsacpe_pred,landsacpe_pred_tmp)
}
close(pb)

landscape_pred_sd <- apply(landsacpe_pred, 1, sd)

#### get spatial object from predictions uncertainties ####
library(dplyr)
landscape_pred_PCoA1_utm <- cbind(landscape_pred_PCoA1_utm, pred_PCoA1_sd = landscape_pred_sd)
#### export spatial predictions with uncertainties ####
st_write(landscape_pred_PCoA1_utm, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_PCoA1_utm.shp")


######################################################################################################
#### Analyze predictions (extrapolated biological attributes from S2) ~ fragmentation indices ####
######################################################################################################
#### load spatial prediction ####
landscape_pred_CWM_SLA_utm <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_CWM_SLA_utm.shp")
landscape_pred_CWM_WD_utm <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_CWM_WD_utm.shp")
landscape_pred_PCoA1_utm <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_PCoA1_utm.shp")

#### add predictions to spatial data ####
poly_cells_30m_ok_landscape_env_S2$preds_CWM_SLA <- landscape_pred_CWM_SLA_utm$pr_CWM_SLA
poly_cells_30m_ok_landscape_env_S2$preds_CWM_SLA_SD <- landscape_pred_CWM_SLA_utm$p_CWM_SLA_
poly_cells_30m_ok_landscape_env_S2$preds_CWM_WD <- landscape_pred_CWM_WD_utm$pr_CWM_WD
poly_cells_30m_ok_landscape_env_S2$preds_CWM_WD_SD <- landscape_pred_CWM_WD_utm$p_CWM_WD_
poly_cells_30m_ok_landscape_env_S2$pred_PCoA1 <- landscape_pred_PCoA1_utm$pr_PCA1
poly_cells_30m_ok_landscape_env_S2$pred_PCoA1_SD <- landscape_pred_PCoA1_utm$p_PCA1_

#### covariation between observed variables and between predictions ####
# observations
cor_obs <- cor( plots_field_spectral[,c("CWM_SLA","CWM_WD", "PCoA_taxo_ok_1")],
                method = "pearson", use = "complete.obs")
cor_obs
# predictions
cor_pred <- cor( data.frame(poly_cells_30m_ok_landscape_env_S2[, c("preds_CWM_SLA", "preds_CWM_WD", "pred_PCoA1")])[,-4],
                 method = "pearson", use = "complete.obs")
cor_pred

#### transform and add aspect in direction #### 
# get aspect in radians
poly_cells_30m_ok_landscape_env_S2$aspect_rad <- poly_cells_30m_ok_landscape_env_S2$aspect/180*pi
# get northness and eastness
poly_cells_30m_ok_landscape_env_S2$northness = cos(poly_cells_30m_ok_landscape_env_S2$aspect_rad)
poly_cells_30m_ok_landscape_env_S2$eastness = sin(poly_cells_30m_ok_landscape_env_S2$aspect_rad)

hist(poly_cells_30m_ok_landscape_env_S2$northness)
hist(poly_cells_30m_ok_landscape_env_S2$eastness)

#### compute perimeter / area ratio for local landscapes ####
poly_cells_30m_ok_landscape_env_S2_df$perimeter.area.ratio_500_centroid <- 
  poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid/poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_500_centroid

#### get dataframe for multivariate model selection  ####
poly_cells_30m_ok_landscape_env_S2_df <- data.frame(poly_cells_30m_ok_landscape_env_S2)
colnames(poly_cells_30m_ok_landscape_env_S2_df)


#### remove outliers from SLA prediction (visual chack shows that negative values are on non forest zones) ####
poly_cells_30m_ok_landscape_env_S2_df <- poly_cells_30m_ok_landscape_env_S2_df[poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA>0,]
#### remove cells centroids a less than 15m from the edge ####
poly_cells_30m_ok_landscape_env_S2_df <- poly_cells_30m_ok_landscape_env_S2_df[poly_cells_30m_ok_landscape_env_S2_df$dist_edge>15,]

#### test relationships ####

plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$log_dist_edge)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$log_dist_edge))

plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge))

plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$elevation)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$elevation))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid))
plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid)
#
plot(poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge)
#
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$total.edge_500_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid))

summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$total.edge_1000_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$edge.density_1000_centroid))


plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_1000_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_1000_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$edge.density_100_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$edge.density_500_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$total.edge_500_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$effective.mesh.size_500_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$effective.mesh.size_1000_centroid)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$effective.mesh.size_500_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$effective.mesh.size_1000_centroid))


plot(poly_cells_30m_ok_landscape_env_S2_df$total.edge_1000_centroid ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge)
plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$total.edge_1000_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_1000_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df$total.edge_1000_centroid ~ poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_1000_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df$total.edge_500_centroid ~ poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_500_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$total.edge_500_centroid)


plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_100_centroid)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_100_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$prop.landscape_100_centroid))


plot(poly_cells_30m_ok_landscape_env_S2_df$perimeter.area.ratio_500_centroid ~ poly_cells_30m_ok_landscape_env_S2_df$dist_edge)
plot(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ log(poly_cells_30m_ok_landscape_env_S2_df$perimeter.area.ratio_500_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df$perimeter.area.ratio_500_centroid))


#### select variables ####
# choose set of explanatory variable

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

# keep this
expl_v_list_list <- list(c( "log_dist_edge", "total.edge_500_centroid", "effective.mesh.size_500_centroid",
                            "elevation","slope",                        
                            "curvature", "twi" ))

expl_v_list_list <- list(c( "log_dist_edge", "prop.landscape_500_centroid",
                            "elevation","slope",                        
                            "curvature", 
                            "northness", "eastness"))


# set of response variables
resp_v_list <- c( "preds_CWM_SLA",  "preds_CWM_SLA_SD",
                  "preds_CWM_WD", "preds_CWM_WD_SD",
                  "pred_PCoA1", "pred_PCoA1_SD")


#### model selection ####

library(MuMIn)
library(parallel)
# total set of variable
data_for_mod_all <- poly_cells_30m_ok_landscape_env_S2_df

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
# write.csv2(list_coef_best, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/rapport/results_models_bio_from_spectral_frag_landscape.csv")



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
levels(plot_coef_list[[1]]$data$term)[levels(plot_coef_list[[1]]$data$term) == "prop.landscape_100_centroid"] <- "Habitat amount"


plot_coef_list[[1]]$labels$title <- "Influence of edge and topography on CWM SLA"
plot_coef_list[[3]]$labels$title <- "Influence of edge and topography on CWM WD"

plot_coef_list
## export png ##
# png("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/plot_coef_pred_SLA.png", width=150, height=70, units = 'mm', res = 300) 
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

data_for_varpart <- poly_cells_30m_ok_landscape_env_S2_df[complete.cases(poly_cells_30m_ok_landscape_env_S2_df[,c("preds_CWM_SLA",
                                                                                                                  "preds_CWM_WD",
                                                                                                                  "pred_PCoA1",
                                                                                                                  "elevation","slope",                        
                                                                                                                  "curvature", "twi",
                                                                                                                  "northness", "eastness",
                                                                                                                  "log_dist_edge")]),]

colnames(data_for_varpart)

Topography <- data.frame(data_for_varpart[,c("elevation","slope",                        
                                             "curvature", "twi",
                                             "northness", "eastness")])

Edge_influence <- data.frame(data_for_varpart[,c("log_dist_edge", "prop.landscape_100_centroid")])


#### preds_sp_rich ####
# relative importance
relaimpo::calc.relimp(data.frame(data_for_varpart[c("preds_CWM_SLA",
                                                    "elevation","slope",                        
                                                    "curvature", "twi",
                                                    "northness", "eastness",
                                                    "log_dist_edge", "prop.landscape_100_centroid")]))

relaimpo::calc.relimp(data.frame(data_for_varpart[c("preds_CWM_WD",
                                                    "elevation","slope",                        
                                                    "curvature", "twi",
                                                    "northness", "eastness",
                                                    "log_dist_edge", "prop.landscape_100_centroid")]))

relaimpo::calc.relimp(data.frame(data_for_varpart[c("pred_PCoA1",
                                                    "elevation","slope",                        
                                                    "curvature", "twi",
                                                    "northness", "eastness",
                                                    "log_dist_edge", "prop.landscape_100_centroid")]))
# variance partitioning with multiple responses 
library("rdacca.hp")
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13800
pred <- data.frame(pred <- data_for_varpart$preds_CWM_SLA)
pred <- data.frame(pred <- data_for_varpart$preds_CWM_WD)
pred <- data.frame(pred <- data_for_varpart$pred_PCoA1)
pred <- data.frame(data_for_varpart[,c("preds_CWM_SLA","preds_CWM_WD", "pred_PCoA1")])
hierpar_pred <- rdacca.hp(pred, list(Edge_influence, Topography), var.part = T)
hierpar_pred_all <- rdacca.hp(pred, cbind(Edge_influence, Topography ), var.part = T)
plot(hierpar_pred_all)
#### plot varpart ####
col_frag <- grDevices::adjustcolor( "forestgreen", alpha.f = 0.2)
col_topo <- grDevices::adjustcolor( "orange", alpha.f = 0.2)
# get custom venn diagram function
source("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/kuebini/showvarparts_custom.R")
Xnames <- c("Edge\ninfluence", "Topography")

hierpar_pred$Var.part <- round(hierpar_pred$Var.part, digits = 2)
labels <- c(paste0("R²=", hierpar_pred$Var.part[1,1], "         "),
            # paste0("\nR²=", ifelse(hierpar_pred$Var.part[3,1]>0, hierpar_pred$Var.part[3,1], 0)),
            paste0("\nR²=", hierpar_pred$Var.part[3,1]),
            paste0("          R²=", hierpar_pred$Var.part[2,1]),
            1-hierpar_pred$Var.part[4,1])

showvarparts_custom(2, labels = labels, Xnames = Xnames, bg = c(col_frag, col_topo), lty = 0)
######################################################################################################
####                                         other approaches                                     ####
######################################################################################################

#### manual LOO validation ####
X <- as.matrix(plots_field_spectral[,c ("S2.mean.B2","S2.mean.B3","S2.mean.B4",
                             "S2.mean.B5","S2.mean.B6","S2.mean.B7","S2.mean.B8",
                             "S2.mean.B8A","S2.mean.B11","S2.mean.B12")])
Y <- plots_field_spectral$CWM_SLA
ncomp <- 5
cvPreds <- matrix(nrow = nrow(X), ncol = ncomp)
for (i in 1:nrow(X)) {
  fit <- simpls.fit(X[-i,], Y[-i], ncomp = ncomp, stripped = TRUE)
  cvPreds[i,] <- (X[i,] - fit$Xmeans) %*% drop(fit$coefficients) +
    fit$Ymeans
}
sqrt(colMeans((cvPreds - Y)^2))

#### model tuning with caret? ####
library(caret)
# test number of components
tuneGrid <- expand.grid(
  ncomp   = seq(1, 10, by = 1)
)
# evaluate model with k-fold corss validation
ctrl <- trainControl(
  method = "cv",
  number = 10,
)
model_caret <- train(
  CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
    S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
    S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
  data = plots_field_spectral,
  method = 'pls',
  preProcess = c("center", "scale"),
  tuneGrid = tuneGrid,
  trControl = ctrl
)
model_caret
plot(model_caret)

plot(varImp(model_caret))
varImp(model_caret)


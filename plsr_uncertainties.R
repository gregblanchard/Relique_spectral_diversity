######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
####                                    PLSR for uncertainties                                    ####
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
#### see https://doi.org/10.1016/j.rse.2018.11.016 for the methods #### 
library(pls)
library(spatial)
library(sp)
library(sf)
library(rgdal)
library(raster)
library(rgeos)
library(RSAGA)
#### load data ####
plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_100sp_formPC1678.rds")
# get fragmentation metrics and mean environmental variables for 30*30m cells with mean S2 bands and mean spectral PCA values
poly_cells_30m_ok_landscape_env_S2 <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/poly_cells_30m_ok_landscape_env_S2.rds")
poly_cells_30m_ok_landscape_env_S2$log_dist_edge <- log(poly_cells_30m_ok_landscape_env_S2$dist_edge)
landscape_data_utm <- poly_cells_30m_ok_landscape_env_S2[,c("mean.B2","mean.B3","mean.B4",                              
                                                            "mean.B5","mean.B6","mean.B7",                              
                                                            "mean.B8","mean.B8A","mean.B11","mean.B12")]

landscape_data <- data.frame(landscape_data_utm)[,c("mean.B2","mean.B3","mean.B4",                              
                                                    "mean.B5","mean.B6","mean.B7",                              
                                                    "mean.B8","mean.B8A","mean.B11","mean.B12")]
colnames(landscape_data) <- paste0("S2.", colnames(landscape_data))
######################################################################################################
####                                           remove plots                                          ####
######################################################################################################

#### this plot is too much outside forest (~open area = 1/3 of the plot) ####
plots_field_spectral_ok <- plots_field_spectral[!plots_field_spectral$locality == "10-10-BIS",]
# save 
# saveRDS(plots_field_spectral_ok, '/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_field_spectral_ok.rds')
plots_field_spectral_ok <- readRDS('/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_field_spectral_ok.rds')

######################################################################################################
####                                          Add taxo data                                        ####
######################################################################################################
#### add PCOA taxo update ####
PCoA_PCs_taxo_CA_update <- readRDS( "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data/PCoA_PCs_taxo_CA_update.rds")

plots_field_spectral_ok$PCoA_taxo_CA_ok_1 <- PCoA_PCs_taxo_CA_update[,1]
plots_field_spectral_ok$PCoA_taxo_CA_ok_2 <- PCoA_PCs_taxo_CA_update[,2]
plots_field_spectral_ok$PCoA_taxo_CA_ok_3 <- PCoA_PCs_taxo_CA_update[,3]

#### add PCOA taxo update ####
nmds <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data/nmds_taxo_CA_update.rds")

cbind(rownames(nmds$points), plots_field_spectral_ok$locality)

plots_field_spectral_ok$NMDS_1 <- nmds$points[,1]
plots_field_spectral_ok$NMDS_2 <- nmds$points[,2]
plots_field_spectral_ok$NMDS_3 <- nmds$points[,3]

plots_field_spectral_ok_env$NMDS_1 <- nmds$points[,1]
plots_field_spectral_ok_env$NMDS_2 <- nmds$points[,2]
plots_field_spectral_ok_env$NMDS_3 <- nmds$points[,3]

######################################################################################################
####                                           test PLSR                                          ####
######################################################################################################
names(plots_field_spectral_ok)
# K-fold CV
model <- plsr(PCoA_taxo_ok_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
pls::R2(model, estimate = "all") 

model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")


summary(model)
plot(RMSEP(model), legendpos = "topright")
pls::R2(model, estimate = "all") 

model <- plsr(CWM_trans_synth_CA_Axis1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="CV",
              segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
plot(model, ncomp = best.dims, asp = 1, line = TRUE)
pls::R2(model, estimate = "all") 


######################################################################################################
####                                         PLSR for SLA                                         ####
######################################################################################################

#### select best model with all data with a cross validation procedure ####

model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="none")
# leave one out 
model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="LOO")
# K-fold CV
model <- plsr(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="CV",
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
  calib_smpl <- sample(nrow(plots_field_spectral_ok),nrow(plots_field_spectral_ok) * cal_prop)
  hist(calib_smpl)
  valid_smpl_tmp <- (1:nrow(plots_field_spectral_ok))[!1:nrow(plots_field_spectral_ok) %in% calib_smpl]
  calib_data_tmp <- plots_field_spectral_ok[calib_smpl,]
  valid_data_tmp <- plots_field_spectral_ok[valid_smpl_tmp,]
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
models_stats_SLA <- c(models_stats_SLA ,models_stats_SLA[1]/(max(plots_field_spectral_ok$CWM_SLA) - min(plots_field_spectral_ok$CWM_SLA))*100)
names(models_stats_SLA)[4] <- "%RMSE"
# get mean predictions and standard deviation in one table
mean_sd_preds_df <- data.frame(aggregate(predictions~sample_data, preds_df, function(x) c(mean = mean(x), sd = sd(x))))
mean_sd_preds_df <- cbind(mean_sd_preds_df[-ncol(mean_sd_preds_df)], mean_sd_preds_df[[ncol(mean_sd_preds_df)]])
mean_sd_preds_df$observed <- plots_field_spectral_ok$CWM_SLA
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
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="none")
# leave one out 
model <- plsr(CWM_WD ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="LOO")
# K-fold CV
model <- plsr(CWM_WD ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="CV",
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
  calib_smpl <- sample(nrow(plots_field_spectral_ok),nrow(plots_field_spectral_ok) * cal_prop)
  hist(calib_smpl)
  valid_smpl_tmp <- (1:nrow(plots_field_spectral_ok))[!1:nrow(plots_field_spectral_ok) %in% calib_smpl]
  calib_data_tmp <- plots_field_spectral_ok[calib_smpl,]
  valid_data_tmp <- plots_field_spectral_ok[valid_smpl_tmp,]
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
models_stats_WD <- c(models_stats_WD ,models_stats_WD[1]/(max(plots_field_spectral_ok$CWM_WD) - min(plots_field_spectral_ok$CWM_WD))*100)
names(models_stats_WD)[4] <- "%RMSE"
# get mean predictions and standard deviation in one table
mean_sd_preds_df <- data.frame(aggregate(predictions~sample_data, preds_df, function(x) c(mean = mean(x), sd = sd(x))))
mean_sd_preds_df <- cbind(mean_sd_preds_df[-ncol(mean_sd_preds_df)], mean_sd_preds_df[[ncol(mean_sd_preds_df)]])
mean_sd_preds_df$observed <- plots_field_spectral_ok$CWM_WD
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

model <- plsr(PCoA_taxo_CA_ok_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="none")
# leave one out 
model <- plsr(PCoA_taxo_CA_ok_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="LOO")
# # K-fold CV
# model <- plsr(PCoA_taxo_CA_ok_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
#                 S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
#                 S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
#               data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="CV",
#               segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
plot(model, ncomp = best.dims, asp = 1, line = TRUE)
pls::R2(model, estimate = "all") 
CV_loo_nobuffer <- pls::R2(model, estimate = "CV") 

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
  calib_smpl <- sample(nrow(plots_field_spectral_ok),nrow(plots_field_spectral_ok) * cal_prop)
  hist(calib_smpl)
  valid_smpl_tmp <- (1:nrow(plots_field_spectral_ok))[!1:nrow(plots_field_spectral_ok) %in% calib_smpl]
  calib_data_tmp <- plots_field_spectral_ok[calib_smpl,]
  valid_data_tmp <- plots_field_spectral_ok[valid_smpl_tmp,]
  model_tmp <- plsr(PCoA_taxo_CA_ok_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                      S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                      S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                    data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
  
  
  preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
  rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
  hist(valid_data_tmp$PCoA_taxo_CA_ok_1)
  plot(valid_data_tmp$PCoA_taxo_CA_ok_1, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
  abline(0, 1, col="red")
  valid_smpl <- cbind(valid_smpl, valid_smpl_tmp)
  valid_data <- cbind(valid_data, valid_data_tmp$PCoA_taxo_CA_ok_1)
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
models_stats_PCoA1 <- c(models_stats_PCoA1 ,models_stats_PCoA1[1]/(max(plots_field_spectral_ok$PCoA_taxo_CA_ok_1) - min(plots_field_spectral_ok$PCoA_taxo_CA_ok_1))*100)
names(models_stats_PCoA1)[4] <- "%RMSE"
# get mean predictions and standard deviation in one table
mean_sd_preds_df <- data.frame(aggregate(predictions~sample_data, preds_df, function(x) c(mean = mean(x), sd = sd(x))))
mean_sd_preds_df <- cbind(mean_sd_preds_df[-ncol(mean_sd_preds_df)], mean_sd_preds_df[[ncol(mean_sd_preds_df)]])
mean_sd_preds_df$observed <- plots_field_spectral_ok$PCoA_taxo_CA_ok_1
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


#### make models for cross validation with spatial LOO with buffer arount each LOO plots ####
# see Ploton et al. 2020 Nat. Com. 
# get spectial object
plots_field_spectral_ok$geo_pt..89
# calibration and validation datasets to compute uncertainties
# for loop for iterations
buffers <- seq(0,4000,500)
valid_smpl_loo_list <- list()
valid_data_loo_list <- list()
preds_loo_list <- list()
rmsep_loo_list <- list()
models_loo_list <- list()
for(j in 1:length(buffers)){
  print(buffers[j])
  buf_tmp <- buffers[j]
  n_pts = nrow(plots_field_spectral_ok)
  valid_smpl_loo<- c()
  valid_data_loo <- c()
  preds_loo <- c()
  rmsep_loo <- c()
  models_loo <- list()
  pb <- txtProgressBar(max = n_pts, style = 3)
  for (i in 1:n_pts){
    setTxtProgressBar(pb, i)
    valid_smpl_tmp  <- i
    # hist(calib_smpl)
    valid_data_tmp <- plots_field_spectral_ok[i,]
    # remove plot within a buffer around validation plot
    valid_plot <- plots_field_spectral_ok$geo_pt..89[i]
    buffer_tmp <- st_buffer(valid_plot, buf_tmp)
    calib_smpl_spatial <- (1:nrow(plots_field_spectral_ok))[lengths(st_intersects(plots_field_spectral_ok$geo_pt..89, buffer_tmp)) == 0]
    if(buf_tmp == 0)     calib_smpl_spatial <- (1:nrow(plots_field_spectral_ok))[-i]
    calib_data_tmp <- plots_field_spectral_ok[calib_smpl_spatial,]
    plot(plots_field_spectral_ok$geo_pt..89)
    plot(buffer_tmp, add = T)
    plot(plots_field_spectral_ok$geo_pt..89[calib_smpl_spatial], lwd = 5, add = T)
    # model
    model_tmp <- plsr(PCoA_taxo_CA_ok_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                        S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                        S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                      data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
    
    preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
    rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
    # hist(valid_data_tmp$PCoA_taxo_CA_ok_1)
    # plot(valid_data_tmp$PCoA_taxo_CA_ok_1, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
    # abline(0, 1, col="red")
    valid_smpl_loo <- cbind(valid_smpl_loo, valid_smpl_tmp)
    valid_data_loo <- cbind(valid_data_loo, valid_data_tmp$PCoA_taxo_CA_ok_1)
    preds_loo <- cbind(preds_loo, preds_tmp)
    rmsep_loo <- c(rmsep_loo, data.frame(rmsep_tmp$val))
    models_loo[[i]] <- model_tmp
  }
  close(pb)
  valid_smpl_loo_list[[j]] <- valid_smpl_loo
  valid_data_loo_list[[j]] <- valid_data_loo
  preds_loo_list[[j]] <- preds_loo
  rmsep_loo_list[[j]] <- rmsep_loo
  models_loo_list[[j]] <- models_loo
}


# get all predictions in one table
preds_loo_pcoa1_taxo_df_all <- list()
models_stats_pcoa1_taxo_loo <- list()
for(j in 1:length(valid_smpl_loo_list)){
  # get prediction vs. observation tables
  preds_loo_df <- data.frame(sample_data = unlist(data.frame(valid_smpl_loo_list[[j]])),
                             observed = unlist(data.frame(valid_data_loo_list[[j]])), 
                             predictions = unlist(data.frame(preds_loo_list[[j]])))
  plot(preds_loo_df$predictions ~ preds_loo_df$observed,main="Test Dataset")
  abline(0,1)
  preds_loo_pcoa1_taxo_df_all[[j]] <- preds_loo_df
  # get statistics from all trhe predictions from iterations
  models_stats_pcoa1_taxo_loo[[j]] <- caret::postResample(preds_loo_df$predictions, preds_loo_df$observed)
  models_stats_pcoa1_taxo_loo[[j]]  <- c(models_stats_pcoa1_taxo_loo[[j]] ,models_stats_pcoa1_taxo_loo[[j]][1]/(max(plots_field_spectral_ok$PCoA_taxo_CA_ok_1) - min(plots_field_spectral_ok$PCoA_taxo_CA_ok_1))*100)
  names(models_stats_pcoa1_taxo_loo[[j]])[4] <- "%RMSE"
  models_stats_pcoa1_taxo_loo[[j]]
  
  # get mean predictions and standard deviation in one table
  mean_sd_preds_loo_df <- data.frame(aggregate(predictions~sample_data, preds_loo_df, function(x) c(mean = mean(x), sd = sd(x))))
  mean_sd_preds_loo_df <- cbind(mean_sd_preds_loo_df[-ncol(mean_sd_preds_loo_df)], mean_sd_preds_loo_df[[ncol(mean_sd_preds_loo_df)]])
  mean_sd_preds_loo_df$observed <- plots_field_spectral_ok$PCoA_taxo_CA_ok_1
  mean_sd_preds_loo_df$plow <- mean_sd_preds_loo_df$mean - mean_sd_preds_loo_df$sd
  mean_sd_preds_loo_df$pup <- mean_sd_preds_loo_df$mean + mean_sd_preds_loo_df$sd
  
  # plot
  plot_CV_pcoa1_taxo_loo <- ggplot(data=mean_sd_preds_loo_df)+aes(x=observed) +
    geom_point(aes(y=mean),color='blue',shape=19) +
    geom_errorbar(aes(ymin=plow,ymax=pup)) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    xlab("Observed") +
    ylab("Predicted") +
    theme_bw()
  plot_CV_pcoa1_taxo_loo
}

# saveRDS(preds_loo_pcoa1_taxo_df_all, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/preds_loo_pcoa1_taxo_df_all.rds")

# saveRDS(models_stats_pcoa1_taxo_loo, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plsr_stats_pcoa1_taxo_bloo.rds")



######################################################################################################
####                                         PLSR for NMDS 1 taxo                                 ####
######################################################################################################
#### see script "updated_pcoa_taxo.R"
#### select best model with all data with a cross validation procedure ####

model <- plsr(NMDS_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="none")
# leave one out 
model <- plsr(NMDS_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="LOO")
# # K-fold CV
# model <- plsr(NMDS_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
#                 S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
#                 S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
#               data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="CV",
#               segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
plot(model, ncomp = best.dims, asp = 1, line = TRUE)
pls::R2(model, estimate = "all") 
CV_loo_nobuffer <- pls::R2(model, estimate = "CV") 

#### get variable importance as the standardized coefficent of predictors ####
library(plsVarSel)
# The regression coefficient represents a predictor variables’ importance in predicting the response
rc <- rev(sort(RC(model, best.dims)))
barplot(rc)
rc_PCoA_1 <- rc/sum(rc)*100
# using caret?
library(caret)
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
landscape_pred_NMDS1_utm <- cbind(landscape_pred_utm, pred_NMDS1 = landsacpe_pred)

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
  calib_smpl <- sample(nrow(plots_field_spectral_ok),nrow(plots_field_spectral_ok) * cal_prop)
  hist(calib_smpl)
  valid_smpl_tmp <- (1:nrow(plots_field_spectral_ok))[!1:nrow(plots_field_spectral_ok) %in% calib_smpl]
  calib_data_tmp <- plots_field_spectral_ok[calib_smpl,]
  valid_data_tmp <- plots_field_spectral_ok[valid_smpl_tmp,]
  model_tmp <- plsr(NMDS_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                      S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                      S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                    data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
  
  
  preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
  rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
  hist(valid_data_tmp$NMDS_1)
  plot(valid_data_tmp$NMDS_1, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
  abline(0, 1, col="red")
  valid_smpl <- cbind(valid_smpl, valid_smpl_tmp)
  valid_data <- cbind(valid_data, valid_data_tmp$NMDS_1)
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
models_stats_NMDS1 <- caret::postResample(preds_df$predictions, preds_df$observed)
models_stats_NMDS1 <- c(models_stats_NMDS1 ,models_stats_NMDS1[1]/(max(plots_field_spectral_ok$NMDS_1) - min(plots_field_spectral_ok$NMDS_1))*100)
names(models_stats_NMDS1)[4] <- "%RMSE"
# get mean predictions and standard deviation in one table
mean_sd_preds_df <- data.frame(aggregate(predictions~sample_data, preds_df, function(x) c(mean = mean(x), sd = sd(x))))
mean_sd_preds_df <- cbind(mean_sd_preds_df[-ncol(mean_sd_preds_df)], mean_sd_preds_df[[ncol(mean_sd_preds_df)]])
mean_sd_preds_df$observed <- plots_field_spectral_ok$NMDS_1
mean_sd_preds_df$plow <- mean_sd_preds_df$mean - mean_sd_preds_df$sd
mean_sd_preds_df$pup <- mean_sd_preds_df$mean + mean_sd_preds_df$sd

# plot
plot_CV_NMDS1 <- ggplot(data=mean_sd_preds_df)+aes(x=observed) +
  geom_point(aes(y=mean),color='blue',shape=19) +
  geom_errorbar(aes(ymin=plow,ymax=pup)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  xlab("Observed") +
  ylab("Predicted") +
  theme_bw()
plot_CV_NMDS1
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
landscape_pred_NMDS1_utm <- cbind(landscape_pred_NMDS1_utm, pred_NMDS1_sd = landscape_pred_sd)
#### export spatial predictions with uncertainties ####
st_write(landscape_pred_NMDS1_utm, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_NMDS1_utm.shp")


#### make models for cross validation with spatial LOO with buffer arount each LOO plots ####
# see Ploton et al. 2020 Nat. Com. 
# get spectial object
plots_field_spectral_ok$geo_pt..89
# calibration and validation datasets to compute uncertainties
# for loop for iterations
buffers <- seq(0,4000,500)
valid_smpl_loo_list <- list()
valid_data_loo_list <- list()
preds_loo_list <- list()
rmsep_loo_list <- list()
models_loo_list <- list()
for(j in 1:length(buffers)){
  print(buffers[j])
  buf_tmp <- buffers[j]
  n_pts = nrow(plots_field_spectral_ok)
  valid_smpl_loo<- c()
  valid_data_loo <- c()
  preds_loo <- c()
  rmsep_loo <- c()
  models_loo <- list()
  pb <- txtProgressBar(max = n_pts, style = 3)
  for (i in 1:n_pts){
    setTxtProgressBar(pb, i)
    valid_smpl_tmp  <- i
    # hist(calib_smpl)
    valid_data_tmp <- plots_field_spectral_ok[i,]
    # remove plot within a buffer around validation plot
    valid_plot <- plots_field_spectral_ok$geo_pt..89[i]
    buffer_tmp <- st_buffer(valid_plot, buf_tmp)
    calib_smpl_spatial <- (1:nrow(plots_field_spectral_ok))[lengths(st_intersects(plots_field_spectral_ok$geo_pt..89, buffer_tmp)) == 0]
    if(buf_tmp == 0)     calib_smpl_spatial <- (1:nrow(plots_field_spectral_ok))[-i]
    calib_data_tmp <- plots_field_spectral_ok[calib_smpl_spatial,]
    plot(plots_field_spectral_ok$geo_pt..89)
    plot(buffer_tmp, add = T)
    plot(plots_field_spectral_ok$geo_pt..89[calib_smpl_spatial], lwd = 5, add = T)
    # model
    model_tmp <- plsr(NMDS_1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                        S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                        S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                      data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
    
    preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
    rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
    # hist(valid_data_tmp$NMDS_1)
    # plot(valid_data_tmp$NMDS_1, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
    # abline(0, 1, col="red")
    valid_smpl_loo <- cbind(valid_smpl_loo, valid_smpl_tmp)
    valid_data_loo <- cbind(valid_data_loo, valid_data_tmp$NMDS_1)
    preds_loo <- cbind(preds_loo, preds_tmp)
    rmsep_loo <- c(rmsep_loo, data.frame(rmsep_tmp$val))
    models_loo[[i]] <- model_tmp
  }
  close(pb)
  valid_smpl_loo_list[[j]] <- valid_smpl_loo
  valid_data_loo_list[[j]] <- valid_data_loo
  preds_loo_list[[j]] <- preds_loo
  rmsep_loo_list[[j]] <- rmsep_loo
  models_loo_list[[j]] <- models_loo
}


# get all predictions in one table
preds_loo_NMDS1_taxo_df_all <- list()
models_stats_NMDS1_taxo_loo <- list()
for(j in 1:length(valid_smpl_loo_list)){
  # get prediction vs. observation tables
  preds_loo_df <- data.frame(sample_data = unlist(data.frame(valid_smpl_loo_list[[j]])),
                             observed = unlist(data.frame(valid_data_loo_list[[j]])), 
                             predictions = unlist(data.frame(preds_loo_list[[j]])))
  plot(preds_loo_df$predictions ~ preds_loo_df$observed,main="Test Dataset")
  abline(0,1)
  preds_loo_NMDS1_taxo_df_all[[j]] <- preds_loo_df
  # get statistics from all trhe predictions from iterations
  models_stats_NMDS1_taxo_loo[[j]] <- caret::postResample(preds_loo_df$predictions, preds_loo_df$observed)
  models_stats_NMDS1_taxo_loo[[j]]  <- c(models_stats_NMDS1_taxo_loo[[j]] ,models_stats_NMDS1_taxo_loo[[j]][1]/(max(plots_field_spectral_ok$NMDS_1) - min(plots_field_spectral_ok$NMDS_1))*100)
  names(models_stats_NMDS1_taxo_loo[[j]])[4] <- "%RMSE"
  models_stats_NMDS1_taxo_loo[[j]]
  
  # get mean predictions and standard deviation in one table
  mean_sd_preds_loo_df <- data.frame(aggregate(predictions~sample_data, preds_loo_df, function(x) c(mean = mean(x), sd = sd(x))))
  mean_sd_preds_loo_df <- cbind(mean_sd_preds_loo_df[-ncol(mean_sd_preds_loo_df)], mean_sd_preds_loo_df[[ncol(mean_sd_preds_loo_df)]])
  mean_sd_preds_loo_df$observed <- plots_field_spectral_ok$NMDS_1
  mean_sd_preds_loo_df$plow <- mean_sd_preds_loo_df$mean - mean_sd_preds_loo_df$sd
  mean_sd_preds_loo_df$pup <- mean_sd_preds_loo_df$mean + mean_sd_preds_loo_df$sd
  
  # plot
  plot_CV_NMDS1_taxo_loo <- ggplot(data=mean_sd_preds_loo_df)+aes(x=observed) +
    geom_point(aes(y=mean),color='blue',shape=19) +
    geom_errorbar(aes(ymin=plow,ymax=pup)) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    xlab("Observed") +
    ylab("Predicted") +
    theme_bw()
  plot_CV_NMDS1_taxo_loo
}

# saveRDS(preds_loo_NMDS1_taxo_df_all, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/preds_loo_NMDS1_taxo_df_all.rds")

# saveRDS(models_stats_NMDS1_taxo_loo, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plsr_stats_NMDS1_taxo_bloo.rds")


######################################################################################################
####                                     PLSR for synthetic traits                                ####
######################################################################################################
#### see script "updated_pcoa_taxo.R"
#### select best model with all data with a cross validation procedure ####

model <- plsr(CWM_trans_synth_CA_Axis1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="none")
# leave one out 
model <- plsr(CWM_trans_synth_CA_Axis1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
              data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="LOO")
# K-fold CV
# model <- plsr(CWM_trans_synth_CA_Axis1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
#                 S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
#                 S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
#               data=plots_field_spectral_ok, scale=TRUE, center = TRUE, validation="CV",
#               segments = 10, segment.type = "random")

summary(model)
plot(RMSEP(model), legendpos = "topright")
cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
plot(model, ncomp = best.dims, asp = 1, line = TRUE)
pls::R2(model, estimate = "all") 

#### get standardized coeficents ####

model
coef(model)

#### get variable importance as the standardized coefficent of predictors ####
library(plsVarSel)
library(caret)
# The regression coefficient represents a predictor variables’ importance in predicting the response
rc <- rev(sort(RC(model, best.dims)))
barplot(rc)
rc_cwm_synth_1 <- rc/sum(rc)*100
# using caret?
var_imp <- varImp(model)
var_imp <- var_imp - min(var_imp)
var_imp <- var_imp / max(var_imp)*100
var_imp <- var_imp[rev(order(var_imp$Overall)),, drop=FALSE]
var_imp
plot(var_imp)
varimp_cwm_synth_1 <- var_imp

#### make spatial predictions at the landscape scale ####
landsacpe_pred <- drop(predict(model, landscape_data, ncomp = best.dims))
#### get spatial object from predictions ####
library(dplyr)
landscape_pred_utm <- landscape_data_utm
landscape_cwm_synth_utm <- cbind(landscape_pred_utm, CWM_trans_synth_CA_Axis1 = landsacpe_pred)

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
  calib_smpl <- sample(nrow(plots_field_spectral_ok),nrow(plots_field_spectral_ok) * cal_prop)
  hist(calib_smpl)
  valid_smpl_tmp <- (1:nrow(plots_field_spectral_ok))[!1:nrow(plots_field_spectral_ok) %in% calib_smpl]
  calib_data_tmp <- plots_field_spectral_ok[calib_smpl,]
  valid_data_tmp <- plots_field_spectral_ok[valid_smpl_tmp,]
  model_tmp <- plsr(CWM_trans_synth_CA_Axis1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                      S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                      S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                    data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
  
  preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
  rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
  hist(valid_data_tmp$CWM_trans_synth_CA_Axis1)
  plot(valid_data_tmp$CWM_trans_synth_CA_Axis1, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
  abline(0, 1, col="red")
  valid_smpl <- cbind(valid_smpl, valid_smpl_tmp)
  valid_data <- cbind(valid_data, valid_data_tmp$CWM_trans_synth_CA_Axis1)
  preds <- cbind(preds, preds_tmp)
  rmsep <- c(rmsep, data.frame(rmsep_tmp$val))
  models[[i]] <- model_tmp
}

# get all predictions in one table
preds_df <- data.frame(sample_data = unlist(data.frame(valid_smpl)),
                       observed = unlist(data.frame(valid_data)), 
                       predictions = unlist(data.frame(preds)))
plot(preds_df$predictions ~ preds_df$observed,main="Test Dataset")
abline(0,1)
# get statistics from all trhe predictions from iterations
models_stats_cwm_synth <- caret::postResample(preds_df$predictions, preds_df$observed)
models_stats_cwm_synth <- c(models_stats_cwm_synth , models_stats_cwm_synth[1]/(max(plots_field_spectral_ok$CWM_trans_synth_CA_Axis1)-min(plots_field_spectral_ok$CWM_trans_synth_CA_Axis1)) *100)
names(models_stats_PCoA1)[4] <- "%RMSE"
models_stats_PCoA1
library(forestmangr)
rmse_per(preds_df, "observed", "predictions", na.rm = TRUE)

# get mean predictions and standard deviation in one table
mean_sd_preds_df <- data.frame(aggregate(predictions~sample_data, preds_df, function(x) c(mean = mean(x), sd = sd(x))))
mean_sd_preds_df <- cbind(mean_sd_preds_df[-ncol(mean_sd_preds_df)], mean_sd_preds_df[[ncol(mean_sd_preds_df)]])
mean_sd_preds_df$observed <- plots_field_spectral_ok$CWM_trans_synth_CA_Axis1
mean_sd_preds_df$plow <- mean_sd_preds_df$mean - mean_sd_preds_df$sd
mean_sd_preds_df$pup <- mean_sd_preds_df$mean + mean_sd_preds_df$sd

# plot
plot_CV_cwm_synth <- ggplot(data=mean_sd_preds_df)+aes(x=observed) +
  geom_point(aes(y=mean),color='blue',shape=19) +
  geom_errorbar(aes(ymin=plow,ymax=pup)) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  xlab("Observed") +
  ylab("Predicted") +
  theme_bw()
plot_CV_cwm_synth


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
landscape_pred_cwm_synth_utm <- cbind(landscape_cwm_synth_utm, CWM_synth_CA_sd = landscape_pred_sd)
#### export spatial predictions with uncertainties ####
st_write(landscape_pred_cwm_synth_utm, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_cwmsynth_utm.shp")


#### make models for cross validation with spatial LOO with buffer arount each LOO plots ####
# see Ploton et al. 2020 Nat. Com. 
# get spectial object
plots_field_spectral_ok$geo_pt..89
# calibration and validation datasets to compute uncertainties
# for loop for iterations
buffers <- seq(0,4000,500)
valid_smpl_loo_list <- list()
valid_data_loo_list <- list()
preds_loo_list <- list()
rmsep_loo_list <- list()
models_loo_list <- list()
for(j in 1:length(buffers)){
  print(buffers[j])
  buf_tmp <- buffers[j]
  n_pts = nrow(plots_field_spectral_ok)
  valid_smpl_loo<- c()
  valid_data_loo <- c()
  preds_loo <- c()
  rmsep_loo <- c()
  models_loo <- list()
  pb <- txtProgressBar(max = n_pts, style = 3)
  for (i in 1:n_pts){
    setTxtProgressBar(pb, i)
    valid_smpl_tmp  <- i
    # hist(calib_smpl)
    valid_data_tmp <- plots_field_spectral_ok[i,]
    # remove plot within a buffer around validation plot
    valid_plot <- plots_field_spectral_ok$geo_pt..89[i]
    buffer_tmp <- st_buffer(valid_plot, buf_tmp)
    calib_smpl_spatial <- (1:nrow(plots_field_spectral_ok))[lengths(st_intersects(plots_field_spectral_ok$geo_pt..89, buffer_tmp)) == 0]
    if(buf_tmp == 0)     calib_smpl_spatial <- (1:nrow(plots_field_spectral_ok))[-i]
    calib_data_tmp <- plots_field_spectral_ok[calib_smpl_spatial,]
    plot(plots_field_spectral_ok$geo_pt..89)
    plot(buffer_tmp, add = T)
    plot(plots_field_spectral_ok$geo_pt..89[calib_smpl_spatial], lwd = 5, add = T)
    # model
    model_tmp <- plsr(CWM_trans_synth_CA_Axis1 ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                        S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                        S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                      data=calib_data_tmp, scale=F, validation="none", ncomp = best.dims)
    
    preds_tmp <- c(drop(predict(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)))
    rmsep_tmp <- RMSEP(model_tmp, newdata = valid_data_tmp, ncomp = best.dims)
    # hist(valid_data_tmp$CWM_trans_synth_CA_Axis1)
    # plot(valid_data_tmp$CWM_trans_synth_CA_Axis1, preds_tmp,main="Test Dataset", xlab="observed", ylab="PLS Predicted")
    # abline(0, 1, col="red")
    valid_smpl_loo <- cbind(valid_smpl_loo, valid_smpl_tmp)
    valid_data_loo <- cbind(valid_data_loo, valid_data_tmp$CWM_trans_synth_CA_Axis1)
    preds_loo <- cbind(preds_loo, preds_tmp)
    rmsep_loo <- c(rmsep_loo, data.frame(rmsep_tmp$val))
    models_loo[[i]] <- model_tmp
  }
  close(pb)
  valid_smpl_loo_list[[j]] <- valid_smpl_loo
  valid_data_loo_list[[j]] <- valid_data_loo
  preds_loo_list[[j]] <- preds_loo
  rmsep_loo_list[[j]] <- rmsep_loo
  models_loo_list[[j]] <- models_loo
}


# get all predictions in one table
preds_loo_cwm_synth_df_all <- list()
models_stats_cwm_synth_loo <- list()
for(j in 1:length(valid_smpl_loo_list)){
    # get prediction vs. observation tables
  preds_loo_df <- data.frame(sample_data = unlist(data.frame(valid_smpl_loo_list[[j]])),
                             observed = unlist(data.frame(valid_data_loo_list[[j]])), 
                             predictions = unlist(data.frame(preds_loo_list[[j]])))
  plot(preds_loo_df$predictions ~ preds_loo_df$observed,main="Test Dataset")
  abline(0,1)
  preds_loo_cwm_synth_df_all[[j]] <- preds_loo_df
  # get statistics from all trhe predictions from iterations
  models_stats_cwm_synth_loo[[j]] <- caret::postResample(preds_loo_df$predictions, preds_loo_df$observed)
  models_stats_cwm_synth_loo[[j]]  <- c(models_stats_cwm_synth_loo[[j]] ,models_stats_cwm_synth_loo[[j]][1]/(max(plots_field_spectral_ok$CWM_trans_synth_CA_Axis1) - min(plots_field_spectral_ok$CWM_trans_synth_CA_Axis1))*100)
  names(models_stats_cwm_synth_loo[[j]])[4] <- "%RMSE"
  models_stats_cwm_synth_loo[[j]]
  
  # get mean predictions and standard deviation in one table
  mean_sd_preds_loo_df <- data.frame(aggregate(predictions~sample_data, preds_loo_df, function(x) c(mean = mean(x), sd = sd(x))))
  mean_sd_preds_loo_df <- cbind(mean_sd_preds_loo_df[-ncol(mean_sd_preds_loo_df)], mean_sd_preds_loo_df[[ncol(mean_sd_preds_loo_df)]])
  mean_sd_preds_loo_df$observed <- plots_field_spectral_ok$CWM_trans_synth_CA_Axis1
  mean_sd_preds_loo_df$plow <- mean_sd_preds_loo_df$mean - mean_sd_preds_loo_df$sd
  mean_sd_preds_loo_df$pup <- mean_sd_preds_loo_df$mean + mean_sd_preds_loo_df$sd
  
  # plot
  plot_CV_cwm_synth_loo <- ggplot(data=mean_sd_preds_loo_df)+aes(x=observed) +
    geom_point(aes(y=mean),color='blue',shape=19) +
    geom_errorbar(aes(ymin=plow,ymax=pup)) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    xlab("Observed") +
    ylab("Predicted") +
    theme_bw()
  plot_CV_cwm_synth_loo
}
# saveRDS(preds_loo_cwm_synth_df_all, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/preds_loo_cwm_synth_df_all.rds")

# saveRDS(models_stats_cwm_synth_loo, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plsr_stats_cwm_synth_bloo.rds")

#### compare R² for taxo and functio ####
# models_stats_pcoa1_taxo_loo <- readRDS( "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plsr_stats_pcoa1_taxo_bloo.rds")
# models_stats_cwm_synth_loo <- readRDS( "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plsr_stats_cwm_synth_bloo.rds")
buffers <- seq(0,4000,500)

models_stats_pcoa1_taxo_loo_df <- data.frame(t(data.frame(models_stats_pcoa1_taxo_loo)))
models_stats_pcoa1_taxo_loo_df$buffer <- buffers

models_stats_cwm_synth_loo_df <- data.frame(t(data.frame(models_stats_cwm_synth_loo)))
models_stats_cwm_synth_loo_df$buffer <- buffers

plot(models_stats_pcoa1_taxo_loo_df$Rsquared ~ models_stats_pcoa1_taxo_loo_df$buffer, col = "darkgreen", ylim = c(0,.7))
lines(models_stats_pcoa1_taxo_loo_df$Rsquared ~ models_stats_pcoa1_taxo_loo_df$buffer, col = "darkgreen")
points(models_stats_cwm_synth_loo_df$Rsquared ~ models_stats_cwm_synth_loo_df$buffer, col = "darkred")
lines(models_stats_cwm_synth_loo_df$Rsquared ~ models_stats_cwm_synth_loo_df$buffer, col = "darkred")

######################################################################################################
#### Analyze predictions (extrapolated biological attributes from S2) ~ fragmentation indices ####
######################################################################################################
#### load spatial prediction ####
landscape_pred_CWM_SLA_utm <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_CWM_SLA_utm.shp")
landscape_pred_CWM_WD_utm <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_CWM_WD_utm.shp")
landscape_pred_PCoA1_utm <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_PCoA1_utm.shp")
landscape_pred_cwm_synth_utm <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_cwmsynth_utm.shp")
landscape_pred_NMDS1_utm <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_NMDS1_utm.shp")

#### lod 30*30 cells with all data from landscape ####
# poly_cells_30m_ok_landscape_env_S2_ok_all <- readRDS( "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/poly_cells_30m_ok_landscape_env_S2_new.rds")

#### add predictions to spatial data ####
poly_cells_30m_ok_landscape_env_S2$preds_CWM_SLA <- landscape_pred_CWM_SLA_utm$pr_CWM_SLA
poly_cells_30m_ok_landscape_env_S2$preds_CWM_SLA_SD <- landscape_pred_CWM_SLA_utm$p_CWM_SLA_
poly_cells_30m_ok_landscape_env_S2$preds_CWM_WD <- landscape_pred_CWM_WD_utm$pr_CWM_WD
poly_cells_30m_ok_landscape_env_S2$preds_CWM_WD_SD <- landscape_pred_CWM_WD_utm$p_CWM_WD_
poly_cells_30m_ok_landscape_env_S2$pred_PCoA1 <- landscape_pred_PCoA1_utm$pr_PCA1
poly_cells_30m_ok_landscape_env_S2$pred_PCoA1_SD <- landscape_pred_PCoA1_utm$p_PCA1_
poly_cells_30m_ok_landscape_env_S2$pred_cwm_synth <- landscape_pred_cwm_synth_utm$CWM___C
poly_cells_30m_ok_landscape_env_S2$pred_cwm_synth_SD <- landscape_pred_cwm_synth_utm$CWM__CA
poly_cells_30m_ok_landscape_env_S2$pred_NMDS1 <- landscape_pred_NMDS1_utm$pr_NMDS1
poly_cells_30m_ok_landscape_env_S2$pred_NMDS1_SD <- landscape_pred_NMDS1_utm$p_NMDS1_

#### transform and add aspect in direction #### 
# get aspect in radians
poly_cells_30m_ok_landscape_env_S2$aspect_rad <- poly_cells_30m_ok_landscape_env_S2$aspect/180*pi
# get northness and eastness
poly_cells_30m_ok_landscape_env_S2$northness = cos(poly_cells_30m_ok_landscape_env_S2$aspect_rad)
poly_cells_30m_ok_landscape_env_S2$eastness = sin(poly_cells_30m_ok_landscape_env_S2$aspect_rad)

hist(poly_cells_30m_ok_landscape_env_S2$northness)
hist(poly_cells_30m_ok_landscape_env_S2$eastness)

#### crop predictions to a defined landscape (exclude litoral forests) ####
landscape_new <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/zone_new/zone_new.shp")
landscape_new <- st_as_sf(landscape_new)
landscape_new <- st_transform(landscape_new, crs(poly_cells_30m_ok_landscape_env_S2))

poly_cells_30m_ok_landscape_env_S2_ok <- st_crop(poly_cells_30m_ok_landscape_env_S2, landscape_new)

#####################################################################################################
####  New edges foret polygons   ####
#####################################################################################################
forest_new <- readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/zone_new/forest_with_roads.shp")
forest_new_utm <- spTransform(forest_new, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
forest_new_utm <- st_as_sf(forest_new_utm)

# check polygon forest validity
forest_new_utm <- st_make_valid(forest_new_utm)
# st_is_valid(forest_new_utm, reason = TRUE)

# make buffer to account for zone edge (no false edge)
zone_buffer <- st_buffer(landscape_new, 500)

# forest_new in the zone
# CP <- as(extent(zone_utm), "SpatialPolygons")
forest_zone_utm <- st_intersection(forest_new_utm, zone_buffer)

forest_zone_utm$id_patch <- 1:length(forest_zone_utm$id)

forest_zone_utm <- as_Spatial(forest_zone_utm)

# patch area
forest_zone_utm$area <- st_area(st_as_sf(forest_zone_utm))
#####################################################################################################
####  FOR ALL landscape : New distance to edge ==> jump to "load" section  ####
#####################################################################################################
### get centroids of cells ####
xy_cells_centroids_utm = st_centroid(poly_cells_30m_ok_landscape_env_S2_ok)
xy_cells_centroids_utm <- as_Spatial(xy_cells_centroids_utm)

# ID patch
patchs_over_cells_utm <- sp::over(xy_cells_centroids_utm,forest_zone_utm)
xy_cells_centroids_utm$id_patch <- patchs_over_cells_utm$id_patch
# patch area
xy_cells_centroids_utm$area <- patchs_over_cells_utm$area

centroids_cells_ok <- xy_cells_centroids_utm
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

poly_cells_30m_ok_landscape_env_S2_ok$dist_edge_new <- dist_edge_centroid
poly_cells_30m_ok_landscape_env_S2_ok$log_dist_edge_new <- log(dist_edge_centroid)

# keep only centroids within zone (not with buffer)
poly_cells_30m_ok_landscape_env_S2_ok <- st_intersection(poly_cells_30m_ok_landscape_env_S2_ok, st_as_sf(zone_utm))

### save ####
# saveRDS(poly_cells_30m_ok_landscape_env_S2_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/poly_cells_30m_ok_landscape_env_S2_new.rds")
#### get all 30*30 cells for the entire landscape with new distane to edge ####
#### load ####

poly_cells_30m_ok_landscape_env_S2_ok_all <- readRDS( "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/poly_cells_30m_ok_landscape_env_S2_new.rds")

#### get same cells for landscape and prediction data ####
poly_cells_30m_ok_landscape_env_S2_ok <- poly_cells_30m_ok_landscape_env_S2_ok[poly_cells_30m_ok_landscape_env_S2_ok$id_cell_poly
                                                                               %in% poly_cells_30m_ok_landscape_env_S2_ok_all$id_cell_poly,]
poly_cells_30m_ok_landscape_env_S2_ok_all <- poly_cells_30m_ok_landscape_env_S2_ok_all[poly_cells_30m_ok_landscape_env_S2_ok_all$id_cell_poly
                                                                               %in% poly_cells_30m_ok_landscape_env_S2_ok$id_cell_poly,]
# verify if the cell correspond
sum(!poly_cells_30m_ok_landscape_env_S2_ok$id_cell_poly ==  poly_cells_30m_ok_landscape_env_S2_ok_all$id_cell_poly)
#### add new distance to edge
poly_cells_30m_ok_landscape_env_S2_ok$dist_edge_new <- poly_cells_30m_ok_landscape_env_S2_ok_all$dist_edge_new
poly_cells_30m_ok_landscape_env_S2_ok$log_dist_edge_new <- log(poly_cells_30m_ok_landscape_env_S2_ok$dist_edge_new)

#### keep only centroids within forest  #### 
poly_cells_30m_ok_landscape_env_S2_ok <- poly_cells_30m_ok_landscape_env_S2_ok_all[!is.na(poly_cells_30m_ok_landscape_env_S2_ok_all$dist_edge_new),]
poly_cells_30m_ok_landscape_env_S2_ok <- poly_cells_30m_ok_landscape_env_S2_ok[!is.na(poly_cells_30m_ok_landscape_env_S2_ok$dist_edge_new),]

#### export map for synthetic traits prediction in new forest fragments ####
landscape_pred_cwm_synth_utm_ok <- poly_cells_30m_ok_landscape_env_S2_ok[,c("geometry", "pred_cwm_synth", "dist_edge_new")]

# st_write(landscape_pred_cwm_synth_utm_ok, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/biovar_preds/landscape_pred_cwm_synth_utm_ok_test.shp")

#####################################################################################################
####  FOR PLOT DATA ONLY: New distance to edge  ####
#####################################################################################################
#### create a spatial object with plots ####
plots_geom_utm_ok <- plots_field_spectral_ok$geo_pt..89

#### get centroids of plots ####
xy_plots_centroids_utm = st_centroid(plots_geom_utm_ok)
xy_plots_centroids_utm <- as_Spatial(xy_plots_centroids_utm)

# ID patch
xy_plots_centroids_utm <- spTransform(xy_plots_centroids_utm, crs(forest_zone_utm))
patchs_over_plots_utm <- sp::over(xy_plots_centroids_utm,forest_zone_utm)
xy_plots_centroids_utm$id_patch <- patchs_over_plots_utm$id_patch
# patch area
xy_plots_centroids_utm$area <- patchs_over_plots_utm$area

centroids_plots_ok <- xy_plots_centroids_utm
#### distance to edge ####
n_pts <- nrow(centroids_plots_ok)
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
      pts_tmp <- centroids_plots_ok[i,]
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
dist_edge_centroid

# add distance to edge 
plots_field_spectral_ok$dist_edge_new <- dist_edge_centroid
plots_field_spectral_ok$log_dist_edge_new <- log(plots_field_spectral_ok$dist_edge_new)

#####################################################################################################
#### (jump to load) FOR ALL landscape : get sol type ####
#####################################################################################################
# geology <- rgdal::readOGR("/home/thesardfou/Documents/GIS/geologie/Geologie_SHP_Lyr_1000000/SurfaceGeologique_1000000.shp")
# geology <- st_as_sf(geology)
# # check polygon forest validity
# geology <- st_make_valid(geology)
# geology <- st_transform(geology, crs(poly_cells_30m_ok_landscape_env_S2_ok))
# ### get centroids of cells ####
# xy_cells_centroids_utm <- st_centroid(poly_cells_30m_ok_landscape_env_S2_ok)
# centroids_cells_ok <- xy_cells_centroids_utm
# centroids_cells_geol <- st_intersection(geology, centroids_cells_ok)
# 
# poly_cells_30m_ok_landscape_env_S2_ok <- poly_cells_30m_ok_landscape_env_S2_ok[poly_cells_30m_ok_landscape_env_S2_ok$id_cell_poly %in% centroids_cells_geol$id_cell_poly,]
# poly_cells_30m_ok_landscape_env_S2_ok$geology <- centroids_cells_geol$lithologie
# 
# table(poly_cells_30m_ok_landscape_env_S2_ok$geology)
# 
# poly_cells_30m_ok_landscape_env_S2_ok$geology[poly_cells_30m_ok_landscape_env_S2_ok$geology %in% c("Dunites", "Gabbros cumulats") ] <- "Péridotites"
# poly_cells_30m_ok_landscape_env_S2_ok$geology[poly_cells_30m_ok_landscape_env_S2_ok$geology %in% c("Fluvio-lacustre") ] <- "Cuirasses"
# poly_cells_30m_ok_landscape_env_S2_ok$geology <- as.factor(poly_cells_30m_ok_landscape_env_S2_ok$geology)
# 
# 
# #### use 1/50000 geology with only iron crust ####
# geology2 <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/zone_new/iron_crust.shp")
# geology2 <- st_as_sf(geology2)
# # check polygon forest validity
# geology2 <- st_make_valid(geology2)
# geology2 <- st_transform(geology2, crs(poly_cells_30m_ok_landscape_env_S2_ok))
# xy_cells_centroids_utm <- st_centroid(poly_cells_30m_ok_landscape_env_S2_ok)
# centroids_cells_ok <- xy_cells_centroids_utm
# centroids_cells_geol2 <- sapply(st_intersects(centroids_cells_ok, geology2),function(x){length(x)>0})
# centroids_cells_geol2 <- ifelse(centroids_cells_geol2, "Iron crust", "Eroded ferritic soils")
# centroids_cells_geol2 <- data.frame(lithologie = centroids_cells_geol2, id_cell_poly = centroids_cells_ok$id_cell_poly)
# poly_cells_30m_ok_landscape_env_S2_ok <- poly_cells_30m_ok_landscape_env_S2_ok[poly_cells_30m_ok_landscape_env_S2_ok$id_cell_poly %in% centroids_cells_geol2$id_cell_poly,]
# poly_cells_30m_ok_landscape_env_S2_ok$geology2 <- centroids_cells_geol2$lithologie
# 
# table(poly_cells_30m_ok_landscape_env_S2_ok$geology2)

#### load ####
poly_cells_30m_ok_landscape_env_S2_ok_env <- poly_cells_30m_ok_landscape_env_S2_ok

# saveRDS(poly_cells_30m_ok_landscape_env_S2_ok_env, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/poly_cells_30m_ok_landscape_env_S2_ok_env.rds")
poly_cells_30m_ok_landscape_env_S2_ok_env <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/fragmentation_data/poly_cells_30m_ok_landscape_env_S2_ok_env.rds")
#####################################################################################################
####  update predictions from plsr if needed (NMDS?) #### 
#####################################################################################################
pred_NMDS1_up <- data.frame(pred_NMDS1 = poly_cells_30m_ok_landscape_env_S2$pred_NMDS1,
                            pred_NMDS1_SD = poly_cells_30m_ok_landscape_env_S2$pred_NMDS1_SD,
                            id_cell_poly = poly_cells_30m_ok_landscape_env_S2$id_cell_poly)

poly_cells_30m_ok_landscape_env_S2_ok_env <- merge(x=poly_cells_30m_ok_landscape_env_S2_ok_env,y=pred_NMDS1_up, 
             by="id_cell_poly", all.x=TRUE)
#####################################################################################################
#### (jump to load) FOR PLOT DATA ONLY: get landscape metrics, environmental variable and sol type ####
#####################################################################################################
# 
# ##### landscape metrics for diferent landscape buffers ##### 
# forest_map_raster_utm <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/forest_raster/forest_raster_GT_utm.tif")
# # buffer zone to get local entire lanscapes around invotory plots (1km)
# zone_utm <- spTransform(zone, CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# zone_buffer <- gBuffer(zone_utm, width = 1000)
# # crop forest raster with zone
# forest_map_raster_zone_utm <- crop(forest_map_raster_utm, zone_buffer)
# 
# plots_geom_utm_ok_sp <- as_Spatial(plots_geom_utm_ok)
# # buffer values 
# buffers = c(50, 100, 250,500,1000)
# # buffers = c(250,500)
# n_pts <- length(plots_geom_utm_ok_sp)
# library(SDMTools)
# library(utils)
# pb <- txtProgressBar(max = n_pts, style = 3)
# local_stats_centroid <- c()
# for(i in 1:n_pts){
#   setTxtProgressBar(pb,i)
#   local_stats_all_tmp <- c()
#   for(j in length(buffers):1)  {
#     points_spl_buffer_tmp <- gBuffer(plots_geom_utm_ok_sp[i], width=buffers[j], byid=TRUE)
#     if(j == length(buffers)){
#       local_habitat =  crop(forest_map_raster_zone_utm, points_spl_buffer_tmp)
#     }else{
#       local_habitat =  crop(local_habitat, points_spl_buffer_tmp)
#     }
#     local_habitat =  mask(local_habitat, points_spl_buffer_tmp)
#     local_stats_tmp = ClassStat(local_habitat)
#     colnames(local_stats_tmp) = paste0( colnames(local_stats_tmp), "_", buffers[j], "_centroid")
#     local_stats_tmp = local_stats_tmp[nrow(local_stats_tmp),]
#     local_stats_all_tmp <- c(local_stats_all_tmp, unlist(local_stats_tmp))
#   }
#   local_stats_centroid <- rbind(local_stats_centroid,local_stats_all_tmp)
# }
#   
# 
# plots_field_spectral_ok <- cbind(plots_field_spectral_ok, local_stats_centroid)
# #### get env data for all plots ####
# library(exactextractr)
# dem <- raster("/home/thesardfou/Documents/GIS/MNT/mnt10_GT.tif")
# slope <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/slope_GT.tif")
# aspect <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/aspect_GT.tif")
# curv <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/curvature_GT.tif")
# twi <- raster("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/twi_GT.tif")
# prec <- raster("/home/thesardfou/Documents/GIS/precipitations/precipitations_utm58s.tif")
# 
# cells_alt <-  exact_extract(dem, plots_geom_utm_ok, "mean")
# cells_slope <-  exact_extract(slope, plots_geom_utm_ok, "mean")
# cells_aspect <-  exact_extract(aspect, plots_geom_utm_ok, "mean")
# cells_curv <-  exact_extract(curv, plots_geom_utm_ok, "mean")
# cells_twi <-  exact_extract(twi, plots_geom_utm_ok, "mean")
# cells_prec <-  exact_extract(prec, plots_geom_utm_ok, "mean")
# 
# plots_field_spectral_ok$elevation <- cells_alt
# plots_field_spectral_ok$slope <- cells_slope
# plots_field_spectral_ok$aspect <- cells_aspect
# plots_field_spectral_ok$curvature <- cells_curv
# plots_field_spectral_ok$twi <- cells_twi
# plots_field_spectral_ok$prec <- cells_prec
# 
# #### soil type ####
# 
# geology <- rgdal::readOGR("/home/thesardfou/Documents/GIS/geologie/Geologie_SHP_Lyr_1000000/SurfaceGeologique_1000000.shp")
# geology <- st_as_sf(geology)
# # check polygon forest validity
# geology <- st_make_valid(geology)
# geology <- st_transform(geology, crs(plots_geom_utm_ok))
# get centroids of cells
# xy_plots_centroids_utm <- st_centroid(plots_geom_utm_ok)
# 
# centroids_plots_ok <- xy_plots_centroids_utm
# 
# centroids_plots_ok_geol <- st_intersection(geology, centroids_plots_ok)
# 
# table(centroids_plots_ok_geol$lithologie)
# 
# plots_field_spectral_ok$geology <- centroids_plots_ok_geol$lithologie
# 
# plots_field_spectral_ok$geology[plots_field_spectral_ok$geology %in% c("Dunites", "Gabbros cumulats") ] <- "Péridotites"
# plots_field_spectral_ok$geology[plots_field_spectral_ok$geology %in% c("Fluvio-lacustre") ] <- "Cuirasses"
# plots_field_spectral_ok$geology <- as.factor(plots_field_spectral_ok$geology)
# plots_field_spectral_ok_env <- plots_field_spectral_ok
### use 1/50000 geology with only iron crust ####
# geology2 <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/zone_new/iron_crust.shp")
# geology2 <- st_as_sf(geology2)
# # check polygon forest validity
# geology2 <- st_make_valid(geology2)
# geology2 <- st_transform(geology2, crs(plots_geom_utm_ok))
# # get centroids of cells
# xy_plots_centroids_utm <- st_centroid(plots_geom_utm_ok)
# 
# centroids_plots_ok <- xy_plots_centroids_utm
# 
# centroids_plots_ok_geol <- sapply(st_intersects(centroids_plots_ok, geology2),function(x){length(x)>0})
# centroids_plots_ok_geol <- ifelse(centroids_plots_ok_geol, "Iron crust", "Eroded ferritic soils")
# table(centroids_plots_ok_geol)
# 
# plots_field_spectral_ok_env$geology2 <- centroids_plots_ok_geol
# save
# saveRDS(plots_field_spectral_ok_env, "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_field_spectral_ok_env.rds")

#### load ####
plots_field_spectral_ok_env <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plot_data_matrices/plots_field_spectral_ok_env.rds")

#####################################################################################################
####  Analyze predictions ~ environment + fragmentation  ####
#####################################################################################################
#### covariation between observed variables and between predictions ####
# observations
cor_obs <- cor( plots_field_spectral_ok_env[,c("CWM_SLA","CWM_WD", "PCoA_taxo_CA_ok_1", "CWM_trans_synth_CA_Axis1")],
                method = "pearson", use = "complete.obs")
cor_obs
# predictions
cor_pred <- cor( data.frame(poly_cells_30m_ok_landscape_env_S2_ok_env[, c("preds_CWM_SLA", "preds_CWM_WD", "pred_PCoA1", "pred_cwm_synth")])[,-5],
                 method = "pearson", use = "complete.obs")
cor_pred <- cor( data.frame(poly_cells_30m_ok_landscape_env_S2_ok_env[, c("preds_CWM_SLA", "preds_CWM_WD", "pred_NMDS1", "pred_cwm_synth")])[,-5],
                 method = "pearson", use = "complete.obs")
cor_pred

cor_pred_predictors <- cor( data.frame(poly_cells_30m_ok_landscape_env_S2_ok_env[, c("log_dist_edge_new",
                                                                                 "elevation","slope",                        
                                                                                 "curvature", "twi",
                                                                                 "northness", "eastness")])[,-8],
                 method = "pearson", use = "complete.obs")
cor_pred_predictors

#### get influence for plots  ####
plots_field_spectral_ok_env$geology2 <- as.factor(plots_field_spectral_ok_env$geology2)
plots_field_spectral_ok_env
summary(lm(plots_field_spectral_ok_env$CWM_trans_synth_CA_Axis1 ~ plots_field_spectral_ok_env$log_dist_edge_new))
plot(plots_field_spectral_ok_env$CWM_trans_synth_CA_Axis1 ~ plots_field_spectral_ok_env$log_dist_edge_new)
plot(plots_field_spectral_ok_env$CWM_trans_synth_CA_Axis1 ~ plots_field_spectral_ok_env$geology2)

#### get dataframe for multivariate model selection  ####
poly_cells_30m_ok_landscape_env_S2_df <- data.frame(poly_cells_30m_ok_landscape_env_S2_ok_env)
colnames(poly_cells_30m_ok_landscape_env_S2_df)

#### compute perimeter / area ratio for local landscapes ####
poly_cells_30m_ok_landscape_env_S2_df$perimeter.area.ratio_50_centroid <- 
  poly_cells_30m_ok_landscape_env_S2_df$edge.density_50_centroid / poly_cells_30m_ok_landscape_env_S2_df$total.area_50_centroid
#### remove outliers from SLA prediction (visual chack shows that negative values are on non forest zones) ####
# poly_cells_30m_ok_landscape_env_S2_df <- poly_cells_30m_ok_landscape_env_S2_df[poly_cells_30m_ok_landscape_env_S2_df$preds_CWM_SLA>0,]

#### remove cells centroids a less than a given distance from the edge ####
 poly_cells_30m_ok_landscape_env_S2_df_test <- poly_cells_30m_ok_landscape_env_S2_df[poly_cells_30m_ok_landscape_env_S2_df$dist_edge_new>15,]

#### Keep only predictions within the observed range ####
# poly_cells_30m_ok_landscape_env_S2_df_test <- 
#   poly_cells_30m_ok_landscape_env_S2_df_test[poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth >
#                                                min(plots_field_spectral_ok_env$CWM_trans_synth_CA_Axis1) & 
#                                                poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth <
#                                                max(plots_field_spectral_ok_env$CWM_trans_synth_CA_Axis1),]

#### test relationships ####
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test$log_dist_edge_new))
hist(poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new)
hist(poly_cells_30m_ok_landscape_env_S2_df_test$log_dist_edge_new)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test$elevation))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test$curvature))

summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test$geology))


poly_cells_30m_ok_landscape_env_S2_df_test_core <- poly_cells_30m_ok_landscape_env_S2_df_test[poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new>100,]
poly_cells_30m_ok_landscape_env_S2_df_test_edge <- poly_cells_30m_ok_landscape_env_S2_df_test[poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new<100,]

summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test_core$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test_core$curvature))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test_core$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test_core$elevation))

# compare with forest withour roads
poly_cells_30m_ok_landscape_env_S2_ok_all
summary(lm(poly_cells_30m_ok_landscape_env_S2_ok_all$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_ok_all$log_dist_edge_new))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test$log_dist_edge_new))
plot(poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test$log_dist_edge_new)

summary(lm(poly_cells_30m_ok_landscape_env_S2_ok_all$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_ok_all$log_dist_edge))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test$log_dist_edge))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth ~ poly_cells_30m_ok_landscape_env_S2_df_test$area))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$log_dist_edge_new))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape.core_100_centroid))
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$log_dist_edge_new)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$log_dist_edge_new))

plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new))

plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$elevation)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$elevation))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$edge.density_500_centroid))
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$edge.density_500_centroid)
#
plot(poly_cells_30m_ok_landscape_env_S2_df_test$edge.density_500_centroid ~ poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge)
#
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$total.edge_500_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$edge.density_500_centroid))

summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$total.edge_1000_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$edge.density_1000_centroid))

plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_1000_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_1000_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$edge.density_100_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$edge.density_500_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$total.edge_500_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$effective.mesh.size_500_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$effective.mesh.size_1000_centroid)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$effective.mesh.size_500_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$effective.mesh.size_1000_centroid))

plot(poly_cells_30m_ok_landscape_env_S2_df_test$total.edge_1000_centroid ~ poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge)
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$total.edge_1000_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_1000_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df_test$total.edge_1000_centroid ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_1000_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df_test$total.edge_500_centroid ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_500_centroid)
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$total.edge_500_centroid)

plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_100_centroid)
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_100_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_100_centroid))

summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$prop.landscape_50_centroid))

plot(poly_cells_30m_ok_landscape_env_S2_df_test$perimeter.area.ratio_500_centroid ~ poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge)
plot(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ log(poly_cells_30m_ok_landscape_env_S2_df_test$perimeter.area.ratio_500_centroid))
summary(lm(poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA ~ poly_cells_30m_ok_landscape_env_S2_df_test$perimeter.area.ratio_500_centroid))

#### select variables ####
# choose set of explanatory variable

expl_v_list_list <- list(c( "log_dist_edge", "prop.landscape_100_centroid", "prop.landscape_250_centroid", "prop.landscape_500_centroid",
                            "elevation","slope","aspect",                        
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


expl_v_list_list <- list(c( "log_dist_edge", "prop.landscape_50_centroid",
                            "elevation","slope",                        
                            "curvature", 
                            "northness", "eastness"))


expl_v_list_list <- list(c( "log_dist_edge_new", 
                            "elevation","slope",                        
                            "curvature", "twi",
                            "geology2",
                            "northness", "eastness"))


# set of response variables
resp_v_list <- c( "preds_CWM_SLA",  "preds_CWM_SLA_SD",
                  "preds_CWM_WD", "preds_CWM_WD_SD",
                  "pred_NMDS1", "pred_NMDS1_SD",
                  "pred_cwm_synth", "pred_cwm_synth_SD")


#### model selection ####

library(MuMIn)
library(parallel)
# total set of variable
data_for_mod_all <- poly_cells_30m_ok_landscape_env_S2_df_test
# no na ?
data_for_mod_all <- data_for_mod_all[!is.na(data_for_mod_all$elevation),]
#### FILTER BY ELEVATION? ####
# data_for_mod_all <- data_for_mod_all[data_for_mod_all$elevation>100,]


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
    if("geology" %in% colnames(data_for_mod_toscale)) data_for_mod_toscale$geology <- as.numeric(data_for_mod_toscale$geology)
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

#### no na  ####

data_for_varpart <- poly_cells_30m_ok_landscape_env_S2_df_test[complete.cases(poly_cells_30m_ok_landscape_env_S2_df_test[,c("preds_CWM_SLA",
                                                                                                                  "preds_CWM_WD",
                                                                                                                  "pred_PCoA1",
                                                                                                                  "elevation","slope",                        
                                                                                                                  "curvature", "twi",
                                                                                                                  "northness", "eastness",
                                                                                                                  "geology",
                                                                                                                  "log_dist_edge_new")]),]

colnames(data_for_varpart)

Topography <- data.frame(data_for_varpart[,c("elevation","slope",                        
                                             "curvature", "twi",
                                             "geology",
                                             "northness", "eastness")])

Edge_influence <- data.frame(log_dist_edge_new = data_for_varpart[,c("log_dist_edge_new")])



#### preds_sp_rich ####
# relative importance
relaimpo::calc.relimp(data.frame(data_for_varpart[c("preds_CWM_SLA",
                                                    "elevation","slope",                        
                                                    "curvature", "twi",
                                                    "northness", "eastness",
                                                    "geology",
                                                    "log_dist_edge_new")]))

relaimpo::calc.relimp(data.frame(data_for_varpart[c("preds_CWM_WD",
                                                    "elevation","slope",                        
                                                    "curvature", "twi",
                                                    "northness", "eastness",
                                                    "geology",
                                                    "log_dist_edge_new")]))

relaimpo::calc.relimp(data.frame(data_for_varpart[c("pred_PCoA1",
                                                    "elevation","slope",                        
                                                    "curvature", "twi",
                                                    "northness", "eastness",
                                                    "geology",
                                                    "log_dist_edge_new")]))

relaimpo::calc.relimp(data.frame(data_for_varpart[c("pred_cwm_synth",
                                                    "elevation","slope",                        
                                                    "curvature", "twi",
                                                    "northness", "eastness",
                                                    "geology",
                                                    "log_dist_edge_new")]), rela = T)
# variance partitioning with multiple responses 
library("rdacca.hp")
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13800
pred <- data.frame(pred <- data_for_varpart$preds_CWM_SLA)
pred <- data.frame(pred <- data_for_varpart$preds_CWM_WD)
pred <- data.frame(pred <- data_for_varpart$pred_PCoA1)
pred <- data.frame(pred <- data_for_varpart$pred_cwm_synth)

# pred <- data.frame(data_for_varpart[,c("preds_CWM_SLA","preds_CWM_WD", "pred_PCoA1", "pred_cwm_synth")])

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
####                                        Plots linear relationships                            ####
######################################################################################################
library(multcompView)
library(ggplot2)
# synthetic trait
data_plot1 <- data.frame(dist_edge = poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new,
                         CWM_PC1 = poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth)

lm1 <- lm(data_plot1$CWM_PC1 ~ log(data_plot1$dist_edge))
smry1 = summary(lm1)
smry1

plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = CWM_PC1)) + 
  geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.3) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = CWM_PC1), formula = y ~ log(x), method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=200, y= 5,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted CWM PC1") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

# synthetic trait with plot data
data_plot1 <- data.frame(dist_edge = poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new,
                         CWM_PC1 = poly_cells_30m_ok_landscape_env_S2_df_test$pred_cwm_synth)

lm1 <- lm(data_plot1$CWM_PC1 ~ log(data_plot1$dist_edge))
smry1 = summary(lm1)
smry1

data_plot2 <- data.frame(dist_edge = plots_field_spectral_ok$dist_edge_new,
                         CWM_PC1 = plots_field_spectral_ok$CWM_trans_synth_CA_Axis1)

lm2 <- lm(data_plot2$CWM_PC1 ~ log(data_plot2$dist_edge))
smry2 = summary(lm2)
smry2


plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = CWM_PC1)) + 
  geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.3) +
  geom_point(data = data_plot2, color='darkgreen', size = 3, alpha =.8) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = CWM_PC1), formula = y ~ log(x), method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  geom_smooth(data = data_plot2, aes(x = dist_edge, y = CWM_PC1), formula = y ~ log(x), method = "lm", color='darkred',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=200, y= 5,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  
  ylab("Predicted CWM PC1") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

# log linear
data_plot2 <- data_plot1
data_plot2$dist_edge <- log(data_plot2$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = log(data_plot1$dist_edge), y = CWM_PC1)) + 
  geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.8) +
  geom_smooth(data = data_plot1, aes(x = log(data_plot1$dist_edge), y = CWM_PC1), formula = y ~ x, method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=5, y= 3,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted CWM PC1") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot
# hist
data_plot1$dist_edge_class <- cut(data_plot1$dist_edge,
                                  breaks=seq(0,500,50))

data_plot1$dist_edge_class <- cut(data_plot1$dist_edge,
                                  breaks=c(0,50,100,500))
library(ggplot2)

ggplot(data = data_plot1, aes(x=as.factor(data_plot1$dist_edge_class), y=CWM_PC1)) +
  geom_boxplot(fill="steelblue") +
  labs( x="Distance to forest edge", y="Predicted CWM PC1") +
  theme_bw()

## tukey test  ##
# analysis of variance
anova <- aov(CWM_PC1 ~ dist_edge_class, data = data_plot1)
summary(anova)
# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)

cld <- multcompLetters4(anova, tukey)
print(cld)
# table with factors and 3rd quantile
Tk <- group_by(data_plot1, dist_edge_class) %>%
  summarise(mean=mean(CWM_PC1), quant = quantile(CWM_PC1, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$dist_edge_class)
Tk$cld <- cld$Letters

print(Tk)
# boxplot
ggplot(data_plot1, aes(dist_edge_class, CWM_PC1)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE, outlier.shape = NA) +
  labs( x="Distance to forest edge", y="Predicted CWM PC1") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = Tk, aes(x = dist_edge_class, y = quant, label = cld), size = 5, vjust=-1, hjust =-1) +
  scale_fill_brewer(palette = "Blues") + 
  coord_cartesian(ylim = quantile(data_plot1$CWM_PC1, c(0.005, 0.9995)))




# SLA
data_plot1 <- data.frame(dist_edge = poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new,
                         CWM_SLA = poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_SLA)

lm1 <- lm(data_plot1$CWM_SLA ~ log(data_plot1$dist_edge))
smry1 = summary(lm1)
smry1

plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = CWM_SLA)) + 
  geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.8) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = CWM_SLA), formula = y ~ log(x), method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=200, y= 14,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted CWM SLA") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot
# log linear
data_plot2 <- data_plot1
data_plot2$dist_edge <- log(data_plot2$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = log(data_plot1$dist_edge), y = CWM_SLA)) + 
  geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.8) +
  geom_smooth(data = data_plot1, aes(x = log(data_plot1$dist_edge), y = CWM_SLA), formula = y ~ x, method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=5, y= 3,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted CWM PC1") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot
# hist
data_plot1$dist_edge_class <- cut(data_plot1$dist_edge,
                                  breaks=seq(0,500,50))

data_plot1$dist_edge_class <- cut(data_plot1$dist_edge,
                                  breaks=c(0,50,100,500))
library(ggplot2)

ggplot(data = data_plot1, aes(x=as.factor(data_plot1$dist_edge_class), y=CWM_SLA)) +
  geom_boxplot(fill="steelblue") +
  labs( x="Distance to forest edge", y="Predicted CWM SLA") +
  theme_bw()

## tukey test  ##
# analysis of variance
anova <- aov(CWM_SLA ~ dist_edge_class, data = data_plot1)
summary(anova)
# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)

cld <- multcompLetters4(anova, tukey)
print(cld)
# table with factors and 3rd quantile
Tk <- group_by(data_plot1, dist_edge_class) %>%
  summarise(mean=mean(CWM_SLA), quant = quantile(CWM_SLA, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$dist_edge_class)
Tk$cld <- cld$Letters

print(Tk)
# boxplot
ggplot(data_plot1, aes(dist_edge_class, CWM_SLA)) + 
  geom_boxplot() +
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  labs( x="Distance to forest edge", y="Predicted CWM SLA") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = Tk, aes(x = dist_edge_class, y = quant, label = cld), size = 3, vjust=-1, hjust =-1) +
  scale_fill_brewer(palette = "Blues")

# WD

data_plot1 <- data.frame(dist_edge = poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new,
                         CWM_WD = poly_cells_30m_ok_landscape_env_S2_df_test$preds_CWM_WD)

lm1 <- lm(data_plot1$CWM_WD ~ log(data_plot1$dist_edge))
smry1 = summary(lm1)
smry1

plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = CWM_WD)) + 
   geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.8) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = CWM_WD), formula = y ~ log(x), method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=200, y= 1,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted CWM WD") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

# no points
plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = CWM_WD)) + 
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = CWM_WD), formula = y ~ log(x), method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=200, y= .77,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted CWM WD") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot


# log linear
data_plot2 <- data_plot1
data_plot2$dist_edge <- log(data_plot2$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = log(data_plot1$dist_edge), y = CWM_WD)) + 
  geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.8) +
  geom_smooth(data = data_plot1, aes(x = log(data_plot1$dist_edge), y = CWM_WD), formula = y ~ x, method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=5, y= .77,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted CWM WD") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot


# hist
data_plot1$dist_edge_class <- cut(data_plot1$dist_edge,
                                  breaks=seq(0,500,50))
data_plot1$dist_edge_class <- cut(data_plot1$dist_edge,
                                  breaks=c(0,50,100,500))
library(ggplot2)

ggplot(data = data_plot1, aes(x=as.factor(data_plot1$dist_edge_class), y=CWM_WD)) +
  geom_boxplot(fill="steelblue") +
  labs( x="Distance to forest edge", y="Predicted CWM WD") +
  theme_bw()
## tukey test  ##
# analysis of variance
anova <- aov(CWM_WD ~ dist_edge_class, data = data_plot1)
summary(anova)
# Tukey's test
tukey <- TukeyHSD(anova)
print(tukey)

cld <- multcompLetters4(anova, tukey)
print(cld)
# table with factors and 3rd quantile
Tk <- group_by(data_plot1, dist_edge_class) %>%
  summarise(mean=mean(CWM_WD), quant = quantile(CWM_WD, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$dist_edge_class)
Tk$cld <- cld$Letters

print(Tk)
# boxplot
ggplot(data_plot1, aes(dist_edge_class, CWM_WD)) + 
  geom_boxplot() +
  geom_boxplot(aes(fill = factor(..middle..)), show.legend = FALSE) +
  labs( x="Distance to forest edge", y="Predicted CWM WD") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = Tk, aes(x = dist_edge_class, y = quant, label = cld), size = 3, vjust=-1, hjust =-1) +
  scale_fill_brewer(palette = "Blues")

# taxo
data_plot1 <- data.frame(dist_edge = poly_cells_30m_ok_landscape_env_S2_df_test$dist_edge_new,
                         CWM_PCoA1 = poly_cells_30m_ok_landscape_env_S2_df_test$pred_PCoA1)

lm1 <- lm(data_plot1$CWM_PCoA1 ~ log(data_plot1$dist_edge))
smry1 = summary(lm1)
smry1

plot <- ggplot(data = data_plot1, aes(x = dist_edge, y = CWM_PCoA1)) + 
  geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.8) +
  geom_smooth(data = data_plot1, aes(x = dist_edge, y = CWM_PCoA1), formula = y ~ log(x), method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=200, y= 0.15,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted PCoA1 taxo") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot


# log linear
data_plot2 <- data_plot1
data_plot2$dist_edge <- log(data_plot2$dist_edge)
plot <- ggplot(data = data_plot1, aes(x = log(data_plot1$dist_edge), y = CWM_PCoA1)) + 
  geom_point(data = data_plot1, color='darkgoldenrod3', size = 1, alpha =.8) +
  geom_smooth(data = data_plot1, aes(x = log(data_plot1$dist_edge), y = CWM_PCoA1), formula = y ~ x, method = "lm", color='navy',
              se=T,fill = "navy", alpha =.2, size=.7) +
  annotate("text", x=5, y= .075,
           label= paste0("R²=", round(smry1$r.squared, digits = 3), gtools::stars.pval(smry1$coefficients[8])),
           color='navy', size = 3.3) +
  ylab("Predicted PCoA1 taxo") + xlab("Distance to forest edge") +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  ) +
  theme_bw()

plot

# hist
data_plot1$dist_edge_class <- cut(data_plot1$dist_edge,
                                 breaks=seq(0,500,50))
library(ggplot2)

ggplot(data = data_plot1, aes(x=as.factor(data_plot1$dist_edge_class), y=CWM_PCoA1)) +
  geom_boxplot(fill="steelblue") +
  labs( x="Distance to forest edge", y="Predicted PCoA1 taxo") +
  theme_bw()

######################################################################################################
####                                         other approaches                                     ####
######################################################################################################

#### manual LOO validation ####
X <- as.matrix(plots_field_spectral_ok[,c ("S2.mean.B2","S2.mean.B3","S2.mean.B4",
                             "S2.mean.B5","S2.mean.B6","S2.mean.B7","S2.mean.B8",
                             "S2.mean.B8A","S2.mean.B11","S2.mean.B12")])
Y <- plots_field_spectral_ok$CWM_SLA
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
  data = plots_field_spectral_ok,
  method = 'pls',
  preProcess = c("center", "scale"),
  tuneGrid = tuneGrid,
  trControl = ctrl
)
model_caret
plot(model_caret)

plot(varImp(model_caret))
varImp(model_caret)


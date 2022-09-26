 

library(raster)
library(sp)
library(sf)
S2 <- stack("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_sudNC_58KFA/58KFA_raw/S2_58KFA_int_raw.tif")

forest <- rgdal::readOGR("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/GIS/forest_map/forest_present_ok.shp")
forest <- spTransform(forest, "+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs")


S2_zone <- crop(S2, extent(forest))
plot(S2_zone)

#### extract sentinel band values for plots ####
merge_all_KB <- readRDS("/home/thesardfou/Documents/projets/Reliques/projet_kuebini/R/data/merge_all_KB.rds")
# do not keep  stations outside forest ! 
merge_all_KB <- merge_all_KB[merge_all_KB$dist_edge>0,]
merge_all_KB <- merge_all_KB[!is.na(merge_all_KB$coordX),]

spdf_merge_all_KB <- SpatialPointsDataFrame(coords = cbind(merge_all_KB$coordX, merge_all_KB$coordY), data = merge_all_KB,
                                            proj4string = CRS("+proj=utm +zone=58 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot_data_sf <- st_as_sf(spdf_merge_all_KB)
plot_areas <- st_buffer( plot_data_sf, dist=11.3, byid=TRUE)

extract_S2 <- extract(S2_zone, plot_areas,fun=mean, na.rm=TRUE )

#### compil data ####
merge_all_KB_S2 <- cbind(merge_all_KB, extract_S2)

#### partil least square regression ####
library(pls)
# data
merge_all_KB_S2_plot <- merge_all_KB_S2[!is.na(merge_all_KB_S2$sp_richness),]
col_S2 <- grep(pattern = "S2_58KFA", colnames(merge_all_KB_S2_plot))
S2_bands <- colnames(merge_all_KB_S2_plot)[col_S2]

# choose variables
colnames(merge_all_KB_S2_plot)
var_resp <- c("SLA", "WD", "sp_richness", "FD.RaoQ")
col_tmp <- which(colnames(merge_all_KB_S2_plot) %in% var_resp)
data_for_mod_S2 <- merge_all_KB_S2_plot[, c(col_tmp ,col_S2)]

# format data
S2 = as.matrix(data_for_mod_S2[,! colnames(data_for_mod_S2) %in% var_resp])
SLA = data_for_mod_S2[,colnames(data_for_mod_S2) == "SLA"]
WD = data_for_mod_S2[,colnames(data_for_mod_S2) == "WD"]
sp_richness = data_for_mod_S2[,colnames(data_for_mod_S2) == "sp_richness"]
FD.RaoQ = data_for_mod_S2[,colnames(data_for_mod_S2) == "FD.RaoQ"]

data_for_mod_S2_ok <- data.frame(I(S2), SLA,WD, sp_richness, FD.RaoQ)

#### SLA ####
model <- plsr(SLA ~ S2, data=data_for_mod_S2_ok, scale=TRUE, validation="CV")
model <- plsr(SLA ~ S2, data=data_for_mod_S2_ok, scale=TRUE, validation="LOO")

summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

plot(RMSEP(model), legendpos = "topright")
### choose the number of components
plot(model, ncomp = 5, asp = 1, line = TRUE) # SLA

R2(model, estimate = "all") 
R2(model, estimate = "train") 




#### WD ####

model <- plsr(WD ~ S2, data=data_for_mod_S2_ok, scale=TRUE, validation="CV")
model <- plsr(WD ~ S2, data=data_for_mod_S2_ok, scale=TRUE, validation="LOO")

summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

plot(RMSEP(model), legendpos = "topright")
### choose the number of components
plot(model, ncomp = 4, asp = 1, line = TRUE) # WD

R2(model, estimate = "all") 
R2(model, estimate = "train") 


#### SP richness ####

model <- plsr(sp_richness ~ S2, data=data_for_mod_S2_ok, scale=TRUE, validation="CV")
model <- plsr(sp_richness ~ S2, data=data_for_mod_S2_ok, scale=TRUE, validation="LOO")

summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

plot(RMSEP(model), legendpos = "topright")
### choose the number of components
plot(model, ncomp = 6, asp = 1, line = TRUE) 

R2(model, estimate = "all") 
R2(model, estimate = "train") 



#### FD ####

model <- plsr(FD.RaoQ ~ S2, data=data_for_mod_S2_ok, scale=TRUE, validation="CV")
model <- plsr(FD.RaoQ ~ S2, data=data_for_mod_S2_ok, scale=TRUE, validation="LOO")

summary(model)
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

plot(RMSEP(model), legendpos = "topright")
### choose the number of components
plot(model, ncomp = 3, asp = 1, line = TRUE) 

R2(model, estimate = "all") 
R2(model, estimate = "train") 




plot(model, plottype = "scores", comps = 1:3)
explvar(model)
plot(model)

plot(model, "loadings", comps = 1:2, legendpos = "topleft")
abline(h = 0)


data(yarn)

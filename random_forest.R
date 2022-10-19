




# Loading package
library(randomForest)
################################################################
#### load data ####
################################################################


plots_field_spectral <- readRDS("/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/NEW/maps/58FA_2021/plots_field_spectral_30x30plots_100sp_formPC1678.rds")

dim(plots_field_spectral)
################################################################
#### Splitting data in train and test data ####
################################################################
# Split into Train and Validation sets
# Training Set : Validation Set = 70 : 30 (random)
set.seed(100)
train <- sample(nrow(plots_field_spectral), 0.7*nrow(plots_field_spectral), replace = FALSE)
TrainSet <- plots_field_spectral[train,]
ValidSet <- plots_field_spectral[-train,]
summary(TrainSet)
summary(ValidSet)

colnames(plots_field_spectral)

################################################################
####  Create a Random Forest model with default parameters #### 
################################################################
model1 <- randomForest(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                         S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                         S2.mean.B8A + S2.mean.B11 + S2.mean.B12,
                       data = TrainSet, ntree = 500,      mtry = 6, importance = TRUE)
model1

# Predicting on train set
predTrain <- predict(model1, TrainSet, type = "class")
confusion_mtx = table(TrainSet$CWM_SLA, predTrain)



################################################################
####  Running random forest regression #### 
################################################################
library(caret)
# we are not going to do any cross-validatin
# and rely on OOB error
trctrl <- trainControl(method = "none")

# we will now train random forest model
rfregFit <- train(CWM_SLA ~ S2.mean.B2 + S2.mean.B3 + S2.mean.B4 +
                    S2.mean.B5 + S2.mean.B6 + S2.mean.B7 + S2.mean.B8 +
                    S2.mean.B8A + S2.mean.B11 + S2.mean.B12, 
                  data = plots_field_spectral, 
                  method = "ranger",
                  trControl=trctrl,
                  # calculate importance
                  importance="permutation", 
                  tuneGrid = data.frame(mtry=6,
                                        min.node.size = 5,
                                        splitrule="variance")
)
https://compgenomr.github.io/book/predicting-continuous-variables-regression-with-machine-learning.html


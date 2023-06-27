rm(list = ls())

library(fastshap)
library(ranger)
library(doParallel)
library(ggbeeswarm)
library(viridis)
library(colorspace)
library(cowplot)
library(tidyverse)
library(caret)

registerDoParallel(12)

#####for Z values#####
###load the composite values 
#compositeVals<-read.csv("biomes_composite_diver_stab.csv") ##this file is too heavy for git
compositeVals<-read.csv("data/raw_data/compositve_values.csv) ##this is a summarized version already subsetting to variables of interest
##transfor for the Grid ID to be the rowname

##match the composite values to the Z-values from all_metrics
all_metrics<-read.csv("data/results/REV_obs_results_200_zvals.csv")  

compositeVals[match(all_metrics$grid_id,compositeVals$Grid),"zPd"]<-all_metrics$Zpd
compositeVals[match(all_metrics$grid_id,compositeVals$Grid),"zMPD"]<-all_metrics$Zmpd
compositeVals[match(all_metrics$grid_id,compositeVals$Grid),"zRao"]<-all_metrics$Zraoq
compositeVals[match(all_metrics$grid_id,compositeVals$Grid),"zfdr"]<-all_metrics$Zfdr


##########################
######for raw values######
##########################

compositeVals[match(all_metrics$grid_id,compositeVals$Grid),"PD"]<-all_metrics$pd
compositeVals[match(all_metrics$grid_id,compositeVals$Grid),"MPD"]<-all_metrics$mpd
compositeVals[match(all_metrics$grid_id,compositeVals$Grid),"Rao"]<-all_metrics$raoq
compositeVals[match(all_metrics$grid_id,compositeVals$Grid),"fdr"]<-all_metrics$fdr



# read in the data and add in an ID field
dt1 <- compositeVals %>%
  mutate(id = 1:nrow(.))
# specify the predictor variables to use

##selection of predicted variable here (must change below in 2 other instances + legend)
##options are fdr, Rao, MPD or PD
coef_det <- function(xtrue, xpred){
  return(1-sum((xtrue-xpred)^2)/sum((xtrue-mean(xtrue))^2))
}
indices<-c("fdr", "Rao", "MPD","PD")
varimp <- data.frame()
for (i in 1:length(indices)){
  evars <- c(
    "CHELSA_BIO_Annual_Mean_Temperature",
    "CHELSA_BIO_Annual_Precipitation",
    "CHELSA_BIO_Max_Temperature_of_Warmest_Month",
    "CHELSA_BIO_Mean_Diurnal_Range",
    "CHELSA_BIO_Isothermality",
    "CHELSA_BIO_Min_Temperature_of_Coldest_Month",
    "CHELSA_BIO_Temperature_Annual_Range",
    "CHELSA_BIO_Precipitation_of_Driest_Month",
    "CHELSA_BIO_Precipitation_of_Warmest_Quarter",
    "CHELSA_exBIO_AridityIndex",
    "CGIAR_PET",
    "EarthEnvTopoMed_Slope",
    "EarthEnvTopoMed_Elevation",
    #"FanEtAl_Depth_to_Water_Table_AnnualMean",
    "SG_Bulk_density_015cm",
    "SG_CEC_015cm",
    "SG_Clay_Content_060cm",
    "SG_H2O_Capacity_015cm",
    "SG_SOC_Content_015cm",
    "WCS_Human_Footprint_2009",
    #	"BesnardEtAl_ForestAge_TC020",
    "CSP_Global_Human_Modification",
    "CHELSAv1_BIO_Precipitation_Seasonality",
    "GranthamEtAl_ForestLandscapeIntegrityIndex",
    "SG_Absolute_depth_to_bedrock", 
    "SG_Coarse_fragments_000cm",
    "SG_Sand_Content_000cm", 
    "SG_Soil_pH_H2O_000cm",
    "StabilityCSI",
    "mrd"
  )
dt <- dt1 %>% select(indices[i], all_of(evars), id) %>% filter(complete.cases(.))
##new names
dt<-rename(dt, c(
"Climate_Annual_Mean_Temperature"="CHELSA_BIO_Annual_Mean_Temperature",
"Climate1_Annual_Precipitation"="CHELSA_BIO_Annual_Precipitation",
"Climate_Max_Temperature_of_Warmest_Month"="CHELSA_BIO_Max_Temperature_of_Warmest_Month",
"Climate_Mean_Diurnal_Range"="CHELSA_BIO_Mean_Diurnal_Range",
"Climate_Isothermality"="CHELSA_BIO_Isothermality",
"Climate_Min_Temperature_of_Coldest_Month"="CHELSA_BIO_Min_Temperature_of_Coldest_Month",
"Climate_Temperature_Annual_Range"="CHELSA_BIO_Temperature_Annual_Range",
"Climate1_Precipitation_of_Driest_Month"="CHELSA_BIO_Precipitation_of_Driest_Month",
"Climate1_Precipitation_of_Warmest_Quarter"="CHELSA_BIO_Precipitation_of_Warmest_Quarter",
"Climate2_AridityIndex"="CHELSA_exBIO_AridityIndex",
"Climate2_PET"="CGIAR_PET",
"Topo_Slope"="EarthEnvTopoMed_Slope",
"Topo_Elevation"="EarthEnvTopoMed_Elevation",
"Human_Human_Footprint_2009"="WCS_Human_Footprint_2009",
#	"BesnardEtAl_ForestAge_TC020",
"Human_Global_Human_Modification"="CSP_Global_Human_Modification",
"Climate1_Precipitation_Seasonality"="CHELSAv1_BIO_Precipitation_Seasonality",
"Human_ForestLandscapeIntegrityIndex"="GranthamEtAl_ForestLandscapeIntegrityIndex",
"Climate_UStabilityCSI"="StabilityCSI",
"A_mrd"="mrd"
)
)
evars <- c(
  "Climate_Annual_Mean_Temperature",
  "Climate1_Annual_Precipitation",
  "Climate_Max_Temperature_of_Warmest_Month",
  "Climate_Mean_Diurnal_Range",
  "Climate_Isothermality",
  "Climate_Min_Temperature_of_Coldest_Month",
  "Climate_Temperature_Annual_Range",
  "Climate1_Precipitation_of_Driest_Month",
  "Climate1_Precipitation_of_Warmest_Quarter",
  "Climate2_AridityIndex",
  "Climate2_PET",
  "Topo_Slope",
  "Topo_Elevation",
  #"FanEtAl_Depth_to_Water_Table_AnnualMean",
  "SG_Bulk_density_015cm",
  "SG_CEC_015cm",
  "SG_Clay_Content_060cm",
  "SG_H2O_Capacity_015cm",
  "SG_SOC_Content_015cm",
  "Human_Human_Footprint_2009",
  #	"BesnardEtAl_ForestAge_TC020",
  "Human_Global_Human_Modification",
  "Climate1_Precipitation_Seasonality",
  "Human_ForestLandscapeIntegrityIndex",
  "SG_Absolute_depth_to_bedrock", 
  "SG_Coarse_fragments_000cm",
  "SG_Sand_Content_000cm", 
  "SG_Soil_pH_H2O_000cm",
  "Climate_UStabilityCSI",
  "A_mrd"
)
set.seed(15)

# number of test vars for the shapley values. should be set to be a reasonable proportion of the full data
ntest <- 1200

# number of sims for shapley values. larger is better (~100), but takes a long time. lower is fine to testing
nsim <- 10

# prediction function for fastshap
pfun <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}


# create the test train split
test_data <- dt %>% sample_n(ntest)
train_data <- dt %>% filter(!id%in%test_data$id)

# Fit the random forest models to the training data
form <- as.formula(paste(indices[i], # name of response
                         "~", # as function of
                        ".")) # all predictors

#r1 <- ranger(form, data = train_data %>% select(indices[i], all_of(evars)), importance = "permutation")
model_caret <- train(form, data = train_data %>% select(indices[i], all_of(evars)),
                     method = "ranger",
                     trControl = trainControl(method="cv", number = 5),
                     importance = 'permutation')
#print(paste("R2 for:",r1$call,":",r1$r.squared,sep=" "))
###get some stats with testing data for manuscript
#preds <- predict(r1, data = test_data)$predictions	
preds <- predict(model_caret$finalModel, data = test_data)$predictions	
results <- test_data %>% mutate(pred = preds) %>% select(indices[i],pred)

#  Estimate the shapley values on the test data
fs1 <- fastshap::explain(model_caret$finalModel, X = test_data %>% select(all_of(evars)) %>% as.matrix(), pred_wrapper = pfun, nsim = nsim, adjust = TRUE)
print(model_caret$finalModel)
print(paste("R2 for:",indices[i],coef_det(results[,indices[i]], results$pred),sep=" "))
print(paste("RMSE for:",indices[i],sqrt(mean((results[,indices[i]] - results$pred)^2)),sep=" "))
# combine the models and merge in the environmental covariates
shap_df <- bind_rows(fs1 %>% as_tibble() %>% mutate(id = test_data$id) %>% gather(var, shap, -id) %>% mutate(type = "shap", axis = 1)) %>% 
	left_join(test_data %>% select(id, all_of(evars)) %>% gather(var, value, -id))


varimp1<- shap_df %>% group_by(var) %>% summarise(feature_imp = mean(abs(shap))) %>% ungroup %>%arrange(var)
varimp1$Index<-rep(indices[i],length(varimp1$var))
varimp1$feature_imp<-varimp1$feature_imp/(sum(varimp1$feature_imp))*100
varimp <- varimp %>% bind_rows(varimp1)
}
#####Create the figure with var imp 
varimp$Index<-as.factor(varimp$Index)
#  mutate(var = factor(var, levels = rev(shap_ord))) %>%  
varimp%>%  mutate(Index = factor(Index,c("PD","fdr","MPD","Rao"))) %>% 
  ggplot(aes(x = Index, y = var,size=feature_imp,col=feature_imp))+
  geom_point()+
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10))+
  coord_equal(ratio = 0.60)
#Create the "correlation" plot [has to be executed after each fastshap analysis)]




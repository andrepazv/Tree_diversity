rm(list = ls())

library(fastshap)
library(ranger)
library(doParallel)
library(ggbeeswarm)
library(viridis)
library(colorspace)
library(cowplot)
library(tidyverse)

registerDoParallel(12)

#####for Z values#####
###load the composite values 
compositeVals<-read.csv("20220718_andrea_grid_sampled.csv")
##transfor for the Grid ID to be the rowname

##match the composite values to the Z-values from all_metrics
all_metrics<-read.csv("obs_null_z_scores_fixed.csv")  #Final simulations from DAn all but FR


metrics<-subset(all_metrics,metric=="pd")
compositeVals[match(metrics$grid_id,compositeVals$Grid),"zPD"]<-metrics$z

metrics<-subset(all_metrics,metric=="mpd")
compositeVals[match(metrics$grid_id,compositeVals$Grid),"zMPD"]<-metrics$z

metrics<-subset(all_metrics,metric=="raoq")
compositeVals[match(metrics$grid_id,compositeVals$Grid),"zRao"]<-metrics$z
all_metrics<-read.csv("obs_null_z_scores_fixFD.csv") #Final simulations from DAn fixing FR
metrics<-subset(all_metrics,metric=="fdr")
compositeVals[match(metrics$grid_id,compositeVals$Grid),"zfdr"]<-metrics$z


##########################
######for raw values######
##########################


###load the composite values 
compositeVals<-read.csv("biomes_composite_diver_stab.csv")
##transfor for the Grid ID to be the rowname

##match the composite values to the values from all_metrics
all_metrics<-read.csv("obs_null_z_scores_fixed.csv")  #Final simulations from DAn all but FR

metrics<-subset(all_metrics,metric=="pd")
compositeVals[match(metrics$grid_id,compositeVals$Grid),"PD"]<-metrics$value

metrics<-subset(all_metrics,metric=="mpd")
compositeVals[match(metrics$grid_id,compositeVals$Grid),"MPD"]<-metrics$value

metrics<-subset(all_metrics,metric=="raoq")
compositeVals[match(metrics$grid_id,compositeVals$Grid),"Rao"]<-metrics$value

all_metrics<-read.csv("obs_null_z_scores_fixFD.csv") #Final simulations from DAn fixing FR

metrics<-subset(all_metrics,metric=="fdr")
compositeVals[match(metrics$grid_id,compositeVals$Grid),"fdr"]<-metrics$value
# read in the data and add in an ID field
dt <- compositeVals %>%
	mutate(id = 1:nrow(.))

# specify the predictor variables to use
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

##selection of predicted variable here (must change below in 2 other instances + legend)
##options are fdr, Rao, MPD or PD
dt <- dt %>% select(MPD, all_of(evars), id) %>% filter(complete.cases(.))
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
r1 <- ranger(MPD~., data = train_data %>% select(MPD, all_of(evars)), importance = "permutation")
print(paste("R2 for:",r1$call,":",r1$r.squared,sep=" "))
#  Estimate the shapley values on the test data
fs1 <- fastshap::explain(r1, X = test_data %>% select(all_of(evars)) %>% as.matrix(), pred_wrapper = pfun, nsim = nsim, .parallel = TRUE, adjust = TRUE)

# combine the models and merge in the environmental covariates
shap_df <- bind_rows(fs1 %>% as_tibble() %>% mutate(id = test_data$id) %>% gather(var, shap, -id) %>% mutate(type = "shap", axis = 1)) %>% 
	left_join(test_data %>% select(id, all_of(evars)) %>% gather(var, value, -id))

# number of variables to plot, by overall importance
#numvars <- 6

# baseline alpha transparency values
#al <- c(0.15, 0.5)

# ylims for beeswarm plot
##change this according to Shapley value range in each prediction
#ylim <- c(-40,40) #fpr MPD
#ylim <- c(-2000,2000) ###for PD
#ylim <- c(-3,3) ###For Rao
#ylim<-c(-0.4,0.4) ##For FDr
# blank plots for top 3 predictors of each axis#
#gg <- vector(mode = "list", length = 1*numvars)

# blank plots for beeswarm plots
#varimp <- vector(mod = "list", length = 1)

# colors and shapes for plotting
#col1 <- "#1A85FF"
#col2 <- "#D41159"
#my_bins <- 15
#my_shape <- 23
#my_cols <- c("red", "darkblue", "orange", "darkslategray3")

# initialize plotting index
#count <- 0

# scale about zero, removing outliers if needed
#symmetric_scale <- function(x, trim = 0.025){
	#x <- x  - median(x)
	#x[x<0] <- x[x<0]/abs(quantile(x, trim))
	#x[x>0] <- x[x>0]/abs(quantile(x, 1-trim))
	#x[x>1] <- 1
	#x[x< (-1)] <- -1
	#return(x)
#}


# plot for each axis/variable of interest. right now just set up for a single variable
#for(j in 1:1){
	
	# reset the count
#	count <- ifelse(count >= numvars*2, 0, count)
	
	# select the axis and clean up some of the variables
#	my_df <- shap_df %>% filter(axis == j) 

	# sort the variables by sum(abs(shap value))
#	shap_ord <- my_df %>% group_by(var) %>% summarise(feature_imp = mean(abs(shap))) %>% ungroup %>% arrange(desc(feature_imp)) %>% select(var) %>% unlist() %>% as.character()
	
	# create the beeswarm plot
#	varimp[[j]] <- my_df %>% group_by(var) %>% mutate(value = symmetric_scale(value)) %>% ungroup %>%
#		mutate(var = factor(var, levels = rev(shap_ord))) %>% 
#		ggplot(aes(x = var, y = shap, col = value))+
#		geom_quasirandom(dodge.width = 0, bandwidth = 0.05, cex = 0.03, groupOnX = TRUE, alpha = 0.95)+coord_flip()+
#		ylab(paste0("Influence on outcome (Shapley value)"))+
#		theme_bw()+
#		theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10))+
#		scale_color_viridis(option = "viridis", limits = c(-1,1), breaks = seq(-1, 1, length = 10), labels = c("Low", rep("", 10-2), "High"))+
#		labs(color = "Feature value")+ylim(ylim)
	
	# get the three strongest variables for plotting the full relationship
#	var_short <- shap_ord[1:numvars]
	
	# cycle through the top 3 and creat the subplot
#	for(i in 1:numvars){
#		# keep track of plotting index
#		count <- count+1
		
		# subset the data to the current variable
#		sub_df <- my_df %>% filter(var == var_short[i])
		
		# get the x range and create the x sequence
#		x_range <- sub_df %>% select(value) %>% distinct() %>% summarize(max = max(value), min = min(value)) %>% unlist %>% as.numeric %>% sort()

		# get symmetrical ylimits
#		ylims <- sub_df %>% filter(value>x_range[1], value<x_range[2]) %>% summarize(max(abs(max(shap)), abs(min(shap)))) %>% unlist()
		
		# scale the alpha transparencies, if needed
#		aluse <- c(1,1)
		
		# create the plot
#		ggtmp <- 
#			ggplot(data = sub_df %>% filter(value>x_range[1], value<x_range[2]), aes(x = value, y = shap))+
#			geom_point(color = "darkblue", alpha = 0.15)+
#			geom_hline(yintercept = 0, linetype = 2, color = "gray50")+
#			theme_bw()+
#			scale_y_continuous(expand = c(0, 0), limits = c(-ylims, ylims), name = paste0("Influence on FD (Shapley value)")) + 
#			scale_color_viridis()+
#			scale_fill_viridis()+
#			scale_alpha_manual(values = aluse)+
#			geom_smooth() + #method = "lm", formula = as.formula(paste0("y~poly(x, ",deg[i],", raw = TRUE)")), color = "black")+
#			scale_x_continuous(expand = c(0, 0), name = var_short[i])  

		# remove the legend
#		gg[[count]] <- ggtmp + theme(legend.position = "none",
#									 plot.margin = unit(c(3,20,3,3), "pt"))
		
#	}
#}


# creat the subplots
#p1 <- plot_grid(plotlist = varimp, ncol = 1, byrow = FALSE) 
#p2 <- plot_grid(plotlist = gg, nrow = 2, byrow = TRUE)

# plot Fig 3
#plot_grid(p1, NULL, p2, ncol = 3, rel_widths = c(0.7, 0.05, 1))


#####Create the figure with var imp 

#Create the "correlation" plot

##for FR

varimp <- shap_df %>% group_by(var) %>% summarise(feature_imp = mean(abs(shap))) %>% ungroup %>%arrange(var)
varimp$Index<-rep("FDr",length(varimp$var))
varimp$feature_imp<-varimp$feature_imp/(sum(varimp$feature_imp))*100
#varimp$feature_imp <- (varimp$feature_imp-min(varimp$feature_imp))/ (max(varimp$feature_imp) -min(varimp$feature_imp))

##For Rao
varimp1 <- shap_df %>% group_by(var) %>% summarise(feature_imp = mean(abs(shap))) %>% ungroup %>%arrange(var)
varimp1$Index<-rep("Rao",length(varimp1$var))
varimp1$feature_imp<-varimp1$feature_imp/(sum(varimp1$feature_imp))*100
#varimp1$feature_imp <- (varimp1$feature_imp-min(varimp1$feature_imp))/ (max(varimp1$feature_imp) -min(varimp1$feature_imp))

##for PD
varimp2 <- shap_df %>% group_by(var) %>% summarise(feature_imp = mean(abs(shap))) %>% ungroup %>%arrange(var)
varimp2$Index<-rep("PD",length(varimp2$var))
varimp2$feature_imp<-varimp2$feature_imp/(sum(varimp2$feature_imp))*100
#varimp2$feature_imp <- (varimp2$feature_imp-min(varimp2$feature_imp))/ (max(varimp2$feature_imp) -min(varimp2$feature_imp))

##for MPD
varimp3 <- shap_df %>% group_by(var) %>% summarise(feature_imp = mean(abs(shap))) %>% ungroup %>%arrange(var)
varimp3$Index<-rep("MPD",length(varimp3$var))
varimp3$feature_imp<-varimp3$feature_imp/(sum(varimp3$feature_imp))*100
#varimp3$feature_imp <- (varimp3$feature_imp-min(varimp3$feature_imp))/ (max(varimp3$feature_imp) -min(varimp3$feature_imp))


varimp<-varimp %>% add_row(varimp1)%>% add_row(varimp2)%>%add_row(varimp3)

varimp$Index<-as.factor(varimp$Index)
#  mutate(var = factor(var, levels = rev(shap_ord))) %>%  
varimp%>%  mutate(Index = factor(Index,c("PD","FDr","MPD","Rao"))) %>% 
 ggplot(aes(x = Index, y = var,size=feature_imp,col=feature_imp))+
    geom_point()+
    scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10))+
  coord_equal(ratio = 0.60)
#  scale_color_viridis(option = "viridis", limits = c(-1,1), breaks = seq(-1, 1, length = 10), labels = c("Low", rep("", 10-2), "High"))+
  #labs(color = "Feature value")+ylim(ylim)

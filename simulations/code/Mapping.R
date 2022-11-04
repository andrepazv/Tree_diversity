####Mapping dimensions of diversity from community matrix


######################################################################
#########MAPPING OBSERVED METRICS AND RESULTS FROM NULL MODELS #######
######################################################################

setwd("simulations/results")
library(raster)
###load a map worth grid numbers or matrix?

###load the matrix with the data and subset by metric 

all_metrics<-read.csv("obs_null_z_scores_fixed.csv")  #Final simulations all but FD indices

#all_metrics<-read.csv("obs_null_z_scores.csv")
#all_metrics<-read.csv("obs_results_no_roots.csv") ##no roots
#all_metrics<-read.csv("obs_results_old_mpd.csv")  ##as mpd
#all_metrics<-read.csv("obs_results_old_no_square.csv")  #no square
#all_metrics<-read.csv("obs_results_no_roots_old.csv")  #mine?
#all_metrics<-read.csv("redo_old_script.csv")  #mine?
#all_metrics<-read.csv("redo_old_script_orig_fixed.csv")  #Dan's redo
#all_metrics<-read.csv("obs_null_z_scores_50vars.csv") #Testing traits imputed with 50 vars


#all_metrics<-all_metrics[all_metrics$raoq>0,]

###Create empty rasters based on richness map
mocklayer <- raster("/SR_all_trees_observed.tif")
mocklayer<-init(mocklayer,"cell")
names(mocklayer) <- "Grid"
r<-mocklayer
values(r)<-NA

###Create all metric rasters for observed, Z-values and average from null models (AV) 

pd_ras<-r
mpd_ras<-r
raoq_ras<-r
SR_ras<-r
fdr_ras<-r
fd_ras<-r

Zpd<-r
Zmpd<-r
Zraoq<-r
Zfdr<-r
Zfd<-r

pdAV<-r
fdAV<-r
raoAV<-r
mpdAV<-r
#2- Assign values to raster: observed and Z vals
metrics<-subset(all_metrics,metric=="pd")
pd_ras[as.numeric(metrics$grid_id)]<-metrics$value
Zpd[as.numeric(metrics$grid_id)]<-metrics$z
pdAV[as.numeric(metrics$grid_id)]<-metrics$null_mean
SR_ras[as.numeric(metrics$grid_id)]<-metrics$n

metrics<-subset(all_metrics,metric=="mpd")
mpd_ras[as.numeric(metrics$grid_id)]<-metrics$value
Zmpd[as.numeric(metrics$grid_id)]<-metrics$z
mpdAV[as.numeric(metrics$grid_id)]<-metrics$null_mean

metrics<-subset(all_metrics,metric=="raoq")
metrics<-metrics[metrics$value>0,]
raoq_ras[as.numeric(metrics$grid_id)]<-metrics$value
Zraoq[as.numeric(metrics$grid_id)]<-metrics$z
raoAV[as.numeric(metrics$grid_id)]<-metrics$null_mean


all_metrics<-read.csv("obs_null_z_scores_fixFD.csv") #Final simulations from DAn fixing FR
metrics<-subset(all_metrics,metric=="fdr")
fdr_ras[as.numeric(metrics$grid_id)]<-metrics$value
Zfdr[as.numeric(metrics$grid_id)]<-metrics$z
fdAV[as.numeric(metrics$grid_id)]<-metrics$null_mean

metrics<-subset(all_metrics,metric=="fd")
fd_ras[as.numeric(metrics$grid_id)]<-metrics$value
Zfd[as.numeric(metrics$grid_id)]<-metrics$z

####WriteRasters to files for maps
writeRaster(pd_ras,"PD_all_trees_observed.tif",format="GTiff")
writeRaster(mpd_ras,"MPD_all_trees_observed.tif",format="GTiff")
writeRaster(SR_ras,"SR_all_trees_observed.tif",format="GTiff")
writeRaster(raoq_ras,"RaoQ_all_trees_observed.tif",format="GTiff")
writeRaster(fdr_ras,"fdr_all_trees_observed.tif",format="GTiff")
writeRaster(fd_ras,"fd_all_trees_observed.tif",format="GTiff")

writeRaster(Zpd,"Zpd_all_trees_observed.tif",format="GTiff")
writeRaster(Zmpd,"Zmpd_all_trees_observed.tif",format="GTiff")
writeRaster(Zraoq,"Zraoq_all_trees_observed.tif",format="GTiff")
writeRaster(Zfdr,"Zfdr_all_trees_observed.tif",format="GTiff")
writeRaster(Zfd,"Zfd_all_trees_observed.tif",format="GTiff")


#####################################
### MAPPING of bootstrap results ###
#####################################


mocklayer <- raster("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Diversity_computation/No_monocots/Trees_nomono_matchPDFD_30_100_PD_Richness.tif")
mocklayer<-init(mocklayer,"cell")
names(mocklayer) <- "Grid"
r<-mocklayer
values(r)<-NA


###Create all metric rasters 

pd_rasB<-r
mpd_rasB<-r
raoq_rasB<-r
SR_rasB<-r
fdr_rasB<-r
fd_rasB<-r



sd_pd_rasB<-r
sd_mpd_rasB<-r
sd_raoq_rasB<-r
sd_SR_rasB<-r
sd_fdr_rasB<-r
sd_fd_rasB<-r

cv_pd_rasB<-r
cv_mpd_rasB<-r
cv_raoq_rasB<-r
cv_SR_rasB<-r
cv_fdr_rasB<-r
cv_fd_rasB<-r

#2- Assign values to raster

boot<-read.csv("bootstrap_aggregated_mean.csv")
names(boot)
pd_rasB[as.numeric(boot$grid_id)]<-boot$pd
mpd_rasB[as.numeric(boot$grid_id)]<-boot$mpd
raoq_rasB[as.numeric(boot$grid_id)]<-boot$raoq
SR_rasB[as.numeric(boot$grid_id)]<-boot$n
fdr_rasB[as.numeric(boot$grid_id)]<-boot$fdr
fd_rasB[as.numeric(boot$grid_id)]<-boot$fd

writeRaster(pd_rasB,"Bootstrap_PD_all_trees.tif",format="GTiff")
writeRaster(mpd_rasB,"Bootstrap_MPD_all_trees.tif",format="GTiff")
writeRaster(SR_rasB,"Bootstrap_SR_all_trees.tif",format="GTiff")
writeRaster(raoq_rasB,"Bootstrap_RaoQ_all_trees.tif",format="GTiff")
writeRaster(fdr_rasB,"Bootstrap_fdr_all_trees.tif",format="GTiff")
writeRaster(fd_rasB,"Bootstrap_fd_all_trees.tif",format="GTiff")

boot<-read.csv("bootstrap_aggregated_sd.csv")


sd_pd_rasB[as.numeric(boot$grid_id)]<-boot$pd
sd_mpd_rasB[as.numeric(boot$grid_id)]<-boot$mpd
sd_raoq_rasB[as.numeric(boot$grid_id)]<-boot$raoq
sd_SR_rasB[as.numeric(boot$grid_id)]<-boot$n
sd_fdr_rasB[as.numeric(boot$grid_id)]<-boot$fdr
sd_fd_rasB[as.numeric(boot$grid_id)]<-boot$fd



writeRaster(sd_pd_rasB,"Bootstrap_SD_PD_all_trees.tif",format="GTiff")
writeRaster(sd_mpd_rasB,"Bootstrap_SD_MPD_all_trees.tif",format="GTiff")
writeRaster(sd_SR_rasB,"Bootstrap_SD_SR_all_trees.tif",format="GTiff")
writeRaster(sd_raoq_rasB,"Bootstrap_SD_RaoQ_all_trees.tif",format="GTiff")
writeRaster(sd_fdr_rasB,"Bootstrap_SD_fdr_all_trees.tif",format="GTiff")
writeRaster(sd_fd_rasB,"Bootstrap_SD_fd_all_trees.tif",format="GTiff")


boot<-read.csv("bootstrap_aggregated_cv.csv")


cv_pd_rasB[as.numeric(boot$grid_id)]<-boot$pd
cv_mpd_rasB[as.numeric(boot$grid_id)]<-boot$mpd
cv_raoq_rasB[as.numeric(boot$grid_id)]<-boot$raoq
cv_SR_rasB[as.numeric(boot$grid_id)]<-boot$n
cv_fdr_rasB[as.numeric(boot$grid_id)]<-boot$fdr
cv_fd_rasB[as.numeric(boot$grid_id)]<-boot$fd



writeRaster(cv_pd_rasB,"Bootstrap_CV_PD_all_trees.tif",format="GTiff")
writeRaster(cv_mpd_rasB,"Bootstrap_CV_MPD_all_trees.tif",format="GTiff")
writeRaster(cv_SR_rasB,"Bootstrap_CV_SR_all_trees.tif",format="GTiff")
writeRaster(cv_raoq_rasB,"Bootstrap_CV_RaoQ_all_trees.tif",format="GTiff")
writeRaster(cv_fdr_rasB,"Bootstrap_CV_fdr_all_trees.tif",format="GTiff")
writeRaster(cv_fd_rasB,"Bootstrap_CV_fd_all_trees.tif",format="GTiff")





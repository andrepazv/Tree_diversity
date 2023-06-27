
#####################################################################################
#########MAPPING OBSERVED METRICS AND RESULTS FROM NULL MODELS AND SIMULATIONS #######
##############################June 2023###############################################
#####################################################################################
####Mapping dimensions of diversity based on results from index calculation, simulations and bootstrapping

library(raster)
setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Manuscript/NatCommsRev/")
###Create empty rasters based on richness map
mocklayer <- raster("SR_all_trees_observed.tif")
mocklayer<-init(mocklayer,"cell")
names(mocklayer) <- "Grid"
r<-mocklayer
values(r)<-NA

############################
######OBERVED DATA##########
###########################

all_metrics<-read.csv("REV_obs_results_200.csv")  #All observed results



###Create all metric rasters for observed
pd_ras<-r
mpd_ras<-r
raoq_ras<-r
SR_ras<-r
fdr_ras<-r
fd_ras<-r


pd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$pd
mpd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$mpd
raoq_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$raoq
SR_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$n
fdr_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$fdr
fd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$fd

writeRaster(pd_ras,"PD_all_trees.tif",format="GTiff")
writeRaster(mpd_ras,"MPD_all_trees.tif",format="GTiff")
writeRaster(SR_ras,"SR_all_trees.tif",format="GTiff")
writeRaster(raoq_ras,"RaoQ_all_trees.tif",format="GTiff")
writeRaster(fdr_ras,"fdr_all_trees.tif",format="GTiff")
writeRaster(fd_ras,"fd_all_trees.tif",format="GTiff")
plot(raoq_ras)
plot(mpd_ras)
############################
######BOOTSTRAP DATA (mean and CV)##########
###########################
all_metrics<-read.csv("REV_downsample_gatti_200.csv") 

pd_rasB<-r
mpd_rasB<-r
raoq_rasB<-r
SR_rasB<-r  
fdr_rasB<-r
fd_rasB<-r

cv_pd_rasB<-r
cv_mpd_rasB<-r
cv_raoq_rasB<-r
cv_fdr_rasB<-r
cv_fd_rasB<-r

boot<-subset(all_metrics,metric=="pd")
pd_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_pd_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean
boot<-subset(all_metrics,metric=="mpd")
mpd_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_mpd_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean
boot<-subset(all_metrics,metric=="raoq")
raoq_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_raoq_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean
boot<-subset(all_metrics,metric=="fdr")
fdr_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_fdr_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean
boot<-subset(all_metrics,metric=="fd")
fd_rasB[as.numeric(boot$grid_id)]<-boot$mean
cv_fd_rasB[as.numeric(boot$grid_id)]<-boot$sd/boot$mean

boot<-subset(all_metrics,metric=="n")
SR_rasB[as.numeric(boot$grid_id)]<-boot$mean

writeRaster(pd_rasB,"Bootstrap_PD_all_trees.tif",format="GTiff")
writeRaster(mpd_rasB,"Bootstrap_MPD_all_trees.tif",format="GTiff")
writeRaster(SR_rasB,"Bootstrap_SR_all_trees.tif",format="GTiff")
writeRaster(raoq_rasB,"Bootstrap_RaoQ_all_trees.tif",format="GTiff")
writeRaster(fdr_rasB,"Bootstrap_fdr_all_trees.tif",format="GTiff")
writeRaster(fd_rasB,"Bootstrap_fd_all_trees.tif",format="GTiff")

writeRaster(cv_pd_rasB,"Bootstrap_CV_PD_all_trees.tif",format="GTiff")
writeRaster(cv_mpd_rasB,"Bootstrap_CV_MPD_all_trees.tif",format="GTiff")
#writeRaster(cv_SR_rasB,"Bootstrap_CV_SR_all_trees.tif",format="GTiff")
writeRaster(cv_raoq_rasB,"Bootstrap_CV_RaoQ_all_trees.tif",format="GTiff")
writeRaster(cv_fdr_rasB,"Bootstrap_CV_fdr_all_trees.tif",format="GTiff")
writeRaster(cv_fd_rasB,"Bootstrap_CV_fd_all_trees.tif",format="GTiff")

plot(cv_pd_rasB)
plot(cv_fdr_rasB)
plot(cv_mpd_rasB)
plot(cv_raoq_rasB)


############################
######NULL MODEL DATA##########
###########################
observed<-read.csv("REV_obs_results_200.csv")
simulations<-read.csv("REV_null_results_200.csv") 
#Z as (obs-mean/sd)

pdSim<-subset(simulations,metric=="pd")
observed <- observed %>% left_join(pdSim, by = c("n"="n")) 
observed$Zpd <- (observed$pd-observed$mean)/observed$sd
observed <- observed %>% select(grid_id,n,pd,mpd,raoq,fd,fdr,Zpd)
mpdSim<-subset(simulations,metric=="mpd")
observed <- observed %>% left_join(mpdSim, by = c("n"="n")) 
observed$Zmpd <- (observed$mpd-observed$mean)/observed$sd
observed <- observed %>% select(grid_id,n,pd,mpd,raoq,fd,fdr,Zpd,Zmpd)

raoSim<-subset(simulations,metric=="raoq")
observed <- observed %>% left_join(raoSim, by = c("n"="n")) 
observed$Zraoq <- (observed$raoq-observed$mean)/observed$sd
observed <- observed %>% select(grid_id,n,pd,mpd,raoq,fd,fdr,Zpd,Zmpd,Zraoq)

fdrSim<-subset(simulations,metric=="fdr")
observed <- observed %>% left_join(fdrSim, by = c("n"="n")) 
observed$Zfdr <- (observed$fdr-observed$mean)/observed$sd
observed <- observed %>% select(grid_id,n,pd,mpd,raoq,fd,fdr,Zpd,Zmpd,Zraoq,Zfdr)

##write file with Z values for Random Forest
write.csv(observed,"REV_obs_results_200_zvals.csv")

Zpd<-r
Zmpd<-r
Zraoq<-r
Zfdr<-r

Zpd[as.numeric(observed$grid_id)]<-observed$Zpd
Zmpd[as.numeric(observed$grid_id)]<-observed$Zmpd
Zraoq[as.numeric(observed$grid_id)]<-observed$Zraoq
Zfdr[as.numeric(observed$grid_id)]<-observed$Zfdr

plot(Zpd)
####WriteRasters to files for maps

writeRaster(Zpd,"Zpd_all_trees_observed.tif",format="GTiff")
writeRaster(Zmpd,"Zmpd_all_trees_observed.tif",format="GTiff")
writeRaster(Zraoq,"Zraoq_all_trees_observed.tif",format="GTiff")
writeRaster(Zfdr,"Zfdr_all_trees_observed.tif",format="GTiff")

###########angiosperms only###############
all_metrics<-read.csv("REV_obs_results_ANGIO_200.csv") #Observed Angio only


###Create all metric rasters for observed
pd_ras<-r
mpd_ras<-r
raoq_ras<-r
SR_ras<-r
fdr_ras<-r
fd_ras<-r


pd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$pd
mpd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$mpd
raoq_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$raoq
SR_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$n
fdr_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$fdr
fd_ras[as.numeric(all_metrics$grid_id)]<-all_metrics$fd

writeRaster(pd_ras,"PD_angio_trees.tif",format="GTiff")
writeRaster(mpd_ras,"MPD_angio_trees.tif",format="GTiff")
writeRaster(SR_ras,"SR_angio_trees.tif",format="GTiff")
writeRaster(raoq_ras,"RaoQ_angio_trees.tif",format="GTiff")
writeRaster(fdr_ras,"fdr_angio_trees.tif",format="GTiff")
writeRaster(fd_ras,"fd_angio_trees.tif",format="GTiff")

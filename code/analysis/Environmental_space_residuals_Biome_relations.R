library(raster)
library(ggplot2)
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/")

###plotting the relation between diversity dimensions
    #Load diversity maps
SR<-raster::raster("Neo_Diversity_2degree_SR.tif")
PD<-raster::raster("Neo_Diversity_2degree_PD.tif")
FD<-raster::raster("Neo_Diversity_2degree_FRic.tif")
srpd_df<-raster::stack(PD,SR)
srpd_df<-raster::as.data.frame(srpd_df)
linear_model_sr_pd<-lm(srpd_df)
plot(srpd_df$Diversity_2degree_SR,srpd_df$Diversity_2degree_PD)
abline(linear_model_sr_pd,col="red")

srfd_df<-raster::stack(FD,SR)
srfd_df<-raster::as.data.frame(srfd_df)
linear_model_sr_fd<-lm(srfd_df)
plot(srfd_df$Diversity_2degree_SR,srfd_df$Diversity_2degree_FRic)
abline(linear_model_sr_fd,col="red")
###comparing biomes? 
wwf_eco<-rgdal::readOGR("~/Dropbox/**PostDoc_ETH/WWF_ecorregions","wwf_terr_ecos")
wwf_eco_Biome<-raster::aggregate(wwf_eco, by = "BIOME" ,dissolve = TRUE)
rgdal::writeOGR(obj=wwf_eco_Biome,dsn=getwd(),layer="WWF_Biomes.shp",driver="ESRI Shapefile")

#select the column of the attribute table that will determine the split of the shp
unique <- as.numeric(na.omit(unique(wwf_eco_Biome$BIOME)))

#create new polygons based on the determined column
for (i in 1:length(unique)) {
tmp <- wwf_eco_Biome[na.omit(wwf_eco_Biome$BIOME== unique[i]), ]
assign(paste("biome_",i,sep=""),tmp)
}
biomes<-ls(pattern="biome_")
for (i in biomes){
  biome<-get(i)
  SR_temp<-raster::crop(SR,biome)
  SR_temp<-raster::mask(SR_temp,biome)
  PD_temp<-raster::crop(PD,biome)
  PD_temp<-raster::mask(PD_temp,biome)
  FD_temp<-raster::crop(FD,biome)
  FD_temp<-raster::mask(FD_temp,biome)  
  srpd_df<-raster::stack(PD_temp,SR_temp)
  srpd_df<-raster::as.data.frame(srpd_df)
 # linear_model_sr_pd<-lm(srpd_df)
  plot(srpd_df$Diversity_2degree_SR,srpd_df$Diversity_2degree_PD,xlim=c(0,6000),ylim=c(0,85000),main=paste("biome",unique(biome$BIOME)))
#  abline(linear_model_sr_pd,col="red")
  lines(loess.smooth(srpd_df$Diversity_2degree_SR,srpd_df$Diversity_2degree_PD), col = "blue")
  srfd_df<-raster::stack(FD_temp,SR_temp)
  srfd_df<-raster::as.data.frame(srfd_df)
#  linear_model_sr_fd<-lm(srfd_df)
  plot(srfd_df$Diversity_2degree_SR,srfd_df$Diversity_2degree_FRic,xlim=c(0,6000),ylim=c(0,25000),main=paste("biome",unique(biome$BIOME)))
  lines(loess.smooth(srfd_df$Diversity_2degree_SR,srfd_df$Diversity_2degree_FRic), col = "blue")
  
  #  abline(linear_model_sr_fd,col="red")
}


###First load or cpmpute residuals maps
###PD res
res_PD<-raster("Neo_Diversity_2degree_SRPDres.tif")
###FD res must build

FD<-raster::raster("Neo_Diversity_2degree_FRic.tif")
SR<-raster::raster("Neo_Diversity_2degree_SR.tif")
srfd_df<-raster::stack(FD,SR)
srfd_df<-raster::as.data.frame(srfd_df)
linear_model_sr_fd<-lm(srfd_df)
linear_residuals_sr_fd<-resid(linear_model_sr_fd)
srfd_df$FD_linear_residuals<-rep(NA,length(srfd_df[,1]))
srfd_df$FD_linear_residuals[as.numeric(names(linear_residuals_sr_fd))]<-linear_residuals_sr_fd[names(linear_residuals_sr_fd)]
linear_residuals_sr_fd_raster<-FD
values(linear_residuals_sr_fd_raster)<-NA
values(linear_residuals_sr_fd_raster)<- srfd_df$FD_linear_residuals
writeRaster(linear_residuals_sr_fd_raster,"Neo_Diversity_2degree_SRFDres.tif",driver="GTiff")
res_FD<-linear_residuals_sr_fd_raster
###reclass the residuals as "+" , "-", "neutral"
m <- c(min(values(res_FD),na.rm=T)-1, -0.5, -1,  -0.5, 0.5, 0,  0.5, max(values(res_FD),na.rm=T)+1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
res_FD_reclass<-reclassify(res_FD,rclmat)
pol_neg <- rasterToPolygons(res_FD_reclass, fun=function(x){x<0},dissolve=T)
pol_pos <- rasterToPolygons(res_FD_reclass, fun=function(x){x>0},dissolve=T)
##polygonize as + -
###Load environmental variables
all_bioclims<-list.files("~/Dropbox//wc2/",full.names = T)
bioclims<-raster::stack(all_bioclims)
#bioclims<-raster::as.data.frame(bioclims)
#bioclims<-na.omit(bioclims)
#pca_bio<-vegan::rda(bioclims,scale=T) 
#plot(pca_bio$CA$u[,1],pca_bio$CA$u[,2])


areas_list<-c("pol_neg","pol_pos")
area_env<-as.data.frame(setNames(replicate(20,numeric(0), simplify = F), letters[1:20]))
names(area_env)<-c(names(bioclims),"area")
area_env<-as.data.frame(area_env)
for (i in areas_list){
  area_dist<-get(i)
  predictors_sp<-crop(bioclims,area_dist)
  predictors_sp<-mask(predictors_sp,area_dist)
  plot(predictors_sp$wc2.1_10m_bio_1,main=i)
  predictors_sp1<-as.data.frame(predictors_sp,xy=T)
  predictors_sp1$area<-i
  area_env<-rbind(area_env, predictors_sp1)
}
area_env1<-na.omit(area_env)
#species_env1$species[species_env1$species=="dominio.shp"]<-"T_AF"
area_env1$area<-as.factor(area_env1$area)
area_env1 <- area_env1[order(area_env1$area),] 


log.ir <- (area_env1[, 3:21])
ir.area <- area_env1[, 22]

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
ir.pca <- prcomp(log.ir,
                 center = TRUE,
                 scale. = TRUE) 

g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 2,var.axes=T,
              group = ir.area, ellipse = TRUE, 
              circle = F,alpha=0.05)
#g<-g+ scale_color_manual(name="",values=c("red","orange", "purple", "green","blue"))   
#g<- g+ scale_shape_manual(name="",values=rep(21,5)) 
#g<-g+  geom_point(aes(colour=ir.species, shape=ir.species), size = 3) 
#g <- g + scale_color_discrete(name = '')
g<- g + scale_color_manual(breaks=c("pol_pos","pol_neg"),values=c("red", "blue"))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g<-g+theme_bw()
print(g)

species_env1_PC<-cbind(area_env1,ir.pca$x[,1:3])

###Can we a) break into more categories
      ### b) create some sort of bootstrap on the map and regression to get a significance of the pixel?
                  ###null model for AF could be a way
      ###Look at whether the regression should be lineal specially in the FD case although plot PD also en do per biome


### traits in region (PCA? cluster ? )
setwd("/Users/andreapaz/Dropbox/**PostDoc_ETH/Trees_PD_FD")

###a-load rasters of species distributions and stack them 
list.files("/Alpha_clean/rasters")
Stack<-raster::stack(all_maps)
species_names<-names(Stack)
###b-Create a trait PCA and plot the two axes with no points

  traits<-arrow::read_feather("ALL_BGCI_species_and_genera_imputed_traits_8_21.feather")
  traits1<-matrix(nrow=length(unique(traits$accepted_bin)),ncol=length(unique(traits$trait_name)),dimnames=list(unique(traits$accepted_bin),unique(traits$trait_name)))   


  for (i in unique(traits$trait_name)){
    test<-subset(traits,trait_name==i)
    traits1[,i]<-test$value
    }
  rownames(traits1)<-gsub(" ","_",rownames(traits1))
  colnames(traits1)<-gsub(" ","_",colnames(traits1))
  traits1<-as.data.frame(traits1)
  traits1<-traits1[,c(1,2,7,8,28,20,26)]
  traits2<-traits1[which(rownames(traits1)%in%all_species),] ##all species should be the names in stack
  ##update species list to match intersection of distribution and traits
  all_species<-rownames(traits2)
  all_species<- all_species[order(all_species)]

  ###Do PCA 
  pca<-vegan::rda(traits2,scale=T) 
  traits3<-as.data.frame(pca$CA$u)
  ###Create the main plot
  plot(traits3$PC1,traits3$PC2,type="n")

###c-Intersect; the rasters to individual regions
  ###Loading WWF info
  wwf_eco<-rgdal::readOGR("~/Dropbox/**PostDoc_ETH/WWF_ecorregions","wwf_terr_ecos")
  wwf_eco_Biome<-raster::aggregate(wwf_eco, by = "BIOME" ,dissolve = TRUE)
  rgdal::writeOGR(obj=wwf_eco_Biome,dsn=getwd(),layer="WWF_Biomes.shp",driver="ESRI Shapefile")

  ###Divide the map into the different biomes
    #select the column of the attribute table that will determine the split of the shp
      unique <- as.numeric(na.omit(unique(wwf_eco_Biome$BIOME)))

    #create new polygons based on the determined column
   for (i in 1:length(unique)) {
    tmp <- wwf_eco_Biome[na.omit(wwf_eco_Biome$BIOME== unique[i]), ]
    assign(paste("biome_",i,sep=""),tmp)
   }
      ##Get all the different polygon names in list
    biomes<-ls(pattern="biome_")
  
    #Generate de intersection between the biomes and communities
    for (i in biomes){
    biome<-get(i)

    Stack_biome<-crop(Stack,biome)
    ##extract trait value for the PC
    traits_biome<-traits3[which(rownames(traits3)%in%names(Stack_biome)),]
    ###d-Plot in x/y (pc1/PC2 space) all values for the region in a color and repeat.
    points(traits_biome$PC1,traits_biome$PC2,col=rainbow(50)[(which(biomes==i))])
    print(paste("biome",i,"is colour",rainbow(which(biomes==i),sep=" ")))
     }


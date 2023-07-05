library(raster)
library(ggplot2)
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot) ### maybe not needed
setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Manuscript/NatCommsRev/")

###plotting the relation between diversity dimensions
    #Load diversity maps
SR<-raster::raster("SR_all_trees.tif")
PD<-raster::raster("PD_all_trees.tif")
FD<-raster::raster("fdr_all_trees.tif")
mpd<-raster("MPD_all_trees.tif")
rao<-raster("RaoQ_all_trees.tif")
    #rename layers for further processing
names(SR) <- "SR"
names(PD) <- "PD"
names(FD) <- "FD"
names(mpd) <- "MPD"
names(rao) <- "RaoQ"
    ##create pairwise stacks for regressions and remove the NA values
srpd_df<-raster::stack(PD,SR)
srpd_df<-raster::as.data.frame(srpd_df)
srpd_df<-na.omit(srpd_df)
srfd_df<-raster::stack(FD,SR)
srfd_df<-raster::as.data.frame(srfd_df)
srfd_df<-na.omit(srfd_df)
fdpd_df<-raster::stack(FD,PD)
fdpd_df<-raster::as.data.frame(fdpd_df)
fdpd_df<-na.omit(fdpd_df)
mpdrao_df<-stack(rao,mpd)
mpdrao_df<-as.data.frame(mpdrao_df)
mpdrao_df<-na.omit(mpdrao_df)
pdmpd_df<-stack(PD,mpd)
pdmpd_df<-as.data.frame(pdmpd_df)
pdmpd_df<-na.omit(pdmpd_df)
fdrao_df<-stack(FD,rao)
fdrao_df<-as.data.frame(fdrao_df)
fdrao_df<-na.omit(fdrao_df)

    ##use this to check the data 
#linear_model_sr_pd<-lm(PD~0+SR,data=srpd_df)
#plot(srpd_df$SR,srpd_df$PD,ylab="PD",xlab="Richness")
#abline(linear_model_sr_pd,col="red")

#linear_model_sr_fd<-lm(FD~0+SR,data=srfd_df)
#plot(srfd_df$SR,srfd_df$FD,xlab="Richness",ylab="FRic")
#abline(linear_model_sr_fd,col="red")

#linear_model_fdpd<-lm(fdpd_df)
#plot(fdpd_df$PD,fdpd_df$FR,xlab="PD",ylab="FRic")
#abline(linear_model_fdpd,col="red")


###To compare biomes first load all biome shapefile
wwf_eco<-rgdal::readOGR("~/Dropbox/**PostDoc_ETH/WWF_ecorregions","wwf_terr_ecos")
wwf_eco_Biome<-raster::aggregate(wwf_eco, by = "BIOME" ,dissolve = TRUE)
rgdal::writeOGR(obj=wwf_eco_Biome,dsn=getwd(),layer="WWF_Biomes.shp",driver="ESRI Shapefile")

#select the column of the attribute table that will determine the split of the shp
unique <- as.numeric(na.omit(unique(wwf_eco_Biome$BIOME)))

#create new polygons based on the determined column
for (i in 1:length(unique)) {
tmp <- wwf_eco_Biome[na.omit(wwf_eco_Biome$BIOME== unique[i]), ]
assign(paste("biome_",unique[i],sep=""),tmp)
}
biomes<-ls(pattern="biome_")
##remove lakes and ice
biomes<-biomes[1:14]


legend<-c("RF","Montane grassland","Tundra","Mediterranean","Deserts",
          "Mangrove","Dry F", "Trop. Coniferous","Temperate mixed",
          "Coniferous", "Boreal", "Trop. Savannas", "Savannas", 
          "Flooded grassland")
selected_biomes<-c(1,9,2,4,7,5,10,12,13,11)
legend<-c("RF","Temperate mixed","Montane grassland","Mediterranean",
         "Dry F","Deserts",
          "Coniferous", "Trop. Savannas", "Savannas","Boreal")
colors_hand1<-c("darkgreen","darkolivegreen3","springgreen3","orange","darkorange2","khaki4","dodgerblue1","khaki","gold2","darkturquoise") 
#colors_hand<-c("darkgreen","darkorange","darkturquoise","darkolivegreen3","gold2","dodgerblue1") #1,5,7,11,10,9
#legend_hand<-c("Rainforest","Dry Forest","Coniferous","Temperate boradleaf and mixed","Desert","Boreal")

#selected_biomes<-c(1,7,10,9,5,11)

#####to get plots of "stars" of relation between measures creates one at a time ######
    ##This is PD vs. FR
    plot(xlim=c(0,27000),ylim=c(0,3.5),fdpd_df$PD,fdpd_df$FD,xlab="PD",ylab="FRic",type="n")
    ###add an average global value 
    points(mean(fdpd_df$PD),mean(fdpd_df$FD),col = "black")
    arrows(mean(fdpd_df$PD)- sd(fdpd_df$PD),  mean(fdpd_df$FD), mean(fdpd_df$PD)+sd(fdpd_df$PD),  mean(fdpd_df$FD), length=0.05, angle=90, code=3,col = "black",lty=2)
    arrows(mean(fdpd_df$PD),  mean(fdpd_df$FD)- sd(fdpd_df$FD), mean(fdpd_df$PD), mean(fdpd_df$FD)+sd(fdpd_df$FD), length=0.05, angle=90, code=3,col = "black",lty=2)
 

    for (i in 1:length(selected_biomes)) {
  biome<-get(biomes[selected_biomes[i]])
    PD_temp<-raster::crop(PD,biome)
    PD_temp<-raster::mask(PD_temp,biome)
    FD_temp<-raster::crop(FD,biome)
    FD_temp<-raster::mask(FD_temp,biome)
    temp_df<-raster::stack(FD_temp,PD_temp)
    temp_df<-raster::as.data.frame(temp_df)
    temp_df<-na.omit(temp_df)
    colnames(temp_df)<-c("y_temp","x_temp")
      sdx = sd(temp_df$x_temp)
      xmean = mean(temp_df$x_temp)
      sdy = sd(temp_df$y_temp)
      ymean = mean(temp_df$y_temp)
    points(xmean,ymean,col = colors_hand1[i])
    arrows(xmean-sdx, ymean, xmean+sdx, ymean, length=0.05, angle=90, code=3,col = colors_hand1[i])
    arrows(xmean, ymean-sdy, xmean, ymean+sdy, length=0.05, angle=90, code=3,col = colors_hand1[i])
}
    ##ADD legend if needed [for manuscript it was added manually to figures]
    #legend(15000, 0.007, legend=legend1,col=colors_hand1, cex=0.65,pch=20)
 
    ##SR vs. PD
    plot(xlim=c(0,1100),ylim=c(0,25000),srpd_df$SR,srpd_df$PD,ylab="PD",xlab="Richness",type="n")
    ###add an average global value 
    points(mean(srpd_df$SR),mean(srpd_df$PD),col = "black")
    arrows(mean(srpd_df$SR)- sd(srpd_df$SR),  mean(srpd_df$PD), mean(srpd_df$SR)+sd(srpd_df$SR),  mean(srpd_df$PD), length=0.05, angle=90, code=3,col = "black",lty=2)
    arrows(mean(srpd_df$SR),  mean(srpd_df$PD)- sd(srpd_df$PD), mean(srpd_df$SR), mean(srpd_df$PD)+sd(srpd_df$PD), length=0.05, angle=90, code=3,col = "black",lty=2)
    for (i in 1:length(selected_biomes)) {
      biome<-get(biomes[selected_biomes[i]])
      SR_temp<-raster::crop(SR,biome)
      SR_temp<-raster::mask(SR_temp,biome)
      PD_temp<-raster::crop(PD,biome)
      PD_temp<-raster::mask(PD_temp,biome)
      temp_df<-raster::stack(PD_temp,SR_temp)
      temp_df<-raster::as.data.frame(temp_df)
      temp_df<-na.omit(temp_df)
      colnames(temp_df)<-c("y_temp","x_temp")
      sdx = sd(temp_df$x_temp)
      xmean = mean(temp_df$x_temp)
      sdy = sd(temp_df$y_temp)
      ymean = mean(temp_df$y_temp)
      points(xmean,ymean,col = colors_hand1[i])
      arrows(xmean-sdx, ymean, xmean+sdx, ymean, length=0.05, angle=90, code=3,col = colors_hand1[i])
      arrows(xmean, ymean-sdy, xmean, ymean+sdy, length=0.05, angle=90, code=3,col = colors_hand1[i])
    } 
    ##ADD legend if needed [for manuscript it was added manually to figures]
    #legend(760, 13000, legend=legend1,col=colors_hand1, cex=0.6,pch=1)

    ##SR vs. FR
    plot(xlim=c(0,1100),ylim=c(0,3.5),srfd_df$SR,srfd_df$FD,ylab="FRic",xlab="Richness",type="n")
    ##add an average global value 
    points(mean(srfd_df$SR),mean(srfd_df$FD),col = "black")
    arrows(mean(srfd_df$SR)- sd(srfd_df$SR),  mean(srfd_df$FD), mean(srfd_df$SR)+sd(srfd_df$SR),  mean(srfd_df$FD), length=0.05, angle=90, code=3,col = "black",lty=2)
    arrows(mean(srfd_df$SR),  mean(srfd_df$FD)- sd(srfd_df$FD), mean(srfd_df$SR), mean(srfd_df$FD)+sd(srfd_df$FD), length=0.05, angle=90, code=3,col = "black",lty=2)
    for (i in 1:length(selected_biomes)) {
      biome<-get(biomes[selected_biomes[i]])
      SR_temp<-raster::crop(SR,biome)
      SR_temp<-raster::mask(SR_temp,biome)
      FD_temp<-raster::crop(FD,biome)
      FD_temp<-raster::mask(FD_temp,biome)
      temp_df<-raster::stack(FD_temp,SR_temp)
      temp_df<-raster::as.data.frame(temp_df)
      temp_df<-na.omit(temp_df)
      colnames(temp_df)<-c("y_temp","x_temp")
      sdx = sd(temp_df$x_temp)
      xmean = mean(temp_df$x_temp)
      sdy = sd(temp_df$y_temp)
      ymean = mean(temp_df$y_temp)
      points(xmean,ymean,col = colors_hand1[i])
      arrows(xmean-sdx, ymean, xmean+sdx, ymean, length=0.05, angle=90, code=3,col = colors_hand1[i])
      arrows(xmean, ymean-sdy, xmean, ymean+sdy, length=0.05, angle=90, code=3,col = colors_hand1[i])
    }   
    ##ADD legend if needed [for manuscript it was added manually to figures]
    #legend(700, 0.007, legend=legend1,col=colors_hand1, cex=0.65,pch=20)
  
    ##PD VS MPD  pdmpd_df
    plot(xlim=c(0,25000),ylim=c(200,450),pdmpd_df$PD,pdmpd_df$MPD,xlab="Phylogenetic richness",ylab="Phylogenetic divergence (MPD)",type="n")
    ###add an average global value 
    points(mean(pdmpd_df$PD),mean(pdmpd_df$MPD),col = "black")
    arrows(mean(pdmpd_df$PD)- sd(pdmpd_df$PD),  mean(pdmpd_df$MPD), mean(pdmpd_df$PD)+sd(pdmpd_df$PD),  mean(pdmpd_df$MPD), length=0.05, angle=90, code=3,col = "black",lty=2)
    arrows(mean(pdmpd_df$PD),  mean(pdmpd_df$MPD)- sd(pdmpd_df$MPD), mean(pdmpd_df$PD), mean(pdmpd_df$MPD)+sd(pdmpd_df$MPD), length=0.05, angle=90, code=3,col = "black",lty=2)
    for (i in 1:length(selected_biomes)) {
      biome<-get(biomes[selected_biomes[i]])
      PD_temp<-raster::crop(PD,biome)
      PD_temp<-raster::mask(PD_temp,biome)
      mpd_temp<-raster::crop(mpd,biome)
      mpd_temp<-raster::mask(mpd_temp,biome)
      temp_df<-raster::stack(mpd_temp,PD_temp)
      temp_df<-raster::as.data.frame(temp_df)
      temp_df<-na.omit(temp_df)
      colnames(temp_df)<-c("y_temp","x_temp")
      sdx = sd(temp_df$x_temp)
      xmean = mean(temp_df$x_temp)
      sdy = sd(temp_df$y_temp)
      ymean = mean(temp_df$y_temp)
      points(xmean,ymean,col = colors_hand1[i])
      arrows(xmean-sdx, ymean, xmean+sdx, ymean, length=0.05, angle=90, code=3,col = colors_hand1[i])
      arrows(xmean, ymean-sdy, xmean, ymean+sdy, length=0.05, angle=90, code=3,col = colors_hand1[i])
      
    }

    ##FR vs. RaoQ  fdrao_df<-stack(FD,rao)
    plot(xlim=c(0.5,3.5),ylim=c(12,25),fdrao_df$FD,fdrao_df$RaoQ,xlab="Functional Richness",ylab="Functional divergence (RaoQ)",type="n")
    ###add an average global value 
    points(mean(fdrao_df$FD),mean(fdrao_df$RaoQ),col = "black")
    arrows(mean(fdrao_df$FD)- sd(fdrao_df$FD),  mean(fdrao_df$RaoQ), mean(fdrao_df$FD)+sd(fdrao_df$FD),  mean(fdrao_df$RaoQ), length=0.05, angle=90, code=3,col = "black",lty=2)
    arrows(mean(fdrao_df$FD),  mean(fdrao_df$RaoQ)- sd(fdrao_df$RaoQ), mean(fdrao_df$FD), mean(fdrao_df$RaoQ)+sd(fdrao_df$RaoQ), length=0.05, angle=90, code=3,col = "black",lty=2)
    for (i in 1:length(selected_biomes)) {
      biome<-get(biomes[selected_biomes[i]])
      FD_temp<-raster::crop(FD,biome)
      FD_temp<-raster::mask(FD_temp,biome)
      rao_temp<-raster::crop(rao,biome)
      rao_temp<-raster::mask(rao_temp,biome)
      temp_df<-raster::stack(rao_temp,FD_temp)
      temp_df<-raster::as.data.frame(temp_df)
      temp_df<-na.omit(temp_df)
      colnames(temp_df)<-c("y_temp","x_temp")
      sdx = sd(temp_df$x_temp)
      xmean = mean(temp_df$x_temp)
      sdy = sd(temp_df$y_temp)
      ymean = mean(temp_df$y_temp)
      points(xmean,ymean,col = colors_hand1[i])
      arrows(xmean-sdx, ymean, xmean+sdx, ymean, length=0.05, angle=90, code=3,col = colors_hand1[i])
      arrows(xmean, ymean-sdy, xmean, ymean+sdy, length=0.05, angle=90, code=3,col = colors_hand1[i])
    }  
    
  
######to get lines and points coloured by biome########
    #for PD vs FR
    ##get the average trend  
      plot(fdpd_df$PD,fdpd_df$FD,xlab="PD",ylab="FRic",type="n",ylim=c(0,4),xlim=c(0,40000))
      loess_points<-loess.smooth(fdpd_df$PD,fdpd_df$FD,evaluation=100,span=1)
      lines(loess_points$x, loess_points$y,lwd=3,lty=2)
      ##add per biome trend to graph
      for (i in 1:length(selected_biomes)) {
      biome<-get(biomes[selected_biomes[i]])
      PD_temp<-raster::crop(PD,biome)
      PD_temp<-raster::mask(PD_temp,biome)
      FD_temp<-raster::crop(FD,biome)
      FD_temp<-raster::mask(FD_temp,biome) 
  
 # rao_temp<-raster::crop(rao,biome)
#  rao_temp<-raster::mask(rao_temp,biome)
 # mpd_temp<-raster::crop(mpd,biome)
#  mpd_temp<-raster::mask(mpd_temp,biome)
 # mpdrao_df<-stack(rao_temp,mpd_temp)
#  mpdrao_df<-as.data.frame(mpdrao_df)
#  points(mpdrao_df$Trees_nomono_matchPDFD_30_100_PD_mpd,mpdrao_df$Functional_raosq_polygon_based_matchedPDFD_30_100_,col = colors_hand[i])
  
 # srpd_df<-raster::stack(PD_temp,SR_temp)
 #srpd_df<-raster::as.data.frame(srpd_df)

  #points(srpd_df$Trees_nomono_30_100_PD_Richness,srpd_df$Trees_nomono_30_100_PD_PD,col = rainbow(24)[i])
  ##get 99% of points
#total_points<-ceiling(length(which(!is.na(srpd_df$Trees_nomono_30_100_PD_PD)))*0.99)
# max_val<-max(sort(srpd_df$Trees_nomono_30_100_PD_PD,na.last=NA)[1:total_points])
 # srpd_df<-srpd_df[which(srpd_df$Trees_nomono_30_100_PD_PD<=max_val),]
  ##get loess for only those 99%
#  loess_points<-loess.smooth(srpd_df$Trees_nomono_30_100_PD_Richness,srpd_df$Trees_nomono_30_100_PD_PD,evaluation=100,span=1)
#  lines(loess_points$x, loess_points$y, col = colors_hand[i],lwd=3)
  #  abline(linear_model_sr_pd,col="red")
#  srfd_df<-raster::stack(FD_temp,SR_temp)
#  srfd_df<-raster::as.data.frame(srfd_df)
 #linear_model_sr_fd<-lm(srfd_df)

##get 99% of points
# total_points<-ceiling(length(which(!is.na(srfd_df$Functional_richness_scaled_1_8_function_polygon_based_matchedPDFD_30_100_)))*0.99)
# max_val<-max(sort(srfd_df$Functional_richness_scaled_1_8_function_polygon_based_matchedPDFD_30_100_,na.last=NA)[1:total_points])
# srfd_df<-srfd_df[which(srfd_df$Functional_richness_scaled_1_8_function_polygon_based_matchedPDFD_30_100_<=max_val),]
 ##get loeass for only those 99%
#loess_points<-loess.smooth(srfd_df$Trees_nomono_30_100_PD_Richness,srfd_df$Functional_richness_scaled_1_8_function_polygon_based_matchedPDFD_30_100_,evaluation=100,span=1)
#  lines(loess_points$x, loess_points$y, col = colors_hand[i],lwd=3)
 
 # points(srfd_df$Trees_nomono_30_100_PD_Richness,srfd_df$Functional_richness_scaled_1_8_function_polygon_based_matchedPDFD_30_100_,col = colors_hand[i])
  
  #  abline(linear_model_sr_fd,col="red")
  ###for dispersion mesaures
  
  ###For FD/PD
   fdpd_df1<-raster::stack(FD_temp,PD_temp)
  fdpd_df1<-raster::as.data.frame(fdpd_df1)
  ##get 99% of points
   total_points<-ceiling(length(which(!is.na(fdpd_df$FD)))*0.99)
  max_val<-max(sort(fdpd_df$FD,na.last=NA)[1:total_points])
  fdpd_df1<-fdpd_df1[which(fdpd_df1$FD<=max_val),]
  ##get loess for only those 99%
  loess_points<-loess.smooth(fdpd_df1$PD,fdpd_df1$FD,evaluation=100,span=1)
    lines(loess_points$x, loess_points$y, col = colors_hand1[i],lwd=3)
  
}


      #for SD vs FR
      ##get the average trend  
      plot(srfd_df$SR,srfd_df$FD,xlab="Species richness",ylab="Functional richness",type="n",ylim=c(0,4))
      loess_points<-loess.smooth(srfd_df$SR,srfd_df$FD,evaluation=100,span=1)
      lines(loess_points$x, loess_points$y,lwd=3,lty=2)
      ##add per biome trend to graph
      for (i in 1:length(selected_biomes)) {
        biome<-get(biomes[selected_biomes[i]])
        SR_temp<-raster::crop(SR,biome)
        SR_temp<-raster::mask(SR_temp,biome)
        FD_temp<-raster::crop(FD,biome)
        FD_temp<-raster::mask(FD_temp,biome) 
       
        srfd_df1<-raster::stack(FD_temp,SR_temp)
        srfd_df1<-raster::as.data.frame(srfd_df1)
        ##get 99% of points
        total_points<-ceiling(length(which(!is.na(srfd_df1$FD)))*0.99)
        max_val<-max(sort(srfd_df1$FD,na.last=NA)[1:total_points])
        srfd_df1<-srfd_df1[which(srfd_df1$FD<=max_val),]
        ##get loess for only those 99%
        loess_points<-loess.smooth(srfd_df1$SR,srfd_df1$FD,evaluation=100,span=1)
        lines(loess_points$x, loess_points$y, col = colors_hand1[i],lwd=3)
        
      }
    
    
      #for SR vs PD
      ##get the average trend  
      plot(srpd_df$SR,srpd_df$PD,xlab="Species richness",ylab="Phylogenetic richness",type="n")
      loess_points<-loess.smooth(srpd_df$SR,srpd_df$PD,evaluation=100,span=1)
      lines(loess_points$x, loess_points$y,lwd=3,lty=2)
      ##add per biome trend to graph
      for (i in 1:length(selected_biomes)) {
        biome<-get(biomes[selected_biomes[i]])
        PD_temp<-raster::crop(PD,biome)
        PD_temp<-raster::mask(PD_temp,biome)
        SR_temp<-raster::crop(SR,biome)
        SR_temp<-raster::mask(SR_temp,biome) 
        
        srpd_df1<-raster::stack(PD_temp,SR_temp)
        srpd_df1<-raster::as.data.frame(srpd_df1)
        ##get 99% of points
        total_points<-ceiling(length(which(!is.na(srpd_df1$PD)))*0.99)
        max_val<-max(sort(srpd_df1$PD,na.last=NA)[1:total_points])
        srpd_df1<-srpd_df1[which(srpd_df1$PD<=max_val),]
        ##get loess for only those 99%
        loess_points<-loess.smooth(srpd_df1$SR,srpd_df1$PD,evaluation=100,span=1)
        lines(loess_points$x, loess_points$y, col = colors_hand1[i],lwd=3)
      }
      
      
  

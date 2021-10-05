### Using Plantago lanceolata as an example of disjunct distribution, same with the White pine
  ###load csv
#points<-read.csv("Dropbox/**PostDoc_ETH/Trees_PD_FD/Alpha_hulls/Plantago_lanceolata_processed_occs.csv")
#i<-"Plantago lanceolata L."
#points<-read.csv("Dropbox/**PostDoc_ETH/Trees_PD_FD/Alpha_hulls/Pinus_strobus_processed_occs.csv")
#  i<-"Pinus strobus L."

#red maple

setwd("/Users/andreapaz/Dropbox/**PostDoc_ETH/Trees_PD_FD")
dir<-getwd()
points<-read.csv("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Alpha_hulls/all_trees_alldb.csv")
##create dataframe to store the alphas for each species
species_names<-unique(points$Species)
alpha_used<-data.frame(NA,NA)
numberused<-data.frame(NA,NA)
colnames(alpha_used)<-c("species","alpha")
colnames(numberused)<-c("species","numberPoly")
for (i in species_names){
  loc<-subset(points,Species==i)
  loc<-loc[!duplicated(loc), ]
  tryCatch({
    ###Add condition of more than 3 points else create point based raster. if
    if(length(loc$species)<=3){
      coordinates(loc)<-~Longitude+Latitude 
      writeOGR(loc,dir,paste("map",i,sep="_"),driver="ESRI Shapefile")
    }
    else
    { alphahull<-getDynamicAlphaHull(loc, fraction = 0.95, partCount = 10, buff = 10000,
                                     initialAlpha = 1, coordHeaders = c('Longitude', 'Latitude'),
                                     clipToCoast = 'terrestrial', proj = "+proj=longlat +datum=WGS84",
                                     alphaIncrement = 1, verbose = FALSE)
    map_poly<-as(alphahull[[1]],"SpatialPolygonsDataFrame")
    alpha_sp<-c(i,alphahull$alpha)
    names(alpha_sp)<-c("species","alpha")
    alpha_used<-rbind(alpha_used,alpha_sp)
    ##Get number of polygons in the map for future flagging
    number<-c(i,length(map_poly@polygons[[1]]@Polygons))
    names(number)<-c("species","numberPoly")
    numberused<-rbind(numberused,number)
    writeOGR(map_poly,dir,paste("map",i,sep="_"),driver="ESRI Shapefile")
    }
    
  },error=function(e){cat("ERROR:",i,conditionMessage(e),"\n")})
}
alpha_used<-alpha_used[2:length(alpha_used$species),]
write.csv(alpha_used,"alpha_per_species.csv")

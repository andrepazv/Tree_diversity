library(raster)

###First step would be to get grid ID for all the plots within the grid of interest 


####Create raster base 
rasterBase<-raster::raster("~/Dropbox/wc2/wc2.1_10m_bio_1.tif") ### This is a 10 min resolution layer (~18km)
##make coarser if needed 
rasterBase<-raster::aggregate(rasterBase,10)
###intersect plots with the raster base to get a list of cells with plots
###plots should be a 2 column data frame , longitude,latitude
plotID<-raster::extract(rasterBase,plotID,cellnumbers=T)
plotID<-plotID[,1]

FocalCellID<-plotID[1]###plot cell ID?
neighbourNumber<-3   ###user det
GettingNeighbours<-function(FocalCellID,neighbourNumber,r=rasterBase){
neighbours<-1:neighbourNumber

FocalCellRow<-rowFromCell(r,FocalCellID)
FocalCellCol<-colFromCell(r,FocalCellID)
neighbourCells<-list()
for(i in 1:neighbourNumber){

neighbourCells[[i]]<-c(cellFromRowCol(r,FocalCellRow-i,FocalCellCol-i),cellFromRowCol(r,FocalCellRow-i,FocalCellCol),
cellFromRowCol(r,FocalCellRow-i,FocalCellCol+i),cellFromRowCol(r,FocalCellRow,FocalCellCol-i),
cellFromRowCol(r,FocalCellRow,FocalCellCol+i),cellFromRowCol(r,FocalCellRow+i,FocalCellCol-i),
cellFromRowCol(r,FocalCellRow+i,FocalCellCol),cellFromRowCol(r,FocalCellRow+i,FocalCellCol+i))
}

#AllCells<-unlist(neighbourCells)
return(neighbourCells)

}
NeighborID<-lapply(plotID,function(x){GettingNeighbours(FocalCellID=x,neighbourNumber=3)})
names(NeighborID)<-plotID


####Get the communities of interest 
##if ran along with the Communities script then use community data else must read the matrix from disk
community_data<-read.csv("CommMatrix.csv")
###Assuming this has a column named Grid from the other script
		###This for a single focus cell then make function
CellsofInterest<-NeighborID[[plotID[1]]]
CommunitiesofInterest<-community_data[which(community_data$Grid==CellsofInterest),]
speciesPool<- names(CommunitiesofInterest)[colSums(CommunitiesofInterest)>=1]



library(raster)

###First step would be to get grid ID for all the plots within the grid of interest 

# set resolution, both for base layer and for community matrix
resolution <- 10

###Setwd to folder 
setwd("~/local_Git/Tree_diversity/data/")

plotID <- read_csv("~/Git/fia_data/data/fia_data/fia_data/OR_PLOT.csv") %>% select(LAT, LON) %>% distinct() %>% rename(Longitude = LON, Latitude = LAT) %>% select(Longitude, Latitude)
#plotID <- tibble(Longitude = c(-119, -110), Latitude = c(43.9, 35.1))


####Create raster base 
rasterBase<-raster::raster(paste0("WorldClim2/wc2.1_",resolution,"m_bio_1.tif")) ### This is a 10 min resolution layer (~18km)
##make coarser if needed 
rasterBase<-raster::aggregate(rasterBase,resolution)

###intersect plots with the raster base to get a list of cells with plots
###plots should be a 2 column data frame , longitude,latitude
plotID<-raster::extract(rasterBase,plotID,cellnumbers=T)
plotID<-plotID[,1]

FocalCellID<-plotID[1]###plot cell ID?
# neighbourNumber<-10   ###user det
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
NeighborID<-lapply(plotID,function(x){GettingNeighbours(FocalCellID=x,neighbourNumber=10)})
names(NeighborID)<-plotID


####Get the communities of interest 
##if ran along with the Communities script then use community data else must read the matrix from disk
community_data<-read.csv(paste0("Community_matrix_res_", resolution,".csv"))
###Assuming this has a column named Grid from the other script
###select number of neighborrs
n<-10
all_neighbors<-1:n
###This for a single focus cell then make function
CellsofInterest<-NeighborID[[plotID[1]]] 
##maybe CellsofInterest<-NeighborID[[as.character(plotID[1])]] 
allCells<-vector()
for(i in 1:n){
	allCells<-c(allCells,CellsofInterest[[i]])
}
allCells<-c(allCells,plotID[1])
CommunitiesofInterest<-community_data[which(community_data$Grid%in%allCells),]
speciesPool<- names(CommunitiesofInterest[-c(1:6)])[colSums(CommunitiesofInterest[-c(1:6)])>=1]



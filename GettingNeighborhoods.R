library(raster)

###First step would be to get grid ID for all the plots within the grid of interest 


####Create raster base 
rasterBase<-raster::raster("~/Dropbox/wc2/wc2.1_10m_bio_1.tif") ### This is a 10 min resolution layer (~18km)
##make coarser if needed 
rasterBase<-raster::aggregate(rasterBase,10)
###intersect plots with the raster base to get a list of cells with plots
###plots should be a 2 column data frame , longitude,latitude
plotID<-raster::extract(rasterBase,plots,cellnumbers=T)
plotID<-plotID$cells

FocalCellID<-plotID[1]###plot cell ID?
neighbourNumber<-3   ###user det
neighbours<-1:neighbourNumber

FocalCellRow<-rowFromCell(FocalCellID)
FocalCellCol<-colFromCell(FocalCellID)
neighborCells<-list()
for(i in 1:neighbourNumber){


 

 


neighborCells[[i]]<-c(cellFromRowCol(FocalCellRow-i,FocalCellCol-i),cellFromRowCol(FocalCellRow-i,FocalCellCol),
cellFromRowCol(FocalCellRow-i,FocalCellCol+i),cellFromRowCol(FocalCellRow,FocalCellCol-i),
cellFromRowCol(FocalCellRow,FocalCellCol+i),cellFromRowCol(FocalCellRow+i,FocalCellCol-i),
cellFromRowCol(FocalCellRow+i,FocalCellCol),cellFromRowCol(FocalCellRow+i,FocalCellCol+i))
}

 cellFromRowCol.
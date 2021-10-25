###October 2021###
library(rgdal)
library(raster)
###Setwd to folder containing the polygons###
setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/")
dir<-getwd()
####Some species have <3 points they are not polygons but points but this would work for both 
### function to rasterize points or polygons
rasterize_species= function (x,mask=selected_mask,r=raster_base) {
  # r<-raster(ncol=300,nrow=400,resolution=resolution,ext=extent(selected_mask))
  # res(r)<-resolution 
  r<-extend(r,selected_mask)
  r<-crop(r,extent(selected_mask))
  r<-mask(r,selected_mask)
  values(r)<-0
  map<-readOGR(dsn=getwd(),layer=x)
  r<-rasterize(map,r,1,update=T,background=0)
  r<-mask(r,selected_mask)
  values<-unique(getValues(r))
  
  if(length(values)==1&&is.na(values)==TRUE){
    
    
  }
  else if (length(values)==2&&values[2]==0){
    
  }
  else {
    writeRaster(r,paste(x,".asc",sep=""))
    return (raster(paste(x,".asc",sep="")))
  }
}

###Load the mask, for world analyses just world_adm and dissolve

selected_mask<-readOGR(dir,"world_adm0") 
selected_mask = aggregate(selected_mask, dissolve = TRUE)

###Get a raster base, can be a worldclim variable for the whole world, it will set the resolution and extent

raster_base<-raster("~/Dropbox/wc2/wc2.1_10m_bio_1.tif") ### use a 10 min WC layer or even coarser?
##make coarser this is toooo fine
raster_base<-raster::aggregate(raster_base,10)
##Get the names of species and distribution file names
setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Alpha_clean/")
distribution_files<-list.files(path=getwd(), pattern= "*.shp$")
##Get files name
distribution_maps<-sub(".shp","",distribution_files) ##really file names
##Get species names
species_names<-sub("map_","",distribution_maps)

##Load all polygons and rasterize them
layers<-lapply(distribution_maps,rasterize_species)

##Stack all maps and give correct species names
Stack_maps<-stack(layers)
names(Stack_maps)<-species_names

########MAP DIVERSITY#########
  ##Taxonomic diversity##
TD<-calc(Stack_maps,sum,na.rm=T) ##taxonomic diversity map
plot(TD)
writeRaster(TD,"richness_trees_10arcmin_agg10.asc")
  
  ##Phylogenetic diversity##
species_names<-gsub("_"," ",species_names)

# sample call, assuming unique_spp_names is in Latin binomial format
phylo_tree <- V.PhyloMaker::phylo.maker(data.frame(species = species_names, genus = stringr::word(species_names, 1), family = rep(NA, length(species_names))))
##note not all trees could be place (~3K not in the tree) 
write.tree(phylo_tree$scenario.3,file="TEST_all_trees_tree.nwk")
###create list of missing trees
species_not_tree<-subset(phylo_tree$species.list,status=="fail to bind")
write.table(species_not_tree,"Species_not_in_tree.txt")
species_names<-gsub(" ","_",species_names)
species_names<-stringr::str_to_title(species_names)

####MMake sure the names match between maps and phylogeny
setdiff(phylo_tree$scenario.3$tip.label, species_names) ##chack var names here
##if it doesnt use this code
###Trim phylogeny to match distribution data (remove non hylids)
pruned.tree<-ape::drop.tip(phylo_tree$scenario.3$tip.label, setdiff(phylo_tree$scenario.3$tip.label, species_names))
##check that all disrtribution data is in tree
test<-as.data.frame(species_names)
rownames(test)<-species_names
check_names<-geiger::name.check(pruned.tree, test, data.names=NULL) #this must be OK



#####Compute and map PD
names_table<-as.data.frame(species_names)
colnames(names_table)<-"Grid"

#Create empty raster of the desired area and resolution to assign pixel numbers
r<-Stack_maps[[1]]
res(r)<-res(Stack_maps) #resolution
#r[is.na(r[])]<-0
#r<-crop(r,selected_mask)
#r<-mask(r,selected_mask)
#r<-aggregate(r,fact=10,fun=max)
grid=r
names(grid)="grid"
grid[1:ncell(grid)]<-1:ncell(grid)
list_grid<-list(grid)
names(list_grid)<-"Grid"
#Change names of models tu species names and add the empty raster in the beggining (pixel number)
names(Stack_maps)<-as.vector(names_table$Grid)
Stack<-stack(List_grid$Grid,Stack_maps) ###Full stack of all maps
#
#Turn maps into dataframe for computation of PD
community_data<-as.data.frame(Stack)
##remove rows that are NA 
community_data1<-community_data[-which(rowSums(!is.na(community_data[,2:ncol(community_data1)]))==0),] ##change number of columns
species_names<-colnames(community_data1)[2:length(community_data1)] #Store species names


#Store 
setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Diversity_computation/")
#write.table(community_data1, file = "communities.txt", append = FALSE,row.names=F,quote=F,sep="\t") #tratar de agregar nombre de mascara
arrow::write_feather(community_data1, "communities.feather")
community_data1[is.na(community_data1)] = 0
#In the community data frame NA must be eliminated done before is this if you load community data and not maps?
#community_data=na.omit(community_data)
#head(community_data)

#Phylogenetic diversity computation 
#computes only Faith�s PD others may be added

pd.result <-picante::pd(community_data1[,2:ncol(community_data1)],pruned.tree,include.root = F) 

#Add the pixel PD value to data frame
community_data1$pd<-pd.result[[1]]

#Write the new matrix to a file to avoid rerunning the script for potential further analyses
#write.table(community_data1, file = "communities_and_pd.txt", append = FALSE,row.names=F,quote=F,sep="\t")
arrow::write_feather(community_data1, "communities_withPD.feather")
#Generate a raster containing PD information per pixel

#1-First generate an empty raster using a base model for resolution and area

values(r)<-0
pd_ras<-r
values(pd_ras)<-NA #Eliminate every value they will be replaced by PD values further down


#2- Assign PD values to raster
pd_ras[community_data1$grid]<-community_data1$pd

#3- Save raster to file 

writeRaster(pd_ras,"Trees_PD.tif",format="GTiff")

#4-Optional plotting map in R 
plot(pd_ras)


#################################################################
#######Functional diversity computation here####################
################################################################
##Read Traits
traits<-arrow::read_feather("ALL_BGCI_species_and_genera_imputed_traits_8_21.feather")
traits1<-matrix(nrow=length(unique(traits$accepted_bin)),ncol=length(unique(traits$trait_name)),dimnames=list(unique(traits$accepted_bin),unique(traits$trait_name)))   
##re-check this part why ?
for (i in unique(traits$trait_name)){
  test<-subset(traits,trait_name==i)
  traits1[,i]<-test$value
}
rownames(traits1)<-gsub(" ","_",rownames(traits1))
colnames(traits1)<-gsub(" ","_",colnames(traits1))
traits1<-as.data.frame(traits1)
traits1<-traits1[,c(1,2,7,8,28,20,26)] ##seelction from Dan's paper
 #remove species not in traits from stack
species_names<-names(Stack_maps)
Stack_maps1<-subset(Stack_maps,intersect(species_names,rownames(traits1)))
species_names<-names(Stack_maps1)
species_names<-species_names[order(species_names)]
##create df of ordered species names 
names_table<-as.data.frame(species_names)
colnames(names_table)<-"Grid"
#Create empty raster of the desired area and resolution to assign pixel numbers
r<-Stack_maps1[[1]]
res(r)<-res(Stack_maps1)
grid=r
names(grid)="grid"
grid[1:ncell(grid)]<-1:ncell(grid)
list_grid<-list(grid)
names(list_grid)<-"Grid"
#Change names of models tu species names and add the empty raster in the beggining (pixel number)
names(Stack_maps1)<-as.vector(names_table$Grid)
Stack<-stack(List_grid$Grid,Stack_maps1) ###Full stack of all maps
#
#Turn maps into dataframe for computation of PD
community_data<-as.data.frame(Stack)
##remove rows that are NA 
community_data1<-community_data[-which(rowSums(!is.na(community_data[,2:ncol(community_data1)]))==0),] ##change number of columns
species_names<-colnames(community_data1)[2:length(community_data1)] #Store species names


#Store 
setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Diversity_computation/")
#write.table(community_data1, file = "communities.txt", append = FALSE,row.names=F,quote=F,sep="\t") #
arrow::write_feather(community_data1, "communities_functional.feather") ###trying to save less heavy
community_data1[is.na(community_data1)] = 0
#In the community data frame NA must be eliminated done before is this if you load community data and not maps?
community_data1=na.omit(community_data1)  #test if this works 
names(community_data1[1:10]) ##check this for final computation, should only be species names)
##Computation of FD here##
FDval<-FD::dbFD(traits1,community_data1,calc.FRic=T,calc.FDiv=F,m=5) ##Test this it might fail at PCoA and would have to create a modified version of function
FDvals<-c(FDval[["FRic"]],FDval[["FDis"]])
#Add the pixel FD values to data frame (2 different columns)
##check structure of FDvals before running
community_data1$fric<-FDvals[,1]
community_data1$fdis<-FDvals[,2]

#Write the new matrix to a file to avoid rerunning the script for potential further analyses
#write.table(community_data1, file = "communities_and_pd.txt", append = FALSE,row.names=F,quote=F,sep="\t")
arrow::write_feather(community_data1, "communities_withFD.feather")


########SPATIALIZE INFORMATION#########
##Generate a raster containing FD information per pixel##

#1-First generate an empty raster using a base model for resolution and area

values(r)<-0
fd_ras<-r
values(fd_ras)<-NA #Eliminate every value they will be replaced by PD values further down
fd_ras1<-fd_ras
#2- Assign FD values to raster
fd_ras[community_data1$grid]<-community_data1$fric
fd_ras1[community_data1$grid]<-community_data1$fdis

#3- Save raster to file 

writeRaster(fd_ras,"Trees_FRic_polygon.tif",format="GTiff")
writeRaster(fd_ras1,"Trees_FDis_polygon.tif",format="GTiff")

#4-Optional plotting map in R 
par(mfrow=c(1,2))
plot(fd_ras)
plot(fd_ras1)

 
    
    if(length(na.omit(species))>1){
      presence<-rep(0,length(all_species))
      presence[which(all_species%in%species)]<-1
      community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",all_species))
      # traits2<-traits_select[which(rownames(traits_select)%in%species),]
      ##replace this by FD::fdisp but would have to replace traits by distance matrix
      FDval<-FD::dbFD(traits_select,community,calc.FRic=T,calc.FDiv=F,m=5)
      FDvals<-c(FDval[["FRic"]],FDval[["FDis"]])
      
    }
    else { FDvals<-c(NA,NA)}
  }
  else { FDvals<-c(NA,NA)}
  return(FDvals)
}














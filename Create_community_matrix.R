###November 2022###
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

raster_base<-raster("~/Dropbox/wc2/wc2.1_10m_bio_1.tif") ### 
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
community_data<-as.data.frame(Stack_maps)

#Store 
setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Diversity_computation/")
write.table(community_data, file = "communities.txt", append = FALSE,row.names=F,quote=F,sep="\t")
arrow::write_feather(community_data, "communities.feather")




##Prepare phylogenetic tree##
species_names<-gsub("_"," ",species_names)

# sample call, assuming unique_spp_names is in Latin binomial format
phylo_tree <- V.PhyloMaker::phylo.maker(data.frame(species = species_names, genus = stringr::word(species_names, 1), family = rep(NA, length(species_names))))
##note not all trees could be place (~3K not in the tree) ###create list of missing trees
species_not_tree<-subset(phylo_tree$species.list,status=="fail to bind")
write.table(species_not_tree,"Species_not_in_tree.txt")
species_names<-gsub(" ","_",species_names)
species_names<-stringr::str_to_title(species_names)

####MMake sure the names match between maps and phylogeny
setdiff(phylo_tree$scenario.3$tip.label, species_names) ##chack var names here
##if it doesnt use this code
###Trim phylogeny to match distribution data (remove non hylids)
pruned.tree<-ape::drop.tip(phylo_tree$scenario.3$tip.label, setdiff(phylo_tree$scenario.3$tip.label, species_names))
##check that all distribution data is in tree
test<-as.data.frame(species_names)
rownames(test)<-species_names
check_names<-geiger::name.check(pruned.tree, test, data.names=NULL) #this must be OK

write.tree(pruned.tree$scenario.3,file="TEST_all_trees_tree.nwk")






















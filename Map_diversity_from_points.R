#https://github.com/jinyizju/S.PhyloMaker
library(dplyr)
install.packages("remotes")
remotes::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
setwd("Dropbox/**PostDoc_ETH/Trees_PD_FD/")
##Call the shapefile with occurrences cleaned one (for now only cleaned of ocean points)
all_trees<-rgdal::readOGR(dsn=".","all_trees_points_Land")

##Get the list of species from the list of names in the occurrences file

all_trees$Species<-gsub("_"," ",all_trees$Species)
all_species<-unique(all_trees$Species)
write.table(all_species,"all_trees_points_neo.txt")
# sample call, assuming unique_spp_names is in Latin binomial format
phylo_tree <- V.PhyloMaker::phylo.maker(data.frame(species = all_species, genus = stringr::word(all_species, 1), family = rep(NA, length(all_species))))
      ##note not all trees could be place (~3K not in the tree)
write.tree(phylo_tree$scenario.3,file="all_trees_tree.nwk")
###create list of missing trees
species_not_tree<-subset(phylo_tree$species.list,status=="fail to bind")
write.table(species_not_tree,"Species_not_in_tree.txt")
all_trees$Species<-gsub(" ","_",all_trees$Species)
all_trees$Species<-stringr::str_to_title(all_trees$Species)
##create maps of diversity measures
    ##first a world empty raster
world<-rgdal::readOGR("world_adm0.shp")
r<-raster::raster(ncols=80000000,nrows=80000000)
res<-2
raster::res(r)<-res
raster::values(r)<-NA
r<-raster::crop(r,world)
raster::crs(r)<-"+proj=longlat +datum=WGS84 +no_defs"
###Define functions to get desired calculations
pdcalc<-function(x,na.rm=T){
  if (na.rm==T){
    if(length(na.omit(x))>1){
    presence<-rep(1,length(na.omit(x)))
    community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",na.omit(x)))
    pdval<-picante::pd(samp=community,tree = phylo_tree$scenario.3,include.root = F)[[1]]
    }
    else{
    pdval<-NA
     }
    }
  else{ 
    if(length(x)>1){
    presence<-rep(1,length(x))
    community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",x))
    pdval<-picante::pd(samp=community,tree = phylo_tree$scenario.3,include.root = F)[[1]]}
    else{
      pdval<-NA
    }}
  return(pdval)
  }

mpdcalc<-function(x,na.rm=T){
  if (na.rm==T){
    if(length(na.omit(x))>1){
      #remove species not in tree from list
      species_not_tree<-geiger::name.check(phylo_tree$scenario.3,na.omit(x),data.names=na.omit(x))
      if(length(species_not_tree)>1){
        species<-setdiff(x,species_not_tree$data_not_tree)
        }
      else {species<-na.omit(x)}
      if(length(na.omit(species))>1){
      presence<-rep(1,length(na.omit(species)))
       #prune tree to only community
      comm_phylo<-keep.tip(phylo_tree$scenario.3, na.omit(species))
      distance_phylo<-ape::cophenetic.phylo(comm_phylo)
      community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",na.omit(species)))
      mpdval<-picante::mpd(samp=community,dis = distance_phylo)
      }
      else {mpdval<-NA}
    }
    else{mpdval<-NA}
  }
  else{ 
    if(length(x)>1){
      #remove species not in tree from list
      species_not_tree<-geiger::name.check(phylo_tree$scenario.3,x,data.names=x)
      if(length(species_not_tree)>1){
        species<-setdiff(x,species_not_tree$data_not_tree)
      }
      else {species<-x}
      if(length(na.omit(species))>1){
        presence<-rep(1,length(na.omit(species)))
        #prune tree to only community
        comm_phylo<-keep.tip(phylo_tree$scenario.3,species)
        distance_phylo<-ape::cophenetic.phylo(comm_phylo)
        community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",species))
        mpdval<-picante::mpd(samp=community,dis = distance_phylo)
      }
      else {mpdval<-NA}
    }
    else{mpdval<-NA}
  }
  return(mpdval)
}

mntdcalc<-function(x,na.rm=T){
  if (na.rm==T){
    if(length(na.omit(x))>1){
      #remove species not in tree from list
      species_not_tree<-geiger::name.check(phylo_tree$scenario.3,na.omit(x),data.names=na.omit(x))
      if(length(species_not_tree)>1){
        species<-setdiff(x,species_not_tree$data_not_tree)
      }
      else {species<-na.omit(x)}
      if(length(na.omit(species))>1){
        presence<-rep(1,length(na.omit(species)))
        #prune tree to only community
        comm_phylo<-keep.tip(phylo_tree$scenario.3, na.omit(species))
        distance_phylo<-ape::cophenetic.phylo(comm_phylo)
        community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",na.omit(species)))
        mntdval<-picante::mntd(samp=community,dis = distance_phylo)
      }
      else { mntdval<-NA}
    }
    else{ mntdval<-NA}
  }
  else{ 
    if(length(x)>1){
      #remove species not in tree from list
      species_not_tree<-geiger::name.check(phylo_tree$scenario.3,x,data.names=x)
      if(length(species_not_tree)>1){
        species<-setdiff(x,species_not_tree$data_not_tree)
      }
      else {species<-x}
      if(length(na.omit(species))>1){
        presence<-rep(1,length(na.omit(species)))
        #prune tree to only community
        comm_phylo<-keep.tip(phylo_tree$scenario.3,species)
        distance_phylo<-ape::cophenetic.phylo(comm_phylo)
        community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",species))
        mpdval<-picante::mpd(samp=community,dis = distance_phylo)
      }
      else { mntdval<-NA}
    }
    else{ mntdval<-NA}
  }
  return( mntdval)
}

richness<-function(x,na.rm=T){
  if (na.rm==T){
  richness<-length(unique(na.omit(x)))
  }
  else {
    richness<-length(unique(x))
    }
  
  return(richness)
}

###TRYING MPD 

##create a stack with all div measures
my_summary <- function(x, na.rm) c(richness = richness(x, na.rm=T), pd = pdcalc(x, na.rm=T),mpd = mpdcalc(x, na.rm=T),mntd = mntdcalc(x,na.rm=T))
#my_summary <- function(x, na.rm) c(mpd = mntdcalc(x, na.rm=T))

diversity <- raster::rasterize(all_trees, r, 'Species', fun=my_summary)
names(diversity)<-c("Richness","PD","mpd","mntd")
###Create residuals for SR/PD
srpd_df<-stack(diversity$PD,diversity$Richness)
srpd_df<-as.data.frame(srpd_df)
linear_model_sr_pd<-lm(srpd_df)
linear_residuals_sr_pd<-resid(linear_model_sr_pd)
srpd_df$PD_linear_residuals<-rep(NA,length(srpd_df[,1]))
srpd_df$PD_linear_residuals[as.numeric(names(linear_residuals_sr_pd))]<-linear_residuals_sr_pd[names(linear_residuals_sr_pd)]
linear_residuals_sr_pd_raster<-raster(diversity$Richness)
values(linear_residuals_sr_pd_raster)<-NA
values(linear_residuals_sr_pd_raster)<- srpd_df$PD_linear_residuals
diversity<-stack(diversity,linear_residuals_sr_pd_raster)
  plot(srpd_df$Richness,srpd_df$PD,xlab="Richness",ylab="Phylogenetic diversity")
  abline(linear_model_sr_pd,col="red")
  names(diversity)<-c("Richness","PD","mpd","mntd","Residuals")
###write all rasters to file
  writeRaster(diversity,c("Neo_Diversity_2degree_SR.tiff","Neo_Diversity_2degree_PD.tiff","Neo_Diversity_2degree_mpd.tiff","Neo_Diversity_2degree_MNTD.tiff","Neo_Diversity_2degree_SRPDres.tiff"),bylayer=T,format="GTiff")

  par(mfrow=c(2,3))
   plot(diversity$Richness,main="Richness")
   plot(diversity$PD,main="PD")
   plot(diversity$mpd,main="MPD")
   plot(diversity$mntd,main="MNTD")
   plot(diversity$Residuals,main="Residuals",zlim=c(-20000,20000))

   
   
####Functions for FD calculations. Trying FDrichness, Fdispersion  and? 
   ##First match species with traits? or within function ? 
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
  #all_species<-rownames(traits1)


 functions<-function(x,all_species=all_species,na.rm=T){
   if(length(na.omit(x))>1){
     #remove species not in traits from list
       species<-intersect(x,rownames(traits1))
       species<-species[order(species)]
   
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
   my_summaryFD <- function(x, na.rm) c(FR= functions(x,all_species=all_species, na.rm=T))
  
   Functionaldiversity1 <- raster::rasterize(all_trees, r, 'Species', fun=my_summaryFD)
   names(Functionaldiversity)<-c("FRic","FDis")
   raster::writeRaster(Functionaldiversity,c("Neo_Diversity_2degree_FRic.tiff","Neo_Diversity_2degree_FDis.tiff"),bylayer=T,format="GTiff")
   
###SESPD must be run with at least 2 communities so cannot be done during rasterize ... leave for alphahull method
sespd<-function(x,na.rm=T){
  if (na.rm==T){
    if(length(na.omit(x))>1){
      presence<-rep(1,length(na.omit(x)))
      community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",na.omit(x)))
      sespd<-picante::ses.pd(samp=community,tree = phylo_tree$scenario.3,include.root = F,null.model="taxa.labels")
      sespd<-sespd[,2]
    }
    else{
      sespd<-NA
    }
  }
  else{ 
    if(length(x)>1){
      presence<-rep(1,length(x))
      community<-matrix(presence,byrow=T,nrow=1,dimnames=list("community",x))
      sespd<-picante::ses.pd(samp=community,tree = phylo_tree$scenario.3,include.root = F,null.model="taxa.labels")
      sespd<-sespd[,2]
    }
    else{
      sespd<-NA
    }}
  return(sespd)
}


#Compute Functional diversity metrics with all communities on the matrix to standarize the PCA!
#all_trees<-sf::st_read("all_trees_points_Land_troptreeNT.shp")

all_trees<-rgdal::readOGR("all_trees_points_Land_troptreeNT.shp")
all_trees<-all_trees[,1]
all_trees$Species<-gsub(" ","_",all_trees$Species)
all_trees$Species<-stringr::str_to_title(all_trees$Species)

#all_trees<-sf::st_as_sf(all_trees)
world<-rgdal::readOGR("world_adm0.shp")
r<-raster::raster(ncols=80000000,nrows=80000000)
res<-2
raster::res(r)<-res
names(r)="Grid"
r[1:raster::ncell(r)]<-1:raster::ncell(r)
r<-raster::crop(r,world)
r<-raster::mask(r,world)
raster::crs(r)<-"+proj=longlat +datum=WGS84 +no_defs"
poly<-raster::rasterToPolygons(r)
#poly<-sf::st_as_sf(poly)
all_trees_grid<-sp::over(poly, all_trees, returnList=T)

#all_species<-rownames(traits1)
createCommunity<-function(x){
  species<-unique(x$Species)
}
all_trees_grid<-lapply(all_trees_grid,createCommunity)
all_species<-unique(unlist(all_trees_grid))
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
traits2<-traits1[which(rownames(traits1)%in%all_species),]
###log transform all traits to be used!!!

all_species<-rownames(traits2)
all_species<- all_species[order(all_species)]

presence<-rep(0,length(all_species))
community<-matrix(rep(presence,length(all_trees_grid)),nrow=length(all_trees_grid),byrow=T,dimnames=list(names(all_trees_grid),all_species))

for (i in 1:length(all_trees_grid)){
  
  species<-all_trees_grid[[i]]
  community[as.character(i),which(all_species%in%species)]<-1
}
####eliminate empty communities
community<-community[as.vector(which(rowSums(community)!=0)),]

###calculate FD indices for communities in matrix

#FDval<-FD::dbFD(traits2,community,calc.FRic=F,calc.FDiv=F) ##this fails out of memory
###Try creating the distance matrix first

#distance<-dist(traits2)
#community<-as.data.frame(community)
#ade4::is.euclid(distance)
#pcoa<-ade4::dudi.pco(distance,scannf=F,nf=5) #fails
#pca<-vegan::rda(traits2,scale=T) 
#distance<-dist(pca$CA$u)
#raoQ<-ade4::divc(community,distance)
#FD::fdisp(distance,community)
#community_transpose <- as.data.frame(t(as.matrix(community)))
#dispersion<-rao_Paz(community_transpose,distance)
##if doing for whole amtrix how s the output stored?
#FDvals<-c(FDval[["FRic"]],FDval[["FDis"]])

#### For an euclidean distance matrix it should be the same doing PCA on df than pcoa on dis so using PCA axes as input for multidimFD should get 
  ####around the problem of memory for pcoa computation
   ###Still must log transform the traits before the PCA
###trnasfor data
###do PCA
pca<-vegan::rda(traits2,scale=T) 
traits3<-pca$CA$u

FD<-multidimFD(traits3, community, check_species_pool=TRUE, verb=TRUE,
                     folder_plot=NULL, nm_asb_plot=NULL, Faxes_plot=NULL, Faxes_nm_plot=NULL, 
                     plot_pool=FALSE, col_bg="grey90", col_sp_pool="grey30", pch_sp_pool="+", cex_sp_pool=1,
                     pch_sp=21, col_sp="#1145F0", transp=50 )
FD_df<-as.data.frame(FD)
FD_df$Grid<-poly$Grid[as.numeric(rownames(FD_df))]
functional_raster<-r
values(functional_raster)<-NA
functional_raster[FD_df$Grid]<-FD_df$FDis
plot(functional_raster)
writeRaster(functional_raster,"Neo_diversity_2degree_FDIS_allcomms.tif",format="GTiff")

##NOW FOR FRic?

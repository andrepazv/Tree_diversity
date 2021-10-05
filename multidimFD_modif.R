###########################################################################################################################################
#
# 'multidimFD': function to compute and illustrate multidimensional functional diversity indices for a set of species assemblages
# For details about indices formulas see Mouillot et al. 2013, Trends in ecology and Evolution (28:167-177) and references therein.
#                 
# INPUTS: 
#       - 'coord': a matrix with coordinates of species (rows) in a multidimensional functional space (columns)
#		    - 'weight': a matrix with weight of species (columns) in a set of assembalges (rows). 
#                   Weight could be any continuous variable:  biomass, relative biomass, number of individuals, percent coverage.
#                   All FD indices are computed using relative weight so that they are not affected by unit (e.g. g or kg for biomass).
#                   Thus, in the special case where 'weight' is filled with only 0/1 (to code for presence/absence), 
#                           then FD indices will be computed assuming all species have the same weight
#       - 'check_species_pool" : a logical value indicating wether two tests are performed to check that all rows and all columns 
#                               of 'weight' have a sum stricly positive (i.e. all species should be present in at least one assemblage 
#                               and all assemblages should have at least one species). Default is TRUE.
#		    - 'verb': a logical value (default=TRUE) indicating whether printing progress of computation process
#
#     NB :  
#           Column names of 'weight' should be the same than row names of 'coord' (i.e. same codes and same order).
#           If 'coord' has no names, default names will be set ('Axis_1', 'Axis_2',...,'Axis_n').
#           'weight' should have row names
#           NA are not allowed in 'coord' or in 'weight'.
#           Only positive numbers are allowed in 'weight'.
#           There should be at least 2 functional axes and 2 species
#           If some assemblages have fewer species than number of axes, FRic and FDiv indices could not be computed (see below).
#
#       - 'folder_plot': a character string for setting the working directory where plots will be saved, default is current working directory.
# 		  - 'nm_asb_plot': a vector with the names of assemblages for which FD indices need to be illustrated. 
#                       Default is NULL (i.e. no plot). To plot FD for all assemblages, set to 'row.names(weight)'.
#       - 'Faxes_plot': a vector with names of axes to plot. Should be of length from 2 to 4. 
#               Default is NULL which means that the four first axes will be plotted.
#               Value of 'Faxes_plot' does not affect the way FD indices are computed (i.e. acccording to all columns of 'coord').
#       - 'Faxes_nm_plot': a vector with titles of axes (default is titles as 'Faxes_plot') 
# 		  - 'plot_pool': a logical value indicating whether all species present in 'coord' need to be illustrated on all plots.
#                      If yes (default value), space filled by the pool of species is in white while background is filled 
#                         with color specified in 'col_bg' (default is light grey), and absent species are plotted with a different symbol (see below).
#       - 'col_bg': a R color name (or hexadecimal code). See above.
#       - "pch_sp_pool", "cex_sp_pool", 'col_sp_pool': Shape, size and color of symbol to plot absent species. See above.
#       - 'pch_sp': a numeric value coding shape of symbol to plot species, as in function 'points' (default is 21, i.e. points).
#       - 'col_sp': a character string for hexadecimal code (e.g. from www.colorpicker.com) for color used for symbols and convex hull.
# 		  - 'transp': a single numeric value indicating the percentage of transparency for convex hull filling. Default is 50%.
#
#
# OUTPUTS: 
#     => a matrix with for each assemblage (row) values for a set of diversity indices (columns):
#       - 'Nb_sp': number of species present
#       - 'Tot_weight': total weight (e.g. biomass, number of individuals)
#       - 'min_f', 'max_f', 'range_f': minimum, maximum and range of values along 'f' functional axis
#       - 'FIde_f': weighted mean position along 'f' functional axis
#       - 'FRic': functional richness (proportion of functional space filled by species present in the assemblage)
#       - 'FDiv': functional divergence (deviation of species weight to the center of the convex hull filled by the species assemblage)
#       - 'FEve': functional evenness (regularity of distribution of species weights in the functional space)
#       - 'FDis': functional dispersion (weighted mean distance to the average position of the species present in the assemblage 
#                                             dividied by half the maximum distance among all the species present in the set of assemblages )
#       - 'FSpe': functional specialization (weighted mean distance to species pool centroid, i.e. average position of all the species 
#                                           present in the set of assemblages, divided by the maximum distance to this centroid)
#       - 'FOri': functional originality (weighted mean distance to nearest species from the species pool
#                                             divided by the maximum distance to the nearest neighbour)
#
#       NB: FRic and FDiv are computed only for assemblages with more species than number of functional axes
#           FEve is computed only for assemblages of at least 3 species
#           Scaling of FDis, FSpe and FOri indices were scaled by maximum value possible given all the species present in the set of assemblages 
#                  so that they range from 0 to 1, they are unitless and easily interpretable (as FRic, FDiv and FEve)
#
#     => for each of the assemblages listed in 'nm_asb_plot':a 6-panels jpeg file named 'AssemblageA_AxisX_AxisY.jpeg' 
#             illustrating FD indices of assemblage 'A' in the functional space defined by axes 'X' And 'Y. 
#             Jpeg file has a resolution of 150dpi and dimensions of 1800x1200 pixels which means a size of around 200ko 
#
#       All plots (i.e. all assemblages and all pairs of axes) have the same axis scale to faithfully represent FD.
#       Species present in the assemblage are plotted with a 'pch_sp' symbol. 
#       Weights of species are illustrated proportionally to symbol area and a legend is displayed in the bottom right corner.
#
#          * top left panel: functional identity and functional dispersion. 
#                 Weighted-mean position of the species on each axis is illustrated with a square and horizontal and vertical dashed lines. 
#                 Distance to this point are shown with segments
# 			   * top middle panel: functional richness.
#					        The colored convex polygon is a projection of the multidimensional convex hull in 2D. 
#                 Filled symbols are species being vertices in the multidimensional space. 
#                 Minimum and maximum values on each axis are illustrated by vertical bars.
# 			   * top right panel: functional divergence.  
#                 The center of gravity of the vertices is show with a diamond and al the distances to it are shown with lines.
# 			   * bottom left panel: functional evenness. 
#                 The minimum spanning tree linking all points in the multidimensional space is shown.
# 			   * bottom middle panel: functional specialization
#                 Center of gravity of all points is shown as a black square and the distances to it for species present are shown with dotted lines.
# 			   * bottom right panel: functional originalty.
#                 Distances to nearest species in the multidimensional functional space are shown with black arrows (basis=focal species, head=nearest neighbour).
#
##########################################################################################################################################


multidimFD<-function(coord, weight, check_species_pool=TRUE, verb=TRUE,
                     folder_plot=NULL, nm_asb_plot=NULL, Faxes_plot=NULL, Faxes_nm_plot=NULL, 
                     plot_pool=FALSE, col_bg="grey90", col_sp_pool="grey30", pch_sp_pool="+", cex_sp_pool=1,
                     pch_sp=21, col_sp="#1145F0", transp=50 ) 
  
{
  
  # library required for indices computation
  require (geometry)
  require(ape)
  
  # saving name of current working directory
  current_wd<-getwd()
  
  ##############################################################################
  # checking inputs
  
  
  # coordinates of species in the functional space
  if( nrow(coord)<2 ) stop(paste(" error: there must be at least 2 species in the dataset"))
  if( ncol(coord)<2 ) stop(paste(" error: there must be at least 2 functional axes"))
  if( is.numeric(coord)==FALSE ) stop(paste(" error: 'coord' is not numeric"))
  if( is.na(sum(coord)) ) stop(paste(" error: NA in 'coord'"))
  if ( is.null(colnames(coord)) ) { colnames(coord)<-paste("Axis", 1:ncol(coord),sep="_") } # if no column names in 'coord' default value
  
  # dominance of species in assemblages
  if( is.matrix(weight)==FALSE ) stop( " 'weight' should be an object of type matrix")
  if( is.numeric(weight)==FALSE ) stop(paste(" error: 'weight' is not numeric"))
  if( is.na(sum(weight)) ) stop(paste(" error: NA in 'weight'"))
  if( min(weight)<0 ) stop(paste(" error: negative value in 'weight'"))
  if(min(apply(weight,1,sum))==0 ) 
    stop(paste(" error: all rows of 'weight' should have a sum stricly positive, i.e. all assemblage must host at least one species"))
  
  # match between datasets
  if( sum(colnames(weight) == row.names(coord))!= nrow(coord) ) stop(paste(" error: 'weight' does not have the same column names than row names of 'coord'"))
  
  # checking graphical parameters
  if( length(pch_sp)!=1 ) stop(paste(" error:'pch_sp' should contain only one value"))
  if( length(col_sp)!=1 ) stop(paste(" error:'col_sp' should contain only one value"))
  if( length(col_bg)!=1 ) stop(paste(" error:'col_bg' should contain only one value"))
  if( length(col_sp_pool)!=1 ) stop(paste(" error:'col_sp_pool' should contain only one value"))
  if( length(pch_sp_pool)!=1 ) stop(paste(" error:'pch_sp_pool' should contain only one value"))
  if( length(cex_sp_pool)!=1 ) stop(paste(" error:'cex_sp_pool' should contain only one value"))
  
  
  # checking species pool
  if (check_species_pool==TRUE)
  {
    if(min(apply(weight,2,sum))==0 ) 
      stop(paste(" error: all columns of 'weight' should have a sum stricly positive, i.e. all species must occur in at least one assemblage"))
  }# end of check species pool
  
  
  ##############################################################################
  # info on study case
  
  # number and names of axes
  nm_axes<-colnames(coord)
  nb_axes<-length(nm_axes)
  
  # number and names of assemblages
  nm_asb<-row.names(weight)
  nb_asb<-length(nm_asb)
  
  # matrix to store results
  indices<- c(paste("min",nm_axes,sep="_"), paste("max",nm_axes,sep="_"), paste("range",nm_axes,sep="_"),paste("FIde",nm_axes,sep="_"),c("FRic","FDiv","FEve","FDis","FSpe", "FOri") )
  FD<-matrix(NA, nb_asb, length(indices), dimnames=list(nm_asb,indices))
  
  ##############################################################################
  # preliminary computation at the species pool level
  
  #######################################
  # originality of each species: distance to nearest neighbour among the global pool of species
  dist_sp<-as.matrix(dist(coord,method="euclidean")) ; dist_sp[which(dist_sp==0)]<-NA
  orig_sp<-apply(dist_sp, 1, min, na.rm=T )
  # identity of Nearest Neighbour
  NN<-dist_sp ; NN<-NN-apply(NN,1,min,na.rm=T) ; NN[which(NN!=0)]<-NA   ; NN[which(NN==0)]<-1
  
  # specialization of each species: distance to centroid of the global pool of species
  centroid_sp<-apply(coord,2,mean) # coordinates of the center of gravity of the vertices (B)
  spec_sp<-apply(coord, 1, function(x) { (sum((x-centroid_sp)^2) )^0.5} )
  
  # convex hull volume of the species pool
  #FRic_pool<-convhulln(coord,"FA")$vol
  
  #######################################
  # setting same graphical parameters for all assemblages
  
  # setting working directory to store jpeg files
  
  ###########################################
  
  # end of preliminary computation
  
  ##############################################################################
  
  # loop on assemblages for computing and plotting functional diversity indices
  for (k in nm_asb)
  {
    
    ###########################################################
    # preparing data
    
    # names, number, weight and coordinates of of species present
    weight_k<-weight[k,]
    nm_sp_k<-row.names(coord)[which(weight_k>0)]
    nb_sp_k<-length(nm_sp_k)
    weight_sp_k<-weight[k,nm_sp_k]
    coord_sp_k<-coord[nm_sp_k,]
    if(nb_sp_k==1) { coord_sp_k<-matrix(coord_sp_k,nrow=1,dimnames=list(nm_sp_k,colnames(coord)) ) } # matrix object
    
    # names of species absent
    nm_sp_absent_k<-names(which(weight[k,]==0))
    
    # total weight
   
    
    #relative weight 
    rel_weight_sp_k<-weight_sp_k/sum(weight_sp_k)
    
    # species richness
    #  FD[k,"Nb_sp"]<-nb_sp_k
    
    ###########################################################
    # computing indices values on each axis
    
    # range of values
    for (z in nm_axes) 
    {
      FD[k,paste("min",z,sep="_")]<-min(coord_sp_k[,z])
      FD[k,paste("max",z,sep="_")]<-max(coord_sp_k[,z])
      FD[k,paste("range",z,sep="_")]<-FD[k,paste("max",z,sep="_")]-FD[k,paste("min",z,sep="_")]
    }# end of z
    
    # abundance-weighted mean values
    FD[k,paste("FIde",nm_axes,sep="_")]<-rel_weight_sp_k%*%coord_sp_k
    # range of values
    
    
    ###########################################################  
    # multivariate indices
    
    
    # indices based on vertices and volume of convex hull, only if more species than number of axes
    
    if (nb_sp_k>nb_axes) {
      
      ########################
      # Functional richness = convex hull volume
      #FD[k,"FRic"]<-round(convhulln(coord_sp_k,"FA")$vol/FRic_pool,6)
      
      
      ########################
      
      
    }# end of if more species than number of axes
    
    ##########################
    # Functional Evenness
    
    
    ##########################
    # Functional Dispersion: abundance-weighted mean distance to abundance-weighted centroid 
    # scaled by maximum value possible given species pool (i.e. the two most distant species have half of total weight)
    dist_centr_sp_k<-apply(coord_sp_k, 1, function(x) { (sum((x-FD[k,paste("FIde",nm_axes,sep="_")])^2) )^0.5} ) # distance to abundance-weighted centroid 
    FD[k,"FDis"]<-(rel_weight_sp_k %*% dist_centr_sp_k) / ( max(dist_sp, na.rm=T) /2 )
    
    
    
    
    ######################################################################################################################
    # End of indices computation
    ######################################################################################################################
    
    
    ###########################################################
    # if graphical output
    
    
    ###########################################################
    
    # printing step achieved
    if (verb==TRUE) print(paste("FD of assemblage '",k,"' computed",sep="") )
    print(k)
  }# end of working on assemblage k
  ###########################################################  
  
  # returning to current working directory
  setwd(current_wd)
  
  # returning results	
  return(FD)	
  
}# end of function multidimFD
########################################################################################################################################
########################################################################################################################################

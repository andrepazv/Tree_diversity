rm(list=ls())
gc()

# library(geometry)
library(doParallel)
# library(vegan)
# library(picante)
# library(GGally)
library(feather)
# library(arrow)
library(tidyverse)

setwd("~/Git/Tree_diversity/simulations/raw_data/")

##remove PD indices from matrix
community_data1<-read.csv("NoMono_matchPD_FD_communities.csv",row.names=1)
# community_data2<-community_data1[,1:43918] ##take out diversity measures if they are in matrix
##take out the monocots... load list of all then get monocots and then only genera
community_data1[is.na(community_data1)] = 0
community_data1=na.omit(community_data1)  
#names(community_data1[1:10]) ##check this for final computation, should only be species names)

#traits<-arrow::read_feather("ALL_BGCI_species_and_genera_imputed_traits_8_21.feather")
traits<-read.csv("Species_trait_values_3_11_22.csv")
traits1<-matrix(nrow=length(unique(traits$accepted_bin)),ncol=length(unique(traits$trait)),dimnames=list(unique(traits$accepted_bin),unique(traits$trait)))   
##re-check this part why ?
for (i in unique(traits$trait)){
	test<-subset(traits,trait==i)
	traits1[,i]<-test$value
}
rownames(traits1)<-gsub(" ","_",rownames(traits1))
colnames(traits1)<-gsub(" ","_",colnames(traits1))
traits1<-as.data.frame(traits1)
traits1<-traits1[,c(1,18,4,13,16,7,6)] ##selection from Dan's paper


##match some more species to traits
matching<-read.csv("matched_trait_names.csv")
matching$raw_name<-gsub(" ","_",matching$raw_name)
matching$trait_name<-gsub(" ","_",matching$trait_name)
matching<-matching[!duplicated(matching$trait_name),]
matching<-matching[!duplicated(matching$raw_name),]
##change routes are different depending on computer
#monocots<-read.csv("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Datasets/Species_trait_values_3_11_22.csv")
monocots<-read.csv("Species_trait_values_3_11_22.csv")

monocots<-subset(monocots,monocot==1)
monocots<-matrix(unique(monocots$accepted_bin)) 
monocots<-sub(" .*", "", monocots)
monocots<-unique(monocots)
test<-list()
for (i in 1:length(monocots)){
	test[[i]]<-which(grepl(monocots[i],colnames(community_data1))) 
}
test<-unlist(test)
community_data2<-community_data1[,-test]
community<-community_data2 
pca<- vegan::rda(traits1,scale=T) 
traits3<-pca$CA$u 
species_names<-colnames(community_data2)

# sort through the species and fix the names
for(i in 1:nrow(matching)){
	rownames(traits3)[which(rownames(traits3) == matching$trait_name[i])]<-matching$raw_name[i]
}
# rownames(traits3)[which(rownames(traits3)%in%matching$trait_name)]<-matching$raw_name


# convert community
community1<-as.data.frame(community)


traits3<-traits3[which(rownames(traits3)%in%names(community1)),]
community1<-community1[,which(names(community1)%in%rownames(traits3))]  
x.dist<-dist(traits3)  #
##function from Dan (see github for documentation)

# renamed them myd and myp since you're subsetting, and this would avoid any issues of subsetting the global variable
pairdist <- function(myd, myp){
	# only works if p is 1/0 or TRUE/FALSE
	rows<-which(myp==1)
	myd <- myd[rows,rows]^2
	n<-sum(myp)
	n2 <- sum(myp)^2
	# get the weighted pairwise dist
	wpd <- myd / n2
	# for rao's q, just sum up
	rq <- sum(wpd)
	# calculate abundance weighted pairwise distance 
	wmpd <- sum(wpd[upper.tri(wpd)]) / ((n*(n-1)/2)*(1/n2))
	# calculate standard unweighted species pairwise distance
	mpd <- mean(myd[upper.tri(myd)])
	return(c("raoq" = rq, "wmpd" = wmpd, "mpd" = mpd))
}


raoq<-list()
distance<-as.matrix(x.dist)

# get the multi species communities
keep_comm <- (apply(community1, 1, sum) > 1)
# get the single communities
single_comms <- rownames(community1)[!keep_comm]
# remove the single communities (tacked on to the final results)
community1 <- community1[keep_comm, ] 

registerDoParallel(72)

# rm(list = setdiff(ls(), c("community1", "distance", "pairdist", "single_comms")))
# gc()

j <- which(rownames(community1) == "11138")


# sort the community matrix -- missing step
community1 <- community1[, rownames(distance)]

# check to make sure the names are matched
if(any(rownames(distance)!=names(community1)) | any(colnames(distance)!=names(community1))){
	stop("Name mismatch")
}

outmat <- foreach(j = 1:nrow(community1), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{
# for(i in 1:nrow(community1)){
	# print(paste(i, "of", nrow(community1)))
	my_vec <- as.numeric(as.vector(community1[j,]))
	res <- pairdist(distance, my_vec)

	return(c("grid_id" = as.numeric(rownames(community1)[j]), "n" = sum(my_vec), "mpd" = as.numeric(res["mpd"]), "raoq" = as.numeric(res["raoq"]), "wmpd" = as.numeric(res["raoq"])))

}

outmat <- outmat %>% as_tibble() %>% bind_rows(tibble(grid_id = as.numeric(single_comms), n = 1, mpd = 0, raoq = 0, wmpd = 0))


write_csv(outmat, "../results/redo_old_script_orig_fixed.csv")




###############
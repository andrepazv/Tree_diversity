#Store 

    ##remove PD indices from matrix
community_data1<-read.csv("NoMono_matchPD_FD_communities.csv",row.names=1)
community_data2<-community_data1[,1:43918] ##take out diversity measures if they are in matrix
##take out the monocots... load list of all then get monocots and then only genera
community_data1[is.na(community_data1)] = 0
community_data1=na.omit(community_data1)  
#names(community_data1[1:10]) ##check this for final computation, should only be species names)


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
rownames(traits3)[which(rownames(traits3)%in%matching$trait_name)]<-matching$raw_name

community1<-as.data.frame(community)
#For RAOQ

  
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

 for(i in 1:nrow(community1)){
  raoq[[i]]<-pairdist(distance, as.numeric(as.vector(community1[i,])))
}
names(raoq)<-rownames(community1)

write.csv(community1)
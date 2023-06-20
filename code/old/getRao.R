pairdist <- function(myd, myp){
  # only works if p is 1/0 or TRUE/FALSE
  rows<-which(myp==1)
  myd <- myd[rows,rows]
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

###then with community being a dataframe of rows for communities and columns for species and distance the pairwise distances between traits but as matrix. 
 for(i in 1:nrow(community)){
raoq[[i]]<-pairdist(distance, as.numeric(as.vector(community[i,])))
}
names(raoq)<-rownames(community)
community$RAO<-NA
raoq1<-lapply(raoq,function(l) l[1]) 
community[names(raoq),"RAO"]<-unlist(raoq1[names(raoq)]) ##this step is necessary to map
rao_raster<-r
values(rao_raster)<-NA
rao_raster[as.numeric(rownames(community))]<-community$RAO 
plot(rao_raster)
       

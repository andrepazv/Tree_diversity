rm(list = ls())
gc()

library(geometry)
library(doParallel)
library(vegan)
library(picante)
library(feather)
library(tidyverse)


# register the parallel cores
registerDoParallel(48)

# raoq formula. assumes equal abundances
pairdist <- function(myd){
	return(sum((myd / nrow(myd))^2))
}


setwd("~/Git/Tree_diversity/simulations/")

comm <- read_csv("results/Community_matrix_for_FD_PD_38k.csv")
dt_mat <- read_csv("results/Trait_matrix_for_FD_PD_38k.csv")
pruned_tree <- read.tree("results/Pruned_tree_for_FD_PD_38k.csv")


# double check that the right files are loaded, all should be zero, otherwise clear the environment
if((sum(!dt_mat$accepted_bin%in%comm$accepted_bin) + 
	sum(!comm$accepted_bin%in%dt_mat$accepted_bin) +
	sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label) +
	sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label) +
	sum(!dt_mat$accepted_bin%in%comm$accepted_bin) +
	sum(!comm$accepted_bin%in%dt_mat$accepted_bin) +
	sum(!unique(comm$accepted_bin)%in%pruned_tree$tip.label) +
	sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label)) >0) {rm(list=ls()); stop("name mismatch")}

# get the plots and shuffle
plts <- sample(unique(comm$grid_id), length(unique(comm$grid_id)), replace = FALSE)

# get no. of traits 
K <- ncol(dt_mat)-1

# get the number of species per plot
nspp <- comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))

# get the unique species
spp <- comm %>% select(accepted_bin) %>% distinct() %>% unlist()


##################
## set params

nrep <- 5
perc_keep <- 0.5 # only keep 50% of species
outfile <- "results/bootstrap_downsample.csv"

####################################
# loop, with each loop in parallel
for(z in 1:nrep){
	outmat_check <- foreach(j = 1:length(plts), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{
		
		set.seed(sample(1:1e9,1))
		# print(j)
		
		# get the focal plot
		my_dt <- comm %>% filter(grid_id == plts[j]) 
		
		# subsample
		my_dt <- my_dt %>% sample_n(round(nrow(my_dt)*perc_keep))
		
		# initialize the variables
		my_dt_sub <- NULL
		pd_cur <- mpd_cur <- pd_sub <- mpd_sub <- my_raoq <- my_fric <- my_dend <- NA
		
		out_res <- c()
		
		if(nrow(my_dt) > 1){
			# subset the data and get the phyometrics
			my_comm <- matrix(1, nrow = 1, ncol = nrow(my_dt)) %>% data.frame() %>% as_tibble() %>% setNames(my_dt$accepted_bin) %>% as.matrix()
			pd_cur <- as.numeric(picante::pd(samp = my_comm, pruned_tree, include.root = FALSE)[1])
			mpd_cur <- picante::mpd(samp = my_comm, dis = cophenetic(ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, my_dt$accepted_bin))))
			
			# # get the trait matrix
			my_tr0 <- dt_mat %>% filter(accepted_bin%in%my_dt$accepted_bin)
			
			# add in noise
			my_tr <- (my_tr0 %>% select(-accepted_bin) %>% as.matrix()) + matrix(runif(nrow(my_tr0)*(ncol(my_tr0)-1), -1e-8,1e-8), nrow = nrow(my_tr0), ncol = ncol(my_tr0)-1)
			
			# add in row names
			rownames(my_tr) <- my_tr0$accepted_bin
			
			# get the pairwise distance between species across all traits
			D <- dist(my_tr) 
			
			# get raoq
			my_raoq <- pairdist(as.matrix(D))
			
			# make sure we have more unique traits than rows
			if(nrow(unique(my_dt)) > ncol(my_tr)){
				my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
			}
			
			# # cluster the traits
			# hc <- hclust(D^2)
			# 
			# # convert to a phylo tree object
			# tr_tree <- as.phylo.hclust(hc)
			# 
			# # get, for example, faith's pd (sum of total branch lengths), for some community abundance matrix
			# dend <- as.numeric(picante::pd(samp = my_comm, tr_tree, include.root = FALSE)[1])
			# 
			# 
			out_res <- c("grid_id" = plts[j], "n" = nrow(my_dt), "pd" = pd_cur, "mpd" = mpd_cur, "raoq" = my_raoq, "fd" = my_fric, "fdr" = my_fric^(1/K)) #, "dend" = dend)
			
		}else{
			my_raoq <- pd_cur <- mpd_cur <- 0
		}
		
		
		return(out_res)
	}
	
	outmat_check <- outmat_check %>% data.frame() %>% as_tibble() #%>% left_join(metrics %>% setNames(paste0("obs_", names(.))) %>% rename(grid_id = obs_grid_id))
	write_csv(outmat_check, outfile, append = file.exists(outfile))
}


rm(outmat_check)

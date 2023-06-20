rm(list = ls())
gc()

library(geometry)
library(doParallel)
library(vegan)
library(picante)
library(feather)
library(tidyverse)


# register the parallel cores
registerDoParallel(12)

# raoq formula. assumes equal abundances
pairdist <- function(myd){
	return(sum((myd / nrow(myd))^2))
}


setwd("~/Git/Tree_diversity/simulations/")

comm0 <- read_csv("results/Community_matrix_for_FD_PD_38k.csv")
dt_mat0 <- read_csv("results/Trait_matrix_for_FD_PD_38k.csv")
pruned_tree0 <- read.tree("results/Pruned_tree_for_FD_PD_38k.csv")

ag <- read_csv("~/Git/dispersion/data/species_names/ANGIO_GYMNO_lookup.csv") %>% 
	mutate(accepted_bin = gsub(" " , "_", accepted_bin)) %>% select(accepted_bin, group)

comm0 <- comm0 %>% left_join(ag)

for(zzz in 1:2){
	
	if(zzz == 1){
		comm <- comm0 %>% filter(group == "Angiosperms")
	}else{
		comm <- comm0 %>% filter(group == "Gymnosperms")	
	}
	
	
	dt_mat <- dt_mat0 %>% filter(accepted_bin%in%comm$accepted_bin)
	pruned_tree <- ape::drop.tip(pruned_tree0, setdiff(pruned_tree0$tip.label, comm$accepted_bin))
	
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
	
	###################################################
	# get the observed fd/pd in parallel
	###################################################
	outmat_check <- foreach(j = 1:length(plts), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{
		
		# get the focal plot
		my_dt <- comm %>% filter(grid_id == plts[j])
		
		# initialize the variables
		my_dt_sub <- NULL
		pd_cur <- mpd_cur <- pd_sub <- mpd_sub <- my_raoq <- my_fric <- NA
		
		out_res <- c()
		
		if(nrow(my_dt)>1){
			
			# subset the data and get the phylometrics
			my_comm <- matrix(1, nrow = 1, ncol = nrow(my_dt)) %>% data.frame() %>% as_tibble() %>% setNames(my_dt$accepted_bin) %>% as.matrix()
			pd_cur <- as.numeric(picante::pd(samp = my_comm, pruned_tree, include.root = FALSE)[1])
			mpd_cur <- picante::mpd(samp = my_comm, dis = cophenetic(ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, my_dt$accepted_bin))))
			
			# get the trait matrix
			my_tr <- dt_mat %>% filter(accepted_bin%in%my_dt$accepted_bin) %>% select(-accepted_bin) %>% as.matrix()
			my_tr <- my_tr + matrix(runif(nrow(my_tr)*ncol(my_tr), -1e-8,1e-8), nrow = nrow(my_tr), ncol = ncol(my_tr))
			
			# get raoq
			my_raoq <- pairdist(as.matrix(dist(my_tr)))
			
			# get convex hull, but only if we have more unique species than traits
			if(nrow(unique(my_dt)) > ncol(my_tr)){
				my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
			}
			
		}else{
			my_raoq <- pd_cur <- mpd_cur <- 0
		}
		
		# save results
		out_res <- c("grid_id" = plts[j], "n" = nrow(my_dt), "pd" = pd_cur, "mpd" = mpd_cur, "raoq" = my_raoq, "fd" = my_fric, "fdr" = my_fric^(1/K))
		
		return(out_res)
	}
	
	
	# save the results
	outmat_check <- outmat_check %>% data.frame() %>% as_tibble() 
	
	# outmat_check %>% gather(metrics, value, -grid_id, -n)
	write_csv(outmat_check %>% gather(metrics, value, -grid_id, -n), paste0("results/obs_results_50vars", ifelse(zzz == 1, "_ANGIO.csv", "_GYMNO.csv")))
}
	
	
	
rm(list = ls())
gc()

library(geometry)
library(doParallel)
library(vegan)
library(picante)
library(feather)
library(tidyverse)


# register the parallel cores
registerDoParallel(round(detectCores()/2))

# # raoq formula. assumes equal abundances
pairdist <- function(myd){
	return(sum((myd / nrow(myd))^2))
}

setwd("~/Git/Tree_diversity/data")

save_suffix <- "_200"


# read in the data
comm0 <- read_csv(paste0("cleaned_data/REV_Community_matrix",save_suffix,".csv"))
dt_mat <- read_csv("cleaned_data/REV_Trait_matrix.csv")
pruned_tree <- read.tree("cleaned_data/REV_Pruned_tree.csv")
ag <- read_csv("raw_data/ANGIO_GYMNO_lookup.csv") %>% 
	mutate(accepted_bin = gsub(" " , "_", accepted_bin)) %>% select(accepted_bin, group)

inv <- read_csv("raw_data/Invasive_glonaf.csv")

print(paste0("Using: ", "REV_Community_matrix",save_suffix,".csv"))

# which null models to do -- 1=all, 2=angio, 3=gymno, 4=no invasives
fit_seq <- c(4)

gg <- 2

for(gg in fit_seq){
	
	if(gg == 2){
		print("fitting ANGIOS")
		comm <- comm0 %>% left_join(ag, by = "accepted_bin") %>% 
			filter(group == "Angiosperms")
		outf <- paste0("results/REV_obs_results_ANGIO", save_suffix, ".csv")
	}else if(gg == 3){
		print("fitting GYMNOS")
		comm <- comm0 %>% left_join(ag, by = "accepted_bin") %>% 
			filter(group == "Gymnosperms")
		outf <- paste0("results/REV_obs_results_GYMNO", save_suffix, ".csv")
	}else if(gg == 4){
		print("removing INVASIVES")
		comm <- comm0 %>% filter(!accepted_bin%in%inv$x)
		outf <- paste0("results/REV_obs_results_INV", save_suffix, ".csv")
	}else{
		print("fitting ALL")
		comm <- comm0
		outf <- paste0("results/REV_obs_results", save_suffix, ".csv")
	}
	
	
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
	write_csv(outmat_check, outf)
}
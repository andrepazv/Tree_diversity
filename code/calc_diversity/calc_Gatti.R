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

# raoq formula. assumes equal abundances
pairdist <- function(myd){
	return(sum((myd / nrow(myd))^2))
}

setwd("~/Git/Tree_diversity/data")

# read in the data
comm <- read_csv("cleaned_data/REV_Community_matrix_CLEANED.csv")
dt_mat <- read_csv("cleaned_data/REV_Trait_matrix.csv")
pruned_tree <- read.tree("cleaned_data/REV_Pruned_tree.csv")


# gatti numbers per contintent
gatti <- rbind(c("Africa", 11875),
			   c("Eurasia", 16264),
			   c("North America", 11131),
			   c("Oceania", 8235),
			   c("South America", 31112)) %>% data.frame() %>% setNames(c("Continent", "n_tot")) %>% as_tibble() %>% mutate(n_tot = as.numeric(n_tot)) %>% 
	mutate(perc = n_tot / sum(n_tot))

# observed per continent
percon <- comm %>% select(accepted_bin, Continent) %>% distinct() %>% group_by(Continent) %>% tally() #%>% mutate(perc = n / sum(n))

# the number of total species to retain all S. American species at the right probability
tot <- floor((percon %>% filter(Continent == "South America") %>% select(n) %>% unlist()) / 
			 	(gatti %>% filter(Continent == "South America") %>% select(perc) %>% unlist()))

# calculate the new number of species per continent, holding S. America fixed
gatti <- gatti %>% mutate(n_new = perc*tot) %>% select(-n_tot) %>% 
	left_join(percon)

# get no. of traits 
K <- ncol(dt_mat)-1

NSIM <- 500

set.seed(10)

# create temp directory
dir.create("gatti_temp", showWarnings = FALSE)

for(i in 1:NSIM){
	
	print(c(i, NSIM))
	
	# randomly downsample per continent to get the desired species
	consp <- comm %>% select(accepted_bin, Continent) %>% distinct() %>% left_join(gatti %>% select(Continent, n_new), by = "Continent") %>% 
		group_by(Continent) %>% sample_n(n_new[1], replace = FALSE) %>% ungroup
	
	# subset the community matrix
	sub_comm <- comm %>% filter(accepted_bin%in%consp$accepted_bin)
	
	# get the number of species per plot
	nspp <- sub_comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))
	
	# get the unique species
	spp <- sub_comm %>% select(accepted_bin) %>% distinct() %>% unlist()
	
	# get the plots and shuffle
	plts <- sample(unique(sub_comm$grid_id), length(unique(sub_comm$grid_id)), replace = FALSE)
	
	###################################################
	# get the observed fd/pd in parallel
	###################################################
	outmat <- foreach(j = 1:length(plts), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{
		
		# get the focal plot
		my_dt <- sub_comm %>% filter(grid_id == plts[j])
		
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
	
	write_csv(outmat %>% data.frame() %>% as_tibble() %>% mutate(rep = i), paste0("gatti_temp/sim",i,".csv"))
}


# read in the files and combine
files <- paste0("gatti_temp/", list.files("gatti_temp"))
data <- files %>% map(read_csv) %>% reduce(rbind)  

# summarize
res <- data %>% select(-rep) %>% gather(metric, value, -grid_id) %>% na.omit() %>% group_by(grid_id, metric) %>% 
	summarize(nsim = length(metric[!is.na(value)]), min = min(value, na.rm=T), max = max(value, na.rm=T), mean = mean(value, na.rm=T), median = median(value, na.rm=T), sd = sd(value, na.rm=T), mad = mad(value, na.rm=T)) %>% 
	ungroup %>% filter(nsim > 0) 

# save the results
write_csv(res, "results/REV_downsample_gatti.csv")

# remove temp directory
unlink("gatti_temp", recursive = TRUE)


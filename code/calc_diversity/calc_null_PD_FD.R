rm(list = ls())
gc()

library(geometry)
library(doParallel)
library(vegan)
library(picante)
library(feather)
library(tidyverse)


# register the parallel cores
nthreads <- round(detectCores()/2)
registerDoParallel(nthreads)

# raoq formula. assumes equal abundances
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
		outf <- paste0("results/REV_null_results_ANGIO", save_suffix, ".csv")
	}else if(gg == 3){
		print("fitting GYMNOS")
		comm <- comm0 %>% left_join(ag, by = "accepted_bin") %>% 
			filter(group == "Gymnosperms")
		outf <- paste0("results/REV_null_results_GYMNO", save_suffix, ".csv")
	}else if(gg == 4){
		print("removing INVASIVES")
		comm <- comm0 %>% filter(!accepted_bin%in%inv$x)
		outf <- paste0("results/REV_null_results_INV", save_suffix, ".csv")		
	}else{
		print("fitting ALL")
		comm <- comm0
		outf <- paste0("results/REV_null_results", save_suffix, ".csv")
	}
	
	# get the plots and shuffle
	plts <- sample(unique(comm$grid_id), length(unique(comm$grid_id)), replace = FALSE)
	
	# get no. of traits 
	K <- ncol(dt_mat)-1
	
	# get the number of species per plot
	nspp <- comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))
	
	# get the unique species
	spp <- comm %>% select(accepted_bin) %>% distinct() %>% unlist()
	
	# get the unique number of species per plot
	unique_n <- nspp %>% select(n) %>% filter(n>1) %>% distinct() %>% unlist() %>% sort() %>% as.numeric()
	
	# set a random seed
	set.seed(sample(1:1e5, 1))
	
	# number of null sums to run 
	targ <- 100
	
	# create a sequence of numbers and shuffly
	nseq <- rep(unique_n, each = targ)
	nseq <- sample(nseq, length(nseq), replace = FALSE)
	
	# create a sequence of bounds to assign to each thread
	bnds <- sort(unique(c(seq(1, length(nseq), by = 100),length(nseq)+1)))
	
	# create temp directory
	dir.create("nullsims_temp", showWarnings = FALSE)
	
	###################################################
	# get the numm fd/pd, in parallel chunks
	###################################################
	outmat <- foreach(i = 1:(length(bnds)-1), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{
		
		set.seed(sample(1:1e5, 1))
		
		# get the range for this chunk
		start <- bnds[i]
		end <- bnds[i+1]-1
		
		outv <- NULL
		for(j in start:end){
			
			# print(j)
			n <- nseq[j]
			
			# initialize the metrics
			my_raoq <- my_raoq_rand <- my_fric <- my_fric_rand <- my_mpd <- my_faiths <- NA
			
			# sample from the unique species
			my_spp <- sample(spp, n, replace = FALSE)
			
			# get phylo metrics
			my_comm <- matrix(1, nrow = 1, ncol = length(my_spp)) %>% data.frame() %>% as_tibble() %>% setNames(my_spp) %>% as.matrix()
			my_mpd <- picante::mpd(samp = my_comm, dis = cophenetic(ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, my_spp))))
			my_faiths <- as.numeric(picante::pd(samp = my_comm, pruned_tree, include.root = FALSE)[1])
			
			# get the trait matrix
			my_dt <- dt_mat %>% filter(accepted_bin%in%my_spp) %>% select(-accepted_bin)
			
			# get the trait matrix
			my_tr <- (my_dt %>% as.matrix()) 
			my_tr <- my_tr + matrix(runif(nrow(my_tr)*ncol(my_tr), -1e-8,1e-8), nrow = nrow(my_tr), ncol = ncol(my_tr))
			
			# get raoq
			my_raoq <- pairdist(as.matrix(dist(my_tr)))
			
			# get functional richness
			if(nrow(unique(my_tr)) > K){
				my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
			}
			
			outv <- rbind(outv, c("n" = n, "raoq" = my_raoq, "mpd" = my_mpd, "pd" = my_faiths, "fd" = my_fric, "fdr" = my_fric^(1/K))) 
			
		}
		
		# write the final time
		if(!is.null(outv)){
			write_csv(outv %>% as_tibble(), paste0("nullsims_temp/sim",i,".csv"))
		}
		
		return(NA)
	}
	
	
	# read in the files and combine
	files <- paste0("nullsims_temp/", list.files("nullsims_temp"))
	data <- files %>% map(read_csv) %>% reduce(rbind)  
	
	res <- data %>% gather(metric, value, -n) %>% group_by(n, metric) %>% summarize(mean = mean(value, na.rm=T), sd = sd(value, na.rm=T), median = median(value, na.rm=T), mad = mad(value, na.rm=T), min = min(value, na.rm=T), max = max(value, na.rm=T)) %>% 
		ungroup %>% 
		mutate(mean = ifelse(is.infinite(mean) | is.nan(mean), NA, mean),
			   min = ifelse(is.infinite(min) | is.nan(min), NA, min),
			   max = ifelse(is.infinite(max) | is.nan(max), NA, max),
			   median = ifelse(is.infinite(median) | is.nan(median), NA, median),
			   sd = ifelse(is.infinite(sd) | is.nan(sd), NA, sd),
			   mad = ifelse(is.infinite(mad) | is.nan(mad), NA, mad))
	
	# save the results
	write_csv(res, outf)
	
	
	unlink("nullsims_temp", recursive = TRUE)
}
	
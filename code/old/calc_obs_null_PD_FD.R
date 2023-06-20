rm(list = ls())
# library(plyr)
library(geometry)
library(doParallel)
library(vegan)
library(picante)
# library(GGally)
library(feather)
# library(arrow)
library(tidyverse)


registerDoParallel(72)


setwd("~/Git/Tree_diversity/simulations/")


# specify the traits
trait_names <- c("Wood_density", 
				 "Specific_leaf_area", 
				 "Leaf_P_per_mass", 
				 "Stem_conduit_diameter", 
				 "Tree_height", 
				 "Seed_dry_mass", 
				 "Bark_thickness", 
				 "Root_depth")

K <- k <- length(trait_names)

# specify monocot orders to make sure none are present
monocots <- c("Acorales", "Alismatales", "Arecales", "Asparagales", "Commelinales", "Dioscoreales", "Liliales", "Pandanales", "Petrosaviales", "Poales", "Zingiberales")
ferns <- c("Lophosoria", "Metaxya", "Sphaeropteris", "Alsophila", "Nephelea", "Trichipteris", "Cyathea" ,"Cnemidaria", "Dicksonia", "Cystodium", "Thyrsopteris", "Culcita", "Cibotium")

# read in the list of angio vs gymno species
ag0 <- read_csv("~/Git/tree_traits/tree_traits_backup_with_names/data/raw_data/ANGIO_GYMNO_lookup.csv") 

# read in the community data
commo <- feather::read_feather("raw_data/NoMono_matchPD_FD_gather.feather") %>% 
	mutate(accepted_bin = trimws(gsub("_", " ", accepted_bin), "both"))

# read in the synonyms
syns <- read_csv("raw_data/matched_trait_names.csv") 

# read in the trait table and subset
dto <- feather::read_feather("raw_data/Estimated_trait_table.feather") %>% 
	mutate(trait = gsub(" ", "_", trait)) %>% 
	filter(trait%in%trait_names) %>% 
	select(accepted_bin, trait, value) %>% 
	spread(trait, value) 

# get the synonyms
new_syns <- commo %>% select(accepted_bin) %>% distinct() %>% filter(!accepted_bin%in%dto$accepted_bin) %>% rename(raw_name = accepted_bin) %>% left_join(syns) %>% 
	filter(raw_name!=trait_name) %>% filter(trait_name%in%dto$accepted_bin) %>% 
	rename(new_name = raw_name, accepted_bin = trait_name) %>% select(-match_type) %>% distinct()

# get the multiples and pick the closest
mults <- new_syns %>% group_by(accepted_bin) %>% tally() %>% filter(n>1)
new_syns <- new_syns %>% filter(!accepted_bin%in%mults$accepted_bin) %>%
	bind_rows(new_syns %>% filter(accepted_bin%in%mults$accepted_bin) %>% rowwise() %>% mutate(d = adist(new_name, accepted_bin)[,1]) %>% ungroup %>% 
	group_by(accepted_bin) %>% summarize(new_name = new_name[d == min(d)][1]) %>% ungroup)

# make sure zero mults
new_syns %>% group_by(accepted_bin) %>% tally() %>% filter(n>1) %>% nrow()

# add in new names 
dt <- dto %>% 
	left_join(new_syns) %>% mutate(accepted_bin = ifelse(!is.na(new_name), new_name, accepted_bin)) %>% select(-new_name) %>%
	mutate(tax_genus = word(accepted_bin, 1)) %>% 
	left_join(ag0 %>% select(accepted_bin, order, group, genus, family)) %>% 
	filter(!order%in%monocots, !genus%in%ferns, !tax_genus%in%ferns, !family%in%c("Osmundaceae"), group%in%c("Angiosperms", "Gymnosperms")) %>% 
	select(-order, -genus, -group, -family, -tax_genus)


# remove ferns/palms
comm <- commo %>%
	mutate(tax_genus = word(accepted_bin, 1)) %>% 
	left_join(ag0 %>% select(accepted_bin, order, group, genus, family)) %>% 
	filter(!order%in%monocots, !genus%in%ferns, !tax_genus%in%ferns, !family%in%c("Osmundaceae"), group%in%c("Angiosperms", "Gymnosperms")) %>% 
	select(grid_id, present, accepted_bin)


# see the mismatches
sum(!dt$accepted_bin%in%comm$accepted_bin)
unique(comm$accepted_bin[!comm$accepted_bin%in%dt$accepted_bin])

# subset the matched species and remove any plots with singletons
comm <- comm %>% filter(accepted_bin%in%dt$accepted_bin) %>% distinct() %>% 
	left_join(comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))) %>% 
	filter(n>1) %>% select(-n)

# subset the matched species, and summarize any duplicates
dt <- dt %>% filter(accepted_bin%in%comm$accepted_bin) %>% distinct() %>% 
	group_by(accepted_bin) %>% summarize_all(.funs = mean)

# double check
sum(!dt$accepted_bin%in%comm$accepted_bin)
unique(comm$accepted_bin[!comm$accepted_bin%in%dt$accepted_bin])

# get the unique species
spp <- comm %>% select(accepted_bin) %>% distinct() %>% unlist()

# get the number fo species per plot
nspp <- comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))

# get the plots with at least one species
plts <- nspp %>% filter(n>1) %>% select(grid_id) %>% distinct() %>% unlist() %>% as.numeric()

# get the PCA of traits
dt_mat <- rda(dt %>% select(-accepted_bin), scale = TRUE)$CA$u %>% data.frame() %>% as_tibble() %>% 
	mutate(accepted_bin = dt$accepted_bin) %>% 
	gather(PC, value, -accepted_bin) %>% 
	group_by(PC) %>% 
	mutate(value = (value - mean(value))/sd(value)) %>% 
	ungroup %>% 
	spread(PC, value)
	
# get the tree and clean the names
pruned_tree <- read.tree("raw_data/no_monocots_tree.nwk")
pruned_tree$tip.label <- gsub("_", " ", pruned_tree$tip.label)

# double check
sum(!dt$accepted_bin%in%pruned_tree$tip.label)
sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label)
sum(!dt_mat$accepted_bin%in%comm$accepted_bin)
sum(!comm$accepted_bin%in%dt_mat$accepted_bin)
sum(!unique(comm$accepted_bin)%in%pruned_tree$tip.label)
sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label)

# get the range of richness values
minr <- min(nspp$n)
maxr <- max(nspp$n)

K <- ncol(dt_mat)-1

pairdist <- function(myd, myp){
	# only works if p is 1/0 or TRUE/FALSE
	rows <- which(myp==1)
	wpd <- (myd[rows,rows] / sum(myp))^2
	rq <- sum(wpd)
	return(rq) 
}

# shuffle the order
plts <- sample(plts, length(plts), replace = FALSE)

# get the observed fd/pd
outmat_check <- foreach(j = 1:length(plts), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{

	print(j)

	# get the focal plot
	my_dt <- comm %>% filter(grid_id == plts[j])

	# initialize the variables
	my_dt_sub <- NULL
	pd_cur <- mpd_cur <- pd_sub <- mpd_sub <- my_raoq <- my_fric <- NA
	if(nrow(my_dt)>1){
		# subset the data and get the phyometrics
		my_comm <- matrix(1, nrow = 1, ncol = nrow(my_dt)) %>% data.frame() %>% as_tibble() %>% setNames(my_dt$accepted_bin) %>% as.matrix()
		pd_cur <- as.numeric(picante::pd(samp = my_comm, pruned_tree, include.root = FALSE)[1])
		mpd_cur <- picante::mpd(samp = my_comm, dis = cophenetic(ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, my_dt$accepted_bin))))

		# get the trait matrix
		my_tr <- (dt_mat %>% filter(accepted_bin%in%my_dt$accepted_bin) %>% select(-accepted_bin) %>% as.matrix())
		my_tr <- my_tr + matrix(runif(nrow(my_tr)*ncol(my_tr), -1e-8,1e-8), nrow = nrow(my_tr), ncol = ncol(my_tr))

		# get raoq
		my_raoq <- pairdist(as.matrix(dist(my_tr)), rep(1, nrow(my_tr)))

		# make sure we have more unique traits than rows
		if(nrow(unique(my_dt)) > ncol(my_tr)){
			my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
		}

	}
	return(c("grid_id" = plts[j], "n" = nrow(my_dt), "pd" = pd_cur, "mpd" = mpd_cur, "raoq" = my_raoq, "fric" = my_fric, "fric_rad" = my_fric^(1/K)))
}



outmat_check <- outmat_check %>% data.frame() %>% as_tibble() #%>% left_join(metrics %>% setNames(paste0("obs_", names(.))) %>% rename(grid_id = obs_grid_id))
write_csv(outmat_check, "results/obs_results.csv")



set.seed(10)

nseq <- nspp$n #seq(2, maxr, by = 1)
nseq <- sample(nseq, length(nseq), replace = FALSE)
neach <- 1

# rm(outmat_check)
# gc()

for(k in 1:50){
	outmat <- foreach(i = 1:length(nseq), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{
		
		print(i)
		n <- nseq[i]
		
		outv <- NULL

		for(j in 1:1){
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
			my_tr <- (my_dt %>% as.matrix()) + matrix(runif(nrow(my_dt)*ncol(my_dt), -1e-8,1e-8), nrow = nrow(my_dt), ncol(my_dt))
			my_tr_rand <- apply(my_dt %>% as.matrix(), 2, function(x) sample(x, n, replace = TRUE)) + matrix(runif(nrow(my_dt)*ncol(my_dt), -1e-8,1e-8), nrow = nrow(my_dt), ncol(my_dt))
			
			# get raoq
			my_raoq <- pairdist(as.matrix(dist(my_tr)), rep(1, nrow(my_tr)))
			my_raoq_rand <- pairdist(as.matrix(dist(my_tr_rand)), rep(1, nrow(my_tr_rand)))
			
			# get functional richness
			if(nrow(unique(my_tr)) > K){
				my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
			}
			if(nrow(unique(my_tr_rand)) > K){
				my_fric_rand <- tryCatch(convhulln(my_tr_rand, "FA")$vol, error = function(e) return(NA))
			}
			
			outv <- rbind(outv, c("n" = n, "raoq" = my_raoq, "raoq_rand" = my_raoq_rand, "mpd" = my_mpd, "pd" = my_faiths, "fric" = my_fric, "fric_rad" = my_fric^(1/K), "fric_rand" = my_fric_rand, "fric_rand_rad" = my_fric_rand^(1/K)))
		}
		
		return(outv)
	}
	
	#
	outmat <- outmat %>% data.frame() %>% as_tibble() %>%
		mutate(fric_vol = fric^K, rep = k)
	#
	write_csv(outmat, "results/null_results_each.csv", append= file.exists("results/null_results_each.csv"))
	
	rm(outmat)
	gc()
}

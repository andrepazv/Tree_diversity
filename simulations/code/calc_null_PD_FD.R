rm(list = ls())
gc()
# library(plyr)
library(geometry)
library(doParallel)
library(vegan)
library(picante)
# library(GGally)
library(feather)
# library(arrow)
library(tidyverse)


setwd("~/Git/Tree_diversity/simulations/")


# specify the traits
trait_names <- c("Wood_density", # 4
				 "Specific_leaf_area", #3117 
				 "Leaf_P_per_mass", # 15 
				 "Stem_conduit_diameter", #281 
				 "Tree_height", #3106
				 "Seed_dry_mass", #26
				 "Bark_thickness",#, #24
				 "Root_depth") #6

K <- k <- length(trait_names)

nthreads <- 72
registerDoParallel(nthreads)

# # specify monocot orders to make sure none are present
monocots <- c("Acorales", "Alismatales", "Arecales", "Asparagales", "Commelinales", "Dioscoreales", "Liliales", "Pandanales", "Petrosaviales", "Poales", "Zingiberales")
ferns <- c("Lophosoria", "Metaxya", "Sphaeropteris", "Alsophila", "Nephelea", "Trichipteris", "Cyathea" ,"Cnemidaria", "Dicksonia", "Cystodium", "Thyrsopteris", "Culcita", "Cibotium")

# read in the list of angio vs gymno species
ag0 <- read_csv("~/Git/tree_traits/tree_traits_backup_with_names/data/raw_data/ANGIO_GYMNO_lookup.csv")

# read in the community data
commo <- feather::read_feather("raw_data/NoMono_matchPD_FD_gather.feather") %>% 
	mutate(accepted_bin = trimws(gsub("_", " ", accepted_bin), "both"))

# read in the synonyms
syns <- read_csv("raw_data/matched_trait_names.csv") %>% 
	filter(!duplicated(trait_name), !duplicated(raw_name))

# # read in the trait table and subset
dt <-
	# feather::read_feather("raw_data/Estimated_trait_table.feather") %>%
	read_csv("~/Git/dispersion/data/combined_datasets/Species_trait_table_50_composite_with_monos.csv") %>% left_join(read_csv("~/Git/dispersion/data/trait_models/Trait_names.csv") %>% select(TraitID, trait)) %>%
	mutate(trait = gsub(" ", "_", trait)) %>%
	filter(trait%in%trait_names) %>%
	select(accepted_bin, trait, value) %>%
	spread(trait, value)


# dto <- read_csv("raw_data/traits_used_FD.csv") 
# names(dto)[1] <- "accepted_bin"
# dto <- dto %>% mutate(accepted_bin = gsub("\\.", " ", gsub("_", " ", accepted_bin))) %>%
# 	mutate(accepted_bin = paste(word(accepted_bin, 1), word(accepted_bin, 2)))
# 
# dt1 <- dt1 %>% setNames(gsub("_", "\\.", names(.))) 
# names(dt1) <- "accepted_bin"
# 
# dto %>% left_join()

# # get the synonyms
# new_syns <- commo %>% select(accepted_bin) %>% distinct() %>% filter(!accepted_bin%in%dto$accepted_bin) %>% rename(raw_name = accepted_bin) %>% left_join(syns) %>%
# 	filter(raw_name!=trait_name) %>% filter(trait_name%in%dto$accepted_bin) %>%
# 	rename(new_name = raw_name, accepted_bin = trait_name) %>% select(-match_type) %>% distinct()
# 
# # get the multiples and pick the closest
# mults <- new_syns %>% group_by(accepted_bin) %>% tally() %>% filter(n>1)
# new_syns <- new_syns %>% filter(!accepted_bin%in%mults$accepted_bin) %>%
# 	bind_rows(new_syns %>% filter(accepted_bin%in%mults$accepted_bin) %>% rowwise() %>% mutate(d = adist(new_name, accepted_bin)[,1]) %>% ungroup %>%
# 			  	group_by(accepted_bin) %>% summarize(new_name = new_name[d == min(d)][1]) %>% ungroup)
# 
# # make sure zero mults
# new_syns %>% group_by(accepted_bin) %>% tally() %>% filter(n>1) %>% nrow()

monocot_genera <- read_csv("raw_data/monocot_genera.csv")

# monocots<-read_csv("Species_trait_values_3_11_22.csv")

# add in new names 
dt <- dto %>% 
	# left_join(new_syns) %>% mutate(accepted_bin = ifelse(!is.na(new_name), new_name, accepted_bin)) %>% select(-new_name) %>%
	mutate(tax_genus = word(accepted_bin, 1)) %>% 
	filter(!tax_genus%in%monocot_genera$tax_genus) %>% 
	select(accepted_bin, tax_genus, all_of(trait_names)) %>% 
	left_join(ag0 %>% select(accepted_bin, order, group, genus, family)) %>%
	filter(!order%in%monocots, !genus%in%ferns, !tax_genus%in%ferns, !family%in%c("Osmundaceae"), group%in%c("Angiosperms", "Gymnosperms")) %>%
	select(accepted_bin, all_of(trait_names)) 



# get the PCA of traits
dt_mat <- rda(dt %>% select(-accepted_bin), scale = TRUE)$CA$u %>% data.frame() %>% as_tibble() %>% 
	mutate(accepted_bin = dt$accepted_bin) %>% 
	gather(PC, value, -accepted_bin) %>%
	group_by(PC) %>%
	mutate(value = (value - mean(value))/sd(value)) %>%
	ungroup %>%
	spread(PC, value)

rm(dt)
rm(dto)

for(i in 1:nrow(syns)){
	dt_mat$accepted_bin[which(dt_mat$accepted_bin == syns$trait_name[i])]<-syns$raw_name[i]
}


# remove ferns/palms
comm <- commo %>%
	# mutate(tax_genus = word(accepted_bin, 1)) %>% 
	# filter(!tax_genus%in%monocot_genera$tax_genus) %>% #%>% 
	# filter(!order%in%monocots, !genus%in%ferns, !tax_genus%in%ferns, !family%in%c("Osmundaceae")) %>%#, group%in%c("Angiosperms", "Gymnosperms")) %>%
	select(grid_id, present, accepted_bin)


# see the mismatches
sum(!dt_mat$accepted_bin%in%comm$accepted_bin)
unique(comm$accepted_bin[!comm$accepted_bin%in%dt_mat$accepted_bin])

# subset the matched species and remove any plots with singletons
comm <- comm %>% filter(accepted_bin%in%dt_mat$accepted_bin) %>% distinct() %>% 
	left_join(comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))) %>% 
	# filter(n>1) %>% 
	select(-n)

# # subset the matched species, and summarize any duplicates
dt_mat <- dt_mat %>% filter(accepted_bin%in%comm$accepted_bin) %>% distinct()# %>%
# 	group_by(accepted_bin) %>% summarize_all(.funs = mean)

# double check
sum(!dt_mat$accepted_bin%in%comm$accepted_bin)
unique(comm$accepted_bin[!comm$accepted_bin%in%dt_mat$accepted_bin])

# get the unique species
spp <- comm %>% select(accepted_bin) %>% distinct() %>% unlist()

# get the number fo species per plot
nspp <- comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))

# get the plots
plts <- nspp %>% select(grid_id) %>% distinct() %>% unlist() %>% as.numeric()

# em <- read_csv("raw_data/empty_pixels.csv") %>% select(Grid) %>% rename(grid_id = Grid)

# plts <- plts %>% filter(grid_id%in%em$grid_id) %>% distinct() %>% unlist() %>% as.numeric()

# get the tree and clean the names
pruned_tree <- read.tree("raw_data/no_monocots_tree.nwk")
pruned_tree$tip.label <- gsub("_", " ", pruned_tree$tip.label)

# double check
sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label)
sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label)
sum(!dt_mat$accepted_bin%in%comm$accepted_bin)
sum(!comm$accepted_bin%in%dt_mat$accepted_bin)
sum(!unique(comm$accepted_bin)%in%pruned_tree$tip.label)
sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label)

# get the range of richness values
minr <- min(nspp$n)
maxr <- max(nspp$n)

K <- ncol(dt_mat)-1

pairdist <- function(myd){
	# only works if p is 1/0 or TRUE/FALSE
	# rows <- which(myp==1)
	return(sum((myd/ nrow(myd))^2))
}

# pairdist2 <- function(myd, myp){
# 	# only works if p is 1/0 or TRUE/FALSE
# 	rows<-which(myp==1)
# 	myd <- myd[rows,rows]^2
# 	n<-sum(myp)
# 	n2 <- sum(myp)^2
# 	# get the weighted pairwise dist
# 	wpd <- myd / n2
# 	# for rao's q, just sum up
# 	rq <- sum(wpd)
# 	# calculate abundance weighted pairwise distance 
# 	wmpd <- sum(wpd[upper.tri(wpd)]) / ((n*(n-1)/2)*(1/n2))
# 	# calculate standard unweighted species pairwise distance
# 	mpd <- mean(myd[upper.tri(myd)])
# 	return(c("raoq" = rq, "wmpd" = wmpd, "mpd" = mpd))
# }


# distance <- as.matrix(dist(dt_mat %>% select(-accepted_bin) %>% as.matrix()))
# rownames(distance) <- colnames(distance) <- dt_mat$accepted_bin

# 
# shuffle the order
plts <- sample(plts, length(plts), replace = FALSE)

# j <- which(plts == 11707)

# comm1 <- read_csv("~/Git/Tree_diversity/simulations/raw_data/comm_matrix_andrea_script.csv")

DO_SUBSAMPLE <- FALSE

nrep <- ifelse(DO_SUBSAMPLE, 50, 1)

# get the observed fd/pd
outmat_check <- foreach(j = 1:length(plts), .inorder = FALSE, .combine = bind_rows, .multicombine = TRUE)%dopar%{

	# print(j)


	# get the focal plot
	my_dt0 <- comm %>% filter(grid_id == plts[j])

	# initialize the variables
	my_dt_sub <- NULL
	pd_cur <- mpd_cur <- pd_sub <- mpd_sub <- my_raoq <- my_fric <- NA

	tmp_res <- c()

	if(nrow(my_dt0)>1){

		for(z in 1:nrep){

			if(DO_SUBSAMPLE){
				my_dt <- my_dt0 %>% sample_n(floor(nrow(my_dt0)/2))
			}else{
				my_dt <- my_dt0
			}

			if(nrow(my_dt) > 1){
				# subset the data and get the phyometrics
				my_comm <- matrix(1, nrow = 1, ncol = nrow(my_dt)) %>% data.frame() %>% as_tibble() %>% setNames(my_dt$accepted_bin) %>% as.matrix()
				pd_cur <- as.numeric(picante::pd(samp = my_comm, pruned_tree, include.root = FALSE)[1])
				mpd_cur <- picante::mpd(samp = my_comm, dis = cophenetic(ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, my_dt$accepted_bin))))
				#
				# # get the trait matrix
				my_tr <- dt_mat %>% filter(accepted_bin%in%my_dt$accepted_bin) %>% select(-accepted_bin) %>% as.matrix()
				my_tr <- my_tr + matrix(runif(nrow(my_tr)*ncol(my_tr), -1e-8,1e-8), nrow = nrow(my_tr), ncol = ncol(my_tr))

				# my_raoq <- pairdist2(distance, rownames(distance) %in% my_dt$accepted_bin)["raoq"]
				# get raoq
				my_raoq <- pairdist(as.matrix(dist(my_tr))) #, rep(1, nrow(my_tr)))

				# # make sure we have more unique traits than rows
				if(nrow(unique(my_dt)) > ncol(my_tr)){
					my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
				}

				tmp_res <- rbind(tmp_res, c("grid_id" = plts[j], "n" = nrow(my_dt), "pd" = pd_cur, "mpd" = mpd_cur, "raoq" = my_raoq, "fd" = my_fric, "fdr" = my_fric^(1/K)))

			}else{
				my_raoq <- pd_cur <- mpd_cur <- 0
			}
		}
	}else{
		my_raoq <- pd_cur <- mpd_cur <- 0
	}

	out_res <- tmp_res %>% data.frame() %>% as_tibble() %>% group_by(grid_id, n) %>% summarize_all(mean) %>% ungroup %>% mutate(fit = "mean") %>%
		bind_rows(tmp_res %>% data.frame() %>% as_tibble() %>% group_by(grid_id, n) %>% summarize_all(sd) %>% ungroup %>% mutate(fit = "sd"))

	return(out_res)
}
#
#
#
outmat_check <- outmat_check %>% data.frame() %>% as_tibble() #%>% left_join(metrics %>% setNames(paste0("obs_", names(.))) %>% rename(grid_id = obs_grid_id))
write_csv(outmat_check, "results/obs_results_50vars.csv")

rm(outmat_check)

# 
# # 
targ <- 100

unique_n <- nspp %>% select(n) %>% distinct() %>% unlist() %>% as.numeric()

# read in the existing data and get the number we still need
if(file.exists("results/null_results_50var.csv")){
	sim <- read_csv("results/null_results_50var.csv")
	sim <- sim %>% mutate(row = 1:nrow(.)+1) %>% filter(!((n>1 & (is.na(raoq)) | is.na(mpd) | is.na(pd)) | n > 2500))
	# write_csv(sim, "results/null_results_fixed.csv")
	ndf <- tibble(n = unique_n) %>% left_join(sim %>% group_by(n) %>% tally()) %>%
		mutate(need = ifelse(is.na(nn), targ, targ - nn)) %>% arrange(desc(need)) %>%
		filter(need > 0)
}else{
	ndf <- tibble(n = unique_n) %>% mutate(nn = 0, need = targ)
}

if(nrow(ndf) == 0){
	stop("NO MORE NEEDED")
}

# set a random seed
set.seed(sample(1:1e5, 1))

# create a sequence of numbers and shuffly
nseq <- rep(ndf$n, ndf$need)
nseq <- sample(nseq, length(nseq), replace = FALSE)

# create a sequence of bounds to assign to each thread
bnds <- sort(unique(c(round(seq(1, length(nseq), length = nthreads - 3)), length(nseq)+1)))

# how often to print
print_each <- 5

outmat <- foreach(i = 1:(length(bnds)-1), .inorder = FALSE, .combine = rbind, .multicombine = TRUE)%dopar%{

	set.seed(sample(1:1e5, 1))
	start <- bnds[i]
	end <- bnds[i+1]-1

	outv <- NULL
	for(j in start:end){
		# print(j)
		n <- nseq[j]

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
		my_tr <- (my_dt %>% as.matrix()) # + matrix(runif(nrow(my_dt)*ncol(my_dt), -1e-8,1e-8), nrow = nrow(my_dt), ncol(my_dt))
		my_tr <- my_tr + matrix(runif(nrow(my_tr)*ncol(my_tr), -1e-8,1e-8), nrow = nrow(my_tr), ncol = ncol(my_tr))

		# my_tr_rand <- apply(my_dt %>% as.matrix(), 2, function(x) sample(x, n, replace = TRUE)) + matrix(runif(nrow(my_dt)*ncol(my_dt), -1e-8,1e-8), nrow = nrow(my_dt), ncol(my_dt))
		#
		# get raoq
		my_raoq <- pairdist(as.matrix(dist(my_tr)))
		# my_raoq_rand <- pairdist(as.matrix(dist(my_tr_rand)), rep(1, nrow(my_tr_rand)))
		#
		# get functional richness
		if(nrow(unique(my_tr)) > K){
			my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
		}
		# if(nrow(unique(my_tr_rand)) > K){
		# 	my_fric_rand <- tryCatch(convhulln(my_tr_rand, "FA")$vol, error = function(e) return(NA))
		# }

		outv <- rbind(outv, c("n" = n, "raoq" = my_raoq, "mpd" = my_mpd, "pd" = my_faiths, "fd" = my_fric, "fdr" = my_fric^(1/K))) #, "fric_rand" = my_fric_rand, "fric_rand_rad" = my_fric_rand^(1/K)))

		if(j%%print_each == 0){
			# print("appendin")
			write_csv(outv %>% as_tibble(), "results/null_results_50var.csv", append= file.exists("results/null_results_50var.csv"))
			outv <- NULL
		}
	}

	# write the final time
	if(!is.null(outv)){
		write_csv(outv %>% as_tibble(), "results/null_results_50var.csv", append= file.exists("results/null_results_50var.csv"))
	}

	return(NA)
}

# 
# 
# # # # ggpairs(outmat, lower = list(continuous = "points"), upper = list(continuous = "points"))
# # # # 
# # # # # 
# obs <- read_csv("results/obs_results_fixed.csv") %>% rename(fd = fric, fdr = fric_rad)
# obs2 <- read_csv("results/obs_results_fixed_zerosfd.csv") %>% rename(fd = fric, fdr = fric_rad)


# # # # # 
# # # # # # 
# # # # # obs %>% ggplot(aes(x = mpd, y = raoq))+geom_point()
# # # # # obs %>% ggplot(aes(x = pd, y = fd))+geom_point()
# # # # # obs %>% ggplot(aes(x = pd, y = fdr))+geom_point()
# # # # # obs %>% ggplot(aes(x = n, y = fdr))+geom_point()
# # # # # obs %>% ggplot(aes(x = n, y = fd))+geom_point()
# # # # # obs %>% ggplot(aes(x = n, y = pd))+geom_point()
# # # # # obs %>% ggplot(aes(x = n, y = raoq))+geom_point()
# # # # # obs %>% ggplot(aes(x = n, y = mpd))+geom_point()
# # # # 
# sim <- read_csv("results/null_results_fixed.csv") %>% rename(fd = fric, fdr = fric_rad)
# sim <- sim %>% filter(!((n>1 & (is.na(raoq)) | is.na(mpd) | is.na(pd)) | n > 2500))
# # # # #
# # sim %>% ggplot(aes(x = mpd, y = raoq))+geom_point()
# # sim %>% ggplot(aes(x = pd, y = fd))+geom_point()
# # sim %>% ggplot(aes(x = pd, y = fdr))+geom_point()
# # sim %>% ggplot(aes(x = n, y = fd))+geom_point()
# # sim %>% ggplot(aes(x = n, y = fd))+geom_point()
# # sim %>% ggplot(aes(x = n, y = pd))+geom_point()
# # sim %>% ggplot(aes(x = n, y = raoq))+geom_point()
# # sim %>% ggplot(aes(x = n, y = mpd))+geom_point()
# # # # 
# # # # 
# sim_means <- sim %>% gather(metric, value, -n) %>% group_by(n, metric) %>% summarize(null_mean = mean(value, na.rm = T), null_sd = sd(value, na.rm=T)) %>% ungroup
# # #
# obsg <- obs %>% gather(metric, value, -n, -grid_id) %>% left_join(sim_means) %>%
# 	mutate(diff = value - null_mean, z = (value - null_mean)/null_sd, prop_diff = (value - null_mean)/null_mean) %>%
# 	arrange(grid_id, n) %>% filter(!is.na(value))
# #
# write_csv(obsg, "results/obs_null_z_scores_fixed.csv")
# # # # 
# # # # 
# obs_z <- obsg %>% select(grid_id, n, metric, z) %>% spread(metric, z)
# # # # 
# # # # 
# # # 
# obs_z %>% ggplot(aes(x = n, y = raoq))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# obs_z %>% ggplot(aes(x = n, y = mpd))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# obs_z %>% ggplot(aes(x = mpd, y = raoq))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 1)
# # # 
# # # 
# obs_z %>% ggplot(aes(x = n, y = fd))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# obs_z %>% ggplot(aes(x = n, y = pd))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# obs_z %>% ggplot(aes(x = pd, y = fd))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# 
# obs_z
# # 
# # # 
# # # 
# # # all <- obs %>% left_join(sim %>% group_by(n) %>% summarize_all(.funs = c(sd, mean), na.rm=T) %>% mutate(fdr2 = fd^(1/8), fdr_rand = fd_rand^(1/8)) %>% 
# # # 						 	setNames(paste0("null_", names(.))) %>% rename(n = null_n))
# # # 
# # # all %>% ggplot(aes(x = null_pd, y = pd))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth(method = "lm", formula = 'y ~ -1 + x')
# # # all %>% ggplot(aes(x = null_fd, y = fd))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth()
# # # all %>% ggplot(aes(x = null_fdr, y = fdr))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth()
# # # all %>% ggplot(aes(x = null_raoq, y = raoq))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth()
# # # all %>% ggplot(aes(x = null_mpd, y = mpd))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth()
# # 
# # # res <- read_csv("For_null_models/null_results5k_new.csv") %>% 
# # # 	gather(trait, value, -n) %>% 
# # # 	group_by(n, trait) %>% 
# # # 	summarize(nsim = length(value[!is.na(value)]), sd = sd(value, na.rm=T), value = mean(value, na.rm = T)) %>% 
# # # 	ungroup %>% 
# # # 	mutate(se = sd / sqrt(nsim)) %>% 
# # # 	filter(!is.na(value)) %>% 
# # # 	rename(null_value = value, null_sd = sd, null_se  = se)
# # # 
# # # 
# # # 
# # # all <- obs %>% left_join(res) %>% 
# # # 	mutate(z = (value - null_value)/null_sd) 
# # # 	
# # # 
# # # ggplot(all %>% filter(trait == "pd"), aes(x = n, y = z))+geom_point()+geom_smooth()
# # # ggplot(all %>% filter(trait == "fric"), aes(x = n, y = z))+geom_point()+geom_smooth()
# # # ggplot(all %>% filter(trait == "raoq"), aes(x = n, y = z))+geom_point()+geom_smooth()
# # # ggplot(all %>% filter(trait == "mpd"), aes(x = n, y = z))+geom_point()+geom_smooth()
# # # 
# # # 
# # # 
# # # all_gath <- all %>% select(grid_id, n, trait, value, null_value) %>% gather(fit, val, -grid_id, -n, -trait) 
# # # 
# # # all_gath %>% filter(trait == "pd") %>% group_by(n, fit) %>% summarize(val = mean(val)) %>% ungroup %>% ggplot(aes(x = n, y = val, color = fit))+geom_point()
# # # all_gath %>% filter(trait == "fric") %>% group_by(n, fit) %>% summarize(val = mean(val)) %>% ungroup %>% ggplot(aes(x = n, y = val, color = fit))+geom_point()
# # # 
# # # all %>% filter(trait%in%c("raoq", "mpd")) %>% select(grid_id, n, trait, z) %>% spread(trait, z) %>% ggplot(aes(x = mpd, y = raoq))+geom_point()+geom_abline(b = 0, a = 1)
# # # all %>% filter(trait%in%c("raoq", "mpd")) %>% select(grid_id, n, trait, value) %>% spread(trait, value) %>% ggplot(aes(x = mpd, y = raoq))+geom_point()+geom_abline(b = 0, a = 1)
# # # all %>% filter(trait%in%c("raoq", "mpd")) %>% select(grid_id, n, trait, null_value) %>% spread(trait, null_value) %>% ggplot(aes(x = mpd, y = raoq))+geom_point()+geom_abline(b = 0, a = 1)
# # # 
# # # 
# # # 
# # # res %>% filter(trait%in%c("raoq", "mpd")) %>% select(n, trait, null_value) %>% spread(trait, null_value) %>% 
# # # 	ggplot(aes(x = mpd, y = raoq))+geom_point()
# # # 
# # # 
# # # ggplot(all %>% filter(n_sub > 7) %>% rowwise() %>% mutate(pd_sub = plyr::round_any(pd_sub, 100)) %>% group_by(pd_sub, type) %>% summarize(fric_vol = mean(fric_vol)), aes(x = pd_sub, y = fric_vol, color = type))+geom_point()
# # # ggplot(all %>% filter(n_sub > 7) %>% rowwise() %>% mutate(pd = plyr::round_any(pd, 100)) %>% group_by(pd, type) %>% summarize(fric = mean(fric)), aes(x = pd, y = fric, color = type))+geom_point()
# # # ggplot(all %>% rowwise() %>% mutate(mpd = plyr::round_any(mpd, 1)) %>% group_by(mpd, type) %>% summarize(raoq = mean(raoq)), aes(x = mpd, y = raoq, color = type))+geom_point()
# # # 
# # # ggplot(all %>% slice(sample(1:nrow(.), 20000)), aes(x = mpd, y = raoq, color = type))+
# # # 	geom_point()
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000)), aes(x = n, y = raoq, color = type))+
# # # 	geom_point()
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000)) %>% filter(type == "obs"), aes(x = n, y = mpd, color = type))+
# # # 	geom_point()
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000)) , aes(x = n, y = pd, color = type))+
# # # 	geom_point()
# # # 
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000) )%>% filter(n_sub > 7) , aes(x = n_sub, y = fric, color = type))+
# # # 	geom_point()
# # # 
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000) )%>% filter(n_sub > 7) %>% mutate(fric_rand = ifelse(type == "obs", fric_vol, fric_rand^7)), aes(x = n_sub, y = fric_rand, color = type))+
# # # 	geom_point()
# # # 
# # # 
# # # all %>%
# # # 
# # # ggplot(all %>% group_by(n, type) %>% summarize(sd = sd(raoq), raoq = mean(raoq)) %>% ungroup, aes(x = n, y = raoq, color = type))+
# # # 	geom_point()
# # # 
# # # ggplot(all %>% group_by(n, type) %>% summarize(sd = sd(mpd), mpd = mean(mpd)) %>% ungroup, aes(x = n, y = mpd, color = type))+
# # # 	geom_point()
# # # 
# # # 
# # # # lower phylogenetic diversity than expected by change
# # # ggplot(all %>% group_by(n, type) %>% summarize(sd = sd(pd), pd = mean(pd)) %>% ungroup, aes(x = n, y = pd, color = type))+
# # # 	geom_point()
# # # 
# # # ggplot(all %>% filter(n_sub > 7) %>% group_by(n_sub, type) %>% summarize(sd = sd(pd), pd = mean(pd)) %>% ungroup, aes(x = n_sub, y = pd, color = type))+
# # # 	geom_point()
# # # 
# # # # maybe comparable equivalent values
# # # ggplot(all %>% filter(n_sub > 7) %>% group_by(n_sub) %>% mutate(n_sub = plyr::round_any(n_sub, 50)) %>% ungroup %>% group_by(n_sub, type) %>% summarize(sd = sd(fric), fric = mean(fric)) %>% ungroup, aes(x = n_sub, y = fric, color = type))+
# # # 	geom_point()
# # # 
# # # 
# # # ggplot(all %>% filter(n_sub > 7) %>% group_by(n_sub) %>% mutate(n_sub = plyr::round_any(n_sub, 50)) %>% ungroup %>% group_by(n_sub, type) %>% summarize(sd = sd(fric_vol), fric_vol = mean(fric_vol)) %>% ungroup, aes(x = n_sub, y = fric_vol, color = type))+
# # # 	geom_point()
# # # 
# # # 				my_spp <- sample(spp, n, replace = FALSE)
# # # 				
# # # 				# get phylo metrics
# # # 				my_comm <- matrix(1, nrow = 1, ncol = length(my_spp)) %>% data.frame() %>% as_tibble() %>% setNames(my_spp) %>% as.matrix()
# # # 				my_mpd <- picante::mpd(samp = my_comm, dis = cophenetic(ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, my_spp))))
# # # 				my_faiths <- as.numeric(picante::pd(samp = my_comm, pruned_tree, include.root = FALSE)[1])
# # # 				
# # # 				# get the trait matrix
# # # 				my_dt <- dt_mat %>% filter(accepted_bin%in%my_spp) %>% select(-accepted_bin)
# # # 				
# # # 				# get the trait matrix
# # # 				my_tr <- (my_dt %>% as.matrix()) + matrix(runif(nrow(my_dt)*ncol(my_dt), -1e-8,1e-8), nrow = nrow(my_dt), ncol(my_dt))
# # # 				my_tr_rand <- apply(my_dt %>% as.matrix(), 2, function(x) sample(x, n, replace = TRUE)) + matrix(runif(nrow(my_dt)*ncol(my_dt), -1e-8,1e-8), nrow = nrow(my_dt), ncol(my_dt))
# # # 				
# # # 				# get raoq
# # # 				my_raoq <- pairdist(as.matrix(dist(my_tr)), rep(1, nrow(my_tr)))
# # # 				my_raoq_rand <- pairdist(as.matrix(dist(my_tr_rand)), rep(1, nrow(my_tr_rand)))
# # # 				
# # # 				# get functional richness
# # # 				if(nrow(unique(my_tr)) > K){
# # # 					my_fric <- tryCatch(convhulln(my_tr, "FA")$vol, error = function(e) return(NA))
# # # 				}
# # # 				if(nrow(unique(my_tr_rand)) > K){
# # # 					my_fric_rand <- tryCatch(convhulln(my_tr_rand, "FA")$vol, error = function(e) return(NA))
# # # 				}
# # # 				
# # # 				outv <- rbind(outv, c("n" = n, "raoq" = my_raoq, "raoq_rand" = my_raoq_rand, "mpd" = my_mpd, "pd" = my_faiths, "fric" = my_fric, "fric_rad" = my_fric^(1/K), "fric_rand" = my_fric_rand, "fric_rand_rad" = my_fric_rand^(1/K)))
# # # 			}
# # # 			
# # # 			return(outv)
# # # 		}
# # # 		
# # # 		#
# # # 		outmat <- outmat %>% data.frame() %>% as_tibble() %>%
# # # 			mutate(fric_vol = fric^K, rep = kk, chunk = j)
# # # 		#
# # # 		write_csv(outmat, "results/null_results_each.csv", append= file.exists("results/null_results_each.csv"))
# # # 		
# # # 		rm(outmat)
# # # 		gc()
# # # 	}
# # # }
# # 
# # 
# # 
# # 
# # 
# # # # ggpairs(outmat, lower = list(continuous = "points"), upper = list(continuous = "points"))
# # # # 
# # # # # 
# # obs <- read_csv("results/obs_results.csv") %>% rename(fd = fric, fdr = fric_rad)
# # # # 
# # # # # 
# # # # obs %>% ggplot(aes(x = mpd, y = raoq))+geom_point()
# # # # obs %>% ggplot(aes(x = pd, y = fd))+geom_point()
# # # # obs %>% ggplot(aes(x = pd, y = fdr))+geom_point()
# # # # obs %>% ggplot(aes(x = n, y = fdr))+geom_point()
# # # # obs %>% ggplot(aes(x = n, y = fd))+geom_point()
# # # # obs %>% ggplot(aes(x = n, y = pd))+geom_point()
# # # # obs %>% ggplot(aes(x = n, y = raoq))+geom_point()
# # # # obs %>% ggplot(aes(x = n, y = mpd))+geom_point()
# # # 
# # sim <- read_csv("results/null_results_new.csv") %>% rename(fd = fric, fdr = fric_rad) %>% 
# # 	select(-raoq_rand)
# # # #
# # sim %>% ggplot(aes(x = mpd, y = raoq))+geom_point()
# # sim %>% ggplot(aes(x = pd, y = fd))+geom_point()
# # sim %>% ggplot(aes(x = pd, y = fdr))+geom_point()
# # sim %>% ggplot(aes(x = n, y = fd))+geom_point()
# # sim %>% ggplot(aes(x = n, y = fd))+geom_point()
# # sim %>% ggplot(aes(x = n, y = pd))+geom_point()
# # sim %>% ggplot(aes(x = n, y = raoq))+geom_point()
# # sim %>% ggplot(aes(x = n, y = mpd))+geom_point()
# # # 
# # # 
# # sim_means <- sim %>% gather(metric, value, -n) %>% group_by(n, metric) %>% summarize(null_mean = mean(value, na.rm = T), null_sd = sd(value, na.rm=T)) %>% ungroup
# # # 
# # obsg <- obs %>% gather(metric, value, -n, -grid_id) %>% left_join(sim_means) %>%
# # 	mutate(diff = value - null_mean, z = (value - null_mean)/null_sd, prop_diff = (value - null_mean)/null_mean) %>%  
# # 	arrange(grid_id, n) %>% filter(!is.na(value))
# # 
# # write_csv(obsg, "results/obs_null_z_scores.csv")
# # # 
# # # 
# # obs_z <- obsg %>% select(grid_id, n, metric, prop_diff) %>% spread(metric, prop_diff)
# # # 
# # # 
# # 
# # obs_z %>% ggplot(aes(x = n, y = raoq))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# # obs_z %>% ggplot(aes(x = n, y = mpd))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# # obs_z %>% filter(n>500) %>% ggplot(aes(x = mpd, y = raoq))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 1)
# # # 
# # # 
# # obs_z %>% ggplot(aes(x = n, y = fd))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# # obs_z %>% ggplot(aes(x = n, y = pd))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# # obs_z %>% filter(n>500) %>% ggplot(aes(x = pd, y = fd))+geom_point()+geom_smooth()+geom_abline(intercept = 0, slope = 0)
# # 
# # # 
# # # 
# # # all <- obs %>% left_join(sim %>% group_by(n) %>% summarize_all(.funs = c(sd, mean), na.rm=T) %>% mutate(fdr2 = fd^(1/8), fdr_rand = fd_rand^(1/8)) %>% 
# # # 						 	setNames(paste0("null_", names(.))) %>% rename(n = null_n))
# # # 
# # # all %>% ggplot(aes(x = null_pd, y = pd))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth(method = "lm", formula = 'y ~ -1 + x')
# # # all %>% ggplot(aes(x = null_fd, y = fd))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth()
# # # all %>% ggplot(aes(x = null_fdr, y = fdr))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth()
# # # all %>% ggplot(aes(x = null_raoq, y = raoq))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth()
# # # all %>% ggplot(aes(x = null_mpd, y = mpd))+geom_point()+geom_abline(intercept = 0, slope = 1) + geom_smooth()
# # 
# # # res <- read_csv("For_null_models/null_results5k_new.csv") %>% 
# # # 	gather(trait, value, -n) %>% 
# # # 	group_by(n, trait) %>% 
# # # 	summarize(nsim = length(value[!is.na(value)]), sd = sd(value, na.rm=T), value = mean(value, na.rm = T)) %>% 
# # # 	ungroup %>% 
# # # 	mutate(se = sd / sqrt(nsim)) %>% 
# # # 	filter(!is.na(value)) %>% 
# # # 	rename(null_value = value, null_sd = sd, null_se  = se)
# # # 
# # # 
# # # 
# # # all <- obs %>% left_join(res) %>% 
# # # 	mutate(z = (value - null_value)/null_sd) 
# # # 	
# # # 
# # # ggplot(all %>% filter(trait == "pd"), aes(x = n, y = z))+geom_point()+geom_smooth()
# # # ggplot(all %>% filter(trait == "fric"), aes(x = n, y = z))+geom_point()+geom_smooth()
# # # ggplot(all %>% filter(trait == "raoq"), aes(x = n, y = z))+geom_point()+geom_smooth()
# # # ggplot(all %>% filter(trait == "mpd"), aes(x = n, y = z))+geom_point()+geom_smooth()
# # # 
# # # 
# # # 
# # # all_gath <- all %>% select(grid_id, n, trait, value, null_value) %>% gather(fit, val, -grid_id, -n, -trait) 
# # # 
# # # all_gath %>% filter(trait == "pd") %>% group_by(n, fit) %>% summarize(val = mean(val)) %>% ungroup %>% ggplot(aes(x = n, y = val, color = fit))+geom_point()
# # # all_gath %>% filter(trait == "fric") %>% group_by(n, fit) %>% summarize(val = mean(val)) %>% ungroup %>% ggplot(aes(x = n, y = val, color = fit))+geom_point()
# # # 
# # # all %>% filter(trait%in%c("raoq", "mpd")) %>% select(grid_id, n, trait, z) %>% spread(trait, z) %>% ggplot(aes(x = mpd, y = raoq))+geom_point()+geom_abline(b = 0, a = 1)
# # # all %>% filter(trait%in%c("raoq", "mpd")) %>% select(grid_id, n, trait, value) %>% spread(trait, value) %>% ggplot(aes(x = mpd, y = raoq))+geom_point()+geom_abline(b = 0, a = 1)
# # # all %>% filter(trait%in%c("raoq", "mpd")) %>% select(grid_id, n, trait, null_value) %>% spread(trait, null_value) %>% ggplot(aes(x = mpd, y = raoq))+geom_point()+geom_abline(b = 0, a = 1)
# # # 
# # # 
# # # 
# # # res %>% filter(trait%in%c("raoq", "mpd")) %>% select(n, trait, null_value) %>% spread(trait, null_value) %>% 
# # # 	ggplot(aes(x = mpd, y = raoq))+geom_point()
# # # 
# # # 
# # # ggplot(all %>% filter(n_sub > 7) %>% rowwise() %>% mutate(pd_sub = plyr::round_any(pd_sub, 100)) %>% group_by(pd_sub, type) %>% summarize(fric_vol = mean(fric_vol)), aes(x = pd_sub, y = fric_vol, color = type))+geom_point()
# # # ggplot(all %>% filter(n_sub > 7) %>% rowwise() %>% mutate(pd = plyr::round_any(pd, 100)) %>% group_by(pd, type) %>% summarize(fric = mean(fric)), aes(x = pd, y = fric, color = type))+geom_point()
# # # ggplot(all %>% rowwise() %>% mutate(mpd = plyr::round_any(mpd, 1)) %>% group_by(mpd, type) %>% summarize(raoq = mean(raoq)), aes(x = mpd, y = raoq, color = type))+geom_point()
# # # 
# # # ggplot(all %>% slice(sample(1:nrow(.), 20000)), aes(x = mpd, y = raoq, color = type))+
# # # 	geom_point()
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000)), aes(x = n, y = raoq, color = type))+
# # # 	geom_point()
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000)) %>% filter(type == "obs"), aes(x = n, y = mpd, color = type))+
# # # 	geom_point()
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000)) , aes(x = n, y = pd, color = type))+
# # # 	geom_point()
# # # 
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000) )%>% filter(n_sub > 7) , aes(x = n_sub, y = fric, color = type))+
# # # 	geom_point()
# # # 
# # # ggplot(all %>% slice(sample(1:nrow(.), 10000) )%>% filter(n_sub > 7) %>% mutate(fric_rand = ifelse(type == "obs", fric_vol, fric_rand^7)), aes(x = n_sub, y = fric_rand, color = type))+
# # # 	geom_point()
# # # 
# # # 
# # # all %>%
# # # 
# # # ggplot(all %>% group_by(n, type) %>% summarize(sd = sd(raoq), raoq = mean(raoq)) %>% ungroup, aes(x = n, y = raoq, color = type))+
# # # 	geom_point()
# # # 
# # # ggplot(all %>% group_by(n, type) %>% summarize(sd = sd(mpd), mpd = mean(mpd)) %>% ungroup, aes(x = n, y = mpd, color = type))+
# # # 	geom_point()
# # # 
# # # 
# # # # lower phylogenetic diversity than expected by change
# # # ggplot(all %>% group_by(n, type) %>% summarize(sd = sd(pd), pd = mean(pd)) %>% ungroup, aes(x = n, y = pd, color = type))+
# # # 	geom_point()
# # # 
# # # ggplot(all %>% filter(n_sub > 7) %>% group_by(n_sub, type) %>% summarize(sd = sd(pd), pd = mean(pd)) %>% ungroup, aes(x = n_sub, y = pd, color = type))+
# # # 	geom_point()
# # # 
# # # # maybe comparable equivalent values
# # # ggplot(all %>% filter(n_sub > 7) %>% group_by(n_sub) %>% mutate(n_sub = plyr::round_any(n_sub, 50)) %>% ungroup %>% group_by(n_sub, type) %>% summarize(sd = sd(fric), fric = mean(fric)) %>% ungroup, aes(x = n_sub, y = fric, color = type))+
# # # 	geom_point()
# # # 
# # # 
# # # ggplot(all %>% filter(n_sub > 7) %>% group_by(n_sub) %>% mutate(n_sub = plyr::round_any(n_sub, 50)) %>% ungroup %>% group_by(n_sub, type) %>% summarize(sd = sd(fric_vol), fric_vol = mean(fric_vol)) %>% ungroup, aes(x = n_sub, y = fric_vol, color = type))+
# # # 	geom_point()
# # # 

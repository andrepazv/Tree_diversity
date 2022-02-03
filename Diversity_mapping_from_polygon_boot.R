rm(list = ls())

###October 2021###
library(rgdal)
library(raster)
library(V.PhyloMaker)
library(ape)
library(feather)
library(doParallel)
library(FD)
library(geosphere)
library(tidyverse)

ncores <- 14

###Setwd to folder containing the polygons###
setwd("~/Git/Tree_diversity/data/")


# replace NA biomes with closest non NA value
fix_na_biomes <- function(cd, nthreads = 14){
	
	# get the poitns with and without na's
	pts <- cd %>% select(x, y, grid, biome) %>% distinct()
	miss_pts <- pts %>% filter(is.na(biome))
	have_pts <- pts %>% filter(!is.na(biome))
	
	new_pts <- tibble()
	registerDoParallel(nthreads)
	# cycle through, in parallel, and get the closest point
	closest_pts <- foreach(i = 1:nrow(miss_pts), .inorder = FALSE, .combine = bind_rows)%dopar%{
		# extract the focal row
		my_pt <- miss_pts %>% slice(i)
		# vector of distances
		val <- rep(0, nrow(have_pts))
		#cycle through non na points and det distance
		for(j in 1:nrow(have_pts)){
			# get the distance between the focal point and new point
			val[j] <- distm(bind_rows(my_pt, have_pts %>% slice(j)) %>% select(x,y))[1,2]
		}
		# return the point with smallest distance
		return(my_pt %>% mutate(biome = have_pts %>% slice(which(val == min(val))[1]) %>% select(biome) %>% unlist))
	}
	
	# replace NAs in the full dataset with biomes, and return
	cd <- cd %>% left_join(closest_pts %>% rename(new_biome = biome), by = c("x", "y", "grid")) %>% mutate(biome = ifelse(is.na(biome), new_biome, biome)) %>% select(-new_biome)
	
	return(cd)
}
# get various pairwise distances
pairdist <- function(d, p){
	# convert to relative abundance
	p <- p/sum(p)
	# get the product matrix of pairwise relative abundances
	pp <- p%*%t(p)
	# weight the distance by probability of picking two species
	wpd <- pp*d
	# for rao's q, just sum up
	rq <- sum(wpd)
	# calculate abundance weighted pairwise distance 
	wmpd <- sum(wpd[upper.tri(wpd)]) / sum(pp[upper.tri(pp)])
	# calculate standard unweighted species pairwise distance
	mpd <- mean(d[upper.tri(d)])
	return(c("raoq" = rq, "wmpd" = wmpd, "mpd" = mpd))
}

# read in the community data
community_data <- read_csv("Community_matrix_res_10.csv")


########MAP DIVERSITY#########
##Taxonomic diversity##
# TD <- calc(Stack_maps, sum, na.rm = T) ##taxonomic diversity map
# plot(TD)
# writeRaster(TD, "richness_trees_10arcmin_agg10.asc")

##Phylogenetic diversity##
# species_names <- gsub("_", " ", species_names)
# 
# ag <- read_csv("~/Git/functional_disp/data/cleaned_data/ANGIO_GYMNO_lookup.csv")
# 
# sp_dat <- tibble(species = species_names) %>% mutate(genus = word(species, 1)) %>% left_join(ag %>% rename(species = accepted_bin)) %>% 
# 	select(species, genus, family) %>% data.frame()
# # sample call, assuming unique_spp_names is in Latin binomial format
# phylo_tree <- V.PhyloMaker::phylo.maker(sp_dat)
# ##note not all trees could be place (~3K not in the tree) 
# # write.tree(phylo_tree$scenario.3, file = "TEST_all_trees_tree.nwk")
# ###create list of missing trees
# species_not_tree <- subset(phylo_tree$species.list, status == "fail to bind")
# # write.table(species_not_tree, "Species_not_in_tree.txt")
# species_names <- gsub(" ", "_", species_names)
# species_names <- stringr::str_to_title(species_names)0.75, 1.25
# 
# ####MMake sure the names match between maps and phylogeny
# setdiff(phylo_tree$scenario.3$tip.label, species_names) ##chack var names here
# ##if it doesnt use this code
# ###Trim phylogeny to match distribution data (remove non hylids)
# pruned.tree <- ape::drop.tip(phylo_tree$scenario.3, setdiff(phylo_tree$scenario.3$tip.label, species_names))
# ##check that all disrtribution data is in tree
# test <- as.data.frame(species_names)
# rownames(test) <- species_names
# check_names <- geiger::name.check(pruned.tree, test, data.names = NULL) #this must be OK
# 
# 
# 
# #####Compute and map PD
# names_table <- as.data.frame(species_names)
# colnames(names_table) <- "Grid"
# 
# #Create empty raster of the desired area and resolution to assign pixel numbers
# r <- Stack_maps[[1]]
# res(r) <- res(Stack_maps) #resolution
# #r[is.na(r[])] <- 0
# #r <- crop(r, selected_mask)
# #r <- mask(r, selected_mask)
# #r <- aggregate(r, fact = 10, fun = max)
# grid = r
# names(grid) = "grid"
# grid[1:ncell(grid)] <- 1:ncell(grid)
# list_grid <- list(grid)
# names(list_grid) <- "Grid"
# #Change names of models tu species names and add the empty raster in the beggining (pixel number)
# names(Stack_maps) <- as.vector(names_table$Grid)
# Stack <- stack(List_grid$Grid, Stack_maps) ###Full stack of all maps
# #
# #Turn maps into dataframe for computation of PD
# community_data <- as.data.frame(Stack)
# ##remove rows that are NA 
# community_data1 <- community_data[-which(rowSums(!is.na(community_data[, 2:ncol(community_data1)])) == 0), ] ##change number of columns
# species_names <- colnames(community_data1)[2:length(community_data1)] #Store species names
# 
# 
# #Store 
# setwd("~/Dropbox/**PostDoc_ETH/Trees_PD_FD/Diversity_computation/")
# #write.table(community_data1, file = "communities.txt", append = FALSE, row.names = F, quote = F, sep = "\t") #tratar de agregar nombre de mascara
# arrow::write_feather(community_data1, "communities.feather")
# community_data1[is.na(community_data1)] = 0
# #In the community data frame NA must be eliminated done before is this if you load community data and not maps?
# #community_data = na.omit(community_data)
# #head(community_data)
# 
# #Phylogenetic diversity computation 
# #computes only Faith�s PD others may be added
# 
# pd.result <- picante::pd(community_data1[, 2:ncol(community_data1)], pruned.tree, include.root = F) 
# 
# #Add the pixel PD value to data frame
# community_data1$pd <- pd.result[[1]]
# 
# #Write the new matrix to a file to avoid rerunning the script for potential further analyses
# #write.table(community_data1, file = "communities_and_pd.txt", append = FALSE, row.names = F, quote = F, sep = "\t")
# arrow::write_feather(community_data1, "communities_withPD.feather")
# #Generate a raster containing PD information per pixel
# 
# #1-First generate an empty raster using a base model for resolution and area
# 
# values(r) <- 0
# pd_ras <- r
# values(pd_ras) <- NA #Eliminate every value they will be replaced by PD values further down
# 
# 
# #2- Assign PD values to raster
# pd_ras[community_data1$grid] <- community_data1$pd
# 
# #3- Save raster to file 
# 
# writeRaster(pd_ras, "Trees_PD.tif", format = "GTiff")
# 
# #4-Optional plotting map in R 
# plot(pd_ras)


#################################################################
#######Functional diversity computation here####################
################################################################
##Read Traits
traits <- read_feather("../ALL_BGCI_species_and_genera_imputed_traits_8_21.feather")
traits1 <- matrix(nrow = length(unique(traits$accepted_bin)), ncol = length(unique(traits$trait_name)), dimnames = list(unique(traits$accepted_bin), unique(traits$trait_name)))	 
##re-check this part why ?
for (i in unique(traits$trait_name)){
	test <- subset(traits, trait_name == i)
	traits1[, i] <- test$value
}
rownames(traits1) <- gsub(" ", "_", rownames(traits1))
colnames(traits1) <- gsub(" ", "_", colnames(traits1))
traits1 <- as.data.frame(traits1)
traits1 <- traits1[, c(1, 2, 7, 8, 28, 20, 26)] ##seelction from Dan's paper
#remove species not in traits from stack
# species_names <- names(Stack_maps)
# Stack_maps1 <- subset(Stack_maps, intersect(species_names, rownames(traits1)))
# species_names <- names(Stack_maps1)
# species_names <- species_names[order(species_names)]
##create df of ordered species names 
names_table <- as.data.frame(species_names)
colnames(names_table) <- "Grid"
#Create empty raster of the desired area and resolution to assign pixel numbers
r <- Stack_maps[[1]]
res(r) <- res(Stack_maps)
grid = r
names(grid) = "grid"
grid[1:ncell(grid)] <- 1:ncell(grid)
list_grid <- list(grid)
names(list_grid) <- "Grid"
#Change names of models tu species names and add the empty raster in the beggining (pixel number)
names(Stack_maps) <- as.vector(names_table$Grid)
Stack <- stack(list_grid$Grid, Stack_maps) ###Full stack of all maps
#
# Turn maps into dataframe for computation of PD
community_data <- as.data.frame(Stack, xy = TRUE)

# remove grids with no species and with NA for biomes
community_data <- community_data %>% as_tibble() %>% gather(metric, value, -grid, -biome, -x, -y) %>% filter(!is.na(value)) %>% filter(value>0) %>% spread(metric, value, fill = 0)

# # visualize NAs
# pd <- community_data %>% select(x, y, biome) %>% distinct() %>% mutate(biome = ifelse(is.na(biome), "red", "black"))
# community_data %>% filter(is.na(biome))
# plot(pd$y~pd$x, col = pd$biome)

# replace NA biomes with closest match
community_data <- fix_na_biomes(community_data)

# visualize biomes
pd <- community_data %>% select(x, y, biome) %>% distinct() 
community_data %>% filter(is.na(biome))
plot(pd$y~pd$x, col = pd$biome)


# gather the data
comm_stack <- community_data %>% gather(species, present, -grid, -biome, -x, -y) %>% mutate(present = ifelse(is.na(present), 0, present))

m <- nrow(comm_stack)
samps12 <- comm_stack %>% select(grid, species, present, biome) %>% filter(present==1) %>% distinct() %>% 
	group_by(species, biome) %>% tally() %>% ungroup %>% filter(n<3) %>% 
	group_by(n, biome) %>% tally() %>% ungroup %>% arrange(biome, n)

samps12 %>% filter(n == 1) %>% rename(q1 = nn) %>% select(-n) %>% left_join(samps12 %>% filter(n == 2) %>% rename(q2 = nn) %>% select(-n)) %>% 
	left_join(comm_stack %>% select(grid, biome) %>% distinct() %>% group_by(biome) %>% tally() %>% rename(m = n)) %>% 
	mutate(j1 = q1*(m - 1)/m, j2 = q1*(2*m-3)/m - (q2*(m-2)^2)/(m*(m-1)))

comm_stack %>% filter(present == 1) %>% select(grid, species, biome) %>% distinct() %>% group_by(grid, biome) %>% tally()

# remove communities with fewer species than traits
comm_stack <- comm_stack %>% filter(grid%in%(comm_stack %>% group_by(grid) %>% summarize(n = sum(present)) %>% ungroup %>% filter(n>(ncol(traits1)+1)) %>% select(grid) %>% unlist))

# remove species with no occurrences
comm_stack <- comm_stack %>% filter(species%in%(comm_stack %>% group_by(species) %>% summarize(n = sum(present)) %>% filter(n>0) %>% select(species) %>% unlist))

# spread the data 
cm <- comm_stack %>% spread(species, present)

# get the species names
species_names <- names(cm %>% select(-grid, -x, -y, -biome))

# log and normalize the traits
traits_norm <- traits1 %>% mutate_all(.funs = log) %>% mutate_all(.funs = function(x) (x - mean(x)) / sd(x))


rm(community_data1, traits1, traits, community_data)

# get the unique species within each biome
spp_biome <- cm %>% select(-x, -y, -grid) %>% gather(species, present, -biome) %>% filter(present == 1) %>% distinct() %>% select(-present)%>% arrange(species, biome) %>% 
	mutate(genus = word(species, 1, sep = "_"))

######################################
### Bootstrapping functional diversity
registerDoParallel(14)

chunk_seq <- sort(unique(seq(1, nrow(cm), by = 14), nrow(cm)))

nboot <- 500
new_spp_perc <- 0
trait_noise_perc <- 0.25
rand_traits <- FALSE

res <- foreach(s = chunk_seq, .inorder = FALSE, .combine = bind_rows, .packages = "FD")%dopar%{
	# get the endpoint
	end <- min(s+25-1, nrow(cm))
	#create output matrix
	all_grids <- cm$grid
	# strip off the grids
	my_commdat <- cm %>% select(-grid)
	# cycle through communities and calculate
	resvec <- NULL
	for(i in s:end){ 
		# reset
		my_fric <- my_fdisp <- base_spp <- no_tr_spp <- base_spp <- no_tr_gen <- miss_n <- have_n <- my_n <- NA
		# get species names
		base_spp <- colnames(my_commdat)[my_commdat[i, ] == 1]
		no_tr_spp <- base_spp[!base_spp%in%rownames(traits_norm)]
		base_spp <- base_spp[!base_spp%in%no_tr_spp]
		no_tr_gen <- table(word(no_tr_spp, 1, sep = "_"))
		# biome
		my_biome <- my_commdat %>% slice(i) %>% select(biome) %>% unlist
		# number of species
		miss_n <- length(no_tr_spp)
		have_n <- length(base_spp)
		my_n <- have_n + miss_n
		# create mock community matrix of 1s
		tmp_comm <- t(rep(1, have_n))
		# set names to current community
		colnames(tmp_comm) <- base_spp
		# subset the traits
		tmp_traits <- traits_norm[base_spp, ]
		# get pairwise distance matrix
		tmp_dist <- dist(tmp_traits)
		# calculate convex hull on unique points, only if enough obs
		if(nrow(tmp_traits) >= (ncol(tmp_traits)+1)){
			tmp_traits <- tmp_traits %>% distinct()
			if(nrow(tmp_traits)>= (ncol(tmp_traits)+1)){		
				my_fric <- geometry::convhulln(tmp_traits, "FA")$vol
			}
		}
		# how many species to add
		new_add <- round(my_n*new_spp_perc)
		# calculate fdisp and convex hull
		resvec <- resvec %>% rbind(c(i, all_grids[i], 0, my_n, have_n, my_fdisp, pairdist(as.matrix(tmp_dist), rep(1/have_n, have_n)), my_fric))
		#get species which are in this biome, in the trait database, and not in the base species
		my_spp_list <- spp_biome %>% filter(biome == my_biome & species%in%rownames(traits_norm) & !species%in%base_spp) %>% select(species) %>% unlist()
		# my_gen_list <- gen_biome %>% filter(biome == my_biome & species%in%rownames(traits_norm)) %>% select(species) %>% unlist()
		
		my_tr <- traits_norm[spp_biome %>% filter(biome == my_biome) %>% filter(species%in%rownames(traits_norm) & !species%in%base_spp) %>% select(species) %>% unlist(),]
		my_tr_gen <- my_tr %>% mutate(genus = word(rownames(.), 1, sep = "_"))
		for(b in 1:nboot){
			# get the traits where we have species 
			tmp_traits <- traits_norm[base_spp, ]
			
			# cycle through all the species with missing trait data and sample at random within genera and biome
			if(any(!is.na(no_tr_spp))){
				tmp1 <- matrix(ncol = ncol(my_tr), nrow = 0)
				for(k in 1:length(no_tr_gen)){
					if(names(no_tr_gen)[k]%in%my_tr_gen$genus){
						# genus is in the biome, sample within genus
						if(rand_traits){
							tmp1 <- tmp1 %>% rbind(matrix(apply(my_tr_gen %>% filter(genus == names(no_tr_gen)[k]) %>% select(-genus), 2, function(x) sample(x, no_tr_gen[k], replace = TRUE)), ncol = ncol(traits_norm)))
							rownames(tmp1) <- paste0("missing", 1:nrow(tmp1))
						}else{
							# sample a species rather than a trait
							tmp1 <- tmp1 %>% rbind(my_tr_gen %>% filter(genus == names(no_tr_gen)[k]) %>% select(-genus) %>% slice(sample(1:nrow(.), min(nrow(.), no_tr_gen[k]), replace = FALSE))%>% as.matrix())
						}
					}else{
						if(rand_traits){
							# genus isn't present in the biome, sample at random across all species in the biome
							tmp1 <- tmp1 %>% rbind(matrix(apply(my_tr, 2, function(x) sample(x, no_tr_gen[k], replace = TRUE)), ncol = ncol(traits_norm)))
							rownames(tmp1) <- paste0("missing", 1:nrow(tmp1))
						}else{
							tmp1 <- tmp1 %>% rbind(my_tr_gen %>% select(-genus) %>% slice(sample(1:nrow(.), min(nrow(.), no_tr_gen[k]), replace = FALSE)) %>% as.matrix())
						}
					}
				}
				colnames(tmp1) <- colnames(traits_norm)
				tmp_traits <- rbind(tmp_traits, tmp1)
			}
			
			# how many species are we adding/subtracting
			adjust_n <- sample(seq(-new_add, new_add, by = 1), 1)
			
			if(adjust_n>0){
				if(rand_traits){
					new_spp_tr <- matrix(apply(my_tr, 2, function(x) sample(x, adjust_n, replace = TRUE)), ncol = ncol(my_tr))
					rownames(new_spp_tr) <- paste0("mock", 1:nrow(new_spp_tr))
				}else{
					new_spp_tr <- my_tr %>% filter(!rownames(my_tr)%in%rownames(tmp_traits)) %>% slice(sample(1:nrow(.), min(nrow(.), adjust_n), replace = FALSE)) %>% as.matrix()
				}
				colnames(new_spp_tr) <- colnames(traits_norm)
				tmp_traits <- rbind(tmp_traits, new_spp_tr)
			}
			if(adjust_n<0){
				tmp_traits <- tmp_traits %>% slice(sample(1:nrow(.), my_n + adjust_n, replace = FALSE))
			}
			
			tmp_comm <- matrix(1, nrow = 1, ncol = nrow(tmp_traits)) %>% data.frame()
			colnames(tmp_comm) <- rownames(tmp_traits)
			# subset the traits
			
			# add a tiny bit of noise in case species were assigned same values
			tmp_traits <- tmp_traits*matrix(runif(nrow(tmp_traits)*ncol(tmp_traits), 1-trait_noise_perc-0.0001, 1+trait_noise_perc+0.0001), nrow = nrow(tmp_traits))
			# get pairwise distance matrix
			tmp_dist <- dist(tmp_traits)
			# or, subset the trait distance matrix and convert to dist object. not possible with too many species
			# tmp_dist <- as.dist(tr_dist[tmp_spp, tmp_spp])
			# my_fdisp <- FD::fdisp(d = tmp_dist, a = tmp_comm)$FDis
			# calculate convex hull on unique points, only if enough obs
			if(nrow(tmp_traits)>= (ncol(tmp_traits)+1)){
				# tmp_traits <- tmp_traits*matrix(runif(nrow(tmp_traits)*ncol(tmp_traits), 0.9999, 1.0001), nrow = nrow(tmp_traits))
				# if(nrow(tmp_traits)>= (ncol(tmp_traits)+1)){		
				my_fric <- geometry::convhulln(tmp_traits, "FA")$vol
				# }
			}
			# calculate fdisp and convex hull
			resvec <- resvec %>% rbind(c(i, all_grids[i], b, my_n, nrow(tmp_traits), my_fdisp, pairdist(as.matrix(tmp_dist), rep(1/nrow(tmp_traits), nrow(tmp_traits))), my_fric))
		}
	}
	# return as a tibble
	return(resvec %>% data.frame() %>% setNames(c("index", "grid", "boot", "nfull", "nspp", "fdisp", "raoq", "wmpd", "mpd", "fric")) %>% as_tibble())
}

# pairs(res %>% select(-grid))

res_comb <- res %>% filter(boot == 0) %>% select(-boot, -fdisp, -wmpd, -raoq, -nspp) %>% 
	gather(metric, value, -index, -grid, -nfull) %>% 
	left_join(res %>% select(-fdisp, -wmpd, -raoq, -nspp, -nfull) %>% gather(metric, value, -index, -grid, -boot) %>% filter(!is.na(value)) %>% 
			  	group_by(grid, metric) %>% 
			  	summarize(mean = mean(value), median = median(value), sd = sd(value), low = quantile(value, 0.025), high = quantile(value, 0.975)) %>% 
			  	ungroup) %>% 
	left_join(res %>% filter(boot > 0) %>% select(index, grid, nspp) %>% group_by(grid, index) %>% 
			  	summarize(n_mean = mean(nspp), n_med = median(nspp), n_low = min(nspp), n_max = max(nspp)) %>% ungroup)


ggplot(res_comb %>% mutate(low1sd = mean - sd, high1sd = mean + sd), aes(x = value, y = median))+geom_point()+facet_wrap(~metric, scale = "free")+geom_abline(intercept = 0, slope = 1)+
	geom_errorbar(aes(ymin = low, ymax = high))


res_comb <- res %>% filter(boot == 0) %>% select(-boot, -fdisp, -wmpd, -raoq, -nspp) %>% 
	gather(metric, value, -index, -grid, -nfull) %>% 
	left_join(res %>% filter(nspp == nfull) %>% select(-fdisp, -wmpd, -raoq, -nspp, -nfull) %>% gather(metric, value, -index, -grid, -boot) %>% filter(!is.na(value)) %>% 
			  	group_by(grid, metric) %>% 
			  	summarize(mean = mean(value), median = median(value), sd = sd(value), low = quantile(value, 0.025), high = quantile(value, 0.975)) %>% 
			  	ungroup) %>% 
	left_join(res %>% filter(nspp == nfull) %>% filter(boot > 0) %>% select(index, grid, nspp) %>% group_by(grid, index) %>% 
			  	summarize(n_mean = mean(nspp), n_med = median(nspp), n_low = min(nspp), n_max = max(nspp)) %>% ungroup)


ggplot(res_comb %>% mutate(low1sd = mean - sd, high1sd = mean + sd), aes(x = value, y = median))+geom_point()+facet_wrap(~metric, scale = "free")+geom_abline(intercept = 0, slope = 1)+
	geom_errorbar(aes(ymin = low, ymax = high))


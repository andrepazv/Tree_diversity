rm(list = ls())
gc()

library(geometry)
library(doParallel)
library(vegan)
library(picante)
library(feather)
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

# number of traits
K <- k <- length(trait_names)

# # specify monocot orders to make sure none are present
monocots <- c("Acorales", "Alismatales", "Arecales", "Asparagales", "Commelinales", "Dioscoreales", "Liliales", "Pandanales", "Petrosaviales", "Poales", "Zingiberales")
ferns <- c("Lophosoria", "Metaxya", "Sphaeropteris", "Alsophila", "Nephelea", "Trichipteris", "Cyathea" ,"Cnemidaria", "Dicksonia", "Cystodium", "Thyrsopteris", "Culcita", "Cibotium")

# read in the list of angio vs gymno species
ag0 <- read_csv("raw_data/ANGIO_GYMNO_lookup.csv")

# read in the community data
commo <- feather::read_feather("raw_data/NoMono_matchPD_FD_gather.feather") %>% 
	mutate(accepted_bin = trimws(gsub("_", " ", accepted_bin), "both"))


# read in the synonyms
syns <- read_csv("raw_data/matched_trait_names.csv") %>% 
	filter(!duplicated(trait_name), !duplicated(raw_name))

tn <- read_csv("raw_data/trait_lookup.csv")

# # read in the trait table and subset
dto <-
	# feather::read_feather("raw_data/Estimated_trait_table.feather") %>%
	read_csv("~/Git/dispersion/data/combined_datasets/Species_trait_table_50_with_monos_AND_OBS_3_23.csv") %>% left_join(tn) %>% 
	filter(!is.na(trait)) %>% 
	mutate(trait = gsub(" ", "_", trait)) %>%
	# filter(trait%in%trait_names) %>%
	select(accepted_bin, trait, value) %>%
	spread(trait, value) 


trait_names <- dto %>% select(-accepted_bin) %>% names(.)

# read in the monocots to exclude
monocot_genera <- read_csv("raw_data/monocot_genera.csv")


# add in new names 
dt <- dto %>% 
	mutate(tax_genus = word(accepted_bin, 1)) %>% 
	# filter(!tax_genus%in%monocot_genera$tax_genus) %>%
	select(accepted_bin, tax_genus, all_of(trait_names)) %>% 
	left_join(ag0 %>% select(accepted_bin, order, group, genus, family)) %>%
	# filter(!order%in%monocots, !genus%in%ferns, !tax_genus%in%ferns, !family%in%c("Osmundaceae"), group%in%c("Angiosperms", "Gymnosperms")) %>%
	select(accepted_bin, all_of(trait_names)) 


# get the PCA of traits
dt_mat <- 
	# rda(dt %>% select(-accepted_bin), scale = TRUE)$CA$u %>% data.frame() %>% as_tibble() %>% 
	# mutate(accepted_bin = dt$accepted_bin) %>% 
	dt %>% 
	gather(PC, value, -accepted_bin) %>%
	# group_by(PC) %>%
	# mutate(value = (value - mean(value))/sd(value)) %>%
	# ungroup %>%
	spread(PC, value)

rm(dt)
rm(dto)

# replace the synonyms
for(i in 1:nrow(syns)){
	dt_mat$accepted_bin[which(dt_mat$accepted_bin == syns$trait_name[i])]<-syns$raw_name[i]
}


# remove ferns/palms
comm <- commo %>%
	select(grid_id, present, accepted_bin)


# see the mismatches
sum(!dt_mat$accepted_bin%in%comm$accepted_bin)
unique(comm$accepted_bin[!comm$accepted_bin%in%dt_mat$accepted_bin])

# subset the matched species and remove any plots with singletons
comm <- comm %>% filter(accepted_bin%in%dt_mat$accepted_bin) %>% distinct() %>% 
	left_join(comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))) %>% 
	select(-n)

# subset the matched species, and summarize any duplicates
dt_mat <- dt_mat %>% filter(accepted_bin%in%comm$accepted_bin) %>% distinct()

# double check
sum(!dt_mat$accepted_bin%in%comm$accepted_bin)
unique(comm$accepted_bin[!comm$accepted_bin%in%dt_mat$accepted_bin])

# get the unique species
spp <- comm %>% select(accepted_bin) %>% distinct() %>% unlist()

# get the number fo species per plot
nspp <- comm %>% select(grid_id, accepted_bin) %>% distinct() %>% group_by(grid_id) %>% tally() %>% ungroup %>% arrange(desc(n))

# get the plots
plts <- nspp %>% select(grid_id) %>% distinct() %>% unlist() %>% as.numeric()

# get the tree and clean the names
pruned_tree <- read.tree("raw_data/no_monocots_tree.nwk")
# pruned_tree <- read.tree("raw_data/all_trees_tree.nwk")

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

# shuffle the order
plts <- sample(plts, length(plts), replace = FALSE)


# # output the community matrix
cd <- comm %>% mutate(accepted_bin = gsub(" ", "_", accepted_bin)) %>% spread(accepted_bin, present, fill = 0) %>%
	rename(Grid = grid_id)

### output the data, replacing species names with underscores
write_csv(comm %>% mutate(accepted_bin = gsub(" ", "_", accepted_bin)), "results/Community_matrix_all_18_with_monos.csv")
write_csv(dt_mat %>% mutate(accepted_bin = gsub(" ", "_", accepted_bin)), "results/Trait_matrix_for_all_18_with_monos.csv")
write.tree(pruned_tree, "results/Pruned_tree_all_18_with_monos.tree")
# write_csv(comm %>% mutate(accepted_bin = gsub(" ", "_", accepted_bin)) %>% spread(accepted_bin, present, fill = 0) %>%
# 		  	rename(Grid = grid_id), "results/Community_matrix_for_FD_PD_38k_spread.csv")

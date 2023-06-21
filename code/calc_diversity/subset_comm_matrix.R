rm(list = ls())
gc()

library(geometry)
library(doParallel)
library(vegan)
library(picante)
library(feather)
library(tidyverse)

# raoq formula. assumes equal abundances
pairdist <- function(myd){
	return(sum((myd / nrow(myd))^2))
}


setwd("~/Git/Tree_diversity/data")

comm0 <- read_csv("cleaned_data/REV_Community_matrix.csv")
dt_mat <- read_csv("cleaned_data/REV_Trait_matrix.csv")
pruned_tree <- read.tree("cleaned_data/REV_Pruned_tree.csv")
dmat <- read_csv("raw_data/distances_plots_obs.csv") 
cont <- read_csv("raw_data/continents_grid.csv") %>% select(Grid, Continent) %>% rename(grid_id = Grid)
ll <- read_csv("raw_data/ll_with_species.csv") %>% left_join(cont)
dcut <- read_csv("cleaned_data/species_cut_ranges.csv") %>% select(accepted_bin, maxd)

# double check that the right files are loaded, all should be zero, otherwise clear the environment
if((sum(!dt_mat$accepted_bin%in%comm0$accepted_bin) + 
	sum(!comm0$accepted_bin%in%dt_mat$accepted_bin) +
	sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label) +
	sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label) +
	sum(!dt_mat$accepted_bin%in%comm0$accepted_bin) +
	sum(!comm0$accepted_bin%in%dt_mat$accepted_bin) +
	sum(!unique(comm0$accepted_bin)%in%pruned_tree$tip.label) +
	sum(!dt_mat$accepted_bin%in%pruned_tree$tip.label)) >0) {rm(list=ls()); stop("name mismatch")}

bad_spp <- c("Rhododendron_watsonii", "Syringa_pinetorum", "Acacia_farnesiana", "Photinia_davidiana", "Paulownia_kawakamii", "Eucalyptus_pulverulenta")

comm <- comm0 %>% left_join(dmat) %>% left_join(dcut) %>% filter(mindist < maxd) %>% 
	filter(!accepted_bin%in%bad_spp)


write_csv(comm, "cleaned_data/REV_Community_matrix_PCUT.csv")


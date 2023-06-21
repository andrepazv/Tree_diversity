library(geometry)
library(doParallel)
library(vegan)
library(picante)
library(feather)
library(tidyverse)

setwd("~/Git/Tree_diversity/data/")
# register the parallel cores
nthreads <- round(detectCores()/2)
registerDoParallel(nthreads)

cm <- read_csv("cleaned_data/REV_Community_matrix_CLEANED.csv")
cm <- read_csv("cleaned_data/REV_Community_matrix.csv")


dmat <- read_csv("raw_data/distances_plots_obs.csv")
dt <- read_csv("results/REV_obs_results.csv")

ll <- read_csv("raw_data/ll_with_species.csv")


dt %>% left_join(ll) %>% arrange(desc(raoq))


dt <- read_csv("/run/media/dsm/Elements/Alpha_hulls_April2022/all_trees_alldb.csv")

spt <- dt %>% mutate(Latitude = round(Latitude), Longitude = round(Longitude)) %>% distinct() %>% 
	group_by(Species) %>% tally() %>% ungroup %>% mutate(accepted_bin = gsub(" ", "_", Species)) %>% 
	rename(nobs = n)


spd <- dmat %>% group_by(accepted_bin) %>% summarize(n = length(mindist), 
													 n_pres = sum(mindist == 0),
													 # n_obs = mean(n_obs),
													 # n500 = sum(mindist > 500), 
													 # n1000 = sum(mindist > 1000), 
													 p500 = mean(mindist > 500), 
													 p1000 = mean(mindist > 1000)) %>% left_join(spt) %>% select(-Species) %>% 
	select(accepted_bin, n, n_pres, nobs, p500, p1000) %>% 
	mutate(maxd = ifelse(p1000 > 0.1 | p500 > 0.5, 100, ifelse(p1000 > 0.05 | p500 > 0.25, 250, ifelse(p1000 > 0.01 | p500 > 0.1, 500, 1000)))) 

table(spd$maxd)

spd %>% filter(accepted_bin == "Sequoia_sempervirens") 

# write_csv(spd, "cleaned_data/species_cut_ranges.csv")

bad <- read_csv("cleaned_data/grid31782_omitted_species.csv")


spd %>% filter(accepted_bin%in%bad$accepted_bin) %>% View()




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


dt <- read_csv("/run/media/dsm/Elements/Alpha_hulls_April2022/all_trees_alldb.csv") %>% 
	mutate(accepted_bin = gsub(" ", "_", Species)) %>% 
	mutate(Latitude = round(Latitude), Longitude = round(Longitude)) %>% 
	select(accepted_bin, Latitude, Longitude) %>% distinct()

spt <- dt %>% group_by(accepted_bin) %>% tally() %>% ungroup %>% rename(nobs = n)


spd <- dmat %>% group_by(accepted_bin) %>% summarize(n = length(mindist), 
													 n_pres = sum(mindist == 0),
													 # n_obs = mean(n_obs),
													 # n500 = sum(mindist > 500), 
													 # n1000 = sum(mindist > 1000), 
													 p500 = mean(mindist > 500), 
													 p1000 = mean(mindist > 1000)) %>% left_join(spt) %>% 
	select(accepted_bin, n, n_pres, nobs, p500, p1000) %>% 
	mutate(maxd = ifelse(p1000 > 0.1 | p500 > 0.5, 100, ifelse(p1000 > 0.05 | p500 > 0.25, 250, ifelse(p1000 > 0.01 | p500 > 0.1, 500, ifelse(p1000 > 0, 1000, 2500))))) 

table(spd$maxd)

spd %>% filter(accepted_bin == "Sequoia_sempervirens") 

spd %>% filter(accepted_bin%in%c("Rhododendron_watsonii", "Syringa_pinetorum", "Acacia_farnesiana", "Photinia_davidiana", "Paulownia_kawakamii", "Eucalyptus_pulverulenta"))

spd %>% filter(p1000 > 0.01) %>% View()

# write_csv(spd, "cleaned_data/species_cut_ranges.csv")

bad <- read_csv("cleaned_data/grid31782_omitted_species.csv")

checkspp <- spd %>% filter(p1000 > 0.01) %>% arrange(desc(nobs)) %>%
	slice(1) %>% select(accepted_bin) %>% unlist()

# checkspp <- "Sequoia_sempervirens"
checkspp <- "Acacia_dealbata"
dmat %>% filter(accepted_bin == checkspp) %>% left_join(ll) %>% 
	mutate(keep = ifelse(mindist <= 200, 0, 1)) %>% 
	# mutate(keep = ifelse(mindist < 100, "100", ifelse(mindist < 200, "200", ifelse(mindist < 500, "500", ifelse(mindist < 1000, "1000", ">1000"))))) %>%
	ggplot(aes(x = Longitude, y = Latitude, color = as.factor(keep))) + geom_point()+
	geom_point(data = dt %>% filter(accepted_bin == checkspp), aes(x = Longitude, y = Latitude), color = "black")




spd %>% filter(accepted_bin%in%bad$accepted_bin) %>% View()




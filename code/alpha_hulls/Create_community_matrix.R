rm(list = ls())

###October 2021###
library(rgdal)
library(raster)
library(doParallel)
library(geosphere)
library(tidyverse)

# read in arguments passed to R, for use on Euler. 
args <- commandArgs(TRUE)
if(length(args)>0){
	for(i in 1:length(args)){eval(parse(text=args[[i]]))}
	print(nthreads)
	print(resolution)
	print(overwrite)
}else{
	# if running locally, just specify values
	nthreads <- 10
	resolution <- 10
	overwrite <- TRUE
}


###Setwd to folder containing the polygons###
setwd("~/local_Git/Tree_diversity/data/")

#######################################################################################
# Functions
######################################################################################

registerDoParallel(nthreads)

####Some species have <3 points they are not polygons but points but this would work for both 
### function to rasterize points or polygons
rasterize_species <- function(x, shape_dir, out_dir, my_mask, my_base, my_res, overwrite) {
	
	file_name <- paste(out_dir,"/", x, "_", my_res, ".asc", sep = "")
	
	if(file.exists(file_name) & !overwrite){
		r <- raster(file_name)
		# return (r)
	}else{
		
		map <- readOGR(dsn = shape_dir, layer = x)
		r <- rasterize(x = map, y = my_base, field = 1, update = T, background = 0)
		r <- mask(r, my_mask)
		r <- crop(r, my_mask)
		names(r) <- gsub("map_", "", x)
		writeRaster(r, file_name, overwrite = TRUE)
		
		# return(r)
	}
	return()
}

# replace NA biomes with closest non NA value
fix_na_biomes <- function(cd, nthreads){
	
	# for debugging
	# cd <- cd %>% mutate(ecoregion = ifelse(grid%in%c(3269, 1945), NA, ecoregion))
	
	# get the points with and without na's
	miss_pts <- cd %>% filter(is.na(biome) | is.na(ecoregion))
	have_pts <- cd %>% filter(!is.na(biome) & !is.na(ecoregion))
	
	#register threads
	registerDoParallel(nthreads)
	
	# make sure we have some missing points
	if(nrow(miss_pts)>0){
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
			# get the closest point
			min_df <- have_pts %>% slice(which(val == min(val))[1])
			# return the point with smallest distance
			return(my_pt %>% mutate(biome = ifelse(is.na(biome), min_df$biome, biome), ecoregion = ifelse(is.na(ecoregion), min_df$ecoregion, ecoregion)))
		}
		
		# replace NAs in the full dataset with biomes, and return
		fixed_cd <- have_pts %>% bind_rows(closest_pts)
	}else{
		return(cd)
	}
	
	return(fixed_cd)
}

# get the mode, used for assigning biomes
Mode <- function(v) {
	uniqv <- unique(v)
	uniqv[which.max(tabulate(match(v, uniqv)))]
}

######################################################################################
######################################################################################


# read in the mask layer
selected_mask <- readOGR("base_rasters/world_adm0/", "world_adm0") 
selected_mask <- aggregate(selected_mask, dissolve = TRUE)

### Get a raster base, can be a worldclim variable for the whole world, it will set the resolution and extent
raster_base <- raster("WorldClim2/wc2.1_10m_bio_1.tif") ### use a 10 min WC layer or even coarser?

## make coarser this is toooo fine, and extend/crop/mask to the selected mask
raster_base <- raster::aggregate(raster_base, resolution)
values(raster_base) <- 0

## Get the names of species and distribution file names
distribution_files <- str_sort(list.files(path = "Alpha_clean/", pattern = "*.shp$"), numeric = TRUE)
distribution_maps <- sub(".shp", "", distribution_files) 

fia_names <- read_csv("~/Git/fia_data/data/metadata/2019 Master Species FGver90_10_01_2019_rev_2_10_2020.csv") %>% 
	mutate(name = str_to_sentence(trimws(paste(Genus, Species), which = "both"))) %>% select(name, Genus, Species) %>% distinct() 

in_fia <- sapply(distribution_maps, function(x) str_to_sentence(trimws(paste(str_split(x, pattern = "_")[[1]][2:3], collapse = " "), which = "both"))%in%fia_names$name)

distribution_maps <- distribution_maps[in_fia]


## rasterize them and write in parallel. This no longer loads the files
mclapply(distribution_maps, 
				   FUN = rasterize_species, 
				   my_mask = selected_mask,
				   my_base = raster_base, 
				   my_res = resolution, 
				   overwrite = overwrite, 
				   shape_dir = "Alpha_clean/", 
				   out_dir = "Alpha_rasters/",
				   mc.cores = nthreads)

# now load the files and stack them
rast_list <- str_sort(list.files(path = "Alpha_rasters/"), numeric = TRUE)
rast_list <- rast_list[grepl(paste0("_", resolution,".asc"), rast_list)]
layers <- stack(paste0("Alpha_rasters/",rast_list))

## read in biomes poly
biomes_poly <- readOGR(dsn = "~/local_Git/Tree_diversity/data/WWF_ecorregions/WWF_Biomes.shp")
biomes_rast <- rasterize(x = biomes_poly, y = raster_base, field = "BIOME", fun = function(x, na.rm) as.numeric(Mode(x[!is.na(x)])))
biomes_rast <- mask(biomes_rast, selected_mask)
biomes_rast <- crop(biomes_rast, selected_mask)
names(biomes_rast) <- "biome"

## read in ecoregion poly
ecoregion_poly <- readOGR(dsn = "~/local_Git/Tree_diversity/data/WWF_ecorregions/wwf_terr_ecos.shp")
eco_ids <- tibble(name = ecoregion_poly@data$ECO_NAME, id = ecoregion_poly@data$ECO_ID) %>% distinct() %>% arrange(id)
ecoregion_rast <- rasterize(x = ecoregion_poly, y = raster_base, field = "ECO_ID", fun = function(x, na.rm) as.numeric(Mode(x[!is.na(x)])))
ecoregion_rast <- mask(ecoregion_rast, selected_mask)
ecoregion_rast <- crop(ecoregion_rast, selected_mask)
names(ecoregion_rast) <- "ecoregion"

## Stack all maps 
Stack_maps <- stack(append(append(layers, biomes_rast), ecoregion_rast))

## get the grid ids and add to the stack as a new layer
mocklayer <- Stack_maps[[1]]
res(mocklayer) <- res(Stack_maps)
names(mocklayer) = "grid"
#mocklayer[1:ncell(mocklayer)] <- 1:ncell(mocklayer)
mocklayer<-init(mocklayer,"cell")
names(mocklayer) <- "Grid"
# list_grid <- list(mocklayer)
Stack <- stack(mocklayer, Stack_maps) 

# Turn maps into dataframe. remove grids with no species
community_data <- as.data.frame(Stack, xy = TRUE) %>% as_tibble() %>% 
	gather(species, present, -Grid, -biome, -ecoregion, -x, -y) %>%
	filter(!is.na(present)) %>%
	filter(present>0) %>%
	spread(species, present, fill = 0) %>% 
	mutate(ecoregion = ifelse(ecoregion<0, NA, ecoregion)) # rock/ice

## visualize NAs, for debugging
# pd <- community_data %>% select(x, y, biome, ecoregion) %>% distinct() %>% 
# 	mutate(biome_col = ifelse(is.na(biome), "red", "black")) %>% 
# 	mutate(eco_col = ifelse(is.na(ecoregion), "red", "black"))
# community_data %>% filter(is.na(ecoregion))
# plot(pd$y~pd$x, col = pd$biome_col)
# plot(pd$y~pd$x, col = pd$eco_col)

# replace NA biomes with closest match
community_data <- fix_na_biomes(cd = community_data, nthreads = nthreads)

### visualize biomes
pd <- community_data %>% select(x, y, biome, ecoregion) %>% distinct()
plot(pd$y~pd$x, col = pd$biome)
plot(pd$y~pd$x, col = pd$ecoregion)

write_csv(community_data, paste0("Community_matrix_res_", resolution, ".csv"))

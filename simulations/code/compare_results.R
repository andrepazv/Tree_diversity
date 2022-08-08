library(tidyverse)

dt1 <- read_csv("~/Git/Tree_diversity/simulations/results/obs_null_z_scores_fixed.csv")
dt2 <- read_csv("~/Git/Tree_diversity/simulations/results/obs_results_subsampled.csv") %>% filter(n>1) %>% 
	rename(fd = fric, fdr = fric_rad) %>% 
	gather(metric, value, -grid_id, -n) %>% select(-n) 


dt1 %>% left_join(dt2 %>% rename(value2 = value)) %>% 
	ggplot(aes(x = value, y = value2))+facet_wrap(~metric, scales = "free")+geom_point(alpha = 0.1)+
	geom_abline(slope = 1, intercept = 0, color = "red")+geom_smooth(method = "lm")

# Functional-fingerprints
#R Code for: Change in avian functional fingerprints of a Neotropical montane forest over 100 years as an indicator of ecosystem integrity. Gomez, C., Tenorio, E., Cadena, C.D. 2021. Conservation Biology.


#1. Functional space and hypervolume estimation. Code adapted from Cooke et al. 2019
####################################################################################

library(pacman)
pacman::p_load(dplyr, plyr, readr, tbible, FD, ade4, cowplot, mice, reshape2, tidyr, ks, hypervolume, alphallhu, purrr, TTR, plotrix, agricolae, psych)

library(factoextra)
library(ggrepel)
library(tibble)

#### Load data

# load: trait data
trait <- readr::read_csv("~/Documents/Funcional San Antonio/TABLA_RASGOS_FINAL.csv")

#Single PCA for all years

#z Transform traits
# function to z-transform data
scale_z <- function(x){
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}

### add scaled trait measures to traits table
trait_z <- trait %>% 
  mutate(log_weight = log(weight)) %>% 
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(TotalCulmen)) %>% 
  mutate(Bill_w_z = scale_z(BillWidth)) %>% 
  mutate(Wing_l_z = scale_z(WingLength)) %>% 
  mutate(Tail_l_z = scale_z(TailLength)) %>% 
  mutate(Tars_z = scale_z(Tarsus)) %>% 
  mutate(HWI_z = scale_z(HWI)) %>% 
  mutate(Hab_br_z = scale_z(Hab_breadth)) %>%
  mutate(Gen_z = scale_z(Gen_length)) %>%
  as.data.frame

#subset scaled measures
trait_z <- trait_z %>% 
  select(binomial,Y1911, Y1959, Y1994, Y2016,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z) 


#table just with columns of morphological traits and species names as rownames
tr_mi_z <- trait_z %>% 
  select(binomial, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z)%>% 
  column_to_rownames(var = "binomial")

#### --------------------------------------------------------------
## PCA ##
#### --------------------------------------------------------------

#Run PCA
pcaTotal <- princomp(tr_mi_z, cor = TRUE, scores = TRUE)
summary(pcaTotal)

#Scores
scoresPCATotal <- as.data.frame(pcaTotal$scores) %>% 
  tibble::rownames_to_column("binomial")

scoresPCATotal <- scoresPCATotal %>% 
  # convert long to wide
  tidyr::gather(key, value, -binomial) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal <- as.data.frame(unclass(pcaTotal$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal <- loadingsPCATotal %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc <- loadingsPCATotal %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Generation length","Habitat breadth","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################

#Subset for species in each period

Y1911 <- trait_z %>% 
  filter(Y1911 == 1) %>% 
  select(binomial, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z) %>% 
  left_join(scoresPCATotal, by = "binomial") 


Y1959 <- trait_z %>% 
  filter(Y1959 == 1) %>% 
  select(binomial, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z) %>%
  left_join(scoresPCATotal, by = "binomial") 


Y1994 <- trait_z %>% 
  filter(Y1994 == 1) %>% 
  select(binomial, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z) %>% 
  left_join(scoresPCATotal, by = "binomial") 

Y2016 <- trait_z %>% 
  filter(Y2016 == 1) %>% 
  select(binomial, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z) %>% 
  left_join(scoresPCATotal, by = "binomial") 


# kernel density estimation for each period
pc_raw_1911 <- Y1911 %>% 
  # extract first two principal components
  dplyr::select(., binomial, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "binomial")

pc_raw_1959 <- Y1959 %>% 
  # extract first two principal components
  dplyr::select(., binomial, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "binomial")

pc_raw_1994 <- Y1994 %>% 
  # extract first two principal components
  dplyr::select(., binomial, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "binomial")

pc_raw_2016 <- Y2016 %>% 
  # extract first two principal components
  dplyr::select(., binomial, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "binomial")

# optimal bandwidth estimation
hpi_1911 <- Hpi(x = pc_raw_1911)
hpi_1959 <- Hpi(x = pc_raw_1959)
hpi_1994 <- Hpi(x = pc_raw_1994)
hpi_2016 <- Hpi(x = pc_raw_2016)

# kernel density estimation 1911   
est_1911 <- kde(x = pc_raw_1911, H = hpi_1911, compute.cont = TRUE)  

den_1911 <- list(est_1911$eval.points[[1]], est_1911$eval.points[[2]], est_1911$estimate)
names(den_1911) <- c("x", "y", "z")
dimnames(den_1911$z) <- list(den_1911$x, den_1911$y)
dcc_1911 <- melt(den_1911$z)

# source: kernel function
## --------------------------------------------------------------
## Name: 2.1.3-kernel-function.R
## Description: Function to calculate kernel density probablities
## Date: October 2017
## Outputs: Function named 'cl'
## Args:
## df = dataframe of kernel density data
## prob = Probabilty level e.g. 0.95 (95% confidence level)
## --------------------------------------------------------------

cl <- function(df, prob) {
  dx <- diff(df$x[1:2])
  dy <- diff(df$y[1:2])
  sz <- sort(df$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

source("2.1.3-kernel-function.R")

# 0.5 probability kernel
cl_50_1911 <- cl(df = den_1911, prob = 0.50)
# 0.95 probability kernel
cl_95_1911 <- cl(df = den_1911, prob = 0.95)
# 0.99 probability kernel
cl_99_1911 <- cl(df = den_1911, prob = 0.99)

# kernel density estimation 1959   
est_1959 <- kde(x = pc_raw_1959, H = hpi_1959, compute.cont = TRUE)  

den_1959 <- list(est_1959$eval.points[[1]], est_1959$eval.points[[2]], est_1959$estimate)
names(den_1959) <- c("x", "y", "z")
dimnames(den_1959$z) <- list(den_1959$x, den_1959$y)
dcc_1959 <- melt(den_1959$z)


# 0.5 probability kernel
cl_50_1959 <- cl(df = den_1959, prob = 0.50)
# 0.95 probability kernel
cl_95_1959 <- cl(df = den_1959, prob = 0.95)
# 0.99 probability kernel
cl_99_1959 <- cl(df = den_1959, prob = 0.99)

# kernel density estimation 1994   
est_1994 <- kde(x = pc_raw_1994, H = hpi_1994, compute.cont = TRUE)  

den_1994 <- list(est_1994$eval.points[[1]], est_1994$eval.points[[2]], est_1994$estimate)
names(den_1994) <- c("x", "y", "z")
dimnames(den_1994$z) <- list(den_1994$x, den_1994$y)
dcc_1994 <- melt(den_1994$z)


# 0.5 probability kernel
cl_50_1994 <- cl(df = den_1994, prob = 0.50)
# 0.95 probability kernel
cl_95_1994 <- cl(df = den_1994, prob = 0.95)
# 0.99 probability kernel
cl_99_1994 <- cl(df = den_1994, prob = 0.99)

# kernel density estimation 2016  
est_2016 <- kde(x = pc_raw_2016, H = hpi_2016, compute.cont = TRUE)  

den_2016 <- list(est_2016$eval.points[[1]], est_2016$eval.points[[2]], est_2016$estimate)
names(den_2016) <- c("x", "y", "z")
dimnames(den_2016$z) <- list(den_2016$x, den_2016$y)
dcc_2016 <- melt(den_2016$z)


# 0.5 probability kernel
cl_50_2016 <- cl(df = den_2016, prob = 0.50)
# 0.95 probability kernel
cl_95_2016 <- cl(df = den_2016, prob = 0.95)
# 0.99 probability kernel
cl_99_2016 <- cl(df = den_2016, prob = 0.99)

# save principal component data
write.csv(Y1911, file = "PCA_1911.csv", row.names = FALSE)
write.csv(Y1959, file = "PCA_1959.csv", row.names = FALSE)
write.csv(Y1994, file = "PCA_1994.csv", row.names = FALSE)
write.csv(Y2016, file = "PCA_2016.csv", row.names = FALSE)

write.csv(scoresPCATotal, file = "PCA_Total.csv")
write.csv(loadingsPCATotal, file = "Loadings_Total.csv", row.names = FALSE)


#### PCA plots
library(ggplot2)

# colour palette
col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)


### Folowing lines indentify species that correspond to extinction,
#recolonization or new additions in different periods
ex <- read.csv("~/Documents/Funcional San Antonio/Extirpated.csv", h = T) #table of extirpated species
ext <- ex %>%
  mutate(binomial = Palacio.et.al.2019) %>% 
  left_join(Y1911, by = 'binomial') ##add scaled trait for extirpated species

extir <- subset(ext, X1911 == 1)
recol<- subset(ext, X2016 == 1)

ne<- read.csv("~/Documents/Funcional San Antonio/NewSpecies.csv", h = T) #table of new species in each period
new <- ne %>% 
  left_join(Y2016, by = "binomial") ##add scaled trait for new species

ex16<- ext <- ex %>%
  mutate(binomial = Palacio.et.al.2019) %>% 
  left_join(Y2016, by = 'binomial')
recol<- subset(ex16, X2016 == 1)

# plot 1911
pca_plot_1911 <- ggplot(dcc_1911, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = Y1911, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_1911, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_1911, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_1911, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Body size (55%)", y = "PC2 - Dispersal ability (13%)") +
  xlim(-5,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 25)) + 
  geom_text(x= 5, y=8, label="1910s", color="black", size = 5)

# display plot
windows()
pca_plot_1911


# plot 1959

ext <- ex %>%
  mutate(binomial = Palacio.et.al.2019) %>% 
  left_join(Y1959, by = 'binomial')

extir <- subset(ext, X1959 == 1)

new <- ne %>% 
  left_join(Y1959, by = "binomial")

ex59<- ext <- ex %>%
  mutate(binomial = Palacio.et.al.2019) %>% 
  left_join(Y1959, by = 'binomial')
recol<- subset(ex59, X1959 == 1)

pca_plot_1959 <- ggplot(dcc_1959, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = Y1959, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = new, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "deepskyblue") +
  # add arrows
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_1959, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_1959, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_1959, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  
  # axis labels - see comp_var
  labs(x = "PC1 - Body size (55%)", y = "PC2 - Dispersal ability (13%)") +
  xlim(-5,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20))+
  geom_text(x= 5, y=8, label="1950s", color="black", size = 5)

# display plot
windows()
pca_plot_1959

# plot 1994

ext <- ex %>%
  mutate(binomial = Palacio.et.al.2019) %>% 
  left_join(Y1994, by = 'binomial')

extir <- subset(ext, X1994 == 1)

new <- ne %>% 
  left_join(Y1994, by = "binomial")

ex94<- ext <- ex %>%
  mutate(binomial = Palacio.et.al.2019) %>% 
  left_join(Y1994, by = 'binomial')
recol<- subset(ex94, X1994 == 1)

pca_plot_1994 <- ggplot(dcc_1994, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = Y1994, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = new, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "deepskyblue") +
  # add arrows
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_1994, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_1994, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_1994, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Body size (55%)", y = "PC2 - Dispersal ability (13%)") +
  xlim(-5,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  geom_text(x= 5, y=8, label="1990s", color="black", size = 5)

# display plot
windows()
pca_plot_1994


# plot 2016

ext <- ex %>%
  mutate(binomial = Palacio.et.al.2019) %>% 
  left_join(Y2016, by = 'binomial')

extir <- subset(ext, X2016 == 1)

new <- ne %>% 
  left_join(Y2016, by = "binomial")

ex16<- ext <- ex %>%
  mutate(binomial = Palacio.et.al.2019) %>% 
  left_join(Y2016, by = 'binomial')
recol<- subset(ex16, X2016 == 1)


pca_plot_2016 <- ggplot(dcc_2016, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = Y2016, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = recol, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.9, colour = "black") +
  geom_point(data = new, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "deepskyblue") +
  # add arrows
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_2016, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_2016, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_2016, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Body size (55%)", y = "PC2 - Dispersal ability (13%)") +
  xlim(-5,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20))+
  geom_text(x=5, y=8, label="2010s", color="black", size = 5)

# display plot
windows()
pca_plot_2016

#################HYPER VOLUMES##################################


#### Years of data (hypervolume overlap) using PCA

# dataframes of trait data (PCA scores)
#Y1911
Y1911A <- Y1911 %>% 
  select(binomial, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5, Comp.6, Comp.7, Comp.8, Comp.9) %>% 
  tibble::column_to_rownames(var = "binomial")

#Y2016
Y2016A <- Y2016 %>% 
  select(binomial, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5, Comp.6, Comp.7, Comp.8, Comp.9) %>% 
  tibble::column_to_rownames(var = "binomial")

# hypervolumes: Takes a long time to run.
set.seed(3)
Y1911_hyper_mi <- hypervolume::hypervolume_svm(Y1911A, name = "1911")
saveRDS(Y1911_hyper_mi, "1911_hyper_mi9.rds")

set.seed(3)
Y2016_hyper_mi <- hypervolume::hypervolume_svm(Y2016A, name = "2016")
saveRDS(Y2016_hyper_mi, "2016_hyper_mi9.rds")

# load: hyper volumes
#Y1911_hyper_mi <- readRDS("1911_hyper_mi9.rds")
#Y2016_hyper_mi <- readRDS("2016_hyper_mi9.rds")


# set hypervolume comparisons 
#1911 vs 2016
set.seed(3)
bm_set_mi_1911vs2016 <- hypervolume::hypervolume_set(Y1911_hyper_mi9, Y2016_hyper_mi9, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_1911vs2016)

# summarise volumes 1911 vs 2016

hypervolume::get_volume(bm_set_mi_1911vs2016)

# Plot 3D hypervolumes 
#1911 vs 2016
bm_set_mi_1911vs2016@HVList$HV1 <- NULL
bm_set_mi_1911vs2016@HVList$HV2 <- NULL
bm_set_mi_1911vs2016@HVList$Union <- NULL

# trait names
names <- c("Size", "Dispersal ability","Habitat breadth", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")

colnames(bm_set_mi_1911vs2016@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_1911vs2016@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_1911vs2016@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_1911vs2016, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_1911vs2016@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_1911vs2016@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique 1911")

uni2 <- data.frame(bm_set_mi_1911vs2016@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique 2016")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume.csv")

colnames(bm_set_mi_1911vs2016@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_1911vs2016@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_1911vs2016@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume.csv", h = T)

pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~Comp.1, y = ~(-Comp.2), z = ~(-Comp.3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'Body size'),
                                   yaxis = list(title = 'Dispersal ability'),
                                   zaxis = list(title = 'Habitat breadth')))

fig


#2. Functional diversity metrics
####################################################################################
#Functional Diversity Metrics with FD
####################################################################################
library(factoextra)
library(FD)
library(ks)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(data.table)
library(picante)
library(tibble)
library(iNEXT)
library(funrar)

TR <- read.csv("~/Documents/Funcional San Antonio/TABLA_RASGOS_FINAL.csv", h = T)

TR <- TR %>% 
  mutate(SD = 1) %>% 
  mutate(log_weight = log(weight))


#z Transform traits
# function to z-transform data
scale_z <- function(x){
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}


trait_z <- TR %>% 
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(TotalCulmen)) %>% 
  mutate(Bill_w_z = scale_z(BillWidth)) %>% 
  mutate(Wing_l_z = scale_z(WingLength)) %>% 
  mutate(Tail_l_z = scale_z(TailLength)) %>% 
  mutate(Tars_z = scale_z(Tarsus)) %>% 
  mutate(HWI_z = scale_z(HWI)) %>% 
  mutate(Hab_br_z = scale_z(Hab_breadth)) %>%
  mutate(Gen_z = scale_z(Gen_length)) %>%
  as.data.frame

tr_mi_z <- trait_z %>% 
  select(binomial,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z)

tr_mi_z1 <- tr_mi_z[order(tr_mi_z$binomial),]

tr_mi_z2 <- data.frame(tr_mi_z1[,-1], row.names=tr_mi_z1[,1])



#Species Year Matrix
dat <- read.csv("SiteSpeciesMatrix.csv", h = T)

dat <- dat[order(dat$binomial),]%>% 
  right_join(tr_mi_z1, by = "binomial") %>% 
  select(binomial, Abundance, Year)

mat <- dat %>% 
  pivot_wider(names_from =  binomial, values_from = Abundance) %>% 
  column_to_rownames(var = "Year") 

mat[is.na(mat)] <- 0


fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = TRUE, calc.FRic = TRUE, m = "max", stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES

#Estimating Standard Effect Sizes with Null model (takes long time to run)


#__(a)_Unconstrained model_####
#____(i)_calculate FD values for random communities----

mat2 <- dat %>% 
  pivot_wider(names_from =  Year, values_from = Abundance) %>% 
  column_to_rownames(var = "binomial")
mat2[is.na(mat2)] <- 0

#mat2 <- as.matrix(mat2)

fdnullall=NULL
for (i in 1:999){  #number of iterations
  randomcom=randomizeMatrix(t(mat2), null.model='richness',iterations = 100)  #Richness maintains species richness
  randomcom <- randomcom[,sort(colnames(randomcom))]
  
  if(0 %in% colSums(randomcom)){ 
    NULL
  } else {
    
    #FRic
    FD0   = dbFD(tr_mi_z2, randomcom,w.abun = TRUE, calc.FRic = TRUE, m = "max", stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE) 
  }
  
  fdcom = as.data.frame(FD0)
  fdnullall= rbind(fdnullall, fdcom)
  
}

write.table(fdnullall,'FDnullmodel.txt')


#2. Boxplots of null model values 

NULL_TAB <- read.table("FDnullmodel.txt")

vec <- rep(c(1911, 1959, 1994, 2016), length.out = 3996)
vec <- data.frame(Year = vec)

dt <- cbind(NULL_TAB, vec) 

#Pvalues = (number of Fnull < Fobs + 1)/(999 + 1)
Y1 <- subset(dt, dt$Year == "1911")
Y2 <- subset(dt, dt$Year == "1959")
Y3 <- subset(dt, dt$Year == "1994")
Y4 <- subset(dt, dt$Year == "2016")
##FRic Pvalues
ob1 <- nrow(Y1[Y1$FRic< 0.59, ]) + 1
ob2 <- nrow(Y2[Y2$FRic< 0.24, ]) + 1
ob3 <- nrow(Y3[Y3$FRic< 0.27, ]) + 1
ob4 <- nrow(Y4[Y4$FRic< 0.45, ]) + 1


P1911 <- ob1/(999+1)# P = 0.193
P1959 <- ob2/(999+1)# P = 0.006
P1994 <- ob3/(999+1) # P = 0.128
P2016 <- ob4/(999+1) # P = 0.021


##FDis P values

##FDis Pvalues
ob1 <- nrow(Y1[Y1$FDis< 2.39, ]) + 1
ob2 <- nrow(Y2[Y2$FDis< 2.19, ]) + 1
ob3 <- nrow(Y3[Y3$FDis< 2.26, ]) + 1
ob4 <- nrow(Y4[Y4$FDis< 2.35, ]) + 1


P1911 <- ob1/(999+1)# P = 0.135
P1959 <- ob2/(999+1)# P = 0.003
P1994 <- ob3/(999+1) # P = 0.045
P2016 <- ob4/(999+1) # P = 0.079

#Boxplots Observed vs Null model
#FRic
windows()
p <- ggplot(dt, aes(x=as.factor(Year), y=FRic)) + 
  geom_boxplot(colour = "gray40", outlier.shape = NA)+ 
  geom_jitter(shape=16, size = 0.1,position=position_jitter(0.1), colour = "gray60") +
  geom_segment(aes(x = 0.55, y = 0.5969007, xend = 1.45, yend = 0.5969007), color = "red", size = 2)+
  geom_segment(aes(x = 1.5, y = 0.2409973, xend = 2.45, yend = 0.2409973), color = "red", size = 2)+
  geom_segment(aes(x = 2.5, y = 0.2754841, xend = 3.45, yend = 0.27548413), color = "red", size = 2)+
  geom_segment(aes(x = 3.5, y = 0.4557709, xend = 4.45, yend = 0.4557709), color = "red", size = 2)+
  geom_text(x=1, y=1, label="p = 0.193")+
  geom_text(x=2, y=1, label="p = 0.006")+
  geom_text(x=3, y=1, label="p = 0.128")+
  geom_text(x=4, y=1, label="p = 0.021")+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Richness")+ scale_x_discrete(labels = c("1910s", "1950s", "1990s", "2010s"))

p



#FDis
windows()
p <- ggplot(dt, aes(x=as.factor(Year), y=FDis)) + 
  geom_boxplot(colour = "gray40", outlier.shape = NA)+
  geom_jitter(shape=16, size = 0.1,position=position_jitter(0.1), colour = "gray60") +
  geom_segment(aes(x = 0.55, y = 2.390626, xend = 1.45, yend = 2.390626), color = "red", size = 2)+
  geom_segment(aes(x = 1.5, y = 2.190313, xend = 2.45, yend = 2.190313), color = "red", size = 2)+
  geom_segment(aes(x = 2.5, y = 2.265837, xend = 3.45, yend = 2.265837), color = "red", size = 2)+
  geom_segment(aes(x = 3.5, y = 2.354284, xend = 4.45, yend = 2.354284), color = "red", size = 2)+
  geom_text(x=1, y=2.82, label="p = 0.135")+
  geom_text(x=2, y=2.82, label="p = 0.003")+
  geom_text(x=3, y=2.82, label="p = 0.045")+
  geom_text(x=4, y=2.82, label="p = 0.079")+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Dispersion")+ scale_x_discrete(labels = c("1910s", "1950s", "1990s", "2010s"))

p

#################

#Estimate SES - Standard Effect Size SES = observed mean - mean null / sd null

dt_Mean <- dt %>% 
  group_by(Year) %>%
  summarise_at(vars(FRic, FDis), funs(mean(., na.rm=TRUE))) %>% 
  transmute(Year = Year, FR_M = FRic, FDis_M = FDis)

dt_SES <- dt %>% 
  group_by(Year) %>%
  summarise_at(vars(FRic, FDis), funs(sd(., na.rm=TRUE))) %>% 
  transmute(Year = Year, FR_SD = FRic, FDis_SD = FDis) %>% 
  left_join(dt_Mean, by = "Year") %>% 
  mutate(FRic =c(0.5969007, 0.2409973, 0.2754841, 0.4557709), FEve = c(0.6782481, 0.6800727, 0.6761855, 0.6850895), FDiv = c(0.6925992, 0.6889875, 0.6887308, 0.7042643 ), FDis = c(2.390626, 2.190313, 2.265837, 2.354284 )) %>% 
  mutate(SES_FRIc = ((FRic - FR_M)/FR_SD),SES_FDis = ((FDis - FDis_M)/FDis_SD))

#Plots SES for FRic and FDis

ses <- dt_SES

windows()
p <- ggplot(ses, aes(x=as.factor(Year), y=SES_FRIc)) + 
  geom_boxplot(colour = "gray40", outlier.shape = NA)+ 
  geom_jitter(shape=16, size = 0.1,position=position_jitter(0.1), colour = "gray60") +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 1.5)+
  
  scale_y_continuous(limit = c(-6, 2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  
  labs(x = "", y = "SES Functional Richness")+ scale_x_discrete(labels = c("1910s", "1950s", "1990s", "2010s"))

p

windows()
p <- ggplot(ses, aes(x=as.factor(Year), y=SES_FDis)) + 
  geom_boxplot(colour = "gray40", outlier.shape = NA)+ 
  geom_jitter(shape=16, size = 0.1,position=position_jitter(0.1), colour = "gray60") +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 1.5)+
  scale_y_continuous(limit = c(-6, 2))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "SES Functional Dispersion")+ scale_x_discrete(labels = c("1910s", "1950s", "1990s", "2010s"))

p


#3. TPD and species group comparisons
#####################################################################################

library(TPD)

###Single PCA with all traits and species regardless of years

tr_mi_z <- trait_z %>% 
  select(binomial,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z) %>% 
  column_to_rownames(var = "binomial")

#Run PCA
pcaTotal <- princomp(tr_mi_z, cor = TRUE, scores = TRUE)
summary(pcaTotal)

#Scores
scoresTotal <- as.data.frame(pcaTotal$scores) %>% 
  tibble::rownames_to_column("binomial") %>% 
  left_join(TR, by = "binomial")

#Community TPDc using PCAs as traits 
#Scores
scoresTotal <- as.data.frame(pcaTotal$scores) %>% 
  tibble::rownames_to_column("binomial") %>% 
  left_join(trait, by = "binomial")

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal %>% 
  mutate(SD = 1) %>%
  select(binomial,Comp.1, Comp.2, Comp.3,id, POP,Hist, Current, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP
ABUN <- TotalPCA %>% 
  select(binomial, Current,POP) %>% 
  pivot_wider(names_from = binomial, values_from = Current) %>% 
  column_to_rownames('POP')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$binomial, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP
ABUN <- TotalPCA %>% 
  select(binomial, Hist,POP) %>% 
  pivot_wider(names_from = binomial, values_from = Hist) %>% 
  column_to_rownames('POP')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$binomial, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP
ABUN <- TotalPCA %>% 
  select(binomial, Hist,POP) %>% 
  pivot_wider(names_from = binomial, values_from = Hist) %>% 
  column_to_rownames('POP')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$binomial, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-smirnov test
library("sm")

A <- subset(TotalPCA1, TotalPCA1$POP == "A")
B <- subset(TotalPCA1, TotalPCA1$POP == "B")
C <- subset(TotalPCA1, TotalPCA1$POP == "C")

#KS test
ks.test(A$Comp.1, B$Comp.1)#P = 0.04

ks.test(A$Comp.2, B$Comp.2)# P = 0.0008

ks.test(A$Comp.3, B$Comp.3)# P = 0.003


##Graphs of TPDs for PC1, 2 and 3

#Multiply PC2 and 3 by -1 to ease interpretation of increasing values.
TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 *-1) %>% 
  mutate(Comp3M = Comp.3 *-1)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (55%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ geom_text(x = 1, y = 0.40, label = " p = 0.04" , size = 6)+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (13%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ geom_text(x = 2.5, y = 0.55, label = " p = 0.0008", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Habitat breadth (11%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none") + geom_text(x = 3, y = 0.62, label = " p = 0.003", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p

#4. Functional Uniqueness and distinctiveness
#####################################################################################

library(funrar)

# load data
trait <- read_csv("trait data.csv")

#Species Year Matrix
dat <- read_csv("~/Documents/Funcional San Antonio/TABLA_RASGOS_FINAL.csv") 

mat <- dat %>% 
  pivot_wider(names_from =  binomial, values_from = Abundance)

mat[is.na(mat)] <- 0

mat <- as.matrix(mat)

# Extirpated and new species

ex <- #filter species with 'current abundance = 0' from dataset

new <- #filter species which were absent from first time period (1910s))

#z Transform traits
# function to z-transform data
scale_z <- function(x){
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}


trait_z <- trait %>% 
  mutate(log_weight = log(weight)) %>% 
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(TotalCulmen)) %>% 
  mutate(Bill_w_z = scale_z(BillWidth)) %>% 
  mutate(Wing_l_z = scale_z(WingLength)) %>% 
  mutate(Tail_l_z = scale_z(TailLength)) %>% 
  mutate(Tars_z = scale_z(Tarsus)) %>% 
  mutate(HWI_z = scale_z(HWI)) %>% 
  mutate(Hab_br_z = scale_z(Hab_breadth)) %>%
  mutate(Gen_z = scale_z(Gen_length)) %>%
  as.data.frame

tra <- trait_z %>% 
  select(binomial,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z,Hab_br_z, Gen_z) %>% 
  column_to_rownames(var = "binomial")

#Compute distance matrix
dist_mat <- compute_dist_matrix(tra, metric = "gower")

# Compute distinctiveness for each species on each site
di_df = distinctiveness_stack(com_df = dat,  # The site x species table
                              sp_col = "binomial",  # Name of the species column
                              com = "Year",  # Name of the community column
                              abund = "Abundance",  # Relative abundances column (facultative)
                              dist_matrix = dist_mat)  # Functional Distance matrix

head(di_df)

Di_Sum <- di_df %>% 
  group_by(binomial) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

Di_Sum <- Di_Sum[order(Di_Sum$Di),] %>% 
  mutate(ID = c(1:233)) %>% 
  right_join(trait, by = "binomial")

### GRaph POP
windows()
p <- ggplot(Di_Sum, aes(x=ID, y=Di, fill = POP)) + 
  geom_bar(stat = "identity")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), breaks=c("Extirpated","New additions","Shared species"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness - Di")

p
#Boxplot  Di POP

windows()
p <- ggplot(Di_Sum, aes(x=POP, y=Di, fill = POP)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  ylim(0.1,0.6)+
  geom_text(x = 2, y = 0.58, label = "p = 0.006")+
  geom_text(x = 1.4, y = 0.53, label = "p = 0.644")+
  geom_text(x = 2.6, y = 0.53, label = "p = 0.542")+
  geom_segment(aes(x = 1, y = 0.56, xend = 3, yend = 0.56))+
  geom_segment(aes(x = 1, y = 0.51, xend = 1.8, yend = 0.51))+
  geom_segment(aes(x = 2.2, y = 0.51, xend = 3, yend = 0.51))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Distinctiveness")

p

#######GRaph F Distinctiveness by period

ext <- ex %>% 
  mutate(binomial = `Palacio et al 2019`) %>% 
  left_join(di_df, by = "binomial")

Y1 <- subset(ext, ext$Year == 1911)

Y2 <- subset(ext, ext$Year == 1959)

Y3 <- subset(ext, ext$Year == 2016)


new <- new %>% left_join(di_df, by = "binomial")
N1 <- subset(new, new$Year == 1959)
N2 <- subset(new, new$Year == 1994)
N3 <- subset(new, new$Year == 2016)

d <- ext %>% 
  group_by(Year) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

windows()
p <- ggplot(di_df, aes(x=as.factor(Year), y=Di)) + 
  geom_violin(trim = FALSE, colour = "gray40")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), colour = "gray40") +
  geom_jitter(data = Y1, aes(x = as.factor(Year), y = Di), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = Y2, aes(x = as.factor(Year), y = Di), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = N1, aes(x = as.factor(Year), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_jitter(data = N2, aes(x = as.factor(Year), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_jitter(data = N3, aes(x = as.factor(Year), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue", position=position_jitter(0.2)) + 
  geom_jitter(data = Y3, aes(x = as.factor(Year), y = Di), size = 3, alpha = 1, colour = "black", position=position_jitter(0.2)) + 
  geom_segment(aes(x = 0.5, y = 0.1619373, xend = 1.5, yend = 0.1619373))+
  geom_segment(aes(x = 1.5, y = 0.1377468, xend = 2.5, yend = 0.1377468))+
  geom_segment(aes(x = 2.5, y = 0.1315327, xend = 3.5, yend = 0.1315327))+
  geom_segment(aes(x = 3.5, y = 0.1541466, xend = 4.5, yend = 0.1541466))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness")+ scale_x_discrete(labels = c("1900s", "1950s", "1990s", "2010s"))

p


#Compute Uniqueness for each species
ui_df = uniqueness_stack(dat, "binomial", dist_mat)

head(ui_df)

### GRaph POP

Ui_Sum <- di_ui %>% 
  group_by(binomial) %>%
  summarise_at(vars(Ui), funs(mean(., na.rm=TRUE)))

Ui_Sum <- Ui_Sum[order(Ui_Sum$Ui),] %>% 
  mutate(ID = c(1:233)) %>% 
  right_join(trait, by = "binomial")

#Boxplot  Ui POP

windows()
p <- ggplot(Ui_Sum, aes(x=POP, y=Ui, fill = POP)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,0.25)+
  geom_text(x = 2, y = 0.25, label = "p = 0.002")+
  geom_text(x = 1.4, y = 0.23, label = "p = 0.434")+
  geom_text(x = 2.6, y = 0.23, label = "p = 0.332")+
  geom_segment(aes(x = 1, y = 0.24, xend = 3, yend = 0.24))+
  geom_segment(aes(x = 1, y = 0.22, xend = 1.8, yend = 0.22))+
  geom_segment(aes(x = 2.2, y = 0.22, xend = 3, yend = 0.22))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Uniqueness")

p

#violin plot

AA <- AA[!AA$Di == "NA",]

u<- AA %>% 
  group_by(Year) %>%
  summarise_at(vars(Ui), funs(mean(., na.rm=TRUE)))


windows()
p <- ggplot(di_ui, aes(x=as.factor(Year), y=Ui)) + 
  geom_violin(trim = FALSE, colour = "gray40")+ 
  geom_point(shape=16, position=position_jitter(0.2), colour = "gray40") +
  geom_jitter(data = Y1, aes(x = as.factor(Year), y = Ui), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = Y2, aes(x = as.factor(Year), y = Ui), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = N1, aes(x = as.factor(Year), y = Ui), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_jitter(data = N2, aes(x = as.factor(Year), y = Ui), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_jitter(data = N3, aes(x = as.factor(Year), y = Ui), size = 3, alpha = 0.5, colour = "cornflowerblue", position=position_jitter(0.2)) + 
  geom_jitter(data = Y3, aes(x = as.factor(Year), y = Ui), size = 3, alpha = 1, colour = "black", position=position_jitter(0.2)) + 
  geom_segment(aes(x = 0.5, y = 0.045, xend = 1.5, yend = 0.045))+
  geom_segment(aes(x = 1.5, y = 0.042, xend = 2.5, yend = 0.042))+
  geom_segment(aes(x = 2.5, y = 0.0451, xend = 3.5, yend = 0.0451))+
  geom_segment(aes(x = 3.5, y = 0.043, xend = 4.5, yend = 0.043))+  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional uniqueness")+ scale_x_discrete(labels = c("1900s", "1950s", "1990s", "2010s"))

p



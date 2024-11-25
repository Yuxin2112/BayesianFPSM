#########################################################################
## This R scrip contains all nessassery code for 
## plotting results of LLE and Cr in each SLAs 


## created 2023-3-31 by Yuxin Huang 
## modify 2024-10-11 by Yuxin Huang (for reviewer's comments)
##########################################################################

## loading pekages

setwd("C:\\Users\\n11117761\\Work\\Data")
library(viridis)
library(lattice)
library(haven)
library(dplyr)
library(ggplot2)
library(spData)
library(spdep)
library(rgdal)
library(sf)
library(terra)
library(rgeos)
library(maptools)
library(sp)
library(scales)
library(gridExtra)
## loading data
lle <- read_dta("C:/Users/n11117761/Work/Data/FinalResult/estimates_dataset/lle_est_gamma05.dta")

cif <- read_dta("C:/Users/n11117761/Work/Data/FinalResult/estimates_dataset/cif_est_gamma05.dta")


#lle <- read_dta("C:/Users/n11117761/Work/Data/WinBugs/further_analysis/model1/lle_est_cage.dta")
## prepare data
# calculate the average estimation in each SLA
lle_avg <- lle %>%
  group_by(grid) %>%
  summarise_at(c("exp","lle_2","lle_50","lle_97","obs_2","obs_50","obs_97"), mean) %>%
  rename(id=grid)


cif_avg <- cif %>%
  group_by(grid) %>%
  summarise_at(c("cr_2","netm_2","cr_50","netm_50","cr_97","netm_97"), mean) %>%
  rename(id=grid)

## loading Shapfile: QLD SLA 2006
map_qld_06  <- readOGR("C:/Users/n11117761/Work/Data/MapInfo file/SLA_QLD_06/SLA_QLD_06.shp"
)

## prepare map data  for plotting
map_df <- fortify(map_qld_06)            # as data frame object(for merging with estimation)
map_df$id <- as.numeric(map_df$id) + 1

map.border <- unionSpatialPolygons(map_qld_06, IDs = rep(1,478)) #map of border (for ggplot) 
map.border <- fortify(map.border)

## creat sub-map if Brisbane

map.proj <- proj4string(map_qld_06)

# User-defined function for creating inset windows
get.inset <- function(x1, x2, y1, y2, map_qld_06, map.proj){
  Inset.window <- as(raster::extent(x1, x2, y1, y2), "SpatialPolygons")
  proj4string(Inset.window) <- map.proj
  map.inset <- gIntersection(map_qld_06, Inset.window, byid = TRUE, drop_lower_td = TRUE, 
                             id = sapply(map_qld_06@polygons, function(x) x@ID))
  return(map.inset)
}

map_brisbane <- get.inset(152.6, 153.6, -28, -27, map_qld_06, map.proj)

map.border.brisbane <- unionSpatialPolygons(map_brisbane, IDs = rep(1, length(map_brisbane)))
map_brisbane_df <- fortify(map_brisbane)
map_brisbane_df$id <- as.numeric(map_brisbane_df$id) + 1

rm(get.inset)

################################ map of median LLE ######################################

## define color and value scales
#Fill.colours.lle <- c("#00004B","#33628DFF","#74ADD1","#ABD9E9","#E0F3F8"
#                    ,"#FFFFCC",
#                    "#FEE090","#FDAE61","#F46D43","#D73027","#A50026")
# not used at the end

#Fill.colours.lle <- c("#d4e6f5", "#abcfeb", "#96c3e6", "#81b7e1", "#276da1", "#1d5178","#0e2b3a")
Fill.colours.viridis <- viridis(9, direction = -1)
#Fill.colours.lle.d <- c("gray100","gray90","gray53","gray30","gray26")

Fill.colours.lle.d <- c("#d4e6f5", "#abcfeb", "#96c3e6", "#81b7e1", "#276da1", "#1d5178","#0e2b3a")

# Values corresponding to the fill colours
########################################################################
Fill.value.lle.50 <- rescale(c(min(lle_avg$lle_50),
                                  quantile(lle_avg$lle_50,name = FALSE,probs = seq(0.05,0.4,length = 4)),
                                  median(lle_avg$lle_50),
                                  quantile(lle_avg$lle_50,name = FALSE,probs = seq(0.65,0.99,length = 4)),
                                  max(lle_avg$lle_50)),
                                from = range(as.numeric(lle_avg$lle_50)))
Fill.value.obs.50 <- rescale(c(min(lle_avg$obs_50),
                               quantile(lle_avg$obs_50,name = FALSE,probs = seq(0.05,0.5,length = 4)),
                               mean(lle_avg$obs_50),
                               quantile(lle_avg$obs_50,name = FALSE,probs = seq(0.65,0.99,length = 4)),
                               max(lle_avg$obs_50)),
                             from = range(as.numeric(lle_avg$obs_50)))

Fill.value.exp <- rescale(c(min(lle_avg$exp),
                               quantile(lle_avg$exp,name = FALSE,probs = seq(0.05,0.5,length = 4)),
                               mean(lle_avg$exp),
                               quantile(lle_avg$exp,name = FALSE,probs = seq(0.65,0.99,length = 4)),
                               max(lle_avg$exp)),
                             from = range(as.numeric(lle_avg$exp)))


#Fill.values.lle.p <- c(min(lle_avg$lle_2),5,10,15,20,25,30,35,max(lle_avg$lle_97))
#Fill.values.lle.p.r <- scales::rescale(Fill.values.lle.p)

#rescale values
#Fill.values.lle.exp.r <- rescale(Fill.values.lle.exp, from = range(as.numeric(lle_avg$exp)))


## merge values with shapfile dataframe
lle_avg$d <- abs(lle_avg$lle_97-lle_avg$lle_2)


lle_grid <- inner_join(map_df,lle_avg,by = "id") 
lle_brisbane <- inner_join(map_brisbane_df, lle_avg, by = "id")

# for credible interval
Fill.value.lle.d <- rescale(c(min(lle_avg$d),
                               quantile(lle_avg$d,name = FALSE,probs = seq(0.05,0.4,length = 4)),
                               median(lle_avg$d),
                               quantile(lle_avg$d,name = FALSE,probs = seq(0.65,0.99,length = 4)),
                               max(lle_avg$d)),
                             from = range(as.numeric(lle_avg$d)))

Fill.values.lle.d <- c(0,1,1.5,2,2.5,3,5,7,8.3)
Fill.values.lle.d.r <- rescale(Fill.values.lle.d, from = range(as.numeric(lle_avg$d)))


## base map
gg.base <- ggplot(data = NULL, aes(x = long, y = lat, group = group)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 15)) +
  coord_map() 

## layer map
gg.lle50 <- gg.base + 
  geom_polygon(data = map.border, fill = "grey30", color = "black", size = 0.8) +
  geom_polygon(data = lle_grid, aes(fill = lle_50), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.viridis, 
                       values = Fill.value.lle.50,
                       breaks=c(min(lle_avg$lle_50),5,10,mean(lle_avg$lle_50),23,31,max(lle_avg$lle_50)),
                       labels=c("lower","5","10","similar","23","31","higher"),
                       limits = range(lle_avg$lle_50)
  ) +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "red", fill = NA,size = 1)

gg.obs50 <- gg.base + 
  geom_polygon(data = map.border, fill = "grey30", color = "black", size = 0.8) +
  geom_polygon(data = lle_grid, aes(fill = obs_50), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.viridis, 
                       values = Fill.value.obs.50,
                       breaks=c(min(lle_avg$obs_50),4,6,mean(lle_avg$obs_50),11,14,max(lle_avg$obs_50)),
                       labels=c("lower","4","6","similar","11","14","higher"),
                       limits = range(lle_avg$obs_50)
  ) +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "red", fill = NA,size = 1)

gg.exp <- gg.base + 
  geom_polygon(data = map.border, fill = "grey30", color = "black", size = 0.8) +
  geom_polygon(data = lle_grid, aes(fill = exp), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.viridis, 
                       values = Fill.value.exp,
                       breaks=c(min(lle_avg$exp),8,16,mean(lle_avg$exp),35,45,max(lle_avg$exp)),
                       labels=c("lower","8","16","similar","35","45","higher"),
                       limits = range(lle_avg$exp)
  ) +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "red", fill = NA,size = 1)

gg.lle_d <- gg.base + 
  geom_polygon(data = map.border, fill = "grey100", color = "black", size = 0.8) +
  geom_polygon(data = lle_grid, aes(fill = d), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.lle.d, 
                       values = Fill.values.lle.d.r,
                       limits = range(lle_avg$d),
                       breaks = c(min(lle_avg$d),0.8,1.8,mean(lle_avg$d),4.5,6.5,max(lle_avg$d)),
                       labels = c("lower",0.8,1.8,"similar",4.5,6.5,"higher")
  ) +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "red", fill = NA,size = 1)


## sub-map of Brisbane
gg.base.inset <- gg.base + guides(fill = FALSE, alpha = FALSE) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.viridis, 
                       values = Fill.value.obs.50,
                       limits = range(lle_avg$obs_50))
gg.inset.brisbane <- gg.base.inset + 
  geom_polygon(data = map.border.brisbane, fill = "grey90", color = "black", size = 0.8) +
  geom_polygon(data = lle_brisbane, aes(fill = obs_50), color = NA)

gg.base.inset.d <- gg.base + guides(fill = FALSE, alpha = FALSE) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.lle.d, 
                       values = Fill.values.lle.d.r,
                       limits = range(lle_avg$d))
gg.inset.brisbane.d <- gg.base.inset.d + 
  geom_polygon(data = map.border.brisbane, fill = "grey90", color = "black", size = 0.8) +
  geom_polygon(data = lle_brisbane, aes(fill = d), color = NA)

## arrange
png(filename = "C:/Users/n11117761/Work/Data/FinalResult/exp_revise.png")
grid.arrange(
  gg.lle_d,
  gg.inset.brisbane50,
  nrow = 1,
  ncol = 2,
  widths = c(10,5)
)
dev.off()

## high res map
tiff(filename = "C:/Users/n11117761/Work/Manuscript/final_review/exp_revise.tiff",
     width = 3200, height = 3200,units = "px",res=800, compression = 'lzw')
grid.arrange(
  gg.exp,
  gg.inset.brisbane,
  nrow = 1,
  ncol = 2,
  widths = c(10,5)
)
dev.off()
tiff(filename = "C:/Users/n11117761/Work/Manuscript/final_review/obs_revise.tiff",
     width = 3200, height = 3200,units = "px",res=800, compression = 'lzw')
grid.arrange(
  gg.obs50,
  gg.inset.brisbane,
  nrow = 1,
  ncol = 2,
  widths = c(10,5)
)
dev.off()

#combine map and CI map in one single plot
lle_full <- grid.arrange(
  gg.lle50,
  gg.inset.brisbane50,
  nrow = 1,
  ncol = 2,
  widths = c(10,5))
lle_d_full <- grid.arrange(
  gg.lle_d,
  gg.inset.brisbane.d,
  nrow = 1,
  ncol = 2,
  widths = c(10,5))

tiff(filename = "C:/Users/n11117761/Work/Manuscript/final_review/lle_review.tiff",
     width = 8000, height = 3200,units = "px",res=800, compression = 'lzw')
grid.arrange(
  lle_full,
  lle_d_full,
  nrow = 1,
  ncol = 2
)
dev.off()
############################### map of Cr #######################################
#Fill.colours.cif <- c("#fffef7", "#ffffc4","#ffffd4","#fee391", "#fec44f","#ff9e33","#d95f0e", "#be3d00", "#993404")

#Fill.colours.cif2 <- c("#fffef7","#fff7c4", "#ffffc4","#ffffd4","#fee391", "#fec44f","#ff9e33","#d95f0e")
#Fill.colours.cif97 <- c("#ffffc4","#fee391", "#fec44f","#ff9e33","#d95f0e", "#be3d00", "#993404")

Fill.values.cif.p <- c(0.12,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.85)
Fill.values.cif.p.r <- scales::rescale(Fill.values.cif.p)

cif_avg$d <- abs(cif_avg$cr_97-cif_avg$cr_2)
cif_avg$netd <- abs(cif_avg$netm_97-cif_avg$netm_2)


cif_grid <- inner_join(map_df,cif_avg,by = "id") 
cif_brisbane <- inner_join(map_brisbane_df, cif_avg, by = "id")


Fill.values.cr.d <- c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.7)
Fill.values.cr.d.r <- rescale(Fill.values.cr.d, from = range(as.numeric(cif_avg$d)))


gg.cif50 <- gg.base + 
  geom_polygon(data = map.border, fill = "grey30", color = "black", size = 0.8) +
  geom_polygon(data = cif_grid, aes(fill = cr_50), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.viridis, 
                       values = Fill.values.cif.p.r,
                       breaks=c(min(cif_avg$cr_50),0.41,0.5,mean(cif_avg$cr_50),0.6,0.65,max(cif_avg$cr_50)),
                       labels=c("lower","0.41","0.5","similar","0.6","0.65","higher"),
                       limits = range(cif_avg$cr_50)
  ) +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "red", fill = NA,size = 1)
gg.netm50 <- gg.base + 
  geom_polygon(data = map.border, fill = "grey30", color = "black", size = 0.8) +
  geom_polygon(data = cif_grid, aes(fill = netm_50), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.viridis, 
                       values = Fill.values.cif.p.r,
                       breaks=c(min(cif_avg$netm_50),0.43,0.53,mean(cif_avg$netm_50),0.67,0.77,max(cif_avg$netm_50)),
                       labels=c("lower","0.43","0.53","similar","0.67","0.77","higher"),
                       limits = range(cif_avg$netm_50)
  ) +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "red", fill = NA,size = 1)

gg.cr_d <- gg.base + 
  geom_polygon(data = map.border, fill = "grey100", color = "black", size = 0.8) +
  geom_polygon(data = cif_grid, aes(fill = d), color = NA) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.lle.d, 
                       values = Fill.values.cr.d.r,
                       limits = range(cif_avg$d),
                       breaks = c(min(cif_avg$d),0.13,mean(cif_avg$d),0.3,0.45,max(cif_avg$d)),
                       labels = c("lower",0.13,"similar",0.3,0.45,"higher")
  ) +
  annotate("rect",
           xmin = 152.6 - 0.01, xmax = 153.6 + 0.01,
           ymin = -28 - 0.01, ymax = -27 + 0.01,
           colour = "red", fill = NA,size = 1)

gg.base.inset.cif50 <- gg.base + guides(fill = FALSE, alpha = FALSE) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.viridis, 
                       values = Fill.values.cif.p.r,
                       limits = range(cif_avg$cr_50))

gg.inset.brisbane50 <- gg.base.inset.cif50 + 
  geom_polygon(data = map.border.brisbane, fill = "grey90", color = "black", size = 0.8) +
  geom_polygon(data = cif_brisbane, aes(fill = cr_50), color = NA)

gg.base.inset.d <- gg.base + guides(fill = FALSE, alpha = FALSE) +
  scale_fill_gradientn("", 
                       colours = Fill.colours.lle.d, 
                       values = Fill.values.cr.d.r,
                       limits = range(cif_avg$d))

gg.inset.brisbane.d <- gg.base.inset.d + 
  geom_polygon(data = map.border.brisbane, fill = "grey90", color = "black", size = 0.8) +
  geom_polygon(data = cif_brisbane, aes(fill = d), color = NA)

png(filename = "C:/Users/n11117761/Work/Data/FinalResult/cif50_R.png")
grid.arrange(
  gg.cr_d,
  gg.inset.brisbane50,
  nrow = 1,
  ncol = 2,
  widths = c(10,5)
)
dev.off()

## high res map
tiff(filename = "C:/Users/n11117761/Work/Data/FinalResult/cr_97.tiff",
     width = 3200, height = 3200,units = "px",res=800, compression = 'lzw')
grid.arrange(
  gg.cif97,
  gg.inset.brisbane97,
  nrow = 1,
  ncol = 2,
  widths = c(10,5)
)
dev.off()

# combine estimation map and CI map into one plot
cr_full <- grid.arrange(
  gg.cif50,
  gg.inset.brisbane50,
  nrow = 1,
  ncol = 2,
  widths = c(10,5))
cr_d_full <- grid.arrange(
  gg.cr_d,
  gg.inset.brisbane.d,
  nrow = 1,
  ncol = 2,
  widths = c(10,5))

tiff(filename = "C:/Users/n11117761/Work/Manuscript/final_review/cr_review.tiff",
     width = 8000, height = 3200,units = "px",res=800, compression = 'lzw')
grid.arrange(
  cr_full,
  cr_d_full,
  nrow = 1,
  ncol = 2
)
dev.off()

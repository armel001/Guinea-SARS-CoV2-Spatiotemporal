rm(list = ls(all=T))


####### VISUALISATION DES INTRODUCTIONS PAR REGION SUR UN DIAGRAM EN BAR####


######################

library(ggplot2)
library(dplyr)
library(cowplot)
library(ggtree)
library(tidyverse)
library(xlsx)
library(readxl)
library(lubridate)
library(ggalluvial)
library(maps)
library(ggmap)

###############################

niceFormat <- function(number) {
  formatC(number, format = 'f', big.mark = ',', digits = 0)
}
theme_USGS_box <- function(base_family = "", ...){
  theme_bw(base_family = base_family, ...) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5),
      axis.ticks.length = unit(-0.05, "in"),
      axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
      axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
      axis.ticks.x = element_blank(),
      #aspect.ratio = 1,
      legend.background = element_rect(color = "black", fill = "white"),
      text = element_text(size = 12, colour = "black"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'pt')
    )}


my_color <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
              '#6a3d9a','#fdbf6f','#ff7f00','#e31a1c',
              "#00868B",'#b15982','#34ace0',"#9ACD32",'#cab2d6',
              "chartreuse1", "#9CAFAA", "yellow3", "#8470FF", 
              "#8B5742", "#9ACD32",'#FAF1E4', "#CD4F39")

metadata <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/Cumul_des_introductions_avec_gps.xlsx")

cumul_intro <- metadata %>% select(Variant, Region) %>% group_by(Region, Variant) %>%
  summarise(count = n())

######
p2 <- cumul_intro %>%
  ggplot(aes(x = count, y = Region, fill = Variant)) +
  geom_col(width = 0.5, color = "grey22") +  # Use geom_col instead of geom_bar for pre-summarized data
  scale_fill_manual(name = 'Variants', values = my_color) +  # Use scale_fill_manual for fill aesthetic
  theme_bw() +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 10),  # Reduce legend text size
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  ) + labs(x = "Number of Introductions", y = "Region")
p2

#### MAP PLOT 

data <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/Cumul_des_introductions_avec_gps.xlsx")

data$Longitude = as.numeric(data$Longitude)
data$Latitude = as.numeric(data$Latitude)

Longitude = as.numeric(data$Longitude)
Latitude = as.numeric(data$Latitude)
str(data)
# Guinea coordinates
guinea_coords <- data.frame(
  Longitude = -9.696645,
  Latitude = 9.945587)

# Get a world map
world_map <- map_data("world")

map <- ggplot() +
  
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = '#f2f2f2', color = "grey1") +
  
  scale_color_manual(values = c(
    'Europe' = '#ffbe0b', 
    'Africa' = '#fb5607', 
    'Asia' = "#ff006e", 
    'North America' = '#8338ec', 
    'South America' = "#3a86ff", 
    'Oceania' = '#219ebc'
  )) +
  
  geom_curve(data = data, aes(x = Longitude, y = Latitude, 
                              xend = guinea_coords$Longitude, yend = guinea_coords$Latitude, 
                              color = Region), linewidth = 0.6) +
  
  geom_point(data = data, aes(x = Longitude, y = Latitude, color = Region), size = 1) +
  
  geom_point(data = guinea_coords, aes(x = Longitude, y = Latitude), color = "red", size = 2) +
  
  # ggtitle("Introductions of Variants of Concern into Guinea from March 2020 to December 2023") +
  
  theme_classic(base_line_size = 0.1) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 12, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    panel.border = element_blank(),
    panel.background = element_blank()
  )  +
  labs(x = "Latitude", y = "Longitude")


map


grid_region_intro <- plot_grid(p2, map, 
                               nrow = 2,
                               labels = c("A", "B"), 
                               label_size = 18) 

ggsave("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/figures/figure4_map-and-region-intro.svg", grid_region_intro, width = 15, height = 20, dpi = 300)





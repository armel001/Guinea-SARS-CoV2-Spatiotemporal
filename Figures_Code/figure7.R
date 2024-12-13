

install.packages("ggalluvial")

library(ggalluvial)
library(ggplot2)
library(dplyr)

data <- read.csv("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/1017-Guinea/annottated_tree_events.csv", sep = ";")
View(data)

my_color <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
              '#6a3d9a','#fdbf6f','#ff7f00','#e31a1c',
              "#00868B",'#b15982','#34ace0',"#9ACD32",'#cab2d6',
              "chartreuse1", "#9CAFAA", "yellow3", "#8470FF", 
              "#8B5742", "#9ACD32",'#FAF1E4', "#CD4F39")

ggplot(data = data,
       aes(axis1 = Origin,   # First variable on the X-axis
           axis2 = Destination, # Second variable on the X-axis
  # Third variable on the X-axis
           )) +
  geom_alluvium(aes(fill = Destination)) +
  scale_color_manual(values = my_color) +
  #geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Origin", "Destination"), 
                 expand = c(0.15, 0.05)) +
  theme_bw() 




# Generate the alluvial plot
flow <- ggplot(data = data,
       aes(axis1 = Origin,   # First variable on the X-axis
           axis2 = Destination)) +  # Second variable on the X-axis
  geom_alluvium(aes(fill = Origin), width = 0.1, alpha = 0.8) +  # Alluvium settings for aesthetics
  geom_stratum(width = 0.2, color = "grey") +   # Add stratum blocks
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 5, vjust = 0.3, hjust = 0.015, color = "black") +   # Text on strata
  scale_fill_manual(name = "Origin", values = my_color) +  # Use fill for colors
  scale_x_discrete(limits = c("Origin", "Destination"), 
                   expand = c(0.15, 0.05)) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),         # Remove Y-axis labels
    axis.ticks.y = element_blank(),        # Remove Y-axis ticks
    axis.line = element_blank(),           # Remove axis lines
    panel.grid = element_blank(),          # Remove grid lines
    plot.title = element_text(hjust = 0.5),# Center the plot title
    text = element_text(size = 14)         # Set text size
  ) 
  #labs(title = "Flow from Origin to Destination", x = "", y = "")

flow


ggsave("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/figures/flow-chart-introduction_guinea.svg", flow, dpi = 300, width = 16, height = 18)

orig <- data %>% group_by(Origin, Destination) %>% summarise(count = n(), .groups = "drop")

table(data$Origin)

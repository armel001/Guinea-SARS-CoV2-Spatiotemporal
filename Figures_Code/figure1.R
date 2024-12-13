rm(list = ls(all=T))

######################

library(ggplot2)
library(dplyr)
library(cowplot)
library(ggtree)
library(tidyverse)
library(xlsx)
library(readxl)
library(lubridate)

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

metadata <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/SARS-CoV-2_Distribution-Guinea/1038-global.metadata.xlsx")

metadata$date <- as.Date(metadata$date, format = "%Y-%m-%d")


data <- metadata
data$years_month <- floor_date(data$date, "month")  # Keep it as date

# Summarize the counts by Lineage and years_month
plot_counts <- data %>%
  dplyr::group_by(Lineage, years_month) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(years_month, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Total = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("Lineage", "Total", dplyr::everything())) %>%
  arrange(desc(Total))

plot_counts_melted <- plot_counts %>%
  dplyr::select(-c('Total')) %>%
  reshape2::melt(id.vars = 'Lineage')

plot_counts_melted$variable <- as.Date(paste0(plot_counts_melted$variable, "-01"), format = "%Y-%m-%d")

str(plot_counts_melted$variable)

######
p2 <- plot_counts_melted %>%
  ggplot(aes(x = variable, y = value, color = Lineage, group = Lineage)) +
  geom_line() +
  geom_point() +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%b") +
  scale_y_continuous(name = 'Noumber of sequence') +
  scale_color_manual(name = 'Variants', values = my_color) +
  theme_bw() +
  theme(legend.position = 'right',
        legend.text = element_text(size = 8),  # Reduce the legend text size
        legend.title = element_text(size = 8),  # Adjust the legend title size
        legend.key.size = unit(0.2, "cm"), 
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt'))
p2

######

p2_area <- plot_counts_melted %>%
  ggplot(aes(x = variable, y = value, fill = Lineage, group = Lineage)) +
  geom_area(position = "stack") +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%b") +
  scale_y_continuous(name = 'Count') +
  scale_fill_manual(name = 'Variants', values = my_color) +
  theme_bw(base_line_size = 0.3)  +
  theme(legend.position = 'right',
        legend.text = element_text(size = 10),  # Reduce the legend text size
        legend.title = element_text(size = 10),  # Adjust the legend title size
        legend.key.size = unit(0.5, "cm"), 
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.title.x = element_text(name = "Sampling date"),
        plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')) +
  labs(x = "Sampling date")
p2_area

ggsave("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/SARS-CoV-2_Distribution-Guinea/variants-distribution-plot_area.svg", p2_area, width = 8, height = 6)


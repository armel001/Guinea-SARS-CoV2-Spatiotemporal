rm(list = ls())

library("dplyr")
library(readr)
library(readxl)
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(ggrepel)
library(ggalluvial)
library(cowplot)
library(xlsx)
library(tidyr)

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






# Origins of Alpha variant 


data1 <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/events_by_variants/alpha_events.xlsx")
data1$EventTime <- as.Date(data1$EventTime)
data1$EventTime <- format(data1$EventTime, "%Y-%b")

# Aggregate the data: Count occurrences of Origin per EventTime
data1 <- data1 %>% 
  select(Variant, Origin, EventTime) %>% 
  group_by(Variant, Origin, EventTime) %>%
  summarise(count = n(), .groups = "drop")  # Use .groups = "drop" to avoid warning


p1 <- ggplot(data1, aes(x = EventTime, fill = Origin)) + 
  geom_bar(position = position_stack(), width = 0.35, color = "grey2") + 
  scale_fill_manual(name = 'Origin', values = my_color) + 
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
  theme_classic(base_line_size = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 9),  
    legend.title = element_text(size = 9, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  labs(title = "Alpha (B.1.1.7)", y = "Count", x = 'Date')

p1




# Origin of BA.1 variant 

data2 <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/events_by_variants/ba1_events.xlsx")
data2$EventTime <- as.Date(data2$EventTime)
data2$EventTime <- format(data2$EventTime, "%Y-%b")

# Aggregate the data: Count occurrences of Origin per EventTime
data2 <- data2 %>% 
  select(Variant, Origin, EventTime) %>% 
  group_by(Variant, Origin, EventTime) %>%
  summarise(count = n(), .groups = "drop")  # Use .groups = "drop" to avoid warning


p2 <- ggplot(data2, aes(x = EventTime, fill = Origin)) + 
  geom_bar(position = position_stack(), width = 0.35, color = "grey2") + 
  scale_fill_manual(name = 'Origin', values = my_color) + 
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  theme_classic(base_line_size = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 9),  
    legend.title = element_text(size = 9, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  labs(title = "BA.1", y = "Count", x = 'Date')

p2



# Origin of BA.2 variant 

data3 <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/events_by_variants/ba2_events.xlsx")
data3$EventTime <- as.Date(data3$EventTime)
data3$EventTime <- format(data3$EventTime, "%Y-%b")

# Aggregate the data: Count occurrences of Origin per EventTime
data3 <- data3 %>% 
  select(Variant, Origin, EventTime) %>% 
  group_by(Variant, Origin, EventTime) %>%
  summarise(count = n(), .groups = "drop")  # Use .groups = "drop" to avoid warning


p3 <- ggplot(data3, aes(x = EventTime, fill = Origin)) + 
  geom_bar(position = position_stack(), width = 0.35, color = "grey2") + 
  scale_fill_manual(name = 'Origin', values = my_color) + 
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  theme_classic(base_line_size = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 9),  
    legend.title = element_text(size = 9, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  labs(title = "BA.2", y = "Count", x = 'Date')

p3


# Origin of Delta variant 

data4 <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/events_by_variants/delta_events.xlsx")
data4$EventTime <- as.Date(data4$EventTime)
data4$EventTime <- format(data4$EventTime, "%Y-%b")

# Aggregate the data: Count occurrences of Origin per EventTime
data4 <- data4 %>% 
  select(Variant, Origin, EventTime) %>% 
  group_by(Variant, Origin, EventTime) %>%
  summarise(count = n(), .groups = "drop")  # Use .groups = "drop" to avoid warning


p4 <- ggplot(data4, aes(x = EventTime, fill = Origin)) + 
  geom_bar(position = position_stack(), width = 0.35, color = "grey2") + 
  scale_fill_manual(name = 'Origin', values = my_color) + 
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 1)) +
  theme_classic(base_line_size = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 9),  
    legend.title = element_text(size = 9, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  labs(title = "Delta (B.1.617.2)", y = "Count", x = 'Date')

p4


# Origin of BA.5 variant 

data5 <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/events_by_variants/ba5_events.xlsx")
data5$EventTime <- as.Date(data5$EventTime)
data5$EventTime <- format(data5$EventTime, "%Y-%b")

# Aggregate the data: Count occurrences of Origin per EventTime
data5 <- data5 %>% 
  select(Variant, Origin, EventTime) %>% 
  group_by(Variant, Origin, EventTime) %>%
  summarise(count = n(), .groups = "drop")  # Use .groups = "drop" to avoid warning


p5 <- ggplot(data5, aes(x = EventTime, fill = Origin)) + 
  geom_bar(position = position_stack(), width = 0.35, color = "grey2") + 
  scale_fill_manual(name = 'Origin', values = my_color) + 
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
  theme_classic(base_line_size = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 9),  
    legend.title = element_text(size = 9, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  labs(title = "BA.5", y = "Count", x = 'Date')

p5


# Origin of Eta variant 

data6 <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/events_by_variants/eta_events.xlsx")
data6$EventTime <- as.Date(data6$EventTime)
data6$EventTime <- format(data6$EventTime, "%Y-%b")

# Aggregate the data: Count occurrences of Origin per EventTime
data6 <- data6 %>% 
  select(Variant, Origin, EventTime) %>% 
  group_by(Origin, EventTime) %>%
  summarise(count = n())  # Use .groups = "drop" to avoid warning


p6 <- ggplot(data6, aes(x = EventTime, fill = Origin)) + 
  geom_bar(position = position_stack(), width = 0.35, color = "grey2") + 
  scale_fill_manual(name = 'Origin', values = my_color) + 
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1)) +
  theme_classic(base_line_size = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 9),  
    legend.title = element_text(size = 9, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  labs(title = "Eta", y = "Count", x = 'Date')

p6


# Origin of XBB variant 

data7 <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/events_by_variants/xbb_events.xlsx")
data7$EventTime <- as.Date(data7$EventTime)
data7$EventTime <- format(data7$EventTime, "%Y-%b")

# Aggregate the data: Count occurrences of Origin per EventTime
data7 <- data7 %>% 
  select(Variant, Origin, EventTime) %>% 
  group_by(Variant, Origin, EventTime) %>%
  summarise(count = n(), .groups = "drop")  # Use .groups = "drop" to avoid warning


p7 <- ggplot(data7, aes(x = EventTime, fill = Origin)) + 
  geom_bar(position = position_stack(), width = 0.35, color = "grey2") + 
  scale_fill_manual(name = 'Origin', values = my_color) + 
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 1)) +
  theme_classic(base_line_size = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 9),  
    legend.title = element_text(size = 9, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  labs(title = "XBB.1*", y = "Count", x = 'Date')

p7


# Origin of BQ.1*

data8 <- read_excel("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/events_by_variants/bq_events.xlsx")
data8$EventTime <- as.Date(data8$EventTime)
data8$EventTime <- format(data8$EventTime, "%Y-%b")

# Aggregate the data: Count occurrences of Origin per EventTime
data8 <- data8 %>% 
  select(Variant, Origin, EventTime) %>% 
  group_by(Variant, Origin, EventTime) %>%
  summarise(count = n(), .groups = "drop")  # Use .groups = "drop" to avoid warning


p8 <- ggplot(data8, aes(x = EventTime, fill = Origin)) + 
  geom_bar(position = position_stack(), width = 0.35, color = "grey2") + 
  scale_fill_manual(name = 'Origin', values = my_color) + 
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
  theme_classic(base_line_size = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 9),  
    legend.title = element_text(size = 9, face = "bold"),  # Bolder legend title
    legend.key.size = unit(0.4, "cm"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title.x = element_blank(),  
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  labs(title = "BQ.1*", y = "Count", x = 'Date')

p8


# Cowplot

g1 <- plot_grid(p1,  p4, p6, p2, p3, p5,  p8, p7, ncol = 2)
g1


ggsave("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/figures/figure5_intro_by_lineages.png", g1, dpi = 300, width = 10, height = 15)


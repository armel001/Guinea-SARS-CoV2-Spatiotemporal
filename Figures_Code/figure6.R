rm(list = ls())

library(svglite)
library(ggplot2)
library(dplyr)
library(janitor)
library("readxl")
library(cowplot)
library(ggthemes)
library(ggrepel)
library(ggalluvial)
library("lubridate")
library(ggtree)
library(ape)
library(treeio)
library(scales) 
library(ggtreeExtra)
library(ggstar)
library(ggnewscale)
library(readr)
library(deeptime)

my_color <- c('#e31a1c','#fdbf6f',"yellow3", '#cab2d6',"chartreuse1", 
              "#8B5742", '#FAF1E4', '#33a02c','#fb9a99',
              '#6a3d9a','#fdbf6f','#ff7f00','#e31a1c',
              "#00868B",'#b15982','#34ace0',"#9ACD32",
              "chartreuse1", "#9CAFAA", "yellow3", 
              "#8B5742", "#9ACD32",'#FAF1E4', "#CD4F39")


####### Alpha Variant #######

tree1 <- read.nexus("/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/alpha.timetree.nexus")

ggtree(tree1, mrsd = "2021-06-22") + theme_tree2()
meta1 <- as.data.frame(tree1$tip.label)

metadata1 <- read_tsv(file = "/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/alpha.subsampled_metadata.tsv")

phylo_meta1 <- metadata1 %>% dplyr::select(c("strain", "date", "region", "country"))

phylo_meta_df1 <- merge(meta1, phylo_meta1, by.x = "tree1$tip.label", by.y = "strain")
row.names(phylo_meta_df1) <- phylo_meta_df1$`tree1$tip.label`
phylo_meta_df1$date <- as.Date(phylo_meta_df1$date, format = "%Y-%m-%d")
phylo_meta_df1$date <- as.Date(cut(phylo_meta_df1$date, breaks = "1 month", format = "%Y-%m-%d", start.on.monday = FALSE))
phylo_meta_df1$date <- format(as.Date(phylo_meta_df1$date), format = "%b, %Y")

phylo_meta_df1 <- phylo_meta_df1 %>%
  mutate(country = if_else(country == "Guinea", "Guinea", "Others"))

options(ignore.negative.edge=TRUE)
ph1 <-ggtree(tree1,  mrsd = "2021-06-22",
            color='grey44',size=0.5) + theme_tree2() 

ph11 <- ph1 %<+% phylo_meta_df1 + 
  geom_tippoint(aes(color = region, 
                    shape = ifelse(country == "Guinea", "Guinea", "Others"),
                    size = ifelse(country == "Guinea", 3, 2),    # Larger size for Guinea
                    stroke = ifelse(country == "Guinea", 0.2, 0.5)), # Thicker stroke for Guinea
                alpha = 1) +
  scale_colour_manual(name = 'Regions', values = my_color) +
  scale_shape_manual(name = "Country", values = c("Guinea" = 17, "Others" = 16)) +  # Custom shapes
  scale_size_identity() +   # Ensure the size mapping is taken from aes()
  theme(legend.position = "right")  +
  labs(title = "Alpha")

ph11

####### Eta #######

tree2 <- read.nexus("/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/eta.timetree.nexus")

ggtree(tree2, mrsd = "2021-07-07") + theme_tree2()
meta2 <- as.data.frame(tree2$tip.label)

metadata2 <- read_tsv(file = "/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/eta.subsampled_metadata.tsv")

phylo_meta2 <- metadata2 %>% dplyr::select(c("strain","date", "region", "country"))

phylo_meta_df2 <- merge(meta2, phylo_meta2, by.x = "tree2$tip.label", by.y = "strain")
row.names(phylo_meta_df2) <- phylo_meta_df2$`tree2$tip.label`
phylo_meta_df2$date <- as.Date(phylo_meta_df2$date, format = "%Y-%m-%d")
phylo_meta_df2$date <- as.Date(cut(phylo_meta_df2$date, breaks = "1 month", format = "%Y-%m-%d", start.on.monday = FALSE))
#phylo_meta_df2$date <- format(as.Date(phylo_meta_df2$date), format = "%b, %Y")

phylo_meta_df2 <- phylo_meta_df2 %>%
  mutate(country = if_else(country == "Guinea", "Guinea", "Others"))

options(ignore.negative.edge=TRUE)

ph2 <-ggtree(tree2,  mrsd = "2021-07-07",
            color='grey44',size=0.5) + theme_tree2() 

ph22 <- ph2 %<+% phylo_meta_df2 + 
  geom_tippoint(aes(color = region, 
                    shape = ifelse(country == "Guinea", "Guinea", "Others"),
                    size = ifelse(country == "Guinea", 4, 2),    # Larger size for Guinea
                    stroke = ifelse(country == "Guinea", 1, 0.5)), # Thicker stroke for Guinea
                alpha = 1) +
  scale_colour_manual(name = 'Regions', values = my_color) +
  scale_shape_manual(name = "Country", values = c("Guinea" = 17, "Others" = 16)) +  # Custom shapes
  scale_size_identity() +   # Ensure the size mapping is taken from aes()
  theme(legend.position = "right")  +
  labs(title = "Eta")

ph22

table(phylo_meta_df2$country)  
summary(phylo_meta_df2)

####### Delta ########

tree3 <- read.nexus("/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/delta.timetree.nexus")

ggtree(tree3, mrsd = "2022-01-04") + theme_tree2()
meta3 <- as.data.frame(tree3$tip.label)

metadata3 <- read_tsv(file = "/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/delta.subsampled_metadata.tsv")

phylo_meta3 <- metadata3 %>% dplyr::select(c("date", "region", "strain", "country"))

phylo_meta3 <- phylo_meta3 %>%
  mutate(country = if_else(country == "Guinea", "Guinea", "Others"))

phylo_meta_df3 <- merge(meta3, phylo_meta3, by.x = "tree3$tip.label", by.y = "strain")
row.names(phylo_meta_df3) <- phylo_meta_df3$`tree3$tip.label`
phylo_meta_df3$date <- as.Date(phylo_meta_df3$date, format = "%Y-%m-%d")
phylo_meta_df3$date <- as.Date(cut(phylo_meta_df3$date, breaks = "1 month", format = "%Y-%m-%d", start.on.monday = FALSE))
phylo_meta_df3$date <- format(as.Date(phylo_meta_df3$date), format = "%b, %Y")

options(ignore.negative.edge=TRUE)
ph3 <-ggtree(tree3,  mrsd = "2022-01-04",
            color='grey44',size=0.5) + theme_tree2() 

ph33 <- ph3 %<+% phylo_meta_df3 + 
  geom_tippoint(aes(color = region, 
                    shape = ifelse(country == "Guinea", "Guinea", "Others"),
                    size = ifelse(country == "Guinea", 3, 2),    # Larger size for Guinea
                    stroke = ifelse(country == "Guinea", 0.2, 0.5)), # Thicker stroke for Guinea
                alpha = 1) +
  scale_colour_manual(name = 'Regions', values = my_color) +
  scale_shape_manual(name = "Country", values = c("Guinea" = 17, "Others" = 16)) +  # Custom shapes
  scale_size_identity() +   # Ensure the size mapping is taken from aes()
  theme(legend.position = "right") +
  labs(title = "Delta and Sublineages")

ph33

########## BA.1 ################

tree4 <- read.nexus("/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/ba1.timetree.nexus")

ggtree(tree4, mrsd = "2022-05-10") + theme_tree2()
meta4 <- as.data.frame(tree4$tip.label)

metadata4 <- read_tsv(file = "/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/ba1.subsampled_metadata.tsv")

phylo_meta4 <- metadata4 %>% dplyr::select(c("date", "region", "strain", "country"))


phylo_meta4 <- phylo_meta4 %>%
  mutate(country = if_else(country == "Guinea", "Guinea", "Others"))

phylo_meta_df4 <- merge(meta4, phylo_meta4, by.x = "tree4$tip.label", by.y = "strain")
row.names(phylo_meta_df4) <- phylo_meta_df4$`tree4$tip.label`
phylo_meta_df4$date <- as.Date(phylo_meta_df4$date, format = "%Y-%m-%d")
phylo_meta_df4$date <- as.Date(cut(phylo_meta_df4$date, breaks = "1 month", format = "%Y-%m-%d", start.on.monday = FALSE))
phylo_meta_df4$date <- format(as.Date(phylo_meta_df4$date), format = "%b, %Y")

options(ignore.negative.edge=TRUE)
ph4 <-ggtree(tree4,  mrsd = "2022-05-10",
            color='grey44',size=0.5) + theme_tree2() 

ph44 <- ph4 %<+% phylo_meta_df4 + 
  geom_tippoint(aes(color = region, 
                    shape = ifelse(country == "Guinea", "Guinea", "Others"),
                    size = ifelse(country == "Guinea", 3, 2),    # Larger size for Guinea
                    stroke = ifelse(country == "Guinea", 0.2, 0.5)), # Thicker stroke for Guinea
                alpha = 1) +
  scale_colour_manual(name = 'Regions', values = my_color) +
  scale_shape_manual(name = "Country", values = c("Guinea" = 17, "Others" = 16)) +  # Custom shapes
  scale_size_identity() +   # Ensure the size mapping is taken from aes()
  theme(legend.position = "right")  +
  labs(title = "BA.1 and Sublineages")

ph44

######### BA.2 ################

tree5 <- read.nexus("/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/ba2.timetree.nexus")

ggtree(tree5, mrsd = "2023-01-04") + theme_tree2()
meta5 <- as.data.frame(tree5$tip.label)

metadata5 <- read_tsv(file = "/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/ba2.subsampled_metadata.tsv")

phylo_meta5 <- metadata5 %>% dplyr::select(c("date", "region", "strain", "country"))

phylo_meta5 <- phylo_meta5 %>%
  mutate(country = if_else(country == "Guinea", "Guinea", "Others"))

phylo_meta_df5 <- merge(meta5, phylo_meta5, by.x = "tree5$tip.label", by.y = "strain")
row.names(phylo_meta_df5) <- phylo_meta_df5$`tree5$tip.label`
phylo_meta_df5$date <- as.Date(phylo_meta_df5$date, format = "%Y-%m-%d")
phylo_meta_df5$date <- as.Date(cut(phylo_meta_df5$date, breaks = "1 month", format = "%Y-%m-%d", start.on.monday = FALSE))
phylo_meta_df5$date <- format(as.Date(phylo_meta_df5$date), format = "%b, %Y")

options(ignore.negative.edge=TRUE)
ph5 <-ggtree(tree5,  mrsd = "2023-01-04",
            color='grey44',size=0.5) + theme_tree2() 

ph55 <- ph5 %<+% phylo_meta_df5 + 
  geom_tippoint(aes(color = region, 
                    shape = ifelse(country == "Guinea", "Guinea", "Others"),
                    size = ifelse(country == "Guinea", 3, 2),    # Larger size for Guinea
                    stroke = ifelse(country == "Guinea", 0.2, 0.5)), # Thicker stroke for Guinea
                alpha = 1) +
  scale_colour_manual(name = 'Regions', values = my_color) +
  scale_shape_manual(name = "Country", values = c("Guinea" = 17, "Others" = 16)) +  # Custom shapes
  scale_size_identity() +   # Ensure the size mapping is taken from aes()
  theme(legend.position = "right")  +
  labs(title = "BA.2 and Sublineages")

ph55

########## BA.5 ###########

tree6 <- read.nexus("/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/ba5.timetree.nexus")

ggtree(tree6, mrsd = "2023-03-04") + theme_tree2()
meta6 <- as.data.frame(tree6$tip.label)

metadata6 <- read_tsv(file = "/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/ba5.subsampled_metadata.tsv")

phylo_meta6 <- metadata6 %>% dplyr::select(c("date", "region", "strain", "country"))

phylo_meta6 <- phylo_meta6 %>%
  mutate(country = if_else(country == "Guinea", "Guinea", "Others"))

phylo_meta_df6 <- merge(meta6, phylo_meta6, by.x = "tree6$tip.label", by.y = "strain")
row.names(phylo_meta_df6) <- phylo_meta_df6$`tree6$tip.label`
phylo_meta_df6$date <- as.Date(phylo_meta_df6$date, format = "%Y-%m-%d")
phylo_meta_df6$date <- as.Date(cut(phylo_meta_df6$date, breaks = "1 month", format = "%Y-%m-%d", start.on.monday = FALSE))
phylo_meta_df6$date <- format(as.Date(phylo_meta_df6$date), format = "%b, %Y")

options(ignore.negative.edge=TRUE)
ph6 <-ggtree(tree6,  mrsd = "2023-03-04",
            color='grey44',size=0.5) + theme_tree2() 

ph66 <- ph6 %<+% phylo_meta_df6 + 
  geom_tippoint(aes(color = region, 
                    shape = ifelse(country == "Guinea", "Guinea", "Others"),
                    size = ifelse(country == "Guinea", 3, 2),    # Larger size for Guinea
                    stroke = ifelse(country == "Guinea", 0.2, 0.5)), # Thicker stroke for Guinea
                alpha = 1) +
  scale_colour_manual(name = 'Regions', values = my_color) +
  scale_shape_manual(name = "Country", values = c("Guinea" = 17, "Others" = 16)) +  # Custom shapes
  scale_size_identity() +   # Ensure the size mapping is taken from aes()
  theme(legend.position = "right")  +
  labs(title = "BA.5 and Sublineages")

ph66

########### BQ.1 #############

tree7 <- read.nexus("/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/bq.timetree.nexus")

ggtree(tree7, mrsd = "2023-11-01") + theme_tree2()
meta7 <- as.data.frame(tree7$tip.label)

metadata7 <- read_tsv(file = "/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/bq.subsampled_metadata.tsv")

phylo_meta7 <- metadata7 %>% dplyr::select(c("date", "region", "strain", "country"))

phylo_meta7 <- phylo_meta7 %>%
  mutate(country = if_else(country == "Guinea", "Guinea", "Others"))

phylo_meta_df7 <- merge(meta7, phylo_meta7, by.x = "tree7$tip.label", by.y = "strain")
row.names(phylo_meta_df7) <- phylo_meta_df7$`tree7$tip.label`
phylo_meta_df7$date <- as.Date(phylo_meta_df7$date, format = "%Y-%m-%d")
phylo_meta_df7$date <- as.Date(cut(phylo_meta_df7$date, breaks = "1 month", format = "%Y-%m-%d", start.on.monday = FALSE))
phylo_meta_df7$date <- format(as.Date(phylo_meta_df7$date), format = "%b, %Y")

options(ignore.negative.edge=TRUE)
ph7 <-ggtree(tree7,  mrsd = "2023-11-01",
            color='grey44',size=0.5) + theme_tree2() 

ph77 <- ph7 %<+% phylo_meta_df7 + 
  geom_tippoint(aes(color = region, 
                    shape = ifelse(country == "Guinea", "Guinea", "Others"),
                    size = ifelse(country == "Guinea", 3, 2),    # Larger size for Guinea
                    stroke = ifelse(country == "Guinea", 0.2, 0.5)), # Thicker stroke for Guinea
                alpha = 1) +
  scale_colour_manual(name = 'Regions', values = my_color) +
  scale_shape_manual(name = "Country", values = c("Guinea" = 17, "Others" = 16)) +  # Custom shapes
  scale_size_identity() +   # Ensure the size mapping is taken from aes()
  theme(legend.position = "right") +
  labs(title = "BQ.1 and Sublineages")

ph77

######## XBB.1 ############

tree8 <- read.nexus("/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/xbb.timetree.nexus")

ggtree(tree8, mrsd = "2023-12-12") + theme_tree2()
meta8 <- as.data.frame(tree8$tip.label)

metadata8 <- read_tsv(file = "/Users/armel001/OneDrive - LECNAM/SARS-CoV-2_Manuscrit/SARS-CoV-2_Phylogenetics_trees/timetree_without_outliers/xbb.subsampled_metadata.tsv")

phylo_meta8 <- metadata8 %>% dplyr::select(c("date", "region", "strain", "country"))

phylo_meta8 <- phylo_meta8 %>%
  mutate(country = if_else(country == "Guinea", "Guinea", "Others"))

phylo_meta_df8 <- merge(meta8, phylo_meta8, by.x = "tree8$tip.label", by.y = "strain")
row.names(phylo_meta_df8) <- phylo_meta_df8$`tree8$tip.label`
phylo_meta_df8$date <- as.Date(phylo_meta_df8$date, format = "%Y-%m-%d")
phylo_meta_df8$date <- as.Date(cut(phylo_meta_df8$date, breaks = "1 month", format = "%Y-%m-%d", start.on.monday = FALSE))
phylo_meta_df8$date <- format(as.Date(phylo_meta_df8$date), format = "%b, %Y")

options(ignore.negative.edge=TRUE)
ph8 <-ggtree(tree8,  mrsd = "2023-12-12",
            color='grey44',size=0.5) + theme_tree2() 

ph88 <- ph8 %<+% phylo_meta_df8 + 
  geom_tippoint(aes(color = region, 
                    shape = ifelse(country == "Guinea", "Guinea", "Others"),
                    size = ifelse(country == "Guinea", 3, 2),    # Larger size for Guinea
                    stroke = ifelse(country == "Guinea", 0.2, 0.5)), # Thicker stroke for Guinea
                alpha = 1) +
  scale_colour_manual(name = 'Regions', values = my_color) +
  scale_shape_manual(name = "Country", values = c("Guinea" = 17, "Others" = 16)) +  # Custom shapes
  scale_size_identity() +   # Ensure the size mapping is taken from aes()
  theme(legend.position = "right") +
  labs(title = "XBB.1 and sublineages")

ph88


grid_plot_tree_1 <- plot_grid(ph11, ph22, ph33, ph44, ncol = 2)

grid_plot_tree_2 <- plot_grid(ph55, ph66, ph77, ph88, ncol = 2)

merge_plot <- plot_grid(grid_plot_tree_1, grid_plot_tree_2, ncol = 2)

  ggsave("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/Phylogenetic-analysis/figures/SubPhyloTree.png", 
         merge_plot, width = 21, height = 15, dpi = 600)



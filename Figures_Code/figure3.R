rm(list = ls())


# Load librairy

library(readr)
library(dplyr)
library(splitstackshape)
library(gtools)
library(cowplot)
library(ggthemes)
library(ggrepel)
library(ggalluvial)
library("lubridate")
library(ggtree)
library(ape)
library(treeio)
library(scales) # to access breaks/formatting functions
library(ggtreeExtra)
library(ggstar)
library(ggnewscale)
#library(pier)
library(conflicted)

# Import data 

data <- read.csv("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/SARS-CoV2_Mutations_analysis/nextclade.csv", sep=";")
view(data)

##Formats
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

mutation <- data %>% select(c("seqName","Lineage", "aaSubstitutions", "substitutions", "Nextclade_pango"))

mutation1 <- cSplit(
  mutation, splitCols = c("substitutions", "aaSubstitutions"),
  sep = ",", direction = "long")

#generate a column with gene names

mutation1$genes <- mutation1$aaSubstitutions
mutation1$genes <- gsub(":.*" , "" ,mutation1$genes, perl = TRUE)

#generate a column with just amino acid mutations

mutation1$aaMutation <- mutation1$aaSubstitutions

mutation1$aaMutation <- gsub(".*:" , "" , mutation1$aaMutation, perl = TRUE) 

plot_counts0 <- mutation1 %>%  dplyr::group_by(aaSubstitutions, Nextclade_pango) %>% dplyr::summarize(count = n()) %>% 
  tidyr::spread(Nextclade_pango, count, fill = 0) %>% dplyr::ungroup() %>% dplyr::mutate(Total = rowSums(.[c(2:ncol(.))])) %>% 
  dplyr::select(c("aaSubstitutions", "Total", dplyr::everything())) %>% arrange(desc(Total))

plot_counts1 <- mutation1 %>%  dplyr::group_by(aaMutation, Nextclade_pango) %>% dplyr::summarize(count = n()) %>% 
  tidyr::spread(Nextclade_pango, count, fill = 0) %>% dplyr::ungroup() %>% dplyr::mutate(Total = rowSums(.[c(2:ncol(.))])) %>% 
  dplyr::select(c("aaMutation", "Total", dplyr::everything())) %>% arrange(desc(Total))


#Convert all aaMutation below frequency of 50 to other aaMutation

plot_counts1 <- plot_counts1 %>% mutate(aaMutation = case_when(Total <= 50 ~ "Others (n<50)",TRUE ~ as.character(aaMutation)))

plot_counts0 <- plot_counts0 %>% mutate(aaSubstitutions = case_when(Total <= 50 ~ "Others (n<50)",TRUE ~ as.character(aaSubstitutions)))

#Summarize the Variants
plot_counts2 <- plot_counts1 %>% dplyr::filter(aaMutation != "Others (n<50)") %>% select(c("aaMutation", "Total")) %>% dplyr::filter(aaMutation != "NA") %>% arrange(desc(Total))
plot_counts01 <- plot_counts0 %>% dplyr::filter(aaSubstitutions != "Others (n<50)") %>% select(c("aaSubstitutions", "Total")) %>% dplyr::filter(aaSubstitutions != "NA") %>% arrange(desc(Total))
plot_counts01 <- plot_counts01[order(as.numeric(gsub("([A-Z])", "", plot_counts01$aaSubstitutions, perl = T)), na.last=FALSE) , ]

#generate a column with just amino acid mutations
plot_counts01$gene <- plot_counts01$aaSubstitutions
plot_counts01$gene <- gsub(":.*","",plot_counts01$gene, perl = TRUE)

plot_counts01 %>% knitr::kable(caption = "Highest Mutations frequency")

#Visualization

library(forcats)
library(ggplot2)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)


custom_colors <- list()
colors_dutch <- c(
  '#2484c1','#86f71a','#cc9fb1','#D980FA','#F79F1F',
  '#EE5A24','#009432','#833471','#1B1464','#0652DD',
  '#006266','#EA2027','#1B1464','#5758BB','#6F1E51'
)
colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)
colors_net <- c(
  "#00CED1", "#90EE90", "#F08080", "#458B74",
  "#CDB5CD", "#8470FF", "#8B5742", "#9ACD32",
  "#8B2500", "#E6E6FA", "#00868B", "#CD4F39")

custom_colors$discrete <- c(colors_dutch, colors_spanish)
custom_colors$nette <- c(colors_net)
custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)

View(custom_colors)

table(plot_counts01$gene)

p1 <- ggplot(plot_counts01, aes(x=aaSubstitutions, y=Total, fill=gene)) + 
  geom_col() + 
  theme_USGS_box() +
  #facet_wrap(~gene, scales = "free_y", ncol = 4) 
  scale_fill_manual(name = 'gene', values = custom_colors$nette) + coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +
  #geom_segment(aes(x=aaSubstitutions, xend=aaSubstitutions,  y=0, yend=Total)) + 
  #theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(t = 20, r = 1, b = 1, l = 1, unit = 'pt')) +
  labs(x = "Amino acid substitutions")

p1 

#p1 + coord_flip()

ggsave("/Users/armel001/Documents/04_Manuscrits/SARS-CoV-2-Phylodynamie-Guinea/SARS-CoV2_Mutations_analysis/Figure3.tiff", p1, width = 14, height = 6)



##Select spike protein and VOCs

#now you can easily filter specific group them
mutation2 <- mutation1 %>% dplyr::filter(genes == "S") 
head(mutation2)

##Now you can extract it and plot it
plot_counts02 <- mutation2 %>%  dplyr::group_by(aaMutation, Lineage) %>% dplyr::summarize(count = n()) %>% 
  tidyr::spread(Lineage, count, fill = 0) %>% dplyr::ungroup() %>% dplyr::mutate(Total = rowSums(.[c(2:ncol(.))])) %>% dplyr::select(c("aaMutation", "Total", dplyr::everything())) %>% arrange(desc(Total))
View(plot_counts02)

#Convert all aaMutation below frequency of 2 to other aaMutation
plot_counts02 <- plot_counts02 %>% mutate(aaMutation = case_when(Total <= 2 ~ "Others (n<2)",TRUE ~ as.character(aaMutation)))
#Summarize the Variants
plot_counts03 <- plot_counts02 %>% filter(aaMutation != "Others (n<2)")
View(plot_counts03)
plot_counts04 <- as.data.frame(plot_counts03)
View(plot_counts04)

row.names(plot_counts04) <- plot_counts04$aaMutation; plot_counts04$aaMutation <- NULL; plot_counts04$Total <- NULL;
plot_counts04[] <- lapply(plot_counts04, function(x) ifelse(x>1, 1, x))

plot_counts04$mutation <- rownames(plot_counts04)
plot_counts04 <- plot_counts04[order(as.numeric(gsub("([A-Z])", "", plot_counts04$mutation, perl = T)), na.last=FALSE) , ]
rownames(plot_counts04) <- plot_counts04$mutation
plot_counts04 <- plot_counts04 %>% dplyr::select(-c("mutation"))
plot_counts04 %>% knitr::kable(caption = "Spike Mutations")
p2 <- pheatmap(plot_counts04,  color = colorRampPalette((brewer.pal(n = 6, name = "Purples")))(100), cluster_rows = F,
               cluster_cols = F, clustering_distance_rows = "correlation",  treeheight_row = 0, treeheight_col = 0)

p2
p2 + coord_flip()

table(plot_counts04)


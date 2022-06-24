#  work in progress

#!/usr/bin/env Rscript

# call from base directory of a project
# usage:
# Rscript --vanilla --quiet src/R/HCPC.R -m M_hierBAPS_FreqTable.txt -y Y_hierBAPS_FreqTable.txt


# Load libraries
library("optparse")
library("rhierbaps")
suppressPackageStartupMessages(library(tidyverse))

 
option_list <- list(

  make_option(c("-m", "--M_FreqTable"),
              type = "character",
              default = NULL,
              help = "M_FreqTable from hierBAPS",
              metavar = "character"),

  make_option(c("-y", "--Y_FreqTable"),
              type = "character",
              default = 2,
              help = "Y_FreqTable from hierBAPS", # nolint
              metavar = "character")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$fasta)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input fasta).n", call. = FALSE)
}


#  Read in variables
m_freq_table <- as.character(opt$M_FreqTable)
y_freq_table <- as.numeric(opt$Y_FreqTable)

m_freq_table <- read_table("output/hierBAPS/Td_100_SB_50_T_5_run_1214366183_M_L3_FreqTable.txt")
y_freq_table <- read_table("output/hierBAPS/Td_100_SB_50_T_5_run_1214366183_Y_L3_FreqTable.txt")

# create necessary dirs where the output will be stored
dir.create("output/HCPC", recursive = T, showWarnings = F)

# extract output name to be the same as input
out_name <- str_split(fasta_path, pattern = "/", simplify = T) %>%
  .[, length(.)] %>%
  str_split(., pattern = ".fasta", simplify = T) %>%
  .[1]

# extract prefix from the input
prfx <- str_split(fasta_path, pattern = "_", simplify = T) %>%
  .[, length(.)] %>%
  str_split(., pattern = ".fasta", simplify = T) %>%
  .[1]




  # Load libraries ---------------------------------------------------------------

library("tidyverse")
library("pheatmap")
library("factoextra")
library("FactoMineR")

# Set variables ----------------------------------------------------------------

tableM_path <- "output/hierBAPS/Td_100_SB_50_T_5_run_1214366183_M_L3_FreqTable.txt"
tableY_path <- "output/hierBAPS/Td_100_SB_50_T_5_run_1214366183_Y_L3_FreqTable.txt"
out_name <- "Td_100_SB_50_T_5_run_1214366183_MY_L3_FreqTable"

# define colors for ggplot -----------------------------------------------------

cols <- c(
  "Blue" = "#67a9cf",
  "Red" = "#ef8a62"
)

# fun for plotting pheatmap in right dimensions ---------------------------------

get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in")) + 1
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in")) + 1
  return(list(height = plot_height, width = plot_width))
}

# read in R --------------------------------------------------------------------

tableM <- read_delim(tableM_path)
tableY <- read_delim(tableY_path)
tableMY <- inner_join(tableM, tableY) %>% drop_na()

df <- bind_cols(
  select_at(tableMY, "POP"),
  select_if(.tbl = tableMY, function(col) is.numeric(col) && sum(col) != 0)
  ) %>% 
  select(-c(Nind_M, Nind_Y)) %>% 
  column_to_rownames(var = "POP")
  

# Pheatmap_Pops  ---------------------------------------------------------------

p_pheatmap_pop <- pheatmap(df,
           clustering_distance_cols = "correlation",
           clustering_distance_rows = "correlation",
           cluster_rows = T, cluster_cols = T, 
           clustering_method = "ward.D2", 
           cellwidth = 5, cellheight = 5, fontsize = 5, 
           treeheight_row = 8, treeheight_col = 8, 
           )

p_pheatmap_pop
  
plot_dims <- get_plot_dims(p_pheatmap_pop)

ggsave(plot = p_pheatmap_pop,
       filename = paste0("MY_pheatmap_", out_name, ".pdf"), 
       path = "out/plots",
       width = plot_dims$width, height = plot_dims$height, 
       units = "in", device = "pdf")


# Weighted CA -----------------------------------------------------------------

weightsCAplot <- c(rep(1/length(grep("^M_", names(df))), times = length(grep("^M_", names(df)))),
                   rep(1/length(grep("^Y_", names(df))), times = length(grep("^Y_", names(df)))))

colforCAplot_pops <- c(rep("Red", length(grep("^M_", names(df)))),
                       rep("Blue", length(grep("^Y_", names(df)))))

res.ca_weights <- df %>% 
  t() %>%  
  CA(graph=FALSE, 
     row.w = weightsCAplot, 
     ncp = 100)

p_CA <- res.ca_weights %>% 
  fviz_ca_biplot(col.row = colforCAplot_pops,
                 col.col = "black", 
                 repel = T,
                 #labelsize = 2, 
                 #pointsize = 0.5, 
                 title="Weighted CA", legend ="none") + 
  scale_color_manual(values = cols)

p_CA

ggsave(plot = p_CA,
       filename = paste0("CA_biplot_", out_name, ".pdf"),
       path = "out/plots",
       width = 20, height = 20, units = "cm", device = "pdf")



# HCPC -------------------------------------------------------------------------

res.hcpc <- res.ca_weights %>% 
  HCPC(graph = F, 
       nb.clust = 2, 
       method="ward", 
       consol = T, 
       description = T, 
       order=T
       )

p_HCPC <- res.hcpc %>% 
  fviz_cluster(res.hcpc, 
               repel = T, show.clust.cent = F,
               #palette = c("#0075dcff","#f0a3ffff"),
               ggtheme = theme_minimal(),
               #labelsize = 25, pointsize = 5,
               main = "HCPC clusters", legend = "none")

p_HCPC

ggsave(plot = p_HCPC,
       filename = paste0("HCPC_plot_", out_name, ".pdf"),
       path = "out/plots",
       width = 20, height = 20, units = "cm", device = "pdf")




# Sex-based gene flow plot -----------------------------------------------------

HCPC_clusters <- tibble(cluster = res.hcpc$data.clust$clust,
                        haplogroup = rownames(res.hcpc$data.clust))

cl_freqs <- df %>% 
  rownames_to_column(var = "POP") %>% 
  pivot_longer(cols = -1, names_to = "haplogroup", values_to = "freq") %>% 
  separate(haplogroup, into = "marker", sep = 1, remove = F) %>% 
  left_join(HCPC_clusters) %>%
  group_by(cluster, POP, marker) %>% 
  summarise(sum = sum(freq)) 

p_SB <- cl_freqs %>% 
  ggplot(aes(sum, POP, colour=marker)) +
  geom_point(aes(shape=marker, color = marker), size=2) + 
  scale_colour_manual(values = c("#e41a1c","#377eb8")) +
  scale_shape_manual(values = c(2,3)) +
  theme_bw() +
  #theme(text = element_text(size=20)) +
  #xlim(0, 1) +
  labs(x = "Frequency", y = "Population") +
  facet_wrap(vars(cluster), ncol = 2)

p_SB

ggsave(plot = p_SB,
       filename = paste0("SB_plot_", out_name, ".pdf"),
       path = "out/plots",
       width = 30, height = 10, units = "cm", device = "pdf")

#!/usr/bin/env Rscript

# call from base directory of a project
# usage:
# Rscript --vanilla --quiet src/R/hierBAPS_clustering.R -f output/fasta/myfastafile.fasta -n 20 -m 3

# Load libraries
library("optparse")
library("rhierbaps")
suppressPackageStartupMessages(library(tidyverse))

 
option_list <- list(

  make_option(c("-f", "--fasta"),
              type = "character",
              default = NULL,
              help = "fasta file e.g. /data/myfasta.fasta",
              metavar = "character"),

  make_option(c("-m", "--max_depth"),
              type = "numeric",
              default = 2,
              help = "Maximum depth of hierarchical search (default = 2) [same as argument max.depth in hierBAPS]", # nolint
              metavar = "numeric"),

  make_option(c("-n", "--n_pops"),
              type = "numeric",
              default = 20,
              help = "Maximum number of populations in the data (same as argument n.pops in hierBAPS)", # nolint
              metavar = "numeric")
);
 
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$fasta)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input fasta).n", call. = FALSE)
}


#  Read in variables
fasta_path <- as.character(opt$fasta)
max_depth <- as.numeric(opt$max_depth)
n_pops <- as.numeric(opt$n_pops)


# create necessary dirs where the output will be stored
dir.create("output/hierBAPS", recursive = T, showWarnings = F)

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

# extract simulated SB from the input
sim_sb <- str_split(fasta_path, pattern = "SB_", simplify = T)[2] %>%
  str_split(pattern = "_", simplify = T) %>%
  .[1] %>%
  as.numeric()

# load fasta as SNP matrix
snp_matrix <- suppressWarnings(load_fasta(fasta_path))

# clustering with hierBAPS
# to run until the algorithm converges to a local optimum
# add n.extra.rounds = Inf
hb <- hierBAPS(snp_matrix,
               max.depth = max_depth,
               n.pops = n_pops,
               quiet = TRUE
)

# make freq tables for each population and each genetic marker
tmp_df <- hb$partition.df %>%
  as_tibble() %>%
  rename_with(~gsub("level ", "L", .x, fixed = TRUE)) %>%
  separate(col = Isolate,
           into = c("POP", "Ind", "SourceP"),
           sep = "_", fill = "right") %>%
  pivot_longer(cols = starts_with("L"),
               names_to = "level",
               values_to = "cluster",
               values_transform = ~ paste0("C", .x)) %>%
  unite(cluster_name,  c(level, cluster), sep = "_", remove = F) %>%
  mutate_at("cluster_name", ~paste0(prfx, "_", .x)) %>%
  mutate(SourceP = coalesce(SourceP, POP)) %>%
  mutate(SourceP = replace(SourceP, SourceP == "SourceP1", "Pop1")) %>%
  mutate(SourceP = replace(SourceP, SourceP == "SourceP2", "Pop2"))


# save freq table for each clustering level
for (x in unique(tmp_df$level)) {
  tmp_df %>%
    filter(level == x) %>%
    count(POP, cluster_name) %>%
    pivot_wider(names_from = cluster_name,
                values_from = n,
                values_fill = 0) %>%
    rowwise() %>%
    mutate(Nind = sum(c_across(starts_with(prfx)))) %>%
    mutate_at(vars(starts_with(prfx)), ~round((. / Nind), 4)) %>%
    rename_with(~gsub("Nind", paste0("Nind_", prfx), .x, fixed = TRUE)) %>%
    rowwise() %>%
    write_delim(file = paste0("output/hierBAPS/", out_name,
                              "_", x, "_", "FreqTable.txt"),
                delim = "\t", append = FALSE)
}


# save table with individuals their source populations and their
# rhierbaps haplogroup clusters
tmp_df %>%
  write_delim(file = paste0("output/hierBAPS/", out_name,
                            "_", "IndHapClust.txt"),
              delim = "\t", append = FALSE)

# save table with frequencies of true_observed and sim_exp per Pop
tmp_df %>%
  select("POP", "SourceP") %>%
  table() %>%
  as_tibble() %>%
  mutate(sim_exp = case_when(
    POP == "Pop1" & SourceP == "Pop1" ~ sum(n) / 3,
    POP == "Pop2" & SourceP == "Pop2" ~ sum(n) / 3,
    POP == "Pop1" & SourceP == "Pop2" ~ 0,
    POP == "Pop2" & SourceP == "Pop1" ~ 0,
    POP == "Pop3" & SourceP == "Pop1" & prfx == "M" ~ 0.5 * (sum(n) / 3),
    POP == "Pop3" & SourceP == "Pop2" & prfx == "M" ~ 0.5 * (sum(n) / 3),
    POP == "Pop3" & SourceP == "Pop1" & prfx == "Y" ~ (sim_sb / 100) * (sum(n) / 3),
    POP == "Pop3" & SourceP == "Pop2" & prfx == "Y" ~ ((100 - sim_sb) / 100) * (sum(n) / 3))) %>%
  rename(true_observed = n) %>%
    write_delim(file = paste0("output/hierBAPS/", out_name,
                            "_", "ObsExp.txt"),
              delim = "\t", append = FALSE)
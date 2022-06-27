# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("factoextra")
library("FactoMineR")


for (file in list.files("output/hierBAPS", pattern = "_M_L", full.names = T)) {

  f <- file

  # extract base name
  base_name <- str_split(f, pattern = "/", simplify = T) %>%
    .[, length(.)] %>%
    str_split(., pattern = "_[MY]_L[0-9]_FreqTable.txt", simplify = T) %>%
    .[1]


  # Set variables based on M file ----------------------------------------------
  tableY_path <- str_replace(f, "_M_", replacement = "_Y_")
  ObsExpM_path <- str_replace(f, "_M_L[0-9]_FreqTable.txt",
                              replacement = "_M_ObsExp.txt")
  ObsExpY_path <- str_replace(ObsExpM_path, "_M_", replacement = "_Y_")

  # create necessary dirs where the output will be stored
  dir.create("output/SB/plots", recursive = T, showWarnings = F)

  # read in R ------------------------------------------------------------------
  tableM <- read_delim(f)
  tableY <- read_delim(tableY_path)
  tableMY <- inner_join(tableM, tableY) %>% drop_na()

  df <- bind_cols(
    select_at(tableMY, "POP"),
    select_if(.tbl = tableMY,
              function(col) is.numeric(col) && sum(col) != 0)) %>%
    select(-c(Nind_M, Nind_Y)) %>%
    column_to_rownames(var = "POP")


  ObsExpM <- read_delim(ObsExpM_path) %>%
              mutate(marker = "M") %>%
              mutate(true_observed = true_observed / 30, sim_exp = sim_exp / 30)

  ObsExpY <- read_delim(ObsExpY_path) %>%
              mutate(marker = "Y") %>%
              mutate(true_observed = true_observed / 30, sim_exp = sim_exp / 30)

  ObsExpMY <- bind_rows(ObsExpM, ObsExpY) %>%
              mutate(SourceP = as.factor(str_replace(SourceP, "Pop", "")))


  # Weighted CA ---------------------------------------------------------------
  weightsCAplot <- c(rep(1 / length(grep("^M_", names(df))),
                        times = length(grep("^M_", names(df)))),
                    rep(1 / length(grep("^Y_", names(df))),
                        times = length(grep("^Y_", names(df)))))

  colforCAplot_pops <- c(rep("Red", length(grep("^M_", names(df)))),
                        rep("Blue", length(grep("^Y_", names(df)))))

  res.ca_weights <- df %>%
    t() %>%
    CA(graph = FALSE, row.w = weightsCAplot, ncp = 100)
    #CA(graph=FALSE, ncp = 100)

  # HCPC -----------------------------------------------------------------------
  res.hcpc <- res.ca_weights %>%
    HCPC(graph = F, nb.clust = 2, method = "ward",
        consol = T, description = T, order = T)

  # Sex-based gene flow plot ---------------------------------------------------
  HCPC_clusters <- tibble(cluster = res.hcpc$data.clust$clust,
                          haplogroup = rownames(res.hcpc$data.clust))

  cl_freqs <- df %>%
    rownames_to_column(var = "POP") %>%
    pivot_longer(cols = -1, names_to = "haplogroup", values_to = "freq") %>%
    separate(haplogroup, into = "marker", sep = 1, remove = F) %>%
    left_join(HCPC_clusters) %>%
    group_by(cluster, POP, marker) %>%
    summarise(hcpc_clust = sum(freq)) %>%
    mutate(cluster = as.factor(cluster))

  tmp <- cl_freqs %>% mutate(POP = case_when(
    POP == "Pop1" ~ 1,
    POP == "Pop2" ~ 2,
    POP == "Pop3" ~ 3)) %>%
    mutate(POP = as.factor(POP)) %>%
    arrange(POP, desc(hcpc_clust))

  # swap clusters to match populations
  if (tmp$cluster[1] != tmp$POP[1]) {
    cl_freqs <- cl_freqs %>% mutate(cluster = case_when(
      cluster == 1 ~ 2,
      cluster == 2 ~ 1)) %>%
      mutate(cluster = as.factor(cluster))
  }

  cl_freqs_ObsExp <- left_join(cl_freqs, ObsExpMY,
    by = c("cluster" = "SourceP", "POP" = "POP", "marker" = "marker")) %>%
    pivot_longer(cols = c("hcpc_clust", "true_observed", "sim_exp"),
                  values_to = "freq", names_to = "feq_type")

  sb_plot <- cl_freqs_ObsExp %>% 
    ggplot(aes(freq, POP, colour=marker)) +
    geom_point(aes(shape=feq_type, color = marker), size=4) + 
    scale_colour_manual(values = c("#e41a1c","#377eb8")) +
    scale_shape_manual(values = c(1,3,4)) +
    theme_bw() +
    labs(x = "Frequency", y = "Population", title= base_name) +
    facet_wrap(vars(cluster), ncol = 2)
  sb_plot

  ggsave(plot = sb_plot,
        filename = paste0(base_name, "SB_plot.png"), 
        path = "output/SB/plots",
        width = 20, height = 10, 
        units = "cm", device = "png")
  }
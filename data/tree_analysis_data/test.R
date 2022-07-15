library(tidyverse)
library(ggplot2)
setwd("D:/Daten/programming_projects/AncestryProject/data/tree_analysis_data")

treedata_Y <- read.csv("treedata_Y.csv", header = FALSE)
treedata_M <- read.csv("treedata_M.csv", header = FALSE)
treedata <- rbind(treedata_M,treedata_Y)
colnames(treedata) <- strsplit(read_lines("headers.txt"), split = ",")[[1]]


ggplot(treedata, aes(x = as.factor(TD), y = G_NO_INFO)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_FALSE)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_UNKNOWN)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = TRUE_P1-G_P1)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = INIT_P1-G_P1)) +
    geom_boxplot()

ggplot(treedata, aes(x = as.factor(TA), y = G_NO_INFO)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TA), y = G_FALSE)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TA), y = G_UNKNOWN)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TA), y = TRUE_P1-G_P1)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TA), y = INIT_P1-G_P1)) +
    geom_boxplot()

ggplot(treedata, aes(x = as.factor(NS), y = G_NO_INFO)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(NS), y = G_FALSE)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(NS), y = G_UNKNOWN)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(NS), y = TRUE_P1-G_P1)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(NS), y = INIT_P1-G_P1)) +
    geom_boxplot()

ggplot(treedata, aes(x = as.factor(TD), y = G_NO_INFO - DM_NO_INFO)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_FALSE - DM_FALSE)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_UNKNOWN - DM_UNKNOWN)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_P1 - DM_P1)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_P1-DM_P1)) +
    geom_boxplot()

ggplot(treedata, aes(x = as.factor(TD), y = G_NO_INFO - MP_NO_INFO)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_FALSE - MP_FALSE)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_UNKNOWN - MP_UNKNOWN)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_P1 - MP_P1)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = G_P1-MP_P1)) +
    geom_boxplot()

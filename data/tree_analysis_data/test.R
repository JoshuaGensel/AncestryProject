library(tidyverse)
library(ggplot2)

treedata <- read_csv("treedata.csv", col_names = FALSE)
colnames(treedata) <- c("ID","SOURCE", "TD", "TA","SB","N_SAMPLES","NOTINFORMATIVE","P1_PROP","P2_PROP","UNKNOWN","FALSE_INF")


ggplot(treedata, aes(x = as.factor(TD), y = NOTINFORMATIVE)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = UNKNOWN)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = FALSE_INF)) +
    geom_boxplot()
ggplot(treedata, aes(x = as.factor(TD), y = SB/100-P1_PROP)) +
    geom_boxplot()




ggplot(treedata, aes(x = TA, y = NOTINFORMATIVE)) +
    geom_point()
ggplot(treedata, aes(x = TA, y = UNKNOWN)) +
    geom_point()
ggplot(treedata, aes(x = TA, y = FALSE_INF)) +
    geom_point()
ggplot(treedata, aes(x = TA, y = SB/100-P1_PROP)) +
    geom_point()

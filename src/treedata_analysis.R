#-------------------------------------------------------------------------------
#reading in data
library(tidyverse)
treedata_M <- read.csv("../data/tree_analysis_data/treedata_M.csv")
colnames(treedata_M) <- strsplit(read_lines("../data/tree_analysis_data/headers.txt"), split = ",")[[1]]
treedata_Y <- read.csv("../data/tree_analysis_data/treedata_Y.csv")
colnames(treedata_Y) <- strsplit(read_lines("../data/tree_analysis_data/headers.txt"), split = ",")[[1]]

treedata_M <- treedata_M %>% mutate(MARKER = "M")
treedata_Y <- treedata_Y %>% mutate(MARKER = "Y")
treedata <- rbind(treedata_M, treedata_Y)
#treedata <- merge(treedata_M,treedata_Y, by = c("ID", "TD", "TA", "NS"), suffixes = c("_M", "_Y"))

#-------------------------------------------------------------------------------
#Exploratory Data Analysis

library(ggplot2)

#Genealogy
#unknown proportion
ggplot(treedata, aes(x = as.factor(NS), y = G_UNKNOWN/NS)) +
    facet_grid(TA~TD, labeller = label_both)+
    geom_boxplot()+ 
    theme_bw()

#false proportion
ggplot(treedata, aes(x = as.factor(NS), y = G_FALSE/NS)) +
    facet_grid(TA~TD, labeller = label_both)+
    geom_boxplot()+ 
    theme_bw()

#difference truth - inferred
ggplot(treedata, aes(x = as.factor(NS), y = abs(TRUE_P1-G_P1))) +
    facet_grid(TA~TD, labeller = label_both)+
    geom_boxplot()+ 
    theme_bw()

#difference initial - inferred
ggplot(treedata, aes(x = as.factor(NS), y = abs(INIT_P1-G_P1))) +
    facet_grid(TA~TD, labeller = label_both)+
    geom_boxplot()+ 
    theme_bw()


#NJ tree
#unknown proportion
ggplot(treedata, aes(x = as.factor(NS), y = NJ_UNKNOWN/NS, col = MARKER)) +
    facet_grid(TA~TD, labeller = label_both)+
    geom_boxplot() + 
    theme_bw()

#false proportion
ggplot(treedata, aes(x = as.factor(NS), y = NJ_FALSE/NS, col = MARKER)) +
    facet_grid(TA~TD, labeller = label_both)+
    geom_boxplot()+ 
    theme_bw()

#difference truth - inferred
ggplot(treedata, aes(x = as.factor(NS), y = abs(TRUE_P1-NJ_P1), col = MARKER)) +
    facet_grid(TA~TD, labeller = label_both)+
    geom_boxplot()+ 
    theme_bw()

#difference initial - inferred
ggplot(treedata, aes(x = as.factor(NS), y = abs(INIT_P1-NJ_P1), col = MARKER)) +
    facet_grid(TA~TD, labeller = label_both)+
    geom_boxplot()+ 
    theme_bw()

#-------------------------------------------------------------------------------
#estimating accuracy

accuracy_df_G <- data.frame()
for (td in c(300,600,1000,3000,10000)) {
    for (ta in c(2,100,200,300,400,500,1000)) {
        for (ns in c(10,50,100)) {
            subdf <- filter(treedata, TD == td & TA == ta & NS == ns)
            prop2.5 = nrow(filter(subdf, abs(TRUE_P1-G_P1)<0.025))/nrow(subdf)
            prop5 = nrow(filter(subdf, abs(TRUE_P1-G_P1)<0.05))/nrow(subdf)
            accuracy_df_G <- rbind(accuracy_df_G, c(td,ta,ns, prop2.5, prop5))
        }
    }
}
colnames(accuracy_df_G) <- c("TD", "TA", "NS", "PROP2.5", "PROP5")


library(reshape2)
G_table <- round(dcast(accuracy_df_G, TA ~ NS + TD, value.var = "PROP5") * 100, 2) %>% mutate(TA = TA/100)

accuracy_df_NJ <- data.frame()
for (td in c(300,600,1000,3000,10000)) {
    for (ta in c(2,100,200,300,400,500,1000)) {
        for (ns in c(10,50,100)) {
            for (marker in c("Y","M")) {
                subdf <- filter(treedata, TD == td & TA == ta & NS == ns & MARKER == marker)
                prop2.5 = nrow(filter(subdf, abs(TRUE_P1-NJ_P1)<0.025))/nrow(subdf)
                prop5 = nrow(filter(subdf, abs(TRUE_P1-NJ_P1)<0.05))/nrow(subdf)
                accuracy_df_NJ <- rbind(accuracy_df_NJ, c(td,ta,ns, marker, prop2.5, prop5))
            }
        }
    }
}
colnames(accuracy_df_NJ) <- c("TD", "TA", "NS", "MARKER", "PROP2.5", "PROP5")
accuracy_df_NJ[,-4] <- as.numeric(unlist(accuracy_df_NJ[,-4]))
accuracy_df_NJ_M <- filter(accuracy_df_NJ, MARKER == "M")
accuracy_df_NJ_Y <- filter(accuracy_df_NJ, MARKER == "Y")

NJ_table_M <- round(dcast(accuracy_df_NJ_M, TA ~ NS + TD, value.var = "PROP5") * 100, 2) %>% mutate(TA = TA/100)
NJ_table_Y <- round(dcast(accuracy_df_NJ_Y, TA ~ NS + TD, value.var = "PROP5") * 100, 2) %>% mutate(TA = TA/100)

#-------------------------------------------------------------------------------
# estimating error due to tree clustering

y = (treedata$G_UNKNOWN + treedata$G_FALSE)/treedata$NS
n = nrow(treedata)
treedata <- mutate(treedata, G_PROP_MISS = (y*(n-1)+0.5)/n)

library(betareg)
fit1 <- betareg(G_PROP_MISS ~ 1, data = treedata)
fit2 <- betareg(G_PROP_MISS ~ TA, data = treedata)
fit3 <- betareg(G_PROP_MISS ~ TD, data = treedata)
fit4 <- betareg(G_PROP_MISS ~ TA + TD, data = treedata)
fit5 <- betareg(G_PROP_MISS ~ TA * TD, data = treedata)
fit6 <- betareg(G_PROP_MISS ~ TA + TD + NS, data = treedata)
fit7 <- betareg(G_PROP_MISS ~ TA + TD * NS, data = treedata)
fit8 <- betareg(G_PROP_MISS ~ TD + TA * NS, data = treedata)


G_AIC_table <- data.frame(fit = 1:8, AIC = c(AIC(fit1),AIC(fit2),AIC(fit3),AIC(fit4),AIC(fit5),AIC(fit6),AIC(fit7),AIC(fit8)))
G_AIC_table[which.min(G_AIC_table$AIC),]

summary(fit7)
par(mfrow=c(2,2))
plot(fit7)
predict_df <- expand.grid(TD = seq(0,10000,10),
                          TA = seq(0,1000,10),
                          NS = unique(treedata$NS))
predict_df$G_PROP_MISS = predict(fit7, predict_df)
ggplot(treedata, aes(x = TA, y = G_PROP_MISS, col = TD)) +
    facet_grid(~NS) +
    geom_line(data = subset(predict_df, TD == min(predict_df$TD))) +
    geom_line(data = subset(predict_df, TD == quantile(predict_df$TD,0.25))) +
    geom_line(data = subset(predict_df, TD == median(predict_df$TD))) +
    geom_line(data = subset(predict_df, TD == quantile(predict_df$TD,.75))) +
    geom_line(data = subset(predict_df, TD == max(predict_df$TD))) +
    geom_count() + 
    theme_bw()

#-------------------------------------------------------------------------------
# estimating error due to tree clustering without subsampling

treedata_M_nosampling <- read.csv("../data/tree_analysis_data/treedata_M_nosampling.csv")
colnames(treedata_M_nosampling) <- strsplit(read_lines("../data/tree_analysis_data/headers.txt"), split = ",")[[1]]
treedata_Y_nosampling <- read.csv("../data/tree_analysis_data/treedata_Y_nosampling.csv")
colnames(treedata_Y_nosampling) <- strsplit(read_lines("../data/tree_analysis_data/headers.txt"), split = ",")[[1]]


treedata_M_nosampling <- treedata_M_nosampling %>% mutate(MARKER = "M")
treedata_Y_nosampling <- treedata_Y_nosampling %>% mutate(MARKER = "Y")
treedata_nosampling <- rbind(treedata_M_nosampling, treedata_Y_nosampling)

false_runs = subset(treedata_nosampling, G_FALSE > 0)
summary(false_runs$TRUE_P1)
summary(false_runs$TD)

unknown_runs = subset(treedata_nosampling, G_UNKNOWN > 0)
summary(unknown_runs$TRUE_P1)
summary(unknown_runs$TD)

sum(treedata_nosampling$G_FALSE > 0 | treedata_nosampling$G_UNKNOWN > 0, na.rm = T)/nrow(treedata_nosampling)
sum(treedata_nosampling$G_FALSE > 0 | treedata_nosampling$G_UNKNOWN > 0, na.rm = T)/nrow(subset(treedata_nosampling, TD <= 200))


#-------------------------------------------------------------------------------
#error using genetic data to construct the tree
ggplot(treedata, aes(x = "Genealogy", y = G_P1)) +
    geom_boxplot() +
    geom_boxplot(aes(x = "NJ tree", y = NJ_P1)) +
    ylab("P1 Estimate") + 
    theme_bw()

par(mfrow=c(1,3))
hist(treedata$G_P1-treedata$NJ_P1)
hist(treedata$G_FALSE-treedata$NJ_FALSE)
hist(treedata$G_UNKNOWN-treedata$NJ_UNKNOWN)

#-------------------------------------------------------------------------------
#sampling error
#probability of sampling given by hypergeometric distribution
#probability of the proportion of samples being within 5% of the proportion of 
#the population can be calculated as follows:
Ne = 5000
nS = 100
p1 = 0.5

phyper((p1+0.05)*nS,Ne*p1,Ne*(1-p1),nS)-phyper((p1-0.05)*nS,Ne*p1,Ne*(1-p1),nS)

dist_param_df <- expand.grid(p = seq(0, 1, .02),
                             n = c(10,50,100,500,1000),
                             N = 5000)
dist_param_df<- mutate(dist_param_df, p_within5 = phyper((p+0.05)*n,N*p,N*(1-p),n)-phyper((p-0.05)*n,N*p,N*(1-p),n))
g <- ggplot(dist_param_df, aes(x = p, y = p_within5, col = as.factor(n))) + 
    theme_bw()
for (i in dist_param_df$n) {
    g <- g + geom_line(data = subset(dist_param_df, n == i))
}
g


dist_param_df2 <- expand.grid(p = seq(0, .5, .1),
                             n = c(10,50,100,500,1000),
                             N = 5000,
                             off= seq(0,.5,0.01))
dist_param_df2 <- mutate(dist_param_df2, p_off = phyper((p+off)*n,N*p,N*(1-p),n)-phyper((p-off)*n,N*p,N*(1-p),n))
pl <- ggplot(dist_param_df, aes(x = off*100, y = p_off, col = p)) +
    xlab("% off") + ylab("P(within off)") +
    facet_grid(N~n, labeller = label_both) + 
    theme_bw()
    
for (i in unique(dist_param_df2$p)) {
    pl <- pl + geom_line(data = subset(dist_param_df2, p == i))
}
pl

#-------------------------------------------------------------------------------
#drift after admixture

G_drift <- treedata %>%
    group_by(TA) %>%
    summarise(G_drift = abs(INIT_P1 - TRUE_P1))
ggplot(G_drift,aes(x = as.factor(TA),y=G_drift)) +
    geom_boxplot() + 
    theme_bw()

accuracy_df <- data.frame()
for (td in c(300,600,1000,3000,10000)) {
    for (ta in c(2,100,200,300,400,500,1000)) {
        for (ns in c(10,50,100)) {
            for (marker in c("Y","M")) {
                subdf <- filter(treedata, TD == td & TA == ta & NS == ns & MARKER == marker)
                prop5 = nrow(filter(subdf, abs(INIT_P1-NJ_P1)<0.05))/nrow(subdf)
                rmse = sqrt(mean((subdf$INIT_P1-subdf$NJ_P1)^2, na.rm = TRUE))
                accuracy_df <- rbind(accuracy_df, c(td,ta,ns, marker, prop5, rmse))
            }
        }
    }
}
colnames(accuracy_df) <- c("TD", "TA", "NS", "MARKER", "PROP5", "RMSE")
accuracy_df[,-4] <- as.numeric(unlist(accuracy_df[,-4]))

acc_table <- round(dcast(accuracy_df, TA ~ MARKER + NS + TD, value.var = "PROP5") * 100, 2) %>% mutate(TA = TA/100)
rmse_table <- round(dcast(accuracy_df, TA ~ MARKER + NS + TD, value.var = "RMSE"), 2)


################################################################################
#Data handling
library(tidyverse)
library(ggplot2)
require(betareg)
setwd("D:/Daten/programming_projects/AncestryProject/data/tree_analysis_data")

treedata_Y <- read.csv("treedata_Y.csv", header = FALSE)
treedata_M <- read.csv("treedata_M.csv", header = FALSE)
treedata <- rbind(treedata_M, treedata_Y)
colnames(treedata) <- strsplit(read_lines("headers.txt"), split = ",")[[1]]
treedata = treedata %>% 
    mutate(G_DIFF = sqrt((TRUE_P1-G_P1)^2))

################################################################################
#Exploratory Analysis
cor(treedata, use = "pairwise.complete")
summary(treedata$G_DIFF)
ggplot(treedata, aes(x = TD, y = G_NO_INFO, col = TA)) +
    facet_grid(~NS) +
    geom_point()
ggplot(treedata, aes(x = TD, y = G_UNKNOWN/NS, col = TA)) +
    facet_grid(~NS) +
    geom_point()
ggplot(treedata, aes(x = TD, y = G_FALSE/NS, col = TA)) +
    facet_grid(~NS) +
    geom_point()
ggplot(treedata, aes(x = TD, y = G_DIFF, col = TA)) +
    facet_grid(~NS) +
    geom_point()
ggplot(treedata, aes(x = TD, y = sqrt((INIT_P1-G_P1)^2), col = TA)) +
    facet_grid(~NS) +
    geom_point()
plot3d(treedata$TD,treedata$TA,treedata$G_DIFF)

################################################################################
#Beta regression for difference

fit1 <- betareg(G_DIFF ~ 1, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit2 <- betareg(G_DIFF ~ TA, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit3 <- betareg(G_DIFF ~ TD, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit4 <- betareg(G_DIFF ~ TA + TD, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit5 <- betareg(G_DIFF ~ TA * TD, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit6 <- betareg(G_DIFF ~ TA + TD + NS, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit7 <- betareg(G_DIFF ~ TA + TD * NS, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit8 <- betareg(G_DIFF ~ TD + TA * NS, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)

G_AIC_table <- data.frame(fit = 1:8, AIC = c(AIC(fit1),AIC(fit2),AIC(fit3),AIC(fit4),AIC(fit5),AIC(fit6),AIC(fit7),AIC(fit8)))
G_AIC_table

summary(fit6)
par(mfrow = c(2,2))
plot(fit6)


G_pred_Df <- expand.grid(TA = seq(1,1000,2),
                         TD = seq(1,1000,2),
                         NS = unique(treedata$NS))
G_pred_Df$G_DIFF = predict(fit6, G_pred_Df)

ggplot(treedata, aes(x = TA, y = G_DIFF, col = TD)) +
    facet_grid(~NS) +
    geom_line(data = subset(G_pred_Df, TD == 1))+
    geom_line(data = subset(G_pred_Df, TD == 251))+
    geom_line(data = subset(G_pred_Df, TD == 501))+
    geom_line(data = subset(G_pred_Df, TD == 751))+
    geom_line(data = subset(G_pred_Df, TD == 999)) +
    geom_point()
rbPal <- colorRampPalette(c('red','blue'))
plot3d(G_pred_Df$TD,G_pred_Df$TA,G_pred_Df$G_DIFF, col = rbPal(50)[G_pred_Df$NS])


################################################################################
#poisson regression for false inferences
treedata <- treedata %>%
    mutate(G_N_FALSE = as.integer(NS*G_FALSE))
summary(treedata$G_N_FALSE)

fit1 <- glm(G_N_FALSE ~ 1, data = treedata, family = poisson)
fit2 <- glm(G_N_FALSE ~ TA, data = treedata, family = poisson)
fit3 <- glm(G_N_FALSE ~ TD, data = treedata, family = poisson)
fit4 <- glm(G_N_FALSE ~ TA + TD, data = treedata, family = poisson)
fit5 <- glm(G_N_FALSE ~ TA * TD, data = treedata, family = poisson)
fit6 <- glm(G_N_FALSE ~ TA + TD + NS, data = treedata, family = poisson)
fit7 <- glm(G_N_FALSE ~ TA + TD * NS, data = treedata, family = poisson)
fit8 <- glm(G_N_FALSE ~ TD + TA * NS, data = treedata, family = poisson)
fit9 <- glm(G_N_FALSE ~ TD * TA * NS, data = treedata, family = poisson)
G_AIC_table <- data.frame(fit = 1:9, AIC = c(AIC(fit1),AIC(fit2),AIC(fit3),AIC(fit4),AIC(fit5),AIC(fit6),AIC(fit7),AIC(fit8),AIC(fit9)))
G_AIC_table

summary(fit9)
par(mfrow = c(2,2))
plot(fit9)
shapiro.test(fit9$residuals)

G_pred_Df$G_N_FALSE = predict(fit9, G_pred_Df)
G_pred_Df$G_FALSE = G_pred_Df$G_N_FALSE/G_pred_Df$NS
ggplot(treedata, aes(x = TA, y = G_FALSE, col = TD)) +
    facet_grid(~NS) +
    geom_line(data = subset(G_pred_Df, TD == 1))+
    geom_line(data = subset(G_pred_Df, TD == 251))+
    geom_line(data = subset(G_pred_Df, TD == 501))+
    geom_line(data = subset(G_pred_Df, TD == 751))+
    geom_line(data = subset(G_pred_Df, TD == 999)) + 
    geom_point()
ggplot(treedata, aes(x = TD, y = G_FALSE, col = TA)) +
    facet_grid(~NS) +
    geom_line(data = subset(G_pred_Df, TA == 1))+
    geom_line(data = subset(G_pred_Df, TA == 251))+
    geom_line(data = subset(G_pred_Df, TA == 501))+
    geom_line(data = subset(G_pred_Df, TA == 751))+
    geom_line(data = subset(G_pred_Df, TA == 999)) +
    geom_point()
plot3d(G_pred_Df$TD,G_pred_Df$TA,G_pred_Df$G_FALSE, col = rbPal(50)[G_pred_Df$NS])

################################################################################
#poisson regression for unknown nodes
treedata <- treedata %>%
    mutate(G_N_UNKNOWN = as.integer(NS*G_UNKNOWN))
summary(treedata$G_N_UNKNOWN)

fit1 <- glm(G_N_UNKNOWN ~ 1, data = treedata, family = poisson)
fit2 <- glm(G_N_UNKNOWN ~ TA, data = treedata, family = poisson)
fit3 <- glm(G_N_UNKNOWN ~ TD, data = treedata, family = poisson)
fit4 <- glm(G_N_UNKNOWN ~ TA + TD, data = treedata, family = poisson)
fit5 <- glm(G_N_UNKNOWN ~ TA * TD, data = treedata, family = poisson)
fit6 <- glm(G_N_UNKNOWN ~ TA + TD + NS, data = treedata, family = poisson)
fit7 <- glm(G_N_UNKNOWN ~ TA + TD * NS, data = treedata, family = poisson)
fit8 <- glm(G_N_UNKNOWN ~ TD + TA * NS, data = treedata, family = poisson)
fit9 <- glm(G_N_UNKNOWN ~ TD * TA * NS, data = treedata, family = poisson)
G_AIC_table <- data.frame(fit = 1:9, AIC = c(AIC(fit1),AIC(fit2),AIC(fit3),AIC(fit4),AIC(fit5),AIC(fit6),AIC(fit7),AIC(fit8),AIC(fit9)))
G_AIC_table

summary(fit9)
par(mfrow = c(2,2))
plot(fit9)
shapiro.test(fit9$residuals)

G_pred_Df$G_N_UNKNOWN = predict(fit9, G_pred_Df)
G_pred_Df$G_UNKNOWN = G_pred_Df$G_N_UNKNOWN/G_pred_Df$NS
ggplot(treedata, aes(x = TA, y = G_UNKNOWN, col = TD)) +
    facet_grid(~NS) +
    geom_line(data = subset(G_pred_Df, TD == 1))+
    geom_line(data = subset(G_pred_Df, TD == 251))+
    geom_line(data = subset(G_pred_Df, TD == 501))+
    geom_line(data = subset(G_pred_Df, TD == 751))+
    geom_line(data = subset(G_pred_Df, TD == 999)) + 
    geom_point()
ggplot(treedata, aes(x = TD, y = G_UNKNOWN, col = TA)) +
    facet_grid(~NS) +
    geom_line(data = subset(G_pred_Df, TA == 1))+
    geom_line(data = subset(G_pred_Df, TA == 251))+
    geom_line(data = subset(G_pred_Df, TA == 501))+
    geom_line(data = subset(G_pred_Df, TA == 751))+
    geom_line(data = subset(G_pred_Df, TA == 999)) +
    geom_point()
plot3d(G_pred_Df$TD,G_pred_Df$TA,G_pred_Df$G_UNKNOWN, col = rbPal(50)[G_pred_Df$NS])

summary(G_pred_Df$G_N_FALSE)
summary(G_pred_Df$G_N_UNKNOWN)
summary(treedata$G_N_FALSE)
summary(treedata$G_N_UNKNOWN)



fit1 <- betareg(G_DIFF ~ TA + TD + G_P1, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit2 <- betareg(G_DIFF ~ TA + TD + NS + G_P1, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit3 <- betareg(G_DIFF ~ TA * TD * NS + G_P1, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)
fit4 <- betareg(G_DIFF ~ TA * TD * NS * G_P1, data = subset(treedata,G_DIFF>0 & G_DIFF<1), weights = NS-NS*G_UNKNOWN)

summary(fit4)

G_AIC_table <- data.frame(fit = 1:4, AIC = c(AIC(fit1),AIC(fit2),AIC(fit3),AIC(fit4)))
G_AIC_table

G_pred_Df <- expand.grid(TA = seq(1,1000,2),
                         TD = seq(1,1000,2),
                         G_P1 = seq(0,1,.1),
                         NS = unique(treedata$NS))
G_pred_Df$G_DIFF = predict(fit4, G_pred_Df)
ggplot(subset(treedata, TD == min(treedata$TD) | TD == max(treedata$TD) ), aes(x = G_P1, y = G_DIFF, col = TA)) +
    facet_grid(TD~NS) +
    geom_line(data = subset(G_pred_Df, TA == 1))+
    geom_line(data = subset(G_pred_Df, TA == 251))+
    geom_line(data = subset(G_pred_Df, TA == 501))+
    geom_line(data = subset(G_pred_Df, TA == 751))+
    geom_line(data = subset(G_pred_Df, TA == 999)) +
    geom_point()
summary(fit4)
ggplot(treedata, aes(x = TA, y = G_DIFF, col = TD)) +
    facet_grid(~NS) +
    geom_line(data = subset(G_pred_Df, TD == 1))+
    geom_line(data = subset(G_pred_Df, TD == 251))+
    geom_line(data = subset(G_pred_Df, TD == 501))+
    geom_line(data = subset(G_pred_Df, TD == 751))+
    geom_line(data = subset(G_pred_Df, TD == 999)) +
    geom_point()
################################################################################
#DM



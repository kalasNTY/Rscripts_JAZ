library(ggplot2)
library(dplyr)
library(titanic)
library(cowplot)
library(mice)
library(missForest)

#test imputation on titanice database from https://appsilon.com/imputation-in-r/

head(titanic_train$Age)

value_imputed <- data.frame(
  original = titanic_train$Age,
  imputed_zero = replace(titanic_train$Age, is.na(titanic_train$Age), 0),
  imputed_mean = replace(titanic_train$Age, is.na(titanic_train$Age), mean(titanic_train$Age, na.rm = TRUE)),
  imputed_median = replace(titanic_train$Age, is.na(titanic_train$Age), median(titanic_train$Age, na.rm = TRUE))
)
value_imputed

h1 <- ggplot(value_imputed, aes(x = original)) +
  geom_histogram(fill = "#ad1538", color = "#000000", position = "identity") +
  ggtitle("Original distribution") +
  theme_classic()
h2 <- ggplot(value_imputed, aes(x = imputed_zero)) +
  geom_histogram(fill = "#15ad4f", color = "#000000", position = "identity") +
  ggtitle("Zero-imputed distribution") +
  theme_classic()
h3 <- ggplot(value_imputed, aes(x = imputed_mean)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("Mean-imputed distribution") +
  theme_classic()
h4 <- ggplot(value_imputed, aes(x = imputed_median)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("Median-imputed distribution") +
  theme_classic()

plot_grid(h1, h2, h3, h4, nrow = 2, ncol = 2)


titanic_numeric <- titanic_train %>%
  select(Survived, Pclass, SibSp, Parch, Age)

md.pattern(titanic_numeric)

mice_imputed <- data.frame(
  original = titanic_train$Age,
  imputed_pmm = complete(mice(titanic_numeric, method = "pmm"))$Age,
  imputed_cart = complete(mice(titanic_numeric, method = "cart"))$Age,
  imputed_lasso = complete(mice(titanic_numeric, method = "lasso.norm"))$Age
)

h1 <- ggplot(mice_imputed, aes(x = original)) +
  geom_histogram(fill = "#ad1538", color = "#000000", position = "identity") +
  ggtitle("Original distribution") +
  theme_classic()
h2 <- ggplot(mice_imputed, aes(x = imputed_pmm)) +
  geom_histogram(fill = "#15ad4f", color = "#000000", position = "identity") +
  ggtitle("pmm-imputed distribution") +
  theme_classic()
h3 <- ggplot(mice_imputed, aes(x = imputed_cart)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("cart-imputed distribution") +
  theme_classic()
h4 <- ggplot(mice_imputed, aes(x = imputed_lasso)) +
  geom_histogram(fill = "#ad8415", color = "#000000", position = "identity") +
  ggtitle("lasso-imputed distribution") +
  theme_classic()

plot_grid(h1, h2, h3, h4, nrow = 2, ncol = 2)



missForest_imputed <- data.frame(
  original = titanic_numeric$Age,
  imputed_missForest = missForest(titanic_numeric)$ximp$Age,
  imputed_rf = complete(mice(titanic_numeric, method = "rf"))$Age
)

h1 <- ggplot(missForest_imputed, aes(x = original)) +
  geom_histogram(fill = "#ad1538", color = "#000000", position = "identity") +
  ggtitle("Original distribution") +
  theme_classic()
h2 <- ggplot(missForest_imputed, aes(x = imputed_missForest)) +
  geom_histogram(fill = "#15ad4f", color = "#000000", position = "identity") +
  ggtitle("missForest-imputed distribution") +
  theme_classic()
h3 <- ggplot(missForest_imputed, aes(x = imputed_rf)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("rf-imputed distribution") +
  theme_classic()
plot_grid(h1, h2, h3, nrow = 1, ncol = 3)

setwd('~/Documents/ARGH/')

dfSPACE <- read.csv("PR736_Karen_20230925_abundances_for_imputation.csv",header = T,sep = ",",na.strings = "")

colnames(dfSPACE) <- c("Accescsion","Description","A_dtag1","A_dtag2","A_dtag3","A_dmso1","A_dmso2","A_dmso3")

dfSPACE_imputed <- data.frame(
  original = dfSPACE[,3:8],
  imputed_missForest = missForest(dfSPACE[,3:8])$ximp,
  imputed_rf = complete(mice(dfSPACE[,3:8], method = "rf"))
)
colnames(dfSPACE_imputed)

dfSPACE_imputed_comb <- data.frame(
  original <- c(dfSPACE_imputed[,1],dfSPACE_imputed[,2],dfSPACE_imputed[,3]),
  imputed_missForest <- c(dfSPACE_imputed[,4],dfSPACE_imputed[,5],dfSPACE_imputed[,6]),
  imputed_rf <- c(dfSPACE_imputed[,7],dfSPACE_imputed[,8],dfSPACE_imputed[,9])
)

h1 <- ggplot(dfSPACE_imputed_comb, aes(x = original)) +
  geom_histogram(fill = "#ad1538", color = "#000000", position = "identity") +
  ggtitle("Original distribution") +
  theme_classic()
h2 <- ggplot(dfSPACE_imputed_comb, aes(x = imputed_missForest)) +
  geom_histogram(fill = "#15ad4f", color = "#000000", position = "identity") +
  ggtitle("missForest-imputed distribution") +
  theme_classic()
h3 <- ggplot(dfSPACE_imputed_comb, aes(x = imputed_rf)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("rf-imputed distribution") +
  theme_classic()
plot_grid(h1, h2, h3, nrow = 1, ncol = 3)

dfSPACE_imputed_mF <- data.frame(
  dfSPACE[,1:2],missForest(dfSPACE[,3:8])$ximp
)
head(dfSPACE_imputed_mF)

dfSPACE_imputed_rf <- data.frame(
  dfSPACE[,1:2],
  complete(mice(dfSPACE[,3:8], method = "rf"))
)
head(dfSPACE_imputed_rf)


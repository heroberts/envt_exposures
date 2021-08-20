## ANALYSIS ----

library(data.table); library(tidyverse)
library(readstata13); library(dplyr)
library(arsenal); library(tableone) 
library(sf); library(car)
library(QuantPsyc); library(fastDummies)
library(psych); library(arm)
library(lmtest); library(Hmisc)
library(corrplot); library(xtable)
library(lubridate); library(e1071)
library(survey)

path_in <- "H:/depression/data"
path_out <- "H:/depression/results"

## Read in data ----
points <- read.csv(file.path(path_in, "final/gps_exposure.csv"), stringsAsFactors = T)
survey.data <- read.csv(file.path(path_in, "final/survey_data.csv"), stringsAsFactors = T)

# formatting
points$timestamp <- as.POSIXct(strptime(points$timestamp, "%Y-%m-%d %H:%M:%S"), tz= "UTC")
points$timestamp <- with_tz(points$timestamp, tzone="Europe/Amsterdam")
survey.data$sex <- factor(survey.data$sex, levels=c("Male", "Female"))
survey.data$income <- factor(survey.data$income, levels=c("Very_low", "Low", "Middle", "High", "Very_high"))
survey.data$educ <- factor(survey.data$educ, levels=c("Low", "Mid", "High"))
survey.data$household <- factor(survey.data$household, levels=c("Couple_with_child", "Couple_without_child","Single_parent", "Other_HHtype"))

## Descriptive statistics (sample) ----

# PHQ9 cronbach's alpha
phq9.items <- survey.data %>%
  dplyr::select(PHQ9_1:PHQ9_9)

psych::alpha(phq9.items)

rm(phq9.items)

# Wilcoxon test of full survey sample PHQ9 score and included sample PHQ9 score
full.sample <- read.dta13("H:/depression/data/temp/survey_clean_linked.dta")

wilcox.test(full.sample$PHQ9_score, survey.data$PHQ9_score, na.rm=T) 

rm(full.sample)

# descriptive table

mycontrols <- tableby.control(
  test = F, total = F,
  numeric.stats = c("meansd", "range"),
  cat.stats = c("countpct"),
  stats.labels = list(meansd = "Mean (SD)", range="Range"),
  control=tableby.control())

surv.desc <- survey.data %>% 
  dplyr::select(c(PHQ9_score, age, sex, empl, educ, marital, household, origin, income,
                  POPDENS16_res50, POPDENS16_res100, FRAG16_res50, FRAG16_res100, DEPRI16_res50, DEPRI16_res100)) 

descs <- tableby( ~ ., data = surv.desc, control = mycontrols, digits.p=3)
summary(descs, title = "", text=T)  
write2word(descs, file.path(path_out, "tab/sample_descriptives.doc")) 

rm(surv.desc)

## Descriptive statistics (environmental exposures) -----

# select residential exposures
res.exp.dat <- survey.data %>% 
  dplyr::select('WE_ID_crypt', 'NDVI2018_res50','LGN18bl_res50','NOIw_res50','p25_res50',
                'NDVI2018_res100','LGN18bl_res100','NOIw_res100','p25_res100')

# change wide to long 
res.exp.dat <- reshape(res.exp.dat, direction='long', 
                       varying=matrix(2:ncol(res.exp.dat), 4), 
                       v.names= c('NDVI2018','LGN18bl','NOIw','p25'), timevar='exposure', idvar='WE_ID_crypt', 
                       times= c('res50','res100'), sep='_') %>%
  dplyr::select(-c(WE_ID_crypt))

# write out to word
descs <- tableby(exposure~., data = res.exp.dat, control = mycontrols, digits.p=3)
summary(descs, title = "", text=T)  
write2word(descs, file.path(path_out, "tab/res_exposure_descriptives.doc"))

rm(res.exp.dat)

# gps exposures descs
descs <- tableby(~ NDVI2018_gps50 + NDVI2018_gps100 + LGN18bl_gps50 + LGN18bl_gps100 + NOIw_gps50 + NOIw_gps100 + p25_gps50 + p25_gps100, data = points, control = mycontrols, digits.p=3)
summary(descs, title = "", text=T)  
write2word(descs, file.path(path_out, "tab/gps_exposure_descriptives.doc"))

rm(descs, mycontrols)

## Calculating unbiased group means (using Croon and van Veldhoven (2007) code) ----

survey.data <- survey.data %>%
  arrange(WE_ID_crypt)

survey.data.dummies <- fastDummies::dummy_cols(survey.data, select_columns = c("sex", "empl", "educ", "marital", "household", "origin", "income"), remove_first_dummy = FALSE)

# code taken directly from paper C+V paper from here

data1 <- na.omit(dplyr::select(survey.data.dummies, age, sex_Female, educ_Mid, educ_High, empl_Unemployed, empl_Non_working, empl_Other, marital_Separated, marital_Widowed,
                               marital_Unmarried, household_Couple_without_child, household_Couple_with_child, household_Single_parent,
                               origin_Western, origin_Non_western, income_Low, income_Middle, income_High, income_Very_high, POPDENS16_res50,
                               FRAG16_res50, DEPRI16_res50, POPDENS16_res100, FRAG16_res100, DEPRI16_res100, PHQ9_score))

data2 <- na.omit(dplyr::select(points, WE_ID_crypt, NDVI2018_gps50, LGN18bl_gps50, NOIw_gps50, p25_gps50, NDVI2018_gps100, LGN18bl_gps100, NOIw_gps100, p25_gps100))

ng <- dim(data1)[1]
nt <- dim(data2)[1]
g <- data2[,1]
ns <- tapply(rep(1,nt),g,sum)
mgroup <- dim(data1)[2]-1
mind <- dim(data2)[2]-1

x <- as.matrix(data2[,1+(1:mind)]) # exposures
z <- as.matrix(data1[,1:mgroup]) # phq
y <- data1[,1+mgroup] 
mux <- apply(x,2,mean)
muz <- apply(z,2,mean)
dz <- z - matrix(rep(muz,ng),ncol=mgroup,byrow=T)
xmean <- matrix(0,ng,mind)

for (i in 1:mind){
  xmean[,i] <- tapply(x[,i],g,mean)
}

vv <- var(cbind(z,xmean))
ind1 <- 1:mgroup
ind2 <- mgroup + (1:mind)
vzz <- vv[ind1,ind1,drop=F]
vzxi <- vv[ind1,ind2,drop=F]
mm <- matrix(rep(c(xmean), rep(ns,mind)), ncol=mind)
d <- x - mm
mse <- t(na.omit(d)) %*% na.omit(d)/(nt-ng)
vu <- mse
d <- mm - matrix(rep(mux,nt), ncol=mind, byrow=T)
msa <- t(d) %*% d/(ng-1)
cc <- (nt - sum(ns^2)/nt) /(ng-1)
vxi <- (msa-mse)/cc
xtilde <- matrix(0, ng, mind)
r2 <- solve(vzz,vzxi)
r1 <- vxi - t(vzxi) %*% r2
id <- diag(mind)

for(i in 1:ng){
  p <- solve(r1 + vu/ns[i], r1)
  q <- r2 %*% (id-p)
  xtilde[i,] <- xmean[i,] %*% p + mux %*% (id-p) + dz[i,] %*% q
  
} 

daf <- data.frame(z,xmean,xtilde,y) %>%
  rename(NDVI2018_gps50 = X1, LGN18bl_gps50 = X2, NOIw_gps50 = X3, p25_gps50 = X4, 
         NDVI2018_gps100 = X5, LGN18bl_gps100 = X6, NOIw_gps100 = X7, p25_gps100 = X8, 
         BLUP.NDVI2018_gps50 = X1.1, BLUP.LGN18bl_gps50 = X2.1, BLUP.NOIw_gps50 = X3.1, BLUP.p25_gps50 = X4.1,
         BLUP.NDVI2018_gps100 = X5.1, BLUP.LGN18bl_gps100 = X6.1, BLUP.NOIw_gps100 = X7.1, BLUP.p25_gps100 = X8.1,
         PHQ9_score = y)

# join daf with survey data to make combined dataset
combined.data <- left_join(survey.data, daf)

rm(d, dz, id, mm, msa, mse, p, q, r1, r2, vu, vv, vxi, vzxi, vzz, x, xmean, xtilde, z, cc, g, i, ind1, ind2, mgroup, mind, mux, muz, ng, ns, nt, y)

## Function for corrected SEs (C+V code) ----

corrected.se <- function(model, model.formula){
  
  tm <- unlist(model.formula[[3]])
  var <- unlist(strsplit(as.character(tm), " "))
  lab <- var[!var=="+"]
  u <- cbind(rep(1,393), as.matrix(combined.data[lab]))
  e <- model$residuals
  p <- solve(t(u) %*% u)
  h <- diag(u %*% p %*% t(u))
  d <- e^2/(1-h)
  v <- p %*% t(u) %*% diag(d) %*% u %*% p
  se <- sqrt(diag(v))
  
  estim <- summary(model)$coefficients[,1:2]
  estim <- cbind(estim, se)
  dimnames(estim)[[2]] <- c("b", "Uncorrected SE", "Corrected SE")
  
  estim
  
}

# C+V code ends here

# Testing differences between residential and mobility exposures ----
options(scipen = 999)

wilcox.test(survey.data$NDVI2018_res50, combined.data$NDVI2018_gps50, na.rm=T, paired=T)
wilcox.test(survey.data$LGN18bl_res50, combined.data$LGN18bl_gps50, na.rm=T, paired=T)
wilcox.test(survey.data$NOIw_res50, combined.data$NOIw_gps50, na.rm=T, paired=T)
wilcox.test(survey.data$p25_res50, combined.data$p25_gps50, na.rm=T, paired=T)

wilcox.test(survey.data$NDVI2018_res100, combined.data$NDVI2018_gps100, na.rm=T, paired=T)
wilcox.test(survey.data$LGN18bl_res100, combined.data$LGN18bl_gps100, na.rm=T, paired=T)
wilcox.test(survey.data$NOIw_res100, combined.data$NOIw_gps100, na.rm=T, paired=T)
wilcox.test(survey.data$p25_res100, combined.data$p25_gps100, na.rm=T, paired=T)


## Plot correlations -------

cormat <- as.matrix(Filter(is.numeric, combined.data %>% 
                             dplyr::select(.,c( "NDVI2018_res50", "NDVI2018_res100", "NDVI2018_gps50",  "NDVI2018_gps100",
                                                "LGN18bl_res50",  "LGN18bl_res100",  "LGN18bl_gps50", "LGN18bl_gps100",
                                                "NOIw_res50",  "NOIw_res100", "NOIw_gps50",  "NOIw_gps100",
                                                "p25_res50", "p25_res100", "p25_gps50", "p25_gps100", "PHQ9_score"))))


colnames(cormat) <- c("NDVI2018_res50", "NDVI2018_res100", "NDVI2018_gps50",  "NDVI2018_gps100",
                      "BLUE_res50",  "BLUE_res100",  "BLUE_gps50", "BLUE_gps100",
                      "NOISE_res50",  "NOISE_res100", "NOISE_gps50",  "NOISE_gps100",
                      "PM2.5_res50", "PM2.5_res100", "PM2.5_gps50", "PM2.5_gps100", "PHQ9")

corrcol2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061"))

corspear <- rcorr(cormat, type="spearman")

corr.coef <- as.data.frame(corspear$r)
corr.coef.table <- xtable(corr.coef, digits=3)

corr.sig <- as.data.frame(corspear$P)
corr.sig.table <- xtable(corr.sig, digits=3)

pdf(file = file.path(path_out, "fig/corr_num.pdf"))

corrplot(corspear$r, method="color", col=rev(corrcol2(200)), type="upper", order="original", addCoef.col = "black",
         tl.col = "black", tl.srt=45, p.mat = corspear$P, sig.level= .05, insig = "blank", diag = F, number.cex=0.6, tl.cex=0.6)

dev.off()

## Residential-based multiple regressions -----

# Function for printing results needed
res_results_func <- function(x){
  print("Unstandardized Results")
  print(summary(x))
  print("Standardized Results")
  print(arm::standardize(x))
  print("AIC Values")
  print(AIC(x))
  print("GVIF Values")
  print(vif(x))
  
}

# 50m buffer
sink(file= file.path(path_out, "tab/multiple_res_50m.csv"), append = F)

model.unadj <- lm(PHQ9_score ~ NDVI2018_res50 + LGN18bl_res50 + NOIw_res50 + p25_res50, data=survey.data)
model.indiv <- update(model.unadj, .~. + age + sex + empl + educ + marital + household + origin + income)
model.area50 <- update(model.indiv, .~. + POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50)
model.area100 <- update(model.indiv, .~. + POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100)

model_list <- list(model.unadj, model.indiv, model.area50, model.area100)

invisible(lapply(model_list, res_results_func))

sink()

# Moderation by sex

sink(file= file.path(path_out, "tab/multiple_res_50m_sex.csv"), append = F)

model.int <- lm(PHQ9_score ~ NDVI2018_res50 + LGN18bl_res50 + NOIw_res50 + p25_res50 + 
                  age + sex + empl + educ + marital + household + origin + income + 
                  POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100 + NDVI2018_res50*sex, data=survey.data)

summary(model.int)
lrtest(model.int, model.area100)

model.int <- lm(PHQ9_score ~ NDVI2018_res50 + LGN18bl_res50 + NOIw_res50 + p25_res50 + 
                  age + sex + empl + educ + marital + household + origin + income + 
                  POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100 + LGN18bl_res50*sex, data=survey.data)

summary(model.int)
lrtest(model.int, model.area100)

model.int <- lm(PHQ9_score ~ NDVI2018_res50 + LGN18bl_res50 + NOIw_res50 + p25_res50 + 
                  age + sex + empl + educ + marital + household + origin + income + 
                  POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100 + NOIw_res50*sex, data=survey.data)

summary(model.int)
lrtest(model.int, model.area100)

model.int <- lm(PHQ9_score ~ NDVI2018_res50 + LGN18bl_res50 + NOIw_res50 + p25_res50 + 
                  age + sex + empl + educ + marital + household + origin + income + 
                  POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100 + p25_res50*sex, data=survey.data)

summary(model.int)
lrtest(model.int, model.area100)

sink()


# 100m buffer

sink(file= file.path(path_out, "tab/multiple_res_100m.csv"), append = F)

model.unadj <- lm(PHQ9_score ~ NDVI2018_res100 + LGN18bl_res100 + NOIw_res100 + p25_res100, data=survey.data)
model.indiv <- update(model.unadj, .~. + age + sex + empl + educ + marital + household + origin + income)
model.area50 <- update(model.indiv, .~. + POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50)
model.area100 <- update(model.indiv, .~. + POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100)

model_list <- list(model.unadj, model.indiv, model.area50, model.area100)

invisible(lapply(model_list, res_results_func))

sink()

# Moderation by sex
sink(file= file.path(path_out, "tab/multiple_res_100m_sex.csv"), append = F)

model.int <- update(model.area100, .~. + NDVI2018_res100*sex, data=survey.data)
summary(model.int)
lrtest(model.int, model.area100)

model.int <- update(model.area100, .~. + LGN18bl_res100*sex, data=survey.data)
summary(model.int)
lrtest(model.int, model.area100)

model.int <- update(model.area100, .~. + NOIw_res100*sex, data=survey.data)
summary(model.int)
lrtest(model.int, model.area100) 

model.int <- update(model.area100, .~. + p25_res100*sex, data=survey.data)
summary(model.int)
lrtest(model.int, model.area100)

sink()


# Mobility-based multiple regressions (unbiased means) ---------

# Function for printing results needed

gps_results_func <- function(x, y){
  print("Unstandardized Results")
  print(summary(x))
  print("Standardized Results")
  print(arm::standardize(x))
  print("Corrected SEs")
  print(corrected.se(x, y))
  print("AIC Value")
  print(AIC(x))
  print("GVIF Values")
  print(vif(x))
  
}

# GPS 50m

sink(file= file.path(path_out, "tab/multiple_gps_50m.csv"), append = F)

# 50m unadjusted
model.unadj <- as.formula(PHQ9_score ~ BLUP.NDVI2018_gps50 + BLUP.LGN18bl_gps50 + BLUP.NOIw_gps50 + BLUP.p25_gps50)

m1 <- lm(PHQ9_score ~ BLUP.NDVI2018_gps50 + BLUP.LGN18bl_gps50 + BLUP.NOIw_gps50 + BLUP.p25_gps50, combined.data)


# 50m + individual variables
model.indiv <- as.formula(PHQ9_score ~ BLUP.NDVI2018_gps50 + BLUP.LGN18bl_gps50 + BLUP.NOIw_gps50 + BLUP.p25_gps50 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high)

m2 <- lm(PHQ9_score ~ BLUP.NDVI2018_gps50 + BLUP.LGN18bl_gps50 + BLUP.NOIw_gps50 + BLUP.p25_gps50 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
           marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
           origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high, combined.data)

# 50m + area 50m variables
model.area50 <- as.formula(PHQ9_score ~ BLUP.NDVI2018_gps50 + BLUP.LGN18bl_gps50 + BLUP.NOIw_gps50 + BLUP.p25_gps50 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                             marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                             origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high + POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50)

m3a <- lm(PHQ9_score ~ BLUP.NDVI2018_gps50 + BLUP.LGN18bl_gps50 + BLUP.NOIw_gps50 + BLUP.p25_gps50 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high + POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50, combined.data)

# 50m + area 100m variables 
model.area100 <- as.formula(PHQ9_score ~ BLUP.NDVI2018_gps50 + BLUP.LGN18bl_gps50 + BLUP.NOIw_gps50 + BLUP.p25_gps50 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                              marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                              origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high + POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100)

m3b <- lm(PHQ9_score ~ BLUP.NDVI2018_gps50 + BLUP.LGN18bl_gps50 + BLUP.NOIw_gps50 + BLUP.p25_gps50 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high + POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100, combined.data)


formula_list <- list(model.unadj, model.indiv, model.area50, model.area100)

model_list <- list(m1, m2, m3a, m3b)

invisible(mapply(gps_results_func, model_list, formula_list))

sink()

# Moderation by sex

model.int <- update(m3b, .~. + BLUP.NDVI2018_gps50*sex_Female, data=combined.data)
summary(model.int)
lrtest(model.int, model.area100)

model.int <- update(m3b, .~. + BLUP.LGN18bl_gps50*sex_Female, data=combined.data)
summary(model.int)
lrtest(model.int, model.area100)

model.int <- update(m3b, .~. + BLUP.NOIw_gps50*sex_Female, data=combined.data)
summary(model.int)
lrtest(model.int, model.area100)

model.int <- update(m3b, .~. + BLUP.p25_gps50*sex_Female, data=combined.data)
summary(model.int)
lrtest(model.int, model.area100)

# GPS 100m
sink(file= file.path(path_out, "tab/multiple_gps_100m.csv"), append = F)

# Unadjusted
model.unadj <- as.formula(PHQ9_score ~ BLUP.NDVI2018_gps100 + BLUP.LGN18bl_gps100 + BLUP.NOIw_gps100 + BLUP.p25_gps100)

m1 <- lm(PHQ9_score ~ BLUP.NDVI2018_gps100 + BLUP.LGN18bl_gps100 + BLUP.NOIw_gps100 + BLUP.p25_gps100, combined.data)

# 100m + individual variables
model.indiv <- as.formula(PHQ9_score ~ BLUP.NDVI2018_gps100 + BLUP.LGN18bl_gps100 + BLUP.NOIw_gps100 + BLUP.p25_gps100 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high)

m2 <- lm(PHQ9_score ~ BLUP.NDVI2018_gps100 + BLUP.LGN18bl_gps100 + BLUP.NOIw_gps100 + BLUP.p25_gps100 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
           marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
           origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high, combined.data)

# 100m + area 50m variables
model.area50 <- as.formula(PHQ9_score ~ BLUP.NDVI2018_gps100 + BLUP.LGN18bl_gps100 + BLUP.NOIw_gps100 + BLUP.p25_gps100 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                             marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                             origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high + POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50)

m3a <- lm(PHQ9_score ~ BLUP.NDVI2018_gps100 + BLUP.LGN18bl_gps100 + BLUP.NOIw_gps100 + BLUP.p25_gps100 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high + POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50, combined.data)

# 100m + area 100m variables 

model.area100 <- as.formula(PHQ9_score ~ BLUP.NDVI2018_gps100 + BLUP.LGN18bl_gps100 + BLUP.NOIw_gps100 + BLUP.p25_gps100 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                              marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                              origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high + POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100)

m3b <- lm(PHQ9_score ~ BLUP.NDVI2018_gps100 + BLUP.LGN18bl_gps100 + BLUP.NOIw_gps100 + BLUP.p25_gps100 + age + sex_Female + educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high + POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100, combined.data)


formula_list <- list(model.unadj, model.indiv, model.area50, model.area100)

model_list <- list(m1, m2, m3a, m3b)

invisible(mapply(gps_results_func, model_list, formula_list))

sink()

# Moderation by sex

model.int <- update(m3b, .~. + BLUP.NDVI2018_gps100*sex_Female, data=combined.data)
summary(model.int)
lrtest(model.int, model.area100)

model.int <- update(m3b, .~. + BLUP.LGN18bl_gps100*sex_Female, data=combined.data)
summary(model.int)
lrtest(model.int, model.area100)

model.int <- update(m3b, .~. + BLUP.NOIw_gps100*sex_Female, data=combined.data)
summary(model.int)
lrtest(model.int, model.area100)

model.int <- update(m3b, .~. + BLUP.p25_gps100*sex_Female, data=combined.data)
summary(model.int)
lrtest(model.int, model.area100)

# Wald tests between fully adjusted residential and mobility-based models -----

combined_50 <- lm(PHQ9_score ~ NDVI2018_res50 + LGN18bl_res50 + NOIw_res50 + p25_res50 + NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50 + age + sex + empl + educ + marital + household + origin + income + POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100, data=combined.data)

regTermTest(combined_50, ~NDVI2018_res50+NDVI_gps50, method ="Wald") 
regTermTest(combined_50, ~LGN18bl_res50+LGN18bl_gps50, method ="Wald")
regTermTest(combined_50, ~NOIw_res50+NOIw_gps50, method ="Wald")
regTermTest(combined_50, ~p25_res50+p25_gps50, method ="Wald")

combined_100 <- lm(PHQ9_score ~ NDVI2018_res100 + LGN18bl_res100 + NOIw_res100 + p25_res100 + NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100 + age + sex + empl + educ + marital + household + origin + income + POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100, data=combined.data)

regTermTest(combined_100, ~NDVI2018_res100+NDVI_gps100, method ="Wald")
regTermTest(combined_100, ~LGN18bl_res100+LGN18bl_gps100, method ="Wald")
regTermTest(combined_100, ~NOIw_res100+NOIw_gps100, method ="Wald")
regTermTest(combined_100, ~p25_res100+p25_gps100, method ="Wald") 


## Mobility-based multiple regression (observed means) -----
# 50m

observed.means <- points %>%
  group_by(WE_ID_crypt) %>% 
  summarise_at(vars(NDVI2018_gps50:p25_gps100), mean)

observed.means <- left_join(survey.data.dummies, observed.means)

sink(file= file.path(path_out, "tab/biased_gps_50m.csv"), append = F)

# unadj

model.unadj <- as.formula(PHQ9_score ~ NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50)

m1 <- lm(PHQ9_score ~ NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50, observed.means)

# indiv

model.indiv <- as.formula(PHQ9_score ~ NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50 + age + sex_Female + 
                            educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high)

m2 <- lm(PHQ9_score ~ NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50 + age + sex_Female + 
           educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
           marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
           origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high, observed.means)

# area 50m

model.area50 <- as.formula(PHQ9_score ~ NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50 + age + sex_Female + 
                             educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                             marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                             origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high +
                             POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50)

m3a <- lm(PHQ9_score ~ NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50 + age + sex_Female + 
            educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high +
            POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50, observed.means)

# area 100m

model.area100 <- as.formula(PHQ9_score ~ NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50 + age + sex_Female + 
                              educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                              marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                              origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high +
                              POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100)

m3b <- lm(PHQ9_score ~ NDVI2018_gps50 + LGN18bl_gps50 + NOIw_gps50 + p25_gps50 + age + sex_Female + 
            educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high +
            POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100, observed.means)

formula_list <- list(model.unadj, model.indiv, model.area50, model.area100)

model_list <- list(m1, m2, m3a, m3b)

invisible(mapply(gps_results_func, model_list, formula_list))

sink()

# 100m 

sink(file= file.path(path_out, "tab/biased_gps_100m.csv"), append = F)

# unadj

model.unadj <- as.formula(PHQ9_score ~ NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100)

m1 <- lm(PHQ9_score ~ NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100, observed.means)

# indiv

model.indiv <- as.formula(PHQ9_score ~ NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100 + age + sex_Female + 
                            educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high)

m2 <- lm(PHQ9_score ~ NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100 + age + sex_Female + 
           educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
           marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
           origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high, observed.means)

# area 50m

model.area50 <- as.formula(PHQ9_score ~ NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100 + age + sex_Female + 
                             educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                             marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                             origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high +
                             POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50)

m3a <- lm(PHQ9_score ~ NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100 + age + sex_Female + 
            educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high +
            POPDENS16_res50 + FRAG16_res50 + DEPRI16_res50, observed.means)

# area 100m

model.area100 <- as.formula(PHQ9_score ~ NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100 + age + sex_Female + 
                              educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
                              marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
                              origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high +
                              POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100)

m3b <- lm(PHQ9_score ~ NDVI2018_gps100 + LGN18bl_gps100 + NOIw_gps100 + p25_gps100 + age + sex_Female + 
            educ_Mid + educ_High + empl_Unemployed + empl_Non_working + empl_Other + marital_Separated + marital_Widowed +
            marital_Unmarried + household_Couple_without_child + household_Couple_with_child + household_Single_parent +
            origin_Western + origin_Non_western + income_Low + income_Middle + income_High + income_Very_high +
            POPDENS16_res100 + FRAG16_res100 + DEPRI16_res100, observed.means)

formula_list <- list(model.unadj, model.indiv, model.area50, model.area100)

model_list <- list(m1, m2, m3a, m3b)

invisible(mapply(gps_results_func, model_list, formula_list))

sink()

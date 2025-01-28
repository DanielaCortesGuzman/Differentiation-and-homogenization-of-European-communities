#######################################################################################################################
## Recovering European river invertebrate communities homogenize or  differentiate depending on anthropogenic stress ##
#######################################################################################################################

library(ggplot2)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(ordbetareg)
library(mgcv)
library(gratia)
library(performance)
library(ggcorrplot)
library(ggExtra)
library(GGally)


dat<-read.csv(file.choose(), header=TRUE, check.names = FALSE) # Data1_communityBasins.csv

#### ---- DATA EXPLORATION

# Define factors

dat$basin<-factor(dat$basin)
dat$country<-factor(dat$country)

# Some basins do not have EQR information: subset data with EQR

dat.eqr<-dat[!is.na(dat$eqr),]

# Visualize land cover change over time

dat %>% ggplot() +
  geom_line(aes(y = urb_mean, x = year, group = basin), color = "lightgray", linewidth = 1) +
  labs(y = "Urban land use", x = "Year") +
  theme_classic() +
  theme(text = element_text(size = 20), legend.position = "none")

dat %>% ggplot() +
  geom_line(aes(y = crop_mean, x = year, group = basin), color = "lightgray", linewidth = 1) +
  labs(y = "Agricultural land use", x = "Year") +
  theme_classic() +
  theme(text = element_text(size = 20), legend.position = "none")

dat %>% ggplot() +
  geom_line(aes(y = forest_mean, x = year, group = basin), color = "lightgray", linewidth = 1) +
  labs(y = "Forest", x = "Year") +
  theme_classic() +
  theme(text = element_text(size = 20), legend.position = "none")

# Spatial variation of land cover and EQR

mean.urban<-dat %>% group_by(basin) %>% summarise(mean_urban = mean(urban_mean))
mean.crop<-dat %>% group_by(basin) %>% summarise(mean_crop = mean(crop_mean))
mean.forest<-dat %>% group_by(basin) %>% summarise(mean_forest = mean(forest_mean))
mean.eqr<-dat.eqr %>% group_by(basin) %>% summarise(mean_eqr = mean(eqr))

dat<-merge(dat, mean.urban, by = c("basin"))
dat<-merge(dat, mean.crop, by = c("basin"))
dat<-merge(dat, mean.forest, by = c("basin"))

dat.eqr<-merge(dat.eqr, mean.eqr, by = c("basin"))
dat.eqr<-merge(dat.eqr, mean.urban, by = c("basin"))
dat.eqr<-merge(dat.eqr, mean.crop, by = c("basin"))
dat.eqr<-merge(dat.eqr, mean.forest, by = c("basin"))

# Predictors' correlation and collinearity

MyVar <- c("temperature_trend", "mean_eqr", "mean_forest", "mean_crop", "mean_urban")
ggpairs(dat.eqr[,MyVar])
corvif(dat.eqr[,MyVar])

# Keep variables with cor<+- 0.5 and VIF<2: temperature_trend, eqr, and mean_urban

### --- Beta diversity distribution by predictors

# taxonomic_bray: taxonomic beta diversity, Bray Curtis index
# biological_bray: biological traits beta diversity, Bray Curtis index
# ecological_bray: ecological traits beta diversity, Bray Curtis index

dat %>% ggplot() +
  geom_point(aes(x = log(mean_urban), y = taxonomic_bray), color = "#008837") +
  geom_smooth(aes(x = log(mean_urban), y = taxonomic_bray), color = "#008837", fill = "#008837", alpha = 0.2) +
  geom_point(aes(x = log(mean_urban), y = biological_bray), color = "#5e3c99") +
  geom_smooth(aes(x = log(mean_urban), y = biological_bray), color = "#5e3c99", fill = "#5e3c99", alpha = 0.2) +
  geom_point(aes(x = log(mean_urban), y = ecological_bray), color = "#e66101") +
  geom_smooth(aes(x = log(mean_urban), y = ecological_bray), color = "#e66101", fill = "#e66101", alpha = 0.2) +
  theme_classic() +
  labs(x = "Proportion of urban areas", y = "β-diversity") +
  theme(text = element_text(size = 20))


dat.eqr %>% ggplot() +
  geom_point(aes(x = mean_eqr, y = taxonomic_bray), color = "#008837") +
  geom_smooth(aes(x = mean_eqr, y = taxonomic_bray), color = "#008837", fill = "#008837", alpha = 0.2) +
  geom_point(aes(x = mean_eqr, y = biological_bray), color = "#5e3c99") +
  geom_smooth(aes(x = mean_eqr, y = biological_bray), color = "#5e3c99", fill = "#5e3c99", alpha = 0.2) +
  geom_point(aes(x = mean_eqr, y = ecological_bray), color = "#e66101") +
  geom_smooth(aes(x = mean_eqr, y = ecological_bray), color = "#e66101", fill = "#e66101", alpha = 0.2) +
  theme_classic() +
  labs(x = "EQR", y = "β-diversity") +
  theme(text = element_text(size = 20))


dat %>% ggplot() +
  geom_point(aes(x = temperature_trend, y = taxonomic_bray), color = "#008837") +
  geom_smooth(aes(x = temperature_trend, y = taxonomic_bray), color = "#008837", fill = "#008837", alpha = 0.2) +
  geom_point(aes(x = temperature_trend, y = biological_bray), color = "#5e3c99") +
  geom_smooth(aes(x = temperature_trend, y = biological_bray), color = "#5e3c99", fill = "#5e3c99", alpha = 0.2) +
  geom_point(aes(x = temperature_trend, y = ecological_bray), color = "#e66101") +
  geom_smooth(aes(x = temperature_trend, y = ecological_bray), color = "#e66101", fill = "#e66101", alpha = 0.2) +
  theme_classic() +
  labs(x = "Temperature slope", y = "β-diversity") +
  theme(text = element_text(size = 20))


#### ----  MODELING

### --- Local trends: are communities recovering?

local.dat<-read.csv(file.choose(), header = TRUE) # Data2_communitySites.csv

local.dat<-local.dat %>% mutate(year_s = scale(year)[,1], fsite = as.factor(site))

year.pred<-seq(from = min(local.dat$year_s), to = max(local.dat$year_s), length = 100)
scale(local.dat$year)# get mean and sd
year.unscaled<-year.pred * 5.274925 + 2011.201

## -- Local EQR (community dissimilarity based on sensitive versus tolerant taxa)

local.dat$beqr<-ifelse(local.dat$eqr_site>=1, 0.9999, local.dat$eqr_site)
range(local.dat$beqr)
local.eqr<-local.dat[!is.na(local.dat$eqr_site), ]
local.dat<-local.dat %>% mutate(fsite = as.factor(site), year_s = scale(year)[,1])

predictors<-dat[, c(1,15,18,19)]
predictors<-predictors[!duplicated(predictors$basin), ]
local.eqr<-merge(local.eqr, predictors, by = c("basin"), all.x = TRUE)

eqr.local<-glmmTMB(beqr ~ year_s + (1|fsite), family = beta_family(link = "logit"), REML = TRUE, data = local.eqr)
summary(eqr.local)

# Predicted values

pred.eqr.local<-predict(eqr.local, newdata = data.frame("year_s"= year.pred, "fsite" = local.eqr$fsite[1]), type = "response", se.fit = TRUE)
pred.eqr.local<-data.frame(pred.eqr.local)
pred.eqr.local$year<-year.unscaled
pred.eqr.local$term<-"allSites"

eqr.basins<-unique(local.eqr$basin)
for(i in 1:length(eqr.basins)){
  eqr.local1<-glmmTMB(beqr ~ year_s + (1|fsite), family = beta_family(link = "logit"), REML = TRUE, data = subset(local.eqr, basin == eqr.basins[i]))
  pred.eqr.local1<-predict(eqr.local1, newdata = data.frame("year_s" = year.pred, "fsite" = subset(local.eqr,  basin == eqr.basins[i])$fsite[1]), type = "response", se.fit = TRUE)
  pred.eqr.local1<-data.frame(pred.eqr.local1)
  pred.eqr.local1$year<-year.pred.unscaled
  pred.eqr.local1$term<-eqr.basins[i]
  pred.eqr.local<-rbind(pred.eqr.local, pred.eqr.local1)
}

pred.eqr.local<-pred.eqr.local %>% relocate(year)
pred.eqr.local<-pred.eqr.local %>% mutate(ci05 = fit - (se.fit * 1.96), ci95 = fit + (se.fit * 1.96))

ggplot() +
  geom_line(aes(x = year, y = fit, group = term), data = filter(pred.eqr.local, !term == "allSites"), color = "gray") +
  geom_line(aes(x = year, y = fit), data = subset(pred.eqr.local, term == "allSites"), color = "#018571", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.eqr.local, term == "allSites"), fill = "#018571", alpha = 0.2, linetype = 0) +
  theme_classic() +
  theme(text = element_text(size = 50)) +
  labs(x = "Year", y = "EQR")

# - EQR trends by predictors

# EQR by urban areas

eqr.URBAN<-glmmTMB(beqr ~ year_s*urban_s + (1|fsite), REML = TRUE, family = beta_family(link = "logit"), data = local.eqr)
summary(eqr.URBAN)

# Predicted values

low.urban.pred<-predict(eqr.URBAN, newdata = data.frame("year_s" = year.pred, "urban_s" = -2, "fsite" = local.eqr$fsite[1]), type = "response", se.fit = TRUE)
low.urban.pred<-data.frame(low.urban.pred)
high.urban.pred<-predict(eqr.URBAN, newdata = data.frame("year_s" = year.pred, "urban_s" = 2, "fsite" = local.eqr$fsite[1]), type = "response", se.fit = TRUE)
high.urban.pred<-data.frame(high.urban.pred)

pred.eqr.urban<-cbind(low.urban.pred, high.urban.pred)
colnames(pred.eqr.urban)<-c("low.fit", "low.se", "high.fit", "high.se")
pred.eqr.urban<-pred.eqr.urban %>% mutate(low.ci05 = low.fit - (1.96 * low.se), low.ci95 = low.fit + (1.96 * low.se),
                                          high.ci05 = high.fit - (1.96 * high.se), high.ci95 = high.fit + (1.96 * high.se))
pred.eqr.urban$year<-year.unscaled
pred.eqr.urban<-pred.eqr.urban %>% relocate(year)

ggplot() +
  ylim(0,1) +
  geom_line(aes(y = low.fit, x = year), data = pred.eqr.urban, color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_line(aes(y = high.fit, x = year), data = pred.eqr.urban, color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_line(aes(x = year, y = fit), data = subset(pred.eqr.local, term == "allSites"), color = "black", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.eqr.local, term == "allSites"), fill = "black", alpha = 0.1, linetype = 0) +
  labs(y = "EQR", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# EQR by average basin EQR

eqr.EQR<-glmmTMB(beqr ~ year_s*eqr_s + (1|fsite), REML = TRUE, family = beta_family(link = "logit"), data = local.eqr)
summary(eqr.EQR)

# Predicted values

low.eqr.pred<-predict(eqr.EQR, newdata = data.frame("year_s" = year.pred, "eqr_s" = -2, "fsite" = local.eqr$fsite[1]), type = "response", se.fit = TRUE)
low.eqr.pred<-data.frame(low.eqr.pred)
high.eqr.pred<-predict(eqr.EQR, newdata = data.frame("year_s" = year.pred, "eqr_s" = 2, "fsite" = local.eqr$fsite[1]), type = "response", se.fit = TRUE)
high.eqr.pred<-data.frame(high.eqr.pred)

pred.eqr.eqr<-cbind(low.eqr.pred, high.eqr.pred)
colnames(pred.eqr.eqr)<-c("low.fit", "low.se", "high.fit", "high.se")
pred.eqr.eqr<-pred.eqr.eqr %>% mutate(low.ci05 = low.fit - (1.96 * low.se), low.ci95 = low.fit + (1.96 * low.se),
                                      high.ci05 = high.fit - (1.96 * high.se), high.ci95 = high.fit + (1.96 * high.se))
pred.eqr.eqr$year<-year.unscaled
pred.eqr.eqr<-pred.eqr.eqr %>% relocate(year)

ggplot() +
  ylim(0,1) +
  geom_line(aes(y = low.fit, x = year), data = pred.eqr.eqr, color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_line(aes(y = high.fit, x = year), data = pred.eqr.eqr, color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_line(aes(x = year, y = fit), data = subset(pred.eqr.local, term == "allSites"), color = "black", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.eqr.local, term == "allSites"), fill = "black", alpha = 0.1, linetype = 0) +
  labs(y = "EQR", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# EQR by temperature trend

eqr.TEMP<-glmmTMB(beqr ~ year_s*temp_s + (1|fsite), REML = TRUE, family = beta_family(link = "logit"), data = local.eqr)
summary(eqr.TEMP)

# Predicted values

low.temp.pred<-predict(eqr.TEMP, newdata = data.frame("year_s" = year.pred, "temp_s" = -2, "fsite" = local.eqr$fsite[1]), type = "response", se.fit = TRUE)
low.temp.pred<-data.frame(low.temp.pred)
high.temp.pred<-predict(eqr.TEMP, newdata = data.frame("year_s" = year.pred, "temp_s" = 2, "fsite" = local.eqr$fsite[1]), type = "response", se.fit = TRUE)
high.temp.pred<-data.frame(high.temp.pred)

pred.eqr.temp<-cbind(low.temp.pred, high.temp.pred)
colnames(pred.eqr.temp)<-c("low.fit", "low.se", "high.fit", "high.se")
pred.eqr.temp<-pred.eqr.temp %>% mutate(low.ci05 = low.fit - (1.96 * low.se), low.ci95 = low.fit + (1.96 * low.se),
                                        high.ci05 = high.fit - (1.96 * high.se), high.ci95 = high.fit + (1.96 * high.se))
pred.eqr.temp$year<-year.unscaled
pred.eqr.temp<-pred.eqr.temp %>% relocate(year)

ggplot() +
  ylim(0,1) +
  geom_line(aes(y = low.fit, x = year), data = pred.eqr.temp, color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_line(aes(y = high.fit, x = year), data = pred.eqr.temp, color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_line(aes(x = year, y = fit), data = subset(pred.eqr.local, term == "allSites"), color = "black", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.eqr.local, term == "allSites"), fill = "black", alpha = 0.1, linetype = 0) +
  labs(y = "EQR", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))



## -- Taxonomic richness

tric.local<-glmmTMB(taxon_richness ~ year_s + (1|fsite), family = poisson(link = "log"), REML = TRUE, data = local.dat)
summary(tric.local)

# Predicted values

pred.tric.local<-predict(tric.local, newdata = data.frame("year_s" = year.pred, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
pred.tric.local<-data.frame(pred.tric.local)
pred.tric.local$year<-year.unscaled
pred.tric.local$term<-"allSites"

for(i in 1:length(basins)){
  tric.local1<-glmmTMB(taxon_richness ~ year_s + (1|fsite), family = poisson(link = "log"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  pred.tric.local1<-predict(tric.local1, newdata=data.frame("year_s" = year.pred, "fsite" = subset(local.dat,  basin == basins[i])$fsite[1]), type = "response", se.fit = TRUE)
  pred.tric.local1<-data.frame(pred.tric.local1)
  pred.tric.local1$year<-year.unscaled
  pred.tric.local1$term<-basins[i]
  pred.tric.local<-rbind(pred.tric.local, pred.tric.local1)
}

pred.tric.local<-pred.tric.local %>% relocate(year)
pred.tric.local<-pred.tric.local %>% mutate(ci05 = fit - (se.fit * 1.96), ci95 = fit + (se.fit * 1.96))

ggplot() +
  geom_line(aes(x = year, y = fit, group = term), data = filter(pred.tric.local, !term == "allSites"), color = "gray") +
  geom_line(aes(x = year, y = fit), data = subset(pred.tric.local, term == "allSites"), color = "#008837", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.tric.local, term == "allSites"), fill = "#008837", alpha = 0.2, linetype = 0) +
  theme_classic() +
  theme(text = element_text(size = 50)) +
  labs(x = "Year", y = "Taxon richness")


# Taxon richness trends by predictors

predictors<-dat[, c(1,15,18,19)]
predictors<-predictors[!duplicated(predictors$basin), ]

local.dat<-merge(local.dat, predictors, by = c("basin"), all.x = TRUE)

# Taxon richness by urban areas

tric.URBAN<-glmmTMB(taxon_richness ~ year_s*urban_s + (1|fsite), REML = TRUE, family = poisson, data = local.dat)
summary(tric.URBAN)

# Predicted values

low.urb.pred<-predict(tric.URBAN, newdata = data.frame("year_s" = year.pred, "urban_s" = -2, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
low.urb.pred<-data.frame(low.urb.pred)
high.urb.pred<-predict(tric.URBAN, newdata = data.frame("year_s" = year.pred, "urban_s" = 2, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
high.urb.pred<-data.frame(high.urb.pred)

pred.tric.urban<-cbind(low.urb.pred, high.urb.pred)
colnames(pred.tric.urban)<-c("low.fit", "low.se", "high.fit", "high.se")
pred.tric.urban<-pred.tric.urban %>% mutate(low.ci05 = low.fit - (1.96 * low.se), low.ci95 = low.fit + (1.96 * low.se),
                                            high.ci05 = high.fit - (1.96 * high.se), high.ci95 = high.fit + (1.96 * high.se))
pred.tric.urban$year<-year.unscaled
pred.tric.urban<-pred.tric.urban %>% relocate(year)

ggplot() +
  ylim(8,31) +
  geom_line(aes(y = low.fit, x = year), data = pred.tric.urban, color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_line(aes(y = high.fit, x = year), data = pred.tric.urban, color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_line(aes(x = year, y = fit), data = subset(pred.tric.local, term == "allSites"), color = "black", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.tric.local, term == "allSites"), fill = "black", alpha = 0.1, linetype = 0) +
  labs(y = "Taxon richness", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# Taxon richness by average basin EQR

tric.EQR<-glmmTMB(taxon_richness ~ year_s*eqr_s + (1|fsite), REML = TRUE, family = poisson, data = local.dat)
summary(tric.EQR)

# Predicted values

low.tric.pred<-predict(tric.EQR, newdata = data.frame("year_s" = year.pred, "eqr_s" = -2, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
low.tric.pred<-data.frame(low.tric.pred)
high.tric.pred<-predict(tric.EQR, newdata = data.frame("year_s" = year.pred, "eqr_s" = 2, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
high.tric.pred<-data.frame(high.tric.pred)

pred.tric.eqr<-cbind(low.tric.pred, high.tric.pred)
colnames(pred.tric.eqr)<-c("low.fit", "low.se", "high.fit", "high.se")
pred.tric.eqr<-pred.tric.eqr %>% mutate(low.ci05 = low.fit - (1.96 * low.se), low.ci95 = low.fit + (1.96 * low.se),
                                        high.ci05 = high.fit - (1.96 * high.se), high.ci95 = high.fit + (1.96 * high.se))
pred.tric.eqr$year<-year.unscaled
pred.tric.eqr<-pred.tric.eqr %>% relocate(year)

ggplot() +
  ylim(8,31) +
  geom_line(aes(y = low.fit, x = year), data = pred.tric.eqr, color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_line(aes(y = high.fit, x = year), data = pred.tric.eqr, color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_line(aes(x = year, y = fit), data = subset(pred.tric.local, term == "allSites"), color = "black", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.tric.local, term == "allSites"), fill = "black", alpha = 0.1, linetype = 0) +
  labs(y = "Taxon richness", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# Taxon richness by temperature trend

tric.TEMP<-glmmTMB(taxon_richness ~ year_s*temp_s + (1|fsite), REML = TRUE, family = poisson, data = local.dat)
summary(tric.TEMP)

# Predicted values

low.temp.pred<-predict(tric.TEMP, newdata = data.frame("year_s" = year.pred, "temp_s" = -2, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
low.temp.pred<-data.frame(low.temp.pred)
high.temp.pred<-predict(tric.TEMP, newdata = data.frame("year_s" = year.pred, "temp_s" = 2, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
high.temp.pred<-data.frame(high.temp.pred)

pred.tric.temp<-cbind(low.temp.pred, high.temp.pred)
colnames(pred.tric.temp)<-c("low.fit", "low.se", "high.fit", "high.se")
pred.tric.temp<-pred.tric.temp %>% mutate(low.ci05 = low.fit - (1.96 * low.se), low.ci95 = low.fit + (1.96 * low.se),
                                          high.ci05 = high.fit - (1.96 * high.se), high.ci95 = high.fit + (1.96 * high.se))
pred.tric.temp$year<-year.unscaled
pred.tric.temp<-pred.tric.temp %>% relocate(year)

ggplot() +
  ylim(8,31) +
  geom_line(aes(y = low.fit, x = year), data = pred.tric.temp, color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_line(aes(y = high.fit, x = year), data = pred.tric.temp, color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_line(aes(x = year, y = fit), data = subset(pred.tric.local, term == "allSites"), color = "black", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.tric.local, term == "allSites"), fill = "black", alpha = 0.1, linetype = 0) +
  labs(y = "Taxon richness", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))



## -- Biological traits richness

bio.local<-glmmTMB(biological_richness ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = local.dat)
summary(bio.local)

# Predicted values

pred.bio.local<-predict(bio.local, newdata = data.frame("year_s"= year.pred, "fsite"=local.dat$fsite[1]), type = "response", se.fit = TRUE)
pred.bio.local<-data.frame(pred.bio.local)
pred.bio.local$year<-year.unscaled
pred.bio.local$term<-"allSites"

for(i in 1:length(basins)){
  bio.local1<-glmmTMB(biological_richness ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  pred.bio.local1<-predict(bio.local1, newdata = data.frame("year_s" = year.pred, "fsite" = subset(local.dat,  basin == basins[i])$fsite[1]), type = "response", se.fit = TRUE)
  pred.bio.local1<-data.frame(pred.bio.local1)
  pred.bio.local1$year<-year.unscaled
  pred.bio.local1$term<-basins[i]
  pred.bio.local<-rbind(pred.bio.local, pred.bio.local1)
}

pred.bio.local<-pred.bio.local %>% relocate(year)
pred.bio.local<-pred.bio.local %>% mutate(ci05 = fit - (se.fit * 1.96), ci95 = fit + (se.fit * 1.96))

ggplot() +
  geom_line(aes(x = year, y = fit, group = term), data = filter(pred.bio.local, !term == "allSites"), color = "gray") +
  geom_line(aes(x = year, y = fit), data = subset(pred.bio.local, term == "allSites"), color = "#5e3c99", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.bio.local, term == "allSites"), fill = "#5e3c99", alpha = 0.2, linetype = 0) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Year", y = "Trait richness")



## -- Ecological traits richness

eric.local<-glmmTMB(ecological_richness ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = local.dat)
summary(eric.local)

# Predicted values

pred.eric.local<-predict(eric.local, newdata = data.frame("year_s"= year.pred, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
pred.eric.local<-data.frame(pred.eric.local)
pred.eric.local$year<-year.unscaled
pred.eric.local$term<-"allSites"

for(i in 1:length(basins)){
  eric.local1<-glmmTMB(ecological_richness ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  pred.eric.local1<-predict(eric.local1, newdata = data.frame("year_s"= year.pred, "fsite"=subset(local.dat,  basin == basins.v[i])$fsite[1]), type = "response", se.fit = TRUE)
  pred.eric.local1<-data.frame(pred.eric.local1)
  pred.eric.local1$year<-year.unscaled
  pred.eric.local1$term<-basins[i]
  pred.eric.local<-rbind(pred.eric.local, pred.eric.local1)
}

pred.eric.local<-pred.eric.local %>% relocate(year)
pred.eric.local<-pred.eric.local %>% mutate(ci05 = fit - (se.fit * 1.96), ci95 = fit + (se.fit * 1.96))

ggplot() +
  geom_line(aes(x = year, y = fit, group = term), data = filter(pred.eric.local, !term == "allSites"), color = "gray") +
  geom_line(aes(x = year, y = fit), data = subset(pred.eric.local, term == "allSites"), color = "#e66101", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.eric.local, term == "allSites"), fill = "#e66101", alpha = 0.2, linetype = 0) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Year", y = "Trait richness")



## -- Taxonomic diversity as Hill number 1 (Shannon entropy)

tax.hill<-glmmTMB(taxonomic_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = local.dat)
summary(tax.hill)

# Predicted values

pred.tax.hill<-predict(tax.hill, newdata = data.frame("year_s" = year.pred, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
pred.tax.hill<-data.frame(pred.tax.hill)
pred.tax.hill$year<-year.unscaled
pred.tax.hill$term<-"allSites"

for(i in 1:length(basins)){
  tax.hill1<-glmmTMB(taxonomic_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  pred.tax.hill1<-predict(tax.hill1, newdata = data.frame("year_s" = year.pred, "fsite" = subset(local.dat,  basin == basins[i])$fsite[1]), type = "response", se.fit = TRUE)
  pred.tax.hill1<-data.frame(pred.tax.hill11)
  pred.tax.hill1$year<-year.unscaled
  pred.tax.hill1$term<-basins[i]
  pred.tax.hill<-rbind(pred.tax.hill1, pred.tax.hill11)
}

pred.tax.hill<-pred.tax.hill %>% relocate(year)
pred.tax.hill<-pred.tax.hill %>% mutate(ci05 = fit - (se.fit * 1.96), ci95 = fit + (se.fit * 1.96))

ggplot() +
  geom_line(aes(x = year, y = fit, group = term), data = filter(pred.tax.hill, !term == "allSites"), color = "gray") +
  geom_line(aes(x = year, y = fit), data = subset(pred.tax.hill, term == "allSites"), color = "#008837", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.tax.hill, term == "allSites"), fill = "#008837", alpha = 0.2, linetype = 0) +
  theme_classic() +
  theme(text = element_text(size = 50)) +
  labs(x = "Year", y = "Diversity")



## -- Biological trait diversity as Hill number 1 (Shannon entropy)

bio.hill<-glmmTMB(biological_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = local.dat)
summary(bio.hill)

# Predicted values

pred.bio.hill<-predict(bio.hill, newdata = data.frame("year_s" = year.pred, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
pred.bio.hill<-data.frame(pred.bio.hill)
pred.bio.hill$year<-year.unscaled
pred.bio.hill$term<-"allSites"

for(i in 1:length(basins)){
  bio.hill1<-glmmTMB(biological_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  pred.bio.hill1<-predict(bio.hill1, newdata = data.frame("year_s" = year.pred, "fsite" = subset(local.dat,  basin == basins[i])$fsite[1]), type = "response", se.fit = TRUE)
  pred.bio.hill1<-data.frame(pred.bio.hill1)
  pred.bio.hill1$year<-year.unscaled
  pred.bio.hill1$term<-basins[i]
  pred.bio.hill<-rbind(pred.bio.hill, pred.bio.hill1)
}

pred.bio.hill<-pred.bio.hill %>% relocate(year)
pred.bio.hill<-pred.bio.hill %>% mutate(ci05 = fit - (se.fit * 1.96), ci95 = fit + (se.fit * 1.96))

ggplot() +
  geom_line(aes(x = year, y = fit, group = term), data = filter(pred.bio.hill, !term == "allSites"), color = "gray") +
  geom_line(aes(x = year, y = fit), data = subset(pred.bio.hill, term == "allSites"), color = "#5e3c99", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.bio.hill, term == "allSites"), fill = "#5e3c99", alpha = 0.2, linetype = 0) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Year", y = "Diversity")



## -- Ecological trait diversity as Hill number 1 (Shannon entropy)

eco.hill<-glmmTMB(ecological_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = local.dat)
summary(eco.hill)

# Predicted values
pred.eco.hill<-predict(eco.hill, newdata = data.frame("year_s" = year.pred, "fsite" = local.dat$fsite[1]), type = "response", se.fit = TRUE)
pred.eco.hill<-data.frame(pred.eco.hill)
pred.eco.hill$year<-year.unscaled
pred.eco.hill$term<-"allSites"

for(i in 1:length(basins)){
  eco.hill1<-glmmTMB(ecological_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  pred.eco.hill1<-predict(eco.hill1, newdata = data.frame("year_s" = year.pred, "fsite" = subset(local.dat,  basin == basins[i])$fsite[1]), type = "response", se.fit = TRUE)
  pred.eco.hill1<-data.frame(pred.eco.hill1)
  pred.eco.hill1$year<-year.unscaled
  pred.eco.hill1$term<-basins[i]
  pred.eco.hill<-rbind(pred.eco.hill, pred.eco.hill1)
}

pred.eco.hill<-pred.eco.hill %>% relocate(year)
pred.eco.hill<-pred.eco.hill %>% mutate(ci05 = fit - (se.fit * 1.96), ci95 = fit + (se.fit * 1.96))

ggplot() +
  geom_line(aes(x = year, y = fit, group = term), data = filter(pred.eco.hill, !term == "allSites"), color = "gray") +
  geom_line(aes(x = year, y = fit), data = subset(pred.eco.hill, term == "allSites"), color = "#e66101", linewidth = 2) +
  geom_ribbon(aes(x = year, ymin = ci05, ymax = ci95), data = subset(pred.eco.hill, term == "allSites"), fill = "#e66101", alpha = 0.2, linetype = 0) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(x = "Year", y = "Diversity")




### --- Beta diversity trends: are communities differentiating?

# Scale predictors

dat<-dat %>% mutate(year_s = scale(year)[,1], temp_s = scale(temperature_trend)[,1], urban_s = scale(mean_urban)[,1])

dat.eqr<-dat.eqr %>% mutate(year_s = scale(year)[,1], temp_s = scale(temperature_trend)[,1], eqr_s = scale(mean_eqr)[,1], urban_s = scale(mean_urban)[,1])

## -- Taxonomic composition

# - Mean trend

tax.mean<-glmmTMB(taxonomic_bray ~ year_s + (1|basin),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(tax.mean)
Esqr <- simulateResiduals(fittedModel = tax.mean, plot = FALSE)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Esqr, testUniformity = TRUE, testOutliers = TRUE, testDispersion = FALSE)

# Predicted values

year.pred<-seq(from=min(dat$year_s), to=max(dat$year_s), length=100)
scale(dat$year)# get mean and sd
year.unscaled<-year.pred * 5.614672 + 2010.45

pred.mean.tax<-predict(tax.mean, newdata=data.frame("year_s" = year.pred, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
pred.mean.tax<-data.frame(pred.mean.tax)
pred.mean.tax$year<-year.unscaled
pred.mean.tax<-pred.mean.tax %>% relocate(year)
pred.mean.tax$index<-"taxonomic"

# - Trends (over time) in beta diversity per basin

tax.year<-glmmTMB(taxonomic_bray ~ basin + basin:year_s,
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(tax.year)

# Predicted values

y.pred.cr01<-predict(tax.year, newdata = data.frame("year_s" = x.pred, "basin" = "CR01"), type = "response", se.fit = TRUE)
tax.trends<-data.frame(y.pred.cr01$fit)
tax.trends$year<-year.unscaled
tax.trends<-tax.trends %>% relocate(year)

basins<-unique(dat$basin)
for (i in 2:length(basins)){
  y.pred<-predict(tax.year, newdata = data.frame("year_s" = x.pred, "basin" = basins[i]), type = "response", se.fit = TRUE)
  tax.trends<-cbind(tax.trends, y.pred$fit)
}

basins.names<-as.vector(basins)
colnames(tax.trends)[2:ncol(tax.trends)]<-c(basins.names)
tax.trends<-tax.trends %>% gather(key = "basin", value = "fit", 2:ncol(tax.trends))
tax.trends$index<-"taxonomic"

# Visualize mean trend and trends per basin

colnames(pred.mean.tax)<-c("year", "fit.mean", "se.mean", "index")
tax<-merge(pred.mean.tax, tax.trends, by = c("year", "index"), all = TRUE)

tax %>% ggplot(aes(x = year, y = fit, group = basin)) +
  ylim (0 ,1) +
  geom_line(linewidth = 1, color = "gray") +
  geom_line(aes(x = year, y = fit.mean), color = "#008837", linewidth = 3, linetype = "longdash") +
  labs(y = "β-diversity", x = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 20))


# - Trends in beta diversity linked to urban areas

tax.urban<-glmmTMB(taxonomic_bray ~ year_s * urban_s + (1|basin),
                   family = beta_family(link = "logit"),
                   REML = TRUE,
                   data = dat)
summary(tax.urban)

# Test spatial autocorrelation

Esqr <- simulateResiduals(fittedModel = tax.urban, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$basin, seed = 123, rotation = "estimated")
testSpatialAutocorrelation(rsims, x = unique(dat$lon_center), y = unique(dat$lat_center), plot = FALSE)

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated")
testTemporalAutocorrelation(rsims, time = unique(dat$year_s), alternative = c("two.sided"), plot = FALSE)

# Taxonomic resolution test

tax.urban.resol<-glmmTMB(tax_bray ~ year_s * urban_s + taxon_resolution + (1|basin),
                         family = beta_family(link = "logit"),
                         REML = TRUE,
                         data = dat)
summary(tax.urban.resol)

# Predicted year effect given values of predictor

urban.pred<-predict(tax.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = mean(dat$urban_s), "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
urban.pred.tax<-data.frame(urban.pred)
urban.pred.low<-predict(tax.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = -2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
urban.pred.high<-predict(tax.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = 2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)

urban.pred.tax<-cbind(urban.pred.tax, urban.pred.low)
urban.pred.tax<-cbind(urban.pred.tax, urban.pred.high)
colnames(urban.pred.tax)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
urban.pred.tax<-urban.pred.tax %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                          low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                          low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
urban.pred.tax$year<-year.unscaled
urban.pred.tax<-urban.pred.tax %>% relocate(year)

urban.pred.tax %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#2c7bb6", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#d7191c", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.1, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# - Trends in beta diversity linked to EQR (average EQR per basin)

tax.eqr<-glmmTMB(taxonomic_bray ~ year_s * eqr_s + (1|basin),
                 family = beta_family(link = "logit"),
                 REML = TRUE,
                 data = dat.eqr)
summary(tax.eqr)

# Test spatial autocorrelation

Esqr <- simulateResiduals(fittedModel = tax.eqr, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat.eqr$basin, seed = 123, rotation = "estimated")
testSpatialAutocorrelation(rsims, x = unique(dat.eqr$lon_center), y = unique(dat.eqr$lat_center), plot = FALSE)

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat.eqr$year_s, seed = 123, rotation = "estimated")
testTemporalAutocorrelation(rsims, time = unique(dat.eqr$year_s), alternative = c("two.sided"), plot = FALSE)

# Taxonomic resolution test

tax.eqr.resol<-glmmTMB(tax_bray ~ year_s * eqr_s + taxon_resolution + (1|basin),
                       family = beta_family(link = "logit"),
                       REML = TRUE,
                       data = dat.eqr)
summary(tax.eqr.resol)

# Predicted year effect given values of predictor

eqr.pred<-predict(tax.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = mean(dat.eqr$eqr_s), "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)
eqr.pred.tax<-data.frame(eqr.pred)
eqr.pred.low<-predict(tax.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = -2, "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)
eqr.pred.high<-predict(tax.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = 2, "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)

eqr.pred.tax<-cbind(eqr.pred.tax, eqr.pred.low)
eqr.pred.tax<-cbind(eqr.pred.tax, eqr.pred.high)
colnames(eqr.pred.tax)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
eqr.pred.tax<-eqr.pred.tax %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                      low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                      low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
eqr.pred.tax$year<-year.unscaled
eqr.pred.tax<-eqr.pred.tax %>% relocate(year)

eqr.pred.tax %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#d7191c", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#2c7bb6", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.1, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# - Trends in beta diversity linked to temperature trends

tax.temp<-glmmTMB(taxonomic_bray ~ year_s * temp_s + (1|basin),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(tax.temp)

# Test spatial autocorrelation

Esqr <- simulateResiduals(fittedModel = tax.temp, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$basin, seed = 123, rotation = "estimated") # recalculate for repeated measures
testSpatialAutocorrelation(rsims, x = unique(dat$lon_center), y = unique(dat$lat_center), plot = FALSE) # coordinates need to be unique

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s), alternative = c("two.sided"), plot = FALSE)

# Test for taxonomic resolution

tax.temp.resol<-glmmTMB(taxonomic_bray ~ year_s * temp_s + taxon_resolution + (1|basin),
                        family = beta_family(link = "logit"),
                        REML = TRUE,
                        data = dat)
summary(tax.temp.resol) # Qualitatively equal

# Predicted year effect given values of predictor

temp.pred<-predict(tax.temp, newdata=data.frame("year_s" = year.pred, "temp_s" = mean(dat$temp_s), "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
temp.pred.tax<-data.frame(temp.pred)
temp.pred.low<-predict(tax.temp, newdata=data.frame("year_s" = year.pred, "temp_s" = -2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
temp.pred.high<-predict(tax.temp, newdata=data.frame("year_s" = year.pred, "temp_s" = 2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)

temp.pred.tax<-cbind(temp.pred.tax, temp.pred.low)
temp.pred.tax<-cbind(temp.pred.tax, temp.pred.high)
colnames(temp.pred.tax)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
temp.pred.tax<-temp.pred.tax %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                        low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                        low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
temp.pred.tax$year<-year.unscaled
temp.pred.tax<-temp.pred.tax %>% relocate(year)

temp.pred.tax %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.1, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))



## -- Biological composition

# - Mean trend

bio.mean<-glmmTMB(biological_bray ~ year_s + (1|basin),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(bio.mean)
Esqr <- simulateResiduals(fittedModel = bio.mean, plot = FALSE)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Esqr, testUniformity = TRUE, testOutliers = TRUE, testDispersion = FALSE)

# Predicted values

pred.mean.bio<-predict(bio.mean, newdata = data.frame("year_s" = year.pred, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
pred.mean.bio<-data.frame(pred.mean.bio)
pred.mean.bio$year<-year.unscaled
pred.mean.bio<-pred.mean.bio %>% relocate(year)
pred.mean.bio$index<-"biological"

# - Trends per basin

bio.year<-glmmTMB(biological_bray ~ basin + basin:year_s,
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(bio.year)

# Predicted trends per basin

y.pred.cr01<-predict(bio.year, newdata = data.frame("year_s" = year.pred, "basin" = "CR01"), type = "response", se.fit = TRUE)
bio.trends<-data.frame(y.pred.cr01$fit)
bio.trends$year<-year.unscaled
bio.trends<-bio.trends %>% relocate(year)

for (i in 2:length(basins)){
  y.pred<-predict(bio.year, newdata = data.frame("year_s" = year.pred, "basin" = basins[i]), type = "response", se.fit = TRUE)
  bio.trends<-cbind(bio.trends, y.pred$fit)
}
colnames(bio.trends)[2:ncol(bio.trends)]<-c(basins.names)
bio.trends<-bio.trends %>% gather(key = "basin", value = "fit", 2:ncol(bio.trends))
bio.trends$index<-"biological"

# Visualize mean trend and by basin

colnames(pred.mean.bio)<-c("year", "fit.mean", "se.mean", "index")
bio<-merge(pred.mean.bio, bio.trends, by = c("year", "index"), all = TRUE)

ggplot(bio, aes(x = year, y = fit, group = basin)) +
  ylim (0 ,1) +
  geom_line(linewidth = 1, color = "gray") +
  geom_line(data = subset(bio, basin == "CR01"), aes(x = year, y = fit.mean), color = "#5e3c99", linewidth = 3)+
  labs(y = "β-diversity", x = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 20))


# - Trends in beta diversity linked to urban areas

bio.urban<-glmmTMB(biological_bray ~ year_s * urban_s + (1|basin),
                   family = beta_family(link = "logit"),
                   REML = TRUE,
                   data = dat)
summary(bio.urban)

# Autocorrelation test

Esqr <- simulateResiduals(fittedModel = bio.urban, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$basin, seed = 123, rotation = "estimated")
testSpatialAutocorrelation(rsims, x = unique(dat$lon_center), y = unique(dat$lat_center), plot = FALSE)

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated")
testTemporalAutocorrelation(rsims, time = unique(dat$year_s), alternative = c("two.sided"), plot = F)

# Taxonomic resolution test

bio.urban.resol<-glmmTMB(biological_bray ~ year_s * urban_s + taxon_resolution + (1|basin),
                         family = beta_family(link = "logit"),
                         REML = TRUE,
                         data = dat)
summary(bio.urban.resol)

# Predicted year effect given values of predictor

urban.pred<-predict(bio.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = mean(dat$urban_s), "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
urban.pred.bio<-data.frame(urban.pred)
urban.pred.low<-predict(bio.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = -2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
urban.pred.high<-predict(bio.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = 2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)

urban.pred.bio<-cbind(urban.pred.bio, urban.pred.low)
urban.pred.bio<-cbind(urban.pred.bio, urban.pred.high)
colnames(urban.pred.bio)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
urban.pred.bio<-urban.pred.bio %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                          low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                          low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
urban.pred.bio$year<-year.unscaled
urban.pred.bio<-urban.pred.bio %>% relocate(year)

urban.pred.bio %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.15, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.15, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# - Trends in beta diversity linked to EQR

bio.eqr<-glmmTMB(biological_bray ~ year_s * eqr_s + (1|basin),
                 family = beta_family(link = "logit"),
                 REML = TRUE,
                 data = dat.eqr)
summary(bio.eqr)

# Test spatial autocorrelation

Esqr <- simulateResiduals(fittedModel = bio.eqr, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat.eqr$basin, seed = 123, rotation = "estimated")
testSpatialAutocorrelation(rsims, x = unique(dat.eqr$lon_center), y = unique(dat.eqr$lat_center), plot = FALSE)

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat.eqr$year_s, seed = 123, rotation = "estimated")
testTemporalAutocorrelation(rsims, time = unique(dat.eqr$year_s), alternative = c("two.sided"), plot = F)

# Taxonomic resolution test

bio.eqr.resol<-glmmTMB(biological_bray ~ year_s * eqr_s + taxon_resolution + (1|basin),
                       family = beta_family(link = "logit"),
                       REML = TRUE,
                       data = dat.eqr)
summary(bio.eqr.resol)

# Predicted year effect given values of predictor

eqr.pred<-predict(bio.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = mean(dat.eqr$eqr_s), "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)
eqr.pred.bio<-data.frame(eqr.pred)
eqr.pred.low<-predict(bio.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = -2, "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)
eqr.pred.high<-predict(bio.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = 2, "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)

eqr.pred.bio<-cbind(eqr.pred.bio, eqr.pred.low)
eqr.pred.bio<-cbind(eqr.pred.bio, eqr.pred.high)
colnames(eqr.pred.bio)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
eqr.pred.bio<-eqr.pred.bio %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                      low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                      low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
eqr.pred.bio$year<-year.unscaled
eqr.pred.bio<-eqr.pred.bio %>% relocate(year)

eqr.pred.bio %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#d7191c", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#2c7bb6", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.1, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# - Trends in beta diversity linked to temperature trends

bio.temp<-glmmTMB(biological_bray ~ year_s * temp_s + (1|basin),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(bio.temp)

# Test spatial autocorrelation

Esqr <- simulateResiduals(fittedModel = bio.temp, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$basin, seed = 123, rotation = "estimated")
testSpatialAutocorrelation(rsims, x = unique(dat$lon_center), y = unique(dat$lat_center), plot = FALSE)

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated")
testTemporalAutocorrelation(rsims, time = unique(dat$year_s), alternative = c("two.sided"), plot = FALSE)

# Taxonomic resolution test

bio.temp.resol<-glmmTMB(biological_bray ~ year_s * temp.slope + taxon_resolution + (1|basin),
                        family = beta_family(link = "logit"),
                        REML = TRUE,
                        data = dat)
summary(bio.temp.resol)

# Predicted year effect given values of predictor

temp.pred<-predict(bio.temp, newdata = data.frame("year_s" = year.pred, "temp_s" = mean(dat$temp_s), "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
temp.pred.bio<-data.frame(temp.pred)
temp.pred.low<-predict(bio.temp, newdata = data.frame("year_s" = year.pred, "temp_s" = -2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
temp.pred.high<-predict(bio.temp, newdata = data.frame("year_s" = year.pred, "temp_s" = 2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)

temp.pred.bio<-cbind(temp.pred.bio, temp.pred.low)
temp.pred.bio<-cbind(temp.pred.bio, temp.pred.high)
colnames(temp.pred.bio)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
temp.pred.bio<-temp.pred.bio %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                        low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                        low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
temp.pred.bio$year<-year.unscaled
temp.pred.bio<-temp.pred.bio %>% relocate(year)

temp.pred.bio %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.1, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))



## -- Ecological composition

# - Mean trend

eco.mean<-glmmTMB(ecological_bray ~ year_s + (1|basin),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(eco.mean)
Esqr <- simulateResiduals(fittedModel = eco.mean, plot = FALSE)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Esqr, testUniformity = TRUE, testOutliers = TRUE, testDispersion = FALSE)

# Predicted mean

pred.mean.eco<-predict(eco.mean, newdata = data.frame("year_s" = year.pred, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
pred.mean.eco<-data.frame(pred.mean.eco)
pred.mean.eco$year<-year.unscaled
pred.mean.eco<-pred.mean.eco %>% relocate(year)
pred.mean.eco$index<-"ecological"

# - Trends per basin

eco.year<-glmmTMB(ecological_bray ~ basin + basin:year_s,
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(eco.year)

# Predicted trend per basin

y.pred.cr01<-predict(eco.year, newdata = data.frame("year_s" = year.pred, "basin" = "CR01"), type = "response", se.fit = TRUE)
eco.trends<-data.frame(y.pred.cr01$fit)
eco.trends$year<-year.unscaled
eco.trends<-eco.trends %>% relocate(year)

for (i in 2:length(basins)){
  y.pred<-predict(eco.year, newdata = data.frame("year_s" = year.pred, "basin" = basins[i]), type = "response", se.fit = TRUE)
  eco.trends<-cbind(eco.trends, y.pred$fit)
}
colnames(eco.trends)[2:ncol(eco.trends)]<-c(basins.names)
eco.trends<-eco.trends %>% gather(key = "basin", value = "fit", 2:ncol(eco.trends))
eco.trends$index<-"ecological"

# Visualize mean trend and by basin

colnames(pred.mean.eco)<-c("year", "fit.mean", "se.mean", "index")
eco<-merge(pred.mean.eco, eco.trends, by = c("year", "index"), all = TRUE)

eco %>% ggplot(aes(x = year, y = fit, group = basin)) +
  ylim (0 ,1) +
  geom_line(linewidth = 1, color = "gray") +
  geom_line(aes(x = year, y = fit.mean), color = "#e66101", linewidth = 3, linetype = "longdash")+
  labs(y = "β-diversity", x = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 20))


# - Trends in beta diversity linked to urban areas

eco.urban<-glmmTMB(ecological_bray ~ year_s * urban_s + (1|basin),
                   family = beta_family(link = "logit"),
                   REML = TRUE,
                   data = dat)
summary(eco.urban)

# Autocorrelation test

Esqr <- simulateResiduals(fittedModel = eco.urban, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$basin, seed = 123, rotation = "estimated")
testSpatialAutocorrelation(rsims, x = unique(dat$lon_center), y = unique(dat$lat_center), plot = FALSE)

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated")
testTemporalAutocorrelation(rsims, time = unique(dat$year_s), alternative = c("two.sided"), plot = FALSE)

# Taxonomic resolution test

eco.urb.resol<-glmmTMB(ecological_bray ~ year_s * urban_s + taxon_resolution + (1|basin),
                       family = beta_family(link = "logit"),
                       REML = TRUE,
                       data = dat)
summary(eco.urb.resol)

# Predicted year effect given values of predictor

urban.pred<-predict(eco.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = mean(dat$urban_s), "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
urban.pred.eco<-data.frame(urban.pred)
urban.pred.low<-predict(eco.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = -2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
urban.pred.high<-predict(eco.urban, newdata = data.frame("year_s" = year.pred, "urban_s" = 2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)

urban.pred.eco<-cbind(urban.pred.eco, urban.pred.low)
urban.pred.eco<-cbind(urban.pred.eco, urban.pred.high)
colnames(urban.pred.eco)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
urban.pred.eco<-urban.pred.eco %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                          low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                          low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
urban.pred.eco$year<-year.unscaled
urban.pred.eco<-urban.pred.eco %>% relocate(year)

urban.pred.eco %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#2c7bb6", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#d7191c", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.1, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# - Trends in beta diversity linked to EQR

eco.eqr<-glmmTMB(ecological_bray ~ year_s * eqr_s + (1|basin),
                 family = beta_family(link = "logit"),
                 REML = TRUE,
                 data = dat.eqr)
summary(eco.eqr)

# Autocorrelation test

Esqr <- simulateResiduals(fittedModel = eco.eqr, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat.eqr$basin, seed = 123, rotation = "estimated")
testSpatialAutocorrelation(rsims, x = unique(dat.eqr$lon_center), y = unique(dat.eqr$lat_center), plot = FALSE)

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat.eqr$year_s, seed = 123, rotation = "estimated")
testTemporalAutocorrelation(rsims, time = unique(dat.eqr$year_s), alternative = c("two.sided"), plot = FALSE)

# Taxonomic resolution test

eco.eqr.resol<-glmmTMB(ecological_bray ~ year_s * eqr_s + taxon_resolution + (1|basin),
                       family = beta_family(link = "logit"),
                       REML = TRUE,
                       data = dat.eqr)
summary(eco.eqr.resol)

# Predicted year effect given values of predictor

eqr.pred<-predict(eco.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = mean(dat.eqr$eqr_s), "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)
eqr.pred.eco<-data.frame(eqr.pred)
eqr.pred.low<-predict(eco.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = -2, "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)
eqr.pred.high<-predict(eco.eqr, newdata = data.frame("year_s" = year.pred, "eqr_s" = 2, "basin" = dat.eqr$basin[1]), type = "response", se.fit = TRUE)

eqr.pred.eco<-cbind(eqr.pred.eco, eqr.pred.low)
eqr.pred.eco<-cbind(eqr.pred.eco, eqr.pred.high)
colnames(eqr.pred.eco)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
eqr.pred.eco<-eqr.pred.eco %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                      low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                      low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
eqr.pred.eco$year<-year.unscaled
eqr.pred.eco<-eqr.pred.eco %>% relocate(year)

eqr.pred.eco %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.1, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


# - Trends in beta diversity linked to temperature trends

eco.temp<-glmmTMB(ecological_bray ~ year_s * temp_s + (1|basin),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(eco.temp)

# Test spatial autocorrelation

Esqr <- simulateResiduals(fittedModel = eco.temp, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$basin, seed = 123, rotation = "estimated")
testSpatialAutocorrelation(rsims, x = unique(dat$lon_center), y = unique(dat$lat_center), plot = FALSE)

# Test temporal autocorrelation

rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated")
testTemporalAutocorrelation(rsims, time = unique(dat$year_s), alternative = c("two.sided"), plot = FALSE)

# Taxonomic resolution test

eco.temp.resol<-glmmTMB(ecological_bray ~ year_s * temp_s + taxon_resolution + (1|basin),
                        family = beta_family(link = "logit"),
                        REML = TRUE,
                        data = dat)
summary(eco.temp.resol)

# Predicted year effect given values of predictor

temp.pred<-predict(eco.temp, newdata = data.frame("year_s" = year.pred, "temp_s" = mean(dat$temp_s), "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
temp.pred.eco<-data.frame(temp.pred)
temp.pred.low<-predict(eco.temp, newdata = data.frame("year_s" = year.pred, "temp_s" = -2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)
temp.pred.high<-predict(eco.temp, newdata = data.frame("year_s" = year.pred, "temp_s" = 2, "basin" = dat$basin[1]), type = "response", se.fit = TRUE)

temp.pred.eco<-cbind(temp.pred.eco, temp.pred.low)
temp.pred.eco<-cbind(temp.pred.eco, temp.pred.high)
colnames(temp.pred.eco)<-c("mean.fit", "mean.se", "fit.low", "se.low", "fit.high", "se.high")
temp.pred.eco<-temp.pred.eco %>% mutate(low.mean = mean.fit - (1.96 * mean.se), up.mean = mean.fit + (1.96 * mean.se),
                                        low.low = fit.low - (1.96 * se.low), up.low = fit.low + (1.96 * se.low),
                                        low.high = fit.high - (1.96 * se.high), up.high = fit.high + (1.96 *se.high))
temp.pred.eco$year<-year.unscaled
temp.pred.eco<-temp.pred.eco %>% relocate(year)

temp.pred.eco %>% ggplot() +
  ylim(0,1) +
  geom_line(aes(y = fit.low, x = year), color = "#2c7bb6", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.low, ymax = up.low), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit.high, x = year), color = "#d7191c", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = low.high, ymax = up.high), alpha = 0.1, fill = "black") +
  labs(y = "β-diversity", x = NULL, title = NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))



### --- Visualize trends per basin and biodiversity facet

trends<-merge(bio.trends, tax.trends, by = c("basin", "year", "index", "fit"), all.x = TRUE, all.y = TRUE)
trends<-merge(trends, eco.trends, by = c("basin", "year", "index", "fit"), all.x = TRUE, all.y = TRUE)

trends %>% ggplot(aes(x = year, y = fit, group = index, color = index)) +
  ylim(0 ,1) +
  geom_line(linewidth = 2) +
  labs(y = NULL, x = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", text = element_text(size = 20)) +
  scale_color_manual(values = c("#e66101", "#5e3c99", "#008837")) +
  facet_wrap(~ basin, ncol = 9) +
  theme(strip.background = element_blank(), strip.text = element_blank())


## -- Trends of beta diversity (trends per basin)

# Extract coefficients from models

tax.slopes<-as.data.frame(coef(summary(tax.year))$cond)
colnames(tax.slopes)<-c("tax_Estimate", "tax_StdError", "tax_zValue", "tax_pValue")
tax.slopes$term<-rownames(tax.slopes)

bio.slopes<-as.data.frame(coef(summary(bio.year))$cond)
colnames(bio.slopes)<-c("bio_Estimate", "bio_StdError", "bio_zValue", "bio_pValue")
bio.slopes$term<-rownames(bio.slopes)

eco.slopes<-as.data.frame(coef(summary(eco.year))$cond)
colnames(eco.slopes)<-c("eco_Estimate", "eco_StdError", "eco_zValue", "eco_pValue")
eco.slopes$term<-rownames(eco.slopes)

beta.slopes<-merge(tax.slopes, bio.slopes, by = c("term"))
beta.slopes<-merge(beta.slopes, eco.slopes, by = c("term"))

beta.slopes<-beta.slopes %>% filter(grepl("year_s", term))
beta.slopes$basin<-basins

beta.slopes<-beta.slopes %>% mutate(tax.low = beta.slopes$tax_Estimate - (beta.slopes$tax_StdError*1.96),
                                    tax.high = beta.slopes$tax_Estimate + (beta.slopes$tax_StdError*1.96),
                                    bio.low = beta.slopes$bio_Estimate - (beta.slopes$bio_StdError*1.96),
                                    bio.high = beta.slopes$bio_Estimate + (beta.slopes$bio_StdError*1.96),
                                    eco.low = beta.slopes$eco_Estimate - (beta.slopes$eco_StdError*1.96),
                                    eco.high = beta.slopes$eco_Estimate + (beta.slopes$eco_StdError*1.96),)

beta.slopes %>% ggplot(aes(x = tax_Estimate, y = reorder(basin, tax_Estimate))) +
  xlim(-2, 2) +
  geom_point(size = 2, color = "#008837") +
  geom_errorbar(aes(xmin = tax.low, xmax = tax.high), color = "#008837") +
  scale_y_discrete(limits = rev) +
  labs(x = "Slope", y = "Basin ID") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  annotate("rect", xmin = c(-2), xmax = c(0), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = c("gray"))

beta.slopes %>% ggplot(aes(x = bio_Estimate, reorder(basin, bio_Estimate))) +
  xlim(-1, 1) +
  geom_point(size = 2, color = "#5e3c99") +
  geom_errorbar(aes(xmin = func.low, xmax = func.high), color = "#5e3c99") +
  scale_y_discrete(limits = rev) +
  labs(x = "Slope", y = "Basin ID") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  annotate("rect", xmin = c(-1), xmax = c(0), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = c("gray"))

beta.slopes %>% ggplot(aes(x = eco_Estimate, reorder(basin, eco_Estimate))) +
  xlim(-1, 1) +
  geom_point(size = 2, color = "#e66101") +
  geom_errorbar(aes(xmin = eco.low, xmax = eco.high), color = "#e66101") +
  scale_y_discrete(limits = rev) +
  labs(x = "Slope", y = "Basin ID") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  annotate("rect", xmin = c(-1), xmax = c(0), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = c("gray"))


### --- Residuals vs basin and data characteristics

dat<-dat %>% mutate(tax_mean_res = residuals(tax.mean, type = c("response")), bio_mean_res = residuals(bio.mean, type = c("response")), eco_mean_res = residuals(eco.mean, type = c("response")))

# 1. Mean distance of sampling sites

dat %>% ggplot(aes(x=log(distance_mean), y = log(tax_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Mean distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

dat %>% ggplot(aes(x=log(distance_mean), y = log(bio_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Mean distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

dat %>% ggplot(aes(x=log(distance_mean), y = log(eco_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Mean distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

# 2. Maximum distance of sampling points

dat %>% ggplot(aes(x=log(distance_max), y = log(tax_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Maximum distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +

dat %>% ggplot(aes(x=log(distance_max), y = log(func_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Maximum distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +

dat %>% ggplot(aes(x=log(distance_max), y = log(eco_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Maximum distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5))

# 3. Number of sampling sites

dat %>% ggplot(aes(x=log(number_sites), y = log(tax_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Number of sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

dat %>% ggplot(aes(x=log(number_sites), y = log(bio_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Number of sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

dat %>% ggplot(aes(x=log(number_sites), y = log(eco_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Number of sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

# 4. Times series length

dat %>% ggplot(aes(x=ts_length, y = log(tax_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Time series length", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

dat %>% ggplot(aes(x=ts_length, y = log(bio_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Time series length", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

dat %>% ggplot(aes(x=ts_length, y = log(eco_mean_res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Time series length", y = "Fitted residuals") +
  theme(text = element_text(size = 15))




### --- Do local trends explain beta diversity trends

## -- Model individual slopes

# - Taxonomic richness

tric.local.mean<-glmmTMB(taxon_richness ~ year_s + (1|fsite), family = poisson(link = "log"), REML = TRUE, data = subset(local.dat, basin == basins[1]))
tric.local.coef<-coef(tric.local.mean)$cond$fsite$year_s[1]
tric.local.coef<-data.frame(tric.local.coef)

for(i in 2:length(basins)){
  tric.local.mean1<-glmmTMB(taxon_richness ~ year_s + (1|fsite), family = poisson(link = "log"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  tric.local.coef1<-coef(tric.local.mean1)$cond$fsite$year_s[1]
  tric.local.coef<-rbind(tric.local.coef, tric.local.coef1)
}

# - Biological traits richness

bio.local.mean<-glmmTMB(biological_richness ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[1]))
bio.local.coef<-coef(bio.local.mean)$cond$fsite$year_s[1]
bio.local.coef<-data.frame(bio.local.coef)

for(i in 2:length(basins)){
  bio.local.mean1<-glmmTMB(biological_richness ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  bio.local.coef1<-coef(bio.local.mean1)$cond$fsite$year_s[1]
  bio.local.coef<-rbind(bio.local.coef, bio.local.coef1)
}

# - Ecological traits richness

eco.local.mean<-glmmTMB(ecological_richness ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[1]))
eco.local.coef<-coef(eco.local.mean)$cond$fsite$year_s[1]
eco.local.coef<-data.frame(eco.local.coef)

for(i in 2:length(basins)){
  eco.local.mean1<-glmmTMB(ecological_richness ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  eco.local.coef1<-coef(eco.local.mean1)$cond$fsite$year_s[1]
  eco.local.coef<-rbind(eco.local.coef, eco.local.coef1)
}

# - Taxonomic diversity

tax.hill.mean<-glmmTMB(taxonomic_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[1]))
tax.hill.coef<-coef(thill.mean)$cond$fsite$year_s[1]
tax.hill.coef<-data.frame(tax.hill.coef)

for(i in 2:length(basins)){
  tax.hill.mean1<-glmmTMB(taxonomic_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  tax.hill.coef1<-coef(tax.hill.mean1)$cond$fsite$year_s[1]
  tax.hill.coef<-rbind(tax.hill.coef, tax.hill.coef1)
}

# - Biological traits diversity

bio.hill.mean<-glmmTMB(biological_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[1]))
bio.hill.coef<-coef(bio.hill.mean)$cond$fsite$year_s[1]
bio.hill.coef<-data.frame(bio.hill.coef)

for(i in 2:length(basins)){
  bio.hill.mean1<-glmmTMB(biological_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  bio.hill.coef1<-coef(bio.hill.mean1)$cond$fsite$year_s[1]
  bio.hill.coef<-rbind(bio.hill.coef, bio.hill.coef1)
}

# - Ecological traits diversity

eco.hill.mean<-glmmTMB(ecological_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[1]))
eco.hill.coef<-coef(eco.hill.mean)$cond$fsite$year_s[1]
eco.hill.coef<-data.frame(eco.hill.coef)

for(i in 2:length(basins)){
  eco.hill.mean1<-glmmTMB(ecological_hill1 ~ year_s + (1|fsite), family = ordbeta(link = "logit"), REML = TRUE, data = subset(local.dat, basin == basins[i]))
  eco.hill.coef1<-coef(eco.hill.mean1)$cond$fsite$year_s[1]
  eco.hill.coef<-rbind(eco.hill.coef, eco.hill.coef1)
}

local.slopes<-data.frame(cbind(tric.local.coef, bio.local.coef, eco.local.coef, tax.hill.coef, bio.hill.coef, eco.hill.coef))
local.slopes$basin<-basins
local.slopes<-merge(local.slopes, beta.slopes[, c(1,3,7,11)], by = c("basin"))

# Linear regression between basin slopes and beta diversity slopes

tax.lm<-lm(tax_Estimate ~ tric.local.coef, data = local.slopes)
summary(tax.lm)
par(mfrow = c(2,2))
plot(tax.lm)

tax.pred<-predict(tax.lm, newdata = data.frame("tric.local.coef" = local.slopes$tric.local.coef), type = "response", se.fit = TRUE)
tax.pred<-data.frame(tax.pred)
colnames(tax.pred)<-c("tax.fit", "tax.se.fit", "tax.df", "tax.Res")
tax.pred$term<-"taxonomic"
tax.pred$tric.local.coef<-local.slopes$tric.local.coef
tax.pred<-tax.pred %>% mutate(ci05 = tax.fit - (1.96* tax.se.fit), ci95 = tax.fit + (1.96* tax.se.fit))


bio.lm<-lm(bio_Estimate ~ bio.local.coef, data = local.slopes)
summary(bio.lm)
par(mfrow = c(2,2))
plot(bio.lm)

bio.pred<-predict(bio.lm, newdata = data.frame("fric.local.coef" = local.slopes$fric.local.coef), type = "response", se.fit = TRUE)
bio.pred<-data.frame(bio.pred)
colnames(bio.pred)<-c("bio.fit", "bio.se.fit", "bio.df", "bio.Res")
bio.pred$term<-"biological"
bio.pred$bio.local.coef<-local.slopes$bio.local.coef
bio.pred<-bio.pred %>% mutate(ci05 = bio.fit - (1.96* bio.se.fit), ci95 = bio.fit + (1.96* bio.se.fit))


eco.lm<-lm(eco_Estimate ~ eco.local.coef, data = local.slopes)
summary(eco.lm)
par(mfrow = c(2,2))
plot(eco.lm)

eco.pred<-predict(eco.lm, newdata = data.frame("fric.local.coef" = local.slopes$fric.local.coef), type = "response", se.fit = TRUE)
eco.pred<-data.frame(eco.pred)
colnames(eco.pred)<-c("eco.fit", "eco.se.fit", "eco.df", "eco.Res")
eco.pred$term<-"ecological"
eco.pred$eco.local.coef<-local.slopes$eco.local.coef
eco.pred<-eco.pred %>% mutate(ci05 = eco.fit - (1.96* eco.se.fit), ci95 = eco.fit + (1.96* eco.se.fit))


ggplot() +
  geom_point(aes(y = tax_Estimate, x = tric.local.coef), data = local.slopes, color = "#008837") +
  geom_line(aes(y = tax.fit, x = tric.local.coef), type = "dashed") +
  geom_ribbon(aes(x = tric.local.coef, ymin = ci05, ymax = ci95), data = tax.pred, fill = "#008837", alpha = 0.1, linetype = 0) +
  geom_point(aes(y = bio_Estimate, x = bio.local.coef), data = local.slopes, color = "#5e3c99") +
  geom_line(aes(y = bio.fit, x = bio.local.coef), data = bio.pred, color = "#5e3c99", linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(x = bio.local.coef, ymin = ci05, ymax = ci95), data = bio.pred, fill = "#5e3c99", alpha = 0.1, linetype = 0) +
  geom_point(aes(y = eco_Estimate, x = eco.local.coef), data = local.slopes, color = "#e66101") +
  geom_line(aes(y = eco.fit, x = eco.local.coef), data = eco.pred, color = "#e66101", linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(x = eco.local.coef, ymin = ci05, ymax = ci95), data = eco.pred, fill = "#e66101", alpha = 0.1, linetype = 0) +
  theme_classic() +
  labs(x = "Average trend in local richness", y = "ꞵ-diversity trend") +
  theme(text = element_text(size = 20))


tax.hill<-lm(tax_Estimate ~ tax.hill.coef, data = local.slopes)
summary(tax.hill)
par(mfrow = c(2,2))
plot(tax.hill)

tax.pred.hill<-predict(tax.hill, newdata = data.frame("tric.hill.coef" = local.slopes$tric.hill.coef), type = "response", se.fit = TRUE)
tax.pred.hill<-data.frame(tax.pred.hill)
colnames(tax.pred.hill)<-c("tax.fit", "tax.se.fit", "tax.df", "tax.Res")
tax.pred.hill$term<-"taxonomic"
tax.pred.hill$tax.hill.coef<-local.slopes$tax.hill.coef
tax.pred.hill<-tax.pred.hill %>% mutate(ci05 = tax.fit - (1.96* tax.se.fit), ci95 = tax.fit + (1.96* tax.se.fit))


bio.hill<-lm(bio_Estimate ~ bio.hill.coef, data = local.slopes)
summary(bio.hill)
par(mfrow = c(2,2))
plot(bio.hill)

bio.pred.hill<-predict(bio.hill, newdata = data.frame("fric.hill.coef" = local.slopes$fric.hill.coef), type = "response", se.fit = TRUE)
bio.pred.hill<-data.frame(bio.pred.hill)
colnames(bio.pred.hill)<-c("bio.fit", "bio.se.fit", "bio.df", "bio.Res")
bio.pred.hill$term<-"biological"
bio.pred.hill$bio.hill.coef<-local.slopes$bio.hill.coef
bio.pred.hill<-bio.pred.hill %>% mutate(ci05 = bio.fit - (1.96* bio.se.fit), ci95 = bio.fit + (1.96* bio.se.fit))


eco.hill<-lm(eco_Estimate ~ eco.hill.coef, data = local.slopes)
summary(eco.hill)
par(mfrow = c(2,2))
plot(eco.hill)

eco.pred.hill<-predict(eco.hill, newdata = data.frame("fric.hill.coef" = local.slopes$fric.hill.coef), type = "response", se.fit = TRUE)
eco.pred.hill<-data.frame(eco.pred.hill)
colnames(eco.pred.hill)<-c("eco.fit", "eco.se.fit", "eco.df", "eco.Res")
eco.pred.hill$term<-"ecological"
eco.pred.hill$eco.hill.coef<-local.slopes$eco.hill.coef
eco.pred.hill<-eco.pred.hill %>% mutate(ci05 = eco.fit - (1.96* eco.se.fit), ci95 = eco.fit + (1.96* eco.se.fit))


ggplot() +
  geom_point(aes(y = tax_Estimate, x = tax.hill.coef), data = local.slopes, color = "#008837") +
  geom_line(aes(y = tax.fit, x = tax.hill.coef), type = "dashed") +
  geom_ribbon(aes(x = tax.hill.coef, ymin = ci05, ymax = ci95), data = tax.pred, fill = "#008837", alpha = 0.1, linetype = 0) +
  geom_point(aes(y = bio_Estimate, x = bio.hill.coef), data = local.slopes, color = "#5e3c99") +
  geom_line(aes(y = bio.fit, x = bio.hill.coef), data = bio.pred, color = "#5e3c99", linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(x = bio.hill.coef, ymin = ci05, ymax = ci95), data = bio.pred, fill = "#5e3c99", alpha = 0.1, linetype = 0) +
  geom_point(aes(y = eco_Estimate, x = eco.hill.coef), data = local.slopes, color = "#e66101") +
  geom_line(aes(y = eco.fit, x = eco.hill.coef), data = eco.pred, color = "#e66101", linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(x = eco.hill.coef, ymin = ci05, ymax = ci95), data = eco.pred, fill = "#e66101", alpha = 0.1, linetype = 0) +
  theme_classic() +
  labs(x = "Average trend in local diversity", y = "ꞵ-diversity trend") +
  theme(text = element_text(size = 20))



### --- Relationship between taxonomic and trait trends

# Chi Square transformation: uplift rare species

chisq.year<-glmmTMB(taxonomic_chiSq ~ basin + basin:year_s,
                    family = beta_family(link = "logit"),
                    REML = TRUE,
                    data = dat)

# Extract coefficients from models

chisq.slopes<-as.data.frame(coef(summary(chisq.year))$cond)
colnames(chisq.slopes)<-c("chisq_Estimate", "chisq_StdError", "chisq_zValue", "chisq_pValue")
chisq.slopes$term<-rownames(chisq.slopes)
chisq.slopes<-chisq.slopes %>% filter(grepl("year_s", term))
beta.slopes<-merge(beta.slopes, chisq.slopes, by = c("term"))

# Taxonomic vs biological traits

tax.bio<-lm(bio_Estimate ~ chisq_Estimate, data = beta.slopes)
summary(tax.bio)

# Test the trend by predictors

beta.slopes<-merge(beta.slopes, predictors, by = c("basin"), all.x = TRUE)

tax.bio.lowurban<-lm(bio_Estimate ~ chisq_Estimate, data = subset(beta.slopes, urban_mean<=0.014))
summary(tax.bio.lowurban)

tax.bio.highurban<-lm(bio_Estimate ~ chisq_Estimate, data = subset(beta.slopes, urban_mean>0.014))
summary(tax.bio.highurban)

tax.bio.loweqr<-lm(bio_Estimate ~ chisq_Estimate, data = subset(beta.slopes, mean_eqr<=0.65))
summary(tax.bio.loweqr)

tax.bio.higheqr<-lm(bio_Estimate ~ chisq_Estimate, data = subset(beta.slopes, mean_eqr>0.65))
summary(tax.bio.higheqr)

tax.bio.lowtemp<-lm(bio_Estimate ~ chisq_Estimate, data = subset(beta.slopes, temperature_trend<=0.03))
summary(tax.bio.lowtemp)

tax.bio.hightemp<-lm(bio_Estimate ~ chisq_Estimate, data = subset(beta.slopes, temperature_trend>0.03))
summary(tax.bio.hightemp)


# Taxonomic vs ecological traits

chisq.eco.lm<-lm(eco_Estimate ~ chisq_Estimate, data = beta.slopes)
summary(chisq.eco.lm)

# Test the trend by predictors

tax.eco.lowurban<-lm(eco_Estimate ~ chisq_Estimate, data = subset(beta.slopes, urban_mean<=0.014))
summary(tax.eco.lowurban)

tax.eco.highurban<-lm(eco_Estimate ~ chisq_Estimate, data = subset(beta.slopes, urban_mean>0.014))
summary(tax.eco.highurban)

tax.eco.loweqr<-lm(eco_Estimate ~ chisq_Estimate, data = subset(beta.slopes, mean_eqr<=0.65))
summary(tax.eco.loweqr)

tax.eco.higheqr<-lm(eco_Estimate ~ chisq_Estimate, data = subset(beta.slopes, mean_eqr>0.65))
summary(tax.eco.higheqr)

tax.eco.lowtemp<-lm(eco_Estimate ~ chisq_Estimate, data = subset(beta.slopes, temperature_trend<=0.03))
summary(tax.eco.lowtemp)

tax.eco.hightemp<-lm(eco_Estimate ~ chisq_Estimate, data = subset(beta.slopes, temperature_trend>0.03))
summary(tax.eco.hightemp)

beta.slopes.long1<-beta.slopes[, c(1,21,7,11)] %>% gather(key = "trait", value = "Estimate", 3:4)

# Slopes all basins

beta.slopes.long1 %>% ggplot(aes(x = chisq_Estimate, y = Estimate, group = trait, color = trait)) +
  xlim(-0.45, 0.3) +
  ylim(-0.5, 0.5) +
  geom_point(size = 2.5) +
  geom_smooth(aes(), method = "lm", linewidth = 2, alpha = 0.2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
  theme_classic() +
  theme(text = element_text(size = 40), legend.position = "none") +
  labs(x = "Trend in taxonomic ꞵ-diversity", y = "Trend in trait ꞵ-diversity") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "grey")


# Slopes by stressor

beta.slopes.long2<-beta.slopes[, c(1,21,7,11,25:27)] %>% gather(key = "trait", value = "Estimate", 3:4)

# Low Urban

beta.slopes.long2 %>% subset(urban_mean<=0.014) %>%
  ggplot(aes(x = chisq_Estimate, y = Estimate, group = trait, color = trait)) +
  xlim(-0.45, 0.3) +
  ylim(-0.5, 0.5) +
  geom_point(size = 2.5) +
  geom_smooth(aes(), method = "lm", linewidth = 2, alpha = 0.2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
  theme_classic() +
  theme(text = element_text(size = 40), legend.position = "none") +
  labs(x = "Trend in taxonomic ꞵ-diversity", y = "Trend in trait ꞵ-diversity") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "grey")

# High Urban

beta.slopes.long2 %>% subset(urban_mean>0.014) %>%
  ggplot(aes(x = chisq_Estimate, y = Estimate, group = trait, color = trait)) +
  xlim(-0.45, 0.3) +
  ylim(-0.5, 0.5) +
  geom_point(size = 2.5) +
  geom_smooth(aes(linetype = trait), method = "lm", linewidth = 2, alpha = 0.2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  theme_classic() +
  theme(text = element_text(size = 40), legend.position = "none") +
  labs(x = "Trend in taxonomic ꞵ-diversity", y = "Trend in trait ꞵ-diversity") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "grey")

# Low EQR

beta.slopes.long2 %>% subset(eqr<=0.65) %>%
  ggplot(aes(x = chisq_Estimate, y = Estimate, group = trait, color = trait)) +
  xlim(-0.45, 0.3) +
  ylim(-0.5, 0.5) +
  geom_point(size = 2.5) +
  geom_smooth(aes(), method = "lm", linewidth = 2, alpha = 0.2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
  theme_classic() +
  theme(text = element_text(size = 40), legend.position = "none") +
  labs(x = "Trend in taxonomic ꞵ-diversity", y = "Trend in trait ꞵ-diversity") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "grey")

# High EQR

beta.slopes.long2 %>% subset(eqr>0.65) %>%
  ggplot(aes(x = chisq_Estimate, y = Estimate, group = trait, color = trait)) +
  xlim(-0.45, 0.3) +
  ylim(-0.5, 0.5) +
  geom_point(size = 2.5) +
  geom_smooth(aes(linetype = trait), method = "lm", linewidth = 2, alpha = 0.2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  theme_classic() +
  theme(text = element_text(size = 40), legend.position = "none") +
  labs(x = "Trend in taxonomic ꞵ-diversity", y = "Trend in trait ꞵ-diversity") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "grey")

# Low temp

beta.slopes.long2 %>% subset(temperature_trend<=0.03) %>%
  ggplot(aes(x = chisq_Estimate, y = Estimate, group = trait, color = trait)) +
  xlim(-0.45, 0.3) +
  ylim(-0.5, 0.5) +
  geom_point(size = 2.5) +
  geom_smooth(aes(), method = "lm", linewidth = 2, alpha = 0.2, linetype = "longdash") +
  scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
  theme_classic() +
  theme(text = element_text(size = 40), legend.position = "none") +
  labs(x = "Trend in taxonomic ꞵ-diversity", y = "Trend in trait ꞵ-diversity") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "grey")

# High temp

beta.slopes.long2 %>% subset(temperature_trend>0.03) %>%
  ggplot(aes(x = chisq_Estimate, y = Estimate, group = trait, color = trait)) +
  xlim(-0.45, 0.3) +
  ylim(-0.5, 0.5) +
  geom_point(size = 2.5) +
  geom_smooth(aes(), method = "lm", linewidth = 2, alpha = 0.2) +
  scale_color_manual(values = c("#d7191c", "#2c7bb6")) +
  theme_classic() +
  theme(text = element_text(size = 40), legend.position = "none") +
  labs(x = "Trend in taxonomic ꞵ-diversity", y = "Trend in trait ꞵ-diversity") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "grey")



#### ---- PARTITION BETA DIVERSITY CHANGE INTO ITS ADDITIVE AND SUBTRACTIVE DYNAMIC COMPONENTS

ecopart<-read.csv(file.choose(), header = TRUE) # Data3_ecopart.csv
ecopart<-merge(ecopart, predictors, by = c("basin"), all.x = TRUE)

allyears<-ecopart[, c(1, 6:14, 33:35)] # year to year change
allyears<-allyears[!is.na(allyears$delta_taxonomic_allYears), ]
tofirst<-ecopart[, c(1, 15:23, 33:35)] # each year to the first year of the time series
tofirst<-tofirst[!is.na(tofirst$delta_taxonomic_toFirst), ]
fromlast<-ecopart[, c(1, 24:32, 33:35)] # last to first year
fromlast<-fromlast[!is.na(fromlast$delta_taxonomic_fromLast), ]

### --- Are additive and subtractive components different from 0? Basins approach

## -- Year to year

tax.add<-glmmTMB(taxa_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = allyears)
tax.coef.add<-data.frame(coef(summary(tax.add))$cond)
tax.coef.add$component<-"additive"
tax.coef.add$index<-"a. taxonomic"

tax.sub<-glmmTMB(taxa_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = allyears)
tax.coef.sub<-data.frame(coef(summary(tax.sub))$cond)
tax.coef.sub$component<-"subtractive"
tax.coef.sub$index<-"a. taxonomic"

bio.add<-glmmTMB(biological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = allyears)
bio.coef.add<-data.frame(coef(summary(bio.add))$cond)
bio.coef.add$component<-"additive"
bio.coef.add$index<-"b. biological"

bio.sub<-glmmTMB(biological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = allyears)
bio.coef.sub<-data.frame(coef(summary(bio.sub))$cond)
bio.coef.sub$component<-"subtractive"
bio.coef.sub$index<-"b. biological"

eco.add<-glmmTMB(ecological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = allyears)
eco.coef.add<-data.frame(coef(summary(eco.add))$cond)
eco.coef.add$component<-"additive"
eco.coef.add$index<-"c. ecological"

eco.sub<-glmmTMB(ecological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = allyears)
eco.coef.sub<-data.frame(coef(summary(eco.sub))$cond)
eco.coef.sub$component<-"subtractive"
eco.coef.sub$index<-"c. ecological"

allyears.coefficients<-rbind(tax.coef.add, tax.coef.sub, bio.coef.add, bio.coef.sub, eco.coef.add, eco.coef.sub)
allyears.coefficients<-allyears.coefficients %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

allyears.coefficients %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)

# - Dynamic components by predictors

# By urban areas

tax.add.lowurban<-glmmTMB(taxa_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean<=0.014))
tax.coef.add.lowurban<-data.frame(coef(summary(tax.add.lowurban))$cond)
tax.coef.add.lowurban$component<-"additive"
tax.coef.add.lowurban$index<-"a. taxonomic"

tax.sub.lowurban<-glmmTMB(taxa_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean<=0.014))
tax.coef.sub.lowurban<-data.frame(coef(summary(tax.sub.lowurban))$cond)
tax.coef.sub.lowurban$component<-"subtractive"
tax.coef.sub.lowurban$index<-"a. taxonomic"


bio.add.lowurban<-glmmTMB(biological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean<=0.014))
bio.coef.add.lowurban<-data.frame(coef(summary(bio.add.lowurban))$cond)
bio.coef.add.lowurban$component<-"additive"
bio.coef.add.lowurban$index<-"b. biological"

bio.sub.lowurban<-glmmTMB(biological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean<=0.014))
bio.coef.sub.lowurban<-data.frame(coef(summary(bio.sub.lowurban))$cond)
bio.coef.sub.lowurban$component<-"subtractive"
bio.coef.sub.lowurban$index<-"b. biological"


eco.add.lowurban<-glmmTMB(ecological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean<=0.014))
eco.coef.add.lowurban<-data.frame(coef(summary(eco.add.lowurban))$cond)
eco.coef.add.lowurban$component<-"additive"
eco.coef.add.lowurban$index<-"c. ecological"

eco.sub.lowurban<-glmmTMB(ecological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears,urban_mean<=0.014))
eco.coef.sub.lowurban<-data.frame(coef(summary(eco.sub.lowurban))$cond)
eco.coef.sub.lowurban$component<-"subtractive"
eco.coef.sub.lowurban$index<-"c. ecological"

allyears.lowurban<-rbind(tax.coef.add.lowurban, tax.coef.sub.lowurban, bio.coef.add.lowurban, bio.coef.sub.lowurban, eco.coef.add.lowurban, eco.coef.sub.lowurban)
allyears.lowurban<-allyears.lowurban %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

allyears.lowurban %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)


tax.add.highurban<-glmmTMB(taxa_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean>0.014))
tax.coef.add.highurban<-data.frame(coef(summary(tax.add.highurban))$cond)
tax.coef.add.highurban$component<-"additive"
tax.coef.add.highurban$index<-"a. taxonomic"

tax.sub.highurban<-glmmTMB(taxa_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean>0.014))
tax.coef.sub.highurban<-data.frame(coef(summary(tax.sub.highurban))$cond)
tax.coef.sub.highurban$component<-"subtractive"
tax.coef.sub.highurban$index<-"a. taxonomic"


bio.add.highurban<-glmmTMB(biological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean>0.014))
bio.coef.add.highurban<-data.frame(coef(summary(bio.add.highurban))$cond)
bio.coef.add.highurban$component<-"additive"
bio.coef.add.highurban$index<-"b. biological"

bio.sub.highurban<-glmmTMB(biological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean>0.014))
bio.coef.sub.highurban<-data.frame(coef(summary(bio.sub.highurban))$cond)
bio.coef.sub.highurban$component<-"subtractive"
bio.coef.sub.highurban$index<-"b. biological"


eco.add.highurban<-glmmTMB(ecological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean>0.014))
eco.coef.add.highurban<-data.frame(coef(summary(eco.add.highurban))$cond)
eco.coef.add.highurban$component<-"additive"
eco.coef.add.highurban$index<-"c. ecological"

eco.sub.highurban<-glmmTMB(ecological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, urban_mean>0.014))
eco.coef.sub.highurban<-data.frame(coef(summary(eco.sub.highurban))$cond)
eco.coef.sub.highurban$component<-"subtractive"
eco.coef.sub.highurban$index<-"c. ecological"

allyears.highurban<-rbind(tax.coef.add.highurban, tax.coef.sub.highurban, bio.coef.add.highurban, bio.coef.sub.highurban, eco.coef.add.highurban, eco.coef.sub.highurban)
allyears.highurban<-allyears.highurban %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

allyears.highurban %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)


# By average EQR

tax.add.loweqr<-glmmTMB(taxa_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr<=0.65))
tax.coef.add.loweqr<-data.frame(coef(summary(tax.add.loweqr))$cond)
tax.coef.add.loweqr$component<-"additive"
tax.coef.add.loweqr$index<-"a. taxonomic"

tax.sub.loweqr<-glmmTMB(taxa_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr<=0.65))
tax.coef.sub.loweqr<-data.frame(coef(summary(tax.sub.loweqr))$cond)
tax.coef.sub.loweqr$component<-"subtractive"
tax.coef.sub.loweqr$index<-"a. taxonomic"


bio.add.loweqr<-glmmTMB(biological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr<=0.65))
bio.coef.add.loweqr<-data.frame(coef(summary(bio.add.loweqr))$cond)
bio.coef.add.loweqr$component<-"additive"
bio.coef.add.loweqr$index<-"b. biological"

bio.sub.loweqr<-glmmTMB(biological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr<=0.65))
bio.coef.sub.loweqr<-data.frame(coef(summary(bio.sub.loweqr))$cond)
bio.coef.sub.loweqr$component<-"subtractive"
bio.coef.sub.loweqr$index<-"b. biological"


eco.add.loweqr<-glmmTMB(ecological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr<=0.65))
eco.coef.add.loweqr<-data.frame(coef(summary(eco.add.loweqr))$cond)
eco.coef.add.loweqr$component<-"additive"
eco.coef.add.loweqr$index<-"c. ecological"

eco.sub.loweqr<-glmmTMB(ecological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr<=0.65))
eco.coef.sub.loweqr<-data.frame(coef(summary(eco.sub.loweqr))$cond)
eco.coef.sub.loweqr$component<-"subtractive"
eco.coef.sub.loweqr$index<-"c. ecological"

allyears.loweqr<-rbind(tax.coef.add.loweqr, tax.coef.sub.loweqr, bio.coef.add.loweqr, bio.coef.sub.loweqr, eco.coef.add.loweqr, eco.coef.sub.loweqr)
allyears.loweqr<-allyears.loweqr %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

allyears.loweqr %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)


tax.add.higheqr<-glmmTMB(taxa_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr>0.65))
tax.coef.add.higheqr<-data.frame(coef(summary(tax.add.higheqr))$cond)
tax.coef.add.higheqr$component<-"additive"
tax.coef.add.higheqr$index<-"a. taxonomic"

tax.sub.higheqr<-glmmTMB(taxa_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr>0.65))
tax.coef.sub.higheqr<-data.frame(coef(summary(tax.sub.higheqr))$cond)
tax.coef.sub.higheqr$component<-"subtractive"
tax.coef.sub.higheqr$index<-"a. taxonomic"


bio.add.higheqr<-glmmTMB(biological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr>0.65))
bio.coef.add.higheqr<-data.frame(coef(summary(bio.add.higheqr))$cond)
bio.coef.add.higheqr$component<-"additive"
bio.coef.add.higheqr$index<-"b. biological"

bio.sub.higheqr<-glmmTMB(biological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr>0.65))
bio.coef.sub.higheqr<-data.frame(coef(summary(bio.sub.higheqr))$cond)
bio.coef.sub.higheqr$component<-"subtractive"
bio.coef.sub.higheqr$index<-"b. biological"


eco.add.higheqr<-glmmTMB(ecological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr>0.65))
eco.coef.add.higheqr<-data.frame(coef(summary(eco.add.higheqr))$cond)
eco.coef.add.higheqr$component<-"additive"
eco.coef.add.higheqr$index<-"c. ecological"

eco.sub.higheqr<-glmmTMB(ecological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, eqr>0.65))
eco.coef.sub.higheqr<-data.frame(coef(summary(eco.sub.higheqr))$cond)
eco.coef.sub.higheqr$component<-"subtractive"
eco.coef.sub.higheqr$index<-"c. ecological"

allyears.higheqr<-rbind(tax.coef.add.higheqr, tax.coef.sub.higheqr, bio.coef.add.higheqr, bio.coef.sub.higheqr, eco.coef.add.higheqr, eco.coef.sub.higheqr)
allyears.higheqr<-allyears.higheqr %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

allyears.higheqr %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)


# By temperature trend

tax.add.lowtemp<-glmmTMB(taxa_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend<=0.03))
tax.coef.add.lowtemp<-data.frame(coef(summary(tax.add.lowtemp))$cond)
tax.coef.add.lowtemp$component<-"additive"
tax.coef.add.lowtemp$index<-"a. taxonomic"

tax.sub.lowtemp<-glmmTMB(taxa_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend<=0.03))
tax.coef.sub.lowtemp<-data.frame(coef(summary(tax.sub.lowtemp))$cond)
tax.coef.sub.lowtemp$component<-"subtractive"
tax.coef.sub.lowtemp$index<-"a. taxonomic"


bio.add.lowtemp<-glmmTMB(biological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend<=0.03))
bio.coef.add.lowtemp<-data.frame(coef(summary(bio.add.lowtemp))$cond)
bio.coef.add.lowtemp$component<-"additive"
bio.coef.add.lowtemp$index<-"b. biological"

bio.sub.lowtemp<-glmmTMB(biological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend<=0.03))
bio.coef.sub.lowtemp<-data.frame(coef(summary(bio.sub.lowtemp))$cond)
bio.coef.sub.lowtemp$component<-"subtractive"
bio.coef.sub.lowtemp$index<-"b. biological"


eco.add.lowtemp<-glmmTMB(ecological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend<=0.03))
eco.coef.add.lowtemp<-data.frame(coef(summary(eco.add.lowtemp))$cond)
eco.coef.add.lowtemp$component<-"additive"
eco.coef.add.lowtemp$index<-"c. ecological"

eco.sub.lowtemp<-glmmTMB(ecological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend<=0.03))
eco.coef.sub.lowtemp<-data.frame(coef(summary(eco.sub.lowtemp))$cond)
eco.coef.sub.lowtemp$component<-"subtractive"
eco.coef.sub.lowtemp$index<-"c. ecological"

allyears.lowtemp<-rbind(tax.coef.add.lowtemp, tax.coef.sub.lowtemp, bio.coef.add.lowtemp, bio.coef.sub.lowtemp, eco.coef.add.lowtemp, eco.coef.sub.lowtemp)
allyears.lowtemp<-allyears.lowtemp %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

allyears.lowtemp %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)


tax.add.hightemp<-glmmTMB(taxa_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend>0.03))
tax.coef.add.hightemp<-data.frame(coef(summary(tax.add.hightemp))$cond)
tax.coef.add.hightemp$component<-"additive"
tax.coef.add.hightemp$index<-"a. taxonomic"

tax.sub.hightemp<-glmmTMB(taxa_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend>0.03))
tax.coef.sub.hightemp<-data.frame(coef(summary(tax.sub.hightemp))$cond)
tax.coef.sub.hightemp$component<-"subtractive"
tax.coef.sub.hightemp$index<-"a. taxonomic"


bio.add.hightemp<-glmmTMB(biological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend>0.03))
bio.coef.add.hightemp<-data.frame(coef(summary(bio.add.hightemp))$cond)
bio.coef.add.hightemp$component<-"additive"
bio.coef.add.hightemp$index<-"b. biological"

bio.sub.hightemp<-glmmTMB(biological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend>0.03))
bio.coef.sub.hightemp<-data.frame(coef(summary(bio.sub.hightemp))$cond)
bio.coef.sub.hightemp$component<-"subtractive"
bio.coef.sub.hightemp$index<-"b. biological"


eco.add.hightemp<-glmmTMB(ecological_addition_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend>0.03))
eco.coef.add.hightemp<-data.frame(coef(summary(eco.add.hightemp))$cond)
eco.coef.add.hightemp$component<-"additive"
eco.coef.add.hightemp$index<-"c. ecological"

eco.sub.hightemp<-glmmTMB(ecological_subtraction_allYears ~ 1 + (1 |basin), family = gaussian, data = subset(allyears, temperature_trend>0.03))
eco.coef.sub.hightemp<-data.frame(coef(summary(eco.sub.hightemp))$cond)
eco.coef.sub.hightemp$component<-"subtractive"
eco.coef.sub.hightemp$index<-"c. ecological"

allyears.hightemp<-rbind(tax.coef.add.hightemp, tax.coef.sub.hightemp, bio.coef.add.hightemp, bio.coef.sub.hightemp, eco.coef.add.hightemp, eco.coef.sub.hightemp)
allyears.hightemp<-allyears.hightemp %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

allyears.hightemp %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)


## -- All years to first

tax.add<-glmmTMB(taxa_addition_toFirst ~ 1 + (1 |basin), family = gaussian, data = tofirst)
tax.coef.add<-data.frame(coef(summary(tax.add))$cond)
tax.coef.add$component<-"additive"
tax.coef.add$index<-"a. taxonomic"

tax.sub<-glmmTMB(taxa_subtraction_toFirst ~ 1 + (1 |basin), family = gaussian, data = tofirst)
tax.coef.sub<-data.frame(coef(summary(tax.sub))$cond)
tax.coef.sub$component<-"subtractive"
tax.coef.sub$index<-"a. taxonomic"

bio.add<-glmmTMB(biological_addition_toFirst ~ 1 + (1 |basin), family = gaussian, data = tofirst)
bio.coef.add<-data.frame(coef(summary(bio.add))$cond)
bio.coef.add$component<-"additive"
bio.coef.add$index<-"b. biological"

bio.sub<-glmmTMB(biological_subtraction_toFirst ~ 1 + (1 |basin), family = gaussian, data = tofirst)
bio.coef.sub<-data.frame(coef(summary(bio.sub))$cond)
bio.coef.sub$component<-"subtractive"
bio.coef.sub$index<-"b. biological"

eco.add<-glmmTMB(ecological_addition_toFirst ~ 1 + (1 |basin), family = gaussian, data = tofirst)
eco.coef.add<-data.frame(coef(summary(eco.add))$cond)
eco.coef.add$component<-"additive"
eco.coef.add$index<-"c. ecological"

eco.sub<-glmmTMB(ecological_subtraction_toFirst ~ 1 + (1 |basin), family = gaussian, data = tofirst)
eco.coef.sub<-data.frame(coef(summary(eco.sub))$cond)
eco.coef.sub$component<-"subtractive"
eco.coef.sub$index<-"c. ecological"

tofirst.coefficients<-rbind(tax.coef.add, tax.coef.sub, bio.coef.add, bio.coef.sub, eco.coef.add, eco.coef.sub)
tofirst.coefficients<-tofirst.coefficients %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

tofirst.coefficients %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)


## -- Last to first year

tax.add<-glmmTMB(taxa_addition_fromLast ~ 1 + (1 |basin), family = gaussian, data = fromlast)
tax.coef.add<-data.frame(coef(summary(tax.add))$cond)
tax.coef.add$component<-"additive"
tax.coef.add$index<-"a. taxonomic"

tax.sub<-glmmTMB(taxa_subtraction_fromLast ~ 1 + (1 |basin), family = gaussian, data = fromlast)
tax.coef.sub<-data.frame(coef(summary(tax.sub))$cond)
tax.coef.sub$component<-"subtractive"
tax.coef.sub$index<-"a. taxonomic"

bio.add<-glmmTMB(biological_addition_fromLast ~ 1 + (1 |basin), family = gaussian, data = fromlast)
bio.coef.add<-data.frame(coef(summary(bio.add))$cond)
bio.coef.add$component<-"additive"
bio.coef.add$index<-"b. biological"

bio.sub<-glmmTMB(biological_subtraction_fromLast ~ 1 + (1 |basin), family = gaussian, data = fromlast)
bio.coef.sub<-data.frame(coef(summary(bio.sub))$cond)
bio.coef.sub$component<-"subtractive"
bio.coef.sub$index<-"b. biological"

eco.add<-glmmTMB(ecological_addition_fromLast ~ 1 + (1 |basin), family = gaussian, data = fromlast)
eco.coef.add<-data.frame(coef(summary(eco.add))$cond)
eco.coef.add$component<-"additive"
eco.coef.add$index<-"c. ecological"

eco.sub<-glmmTMB(ecological_subtraction_fromLast ~ 1 + (1 |basin), family = gaussian, data = fromlast)
eco.coef.sub<-data.frame(coef(summary(eco.sub))$cond)
eco.coef.sub$component<-"subtractive"
eco.coef.sub$index<-"c. ecological"

fromlast.coefficients<-rbind(tax.coef.add, tax.coef.sub, bio.coef.add, bio.coef.sub, eco.coef.add, eco.coef.sub)
fromlast.coefficients<-fromlast.coefficients %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))

fromlast.coefficients %>% ggplot(aes(x = Estimate, y = index, group = component)) +
  geom_point(aes(color = component), size = 5) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = component), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Ecological traits", "Biological traits", "Taxonomic")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "#d7191c")) +
  labs(y = NULL)


### --- Are additive and subtractive components different from 0? Traits approach

## -- Biological traits

biotraits<-read.csv(file.choose(), header = TRUE) # Data4_ecopart_biotraits.csv
bio.traits<-unique(biotraits$biotrait)

bio.add<-glmmTMB(additive ~ 1 + (1|basin), family = gaussian, data = subset(biotraits, biotrait == bio.traits[1]))
bio.add.coef<-data.frame(coef(summary(bio.add))$cond)

for (i in 2:length(bio.traits)){
  bio.add1<-glmmTMB(additive ~ 1 + (1|basin), family = gaussian, data = subset(biotraits, biotrait == bio.traits[i]))
  coef<-data.frame(coef(summary(bio.add1))$cond)
  bio.add.coef<-rbind(bio.add.coef, coef)
}
bio.add.coef$biotrait<-traits
bio.add.coef$component<-"additions"


bio.sub<-glmmTMB(subtractive ~ 1 + (1|basin), family = gaussian, data = subset(biotraits, biotrait == bio.traits[1]))
bio.sub.coef<-data.frame(coef(summary(bio.sub))$cond)

for (i in 2:length(bio.traits)){
  bio.sub1<-glmmTMB(subtractive ~ 1 + (1|basin), family = gaussian, data = subset(biotraits, biotrait == bio.traits[i]))
  coef<-data.frame(coef(summary(bio.sub1))$cond)
  bio.sub.coef<-rbind(bio.sub.coef, coef)
}
bio.sub.coef$biotrait<-traits
bio.sub.coef$component<-"subtractions"

bio.coef<-rbind(bio.add.coef, bio.sub.coef)
bio.coef<-bio.coef %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))
bio.coef$significance<-ifelse(bio.coef$Pr...z..<0.05, bio.coef$component, "not significant")
bio.coef$biotrait<-factor(bio.coef$biotrait,
                          levels = c("tachaqua_egg", "tachaqua_larva", "tachaqua_nymph", "tachaqua_adult",
                                     "tachdisp_aeract", "tachdisp_aerpass", "tachdisp_aquaact", "tachdisp_aquapass",
                                     "tachrepcyc_less1", "tachrepcyc_one", "tachrepcyc_gr1",
                                     "feeding_gra", "feeding_shr", "feeding_gat", "feeding_aff", "feeding_pff",
                                     "feeding_pre", "feeding_min", "feeding_xyl", "feeding_par", "feeding_oth",
                                     "tachresist_coc", "tachresist_did", "tachresist_egg", "tachresist_hou", "tachresist_non",
                                     "tachsize_less0.25", "tachsize_0.25to0.5", "tachsize_0.5to1", "tachsize_1to2",
                                     "tachsize_2to4", "tachsize_4to8", "tachsize_gr8"))


bio.coef %>% ggplot(aes(x = Estimate, y = biotrait, group = significance)) +
  geom_point(aes(color = significance), size = 3) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = significance), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c(">8 cm", "4-8 cm", "2-4 cm", "1-2 cm", "0.5-1 cm", "0.25-0.5 cm", "<0.25 cm",
                                            "None", "Housing", "Dormant eggs", "Diapause", "Cocoons",
                                            "Other", "Parasites", "Xylophagous", "Miners", "Predators",
                                            "Passive filter feeders", "Active filter feeders", "Gatherers", "Shedders", "Grazers",
                                            ">1 per year", "1 per year", "<1 per year",
                                            "Aquatic passive", "Aquatic active", "Aerial passive", "Aerial active",
                                            "Adult", "Nymph", "Larva", "Egg")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 15), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "gray", "#d7191c")) +
  labs(y = NULL, )


## -- Ecological traits

ecotraits<-read.csv(file.choose(), header = TRUE) # Data5_ecopart_ecotraits.csv
eco.traits<-unique(ecotraits$ecotrait)

eco.add<-glmmTMB(additive ~ 1 + (1|basin), family = gaussian, data = subset(ecotraits, ecotrait == eco.traits[1]))
eco.add.coef<-data.frame(coef(summary(eco.add))$cond)

for (i in 2:length(eco.traits)){
  eco.add1<-glmmTMB(additive ~ 1 + (1|basin), family = gaussian, data = subset(ecotraits, ecotrait == eco.traits[i]))
  coef<-data.frame(coef(summary(eco.add1))$cond)
  eco.add.coef<-rbind(eco.add.coef, coef)
}
eco.add.coef$ecotrait<-traits
eco.add.coef$component<-"additions"


eco.sub<-glmmTMB(subtractive ~ 1 + (1|basin), family = gaussian, data = subset(ecotraits, ecotrait == eco.traits[1]))
eco.sub.coef<-data.frame(coef(summary(eco.sub))$cond)

for (i in 2:length(eco.traits)){
  eco.sub1<-glmmTMB(subtractive ~ 1 + (1|basin), family = gaussian, data = subset(ecotraits, ecotrait == eco.traits[i]))
  coef<-data.frame(coef(summary(eco.sub1))$cond)
  eco.sub.coef<-rbind(eco.sub.coef, coef)
}
eco.sub.coef$ecotrait<-eco.traits
eco.sub.coef$component<-"subtractions"

eco.coef<-rbind(eco.add.coef, eco.sub.coef)
eco.coef<-eco.coef %>% mutate(CI05 = Estimate - (Std..Error*1.96), CI95 = Estimate + (Std..Error*1.96))
eco.coef$significance<-ifelse(eco.coef$Pr...z..<0.05, eco.coef$component, "not significant")

eco.coef$ecotrait<-factor(eco.coef$ecotrait,
                          levels = c("tachtemp_eut", "tachtemp_psy", "tachtemp_the",
                                     "microhab_arg", "microhab_pel", "microhab_psa", "microhab_aka", "microhab_lit", "microhab_phy",
                                     "microhab_pom", "microhab_oth",
                                     "tachsap_x", "tachsap_o", "tachsap_b", "tachsap_a", "tachsap_p"))

eco.coef %>% ggplot(aes(x = Estimate, y = ecotrait, group = significance)) +
  xlim(-0.0015, 0.0015) +
  geom_point(aes(color = significance), size = 3) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = significance), linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Polysaprobic", "α-Mesosaprobic", "ꞵ-Mesosaprobic", "Oligosaprobic", "Xenosaprobic",
                                            "Other", "POM", "Phytal", "Lithral", "Sense", "Psammal", "Pelal", "Argyllal",
                                            "Thermophilic", "Psychrophilic", "Eurythermal")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 15), legend.position = "none") +
  scale_color_manual(values = c("#2c7bb6", "gray", "#d7191c")) +
  labs(y = NULL)



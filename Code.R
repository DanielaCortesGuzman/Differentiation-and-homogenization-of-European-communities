######################################################################################################
####################  Homogenization or differentiation of communities  ####################
######################################################################################################

library(ggplot2)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(ordbetareg)
library(performance)
library(ggExtra)
library(GGally)
library(marginaleffects)
library(ggbreak)
library(ggeffects)
library(scales)


## -- ꞵ-diversity data (basin level)
dat<-read.csv(file.choose(), header=TRUE, check.names = FALSE) # ꞵ-diversity by basin and year.csv

## -- Alpha diversity data (site level)
alpha<-read.csv(file.choose(), header=TRUE, check.names = FALSE)  # Diversity indices by site and year.csv

## -- Environmental data (site and basin level)
env<-read.csv(file.choose(), header=TRUE, check.names = FALSE) # Environmental data per site and year.csv


# Standardized initial predictors and define factors

dat<-dat %>% mutate(year_s = scale(year)[,1],
                    basin_f = as.factor(basin))
basins<-unique(dat$basin)

alpha<-alpha %>% mutate(year_s = scale(year)[,1],
                        basin_f = as.factor(basin),
                        site_f = as.factor(site))
sites<-unique(alpha$site)

env<-env %>% mutate(year_s = scale(year)[,1],
                    basin_f = as.factor(basin),
                    site_f = as.factor(site))


### --- Temperature trend per basin

temp.change<-env[!duplicated(env$basin), c("basin", "temp.change")]
dat<-merge(dat, temp.change, by = c("basin"), all.x = TRUE)


### --- Land cover PCA per basin
env$total.lc<-env$tree.full + env$crop.full + env$urb.full
mean(env$total.lc, na.rm = TRUE)
median(env$total.lc, na.rm = TRUE)

land.mean<-env %>% group_by(basin) %>% summarise(mean.tree = mean(tree.full, na.rm = TRUE),
                                                 mean.crop = mean(crop.full, na.rm = TRUE),
                                                 mean.urb = mean(urb.full, na.rm = TRUE), .groups = "drop")

cor(land.mean[, c("mean.tree", "mean.crop", "mean.urb")])
pca<-prcomp(land.mean[, c("mean.tree", "mean.crop", "mean.urb")], scale. = TRUE) # scale true because variances are too different
summary(pca)
pca$rotation # contribution of each type to the axes
# PC1: more forests, less crops and urban; PC2: more urban, less crops

land.pca<-data.frame(pca$x)
land.pca$basin<-basins

dat<-merge(dat, land.pca, by = c("basin"), all.x = TRUE)


### --- Predictors' correlation and collinearity

MyVar <- c("eqr.mean", "temp.change", "PC1", "PC2")
ggpairs(dat[,MyVar])
corvif(dat[,MyVar])

# Keep variables with cor<+- 0.5 and VIF<2



### --- Mean environmental variables per basin

land.mean<-env %>% group_by(basin) %>% summarise(urb.mean = mean(urb.full, na.rm = TRUE),
                                                 crop.mean = mean(crop.full, na.rm = TRUE),
                                                 tree.mean = mean(tree.full, na.rm = TRUE),
                                                 eqr.mean = mean(eqr, na.rm = TRUE), .groups = "drop")

dat<-merge(dat, land.mean, by = c("basin"), all.x = TRUE)



### --- Predictors' distributions

land.mean %>% ggplot(aes(x = eqr.mean)) +
  geom_histogram(color = "#C49EC4", fill = "#C49EC4") +
  labs(x = "Mean EQR", y = "Count") +
  theme_bw() +
  theme(text = element_text(size = 30))

env.temp<-env[!duplicated(env$basin), c("basin", "temp.change")]
env.temp %>% ggplot(aes(x = temp.change)) +
  geom_histogram(color = "#B22222", fill = "#B22222") +
  labs(x = "Temperature slope", y = "Count") +
  theme_bw() +
  theme(text = element_text(size = 30))

land.mean %>% ggplot(aes(x = tree.mean)) +
  geom_histogram(color = "#A1B56C", fill = "#A1B56C") +
  labs(x = "Proportion of forested areas", y = "Count") +
  theme_bw() +
  theme(text = element_text(size = 30))

land.mean %>% ggplot(aes(x = crop.mean)) +
  geom_histogram(color = "#5FB0A9", fill = "#5FB0A9") +
  labs(x = "Proportion of croplands", y = "Count") +
  theme_bw() +
  theme(text = element_text(size = 30))

land.mean %>% ggplot(aes(x = urb.mean)) +
  geom_histogram(color = "#7E587E", fill = "#7E587E") +
  labs(x = "Proportion of urban areas", y = "Count") +
  theme_bw() +
  theme(text = element_text(size = 30))




##### ----- Evidence of recovery

### --- Taxonomic richness trends per site

## -- Mean trend over time

richness.mean<-glmmTMB(tax.richness ~ year_s + (year_s|basin_f) + (1 | site_f),
                       family = nbinom2,
                       REML = TRUE,
                       data = alpha)
summary(richness.mean)

tax.basins<-coef(richness.mean)$cond$basin_f

example.sites<-alpha %>% group_by(basin_f) %>% summarise(site_f = first(site_f))
pred.tax.richness<-expand.grid(year_s = seq(min(alpha$year_s), max(alpha$year_s), length.out = 100), basin_f = unique(alpha$basin_f)) %>%
  left_join(example.sites, by = "basin_f")
scale(alpha$year)
pred.tax.richness$year<-pred.tax.richness$year_s * 5.581855 + 2010.928
pred.tax.richness$fit<-predict(richness.mean, newdata = pred.tax.richness, type = "response", allow.new.levels = TRUE)

tax.fixed<-data.frame(year_s = seq(min(pred.tax.richness$year_s), max(pred.tax.richness$year_s), length.out = 100))
tax.fixed$fit <- predict(richness.mean, newdata = tax.fixed, type = "response", re.form = NA)
tax.fixed$year<-year.unscaled

ggplot() +
  geom_line(data = pred.tax.richness, aes(x = year, y = fit, group = basin_f), linewidth = 1, color = "gray") +
  geom_line(data = tax.fixed, aes(x = year, y = fit), color = "#008837", linewidth = 3, linetype = "solid")+
  labs(y = "Richness", x = "Year") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))



### --- Biological trait richness trends per site

## -- Trends (over time) per basin

bio.richness.mean<-glmmTMB(bio.richness ~ year_s + (year_s|basin_f) + (1 | site_f),
                           family = beta_family(link = "logit"),
                           REML = TRUE,
                           data = alpha)
summary(bio.richness.mean)

bio.basins<-coef(bio.richness.mean)$cond$basin_f

pred.bio.richness<-expand.grid(year_s = seq(min(alpha$year_s), max(alpha$year_s), length.out = 100), basin_f = unique(alpha$basin_f)) %>%
  left_join(example.sites, by = "basin_f")
pred.bio.richness$year<-pred.bio.richness$year_s * 5.581855 + 2010.928
pred.bio.richness$fit<-predict(bio.richness.mean, newdata = pred.bio.richness, type = "response", allow.new.levels = TRUE)

bio.fixed<-data.frame(year_s = seq(min(pred.bio.richness$year_s), max(pred.bio.richness$year_s), length.out = 100))
bio.fixed$fit <- predict(bio.richness.mean, newdata = bio.fixed, type = "response", re.form = NA)
bio.fixed$year<-year.unscaled

ggplot() +
  geom_line(data = pred.bio.richness, aes(x = year, y = fit, group = basin_f), linewidth = 1, color = "gray") +
  geom_line(data = bio.fixed, aes(x = year, y = fit), color = "#5e3c99", linewidth = 3, linetype = "solid")+
  labs(y = "Richness", x = "Year") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))



### --- Ecological trait richness trends per site

## -- Trends (over time) per basin

eco.richness.mean<-glmmTMB(eco.richness ~ year_s + (year_s|basin_f) + (1 | site_f),
                           family = beta_family(link = "logit"),
                           REML = TRUE,
                           data = alpha)
summary(eco.richness.mean)

eco.basins<-coef(eco.richness.mean)$cond$basin_f

pred.eco.richness<-expand.grid(year_s = seq(min(alpha$year_s), max(alpha$year_s), length.out = 100), basin_f = unique(alpha$basin_f)) %>%
  left_join(example.sites, by = "basin_f")
pred.eco.richness$year<-pred.eco.richness$year_s * 5.581855 + 2010.928
pred.eco.richness$fit<-predict(eco.richness.mean, newdata = pred.eco.richness, type = "response", allow.new.levels = TRUE)

eco.fixed<-data.frame(year_s = seq(min(pred.eco.richness$year_s), max(pred.eco.richness$year_s), length.out = 100))
eco.fixed$fit <- predict(eco.richness.mean, newdata = eco.fixed, type = "response", re.form = NA)
eco.fixed$year<-year.unscaled

ggplot() +
  geom_line(data = pred.eco.richness, aes(x = year, y = fit, group = basin_f), linewidth = 1, color = "gray") +
  geom_line(data = eco.fixed, aes(x = year, y = fit), color = "#e66101", linewidth = 3, linetype = "solid")+
  labs(y = "Richness", x = "Year") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


df1<-tax.basins %>% rename(tax.slopes = year_s, tax.intercept = "(Intercept)") %>% mutate(basin = rownames(tax.basins))
df2<-bio.basins %>% rename(bio.slopes = year_s, bio.intercept = "(Intercept)") %>% mutate(basin = rownames(bio.basins))
df3<-eco.basins %>% rename(eco.slopes = year_s, eco.intercept = "(Intercept)") %>% mutate(basin = rownames(eco.basins))

basin.summary<-reduce(list(df1, df2, df3), full_join, by = "basin")

basin.summary<-merge(basin.summary, land.mean[, c("basin", "eqr.mean")], by = c("basin"), all.x = TRUE)
basin.summary<-merge(basin.summary, temp.change, by = c("basin"), all.x = TRUE)
basin.summary<-merge(basin.summary, land.pca[, c("basin", "PC1", "PC2")], by = c("basin"), all.x = TRUE)

basin.summary.long<-basin.summary %>% pivot_longer(cols = ends_with("slopes"),
                                                   names_to = "variable",
                                                   values_to = "slopes")

basin.summary.long %>% ggplot(aes(y = slopes, x = eqr.mean, color = variable)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("#5e3c99", "#e66101", "#008837")) +
  labs(y = "Richness trend", x = "Mean EQR") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  scale_x_reverse()

basin.summary.long %>% ggplot(aes(y = slopes, x = temp.change, color = variable)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("#5e3c99", "#e66101", "#008837")) +
  labs(y = "Richness trend", x = "Temperature trend") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

basin.summary.long %>% ggplot(aes(y = slopes, x = PC2, color = variable)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("#5e3c99", "#e66101", "#008837")) +
  labs(y = "Richness trend", x = "Land cover") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))



##### -----  Modeling ꞵ-diversity trends

# Scale predictors

dat<-dat %>% mutate(eqr.mean = scale(eqr.mean)[,1],
                    temp.change = scale(temp.change)[,1]) # PC1 and PC2 don't need to be scaled again


# Obtain initial ꞵ-diversity values

beta.i<-dat %>% group_by(basin) %>% arrange(year, .by_group = TRUE) %>% slice(1) %>%  # Take the first row per group
  select(basin, year, tax_bray, bio_bray, eco_bray)
colnames(beta.i)[c(3:5)]<-c("taxi", "bioi", "ecoi")

dat<-merge(dat, beta.i[, c("basin", "taxi", "bioi", "ecoi")], by = c("basin"), all.x = TRUE)
dat<-dat %>% mutate(taxi = scale(taxi)[,1], bioi = scale(bioi)[,1], ecoi = scale(ecoi)[,1],)



### --- Taxonomic composition

## -- Mean trend

tax.mean<-glmmTMB(tax_bray ~ year_s + taxi + (1|basin_f),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(tax.mean)

Esqr <- simulateResiduals(fittedModel = tax.mean, plot = FALSE)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Esqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = FALSE)
plotResiduals(Esqr, quantreg = TRUE, smoothScatter = FALSE)
plotResiduals(Esqr, form = dat$year_s, xlab = "Rank-transformed X")

# Predicted values

pred.mean.tax<-predict(tax.mean, newdata = data.frame("year_s" = year.pred, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
pred.mean.tax<-data.frame(pred.mean.tax)
pred.mean.tax$year<-year.unscaled
pred.mean.tax<-pred.mean.tax %>% relocate(year)
pred.mean.tax$index<-"taxonomic"


## -- Trends (over time) in homogenization per basin

tax.basin<-glmmTMB(tax_bray ~ basin_f + basin_f:year_s,
                   family = beta_family(link = "logit"),
                   REML = TRUE,
                   data = dat)
summary(tax.basin)

# Predicted values

y.pred.be01<-predict(tax.basin, newdata = data.frame("year_s" = year.pred, "basin_f" = "BE01"), type = "response", se.fit = TRUE)
tax.trends<-data.frame(y.pred.be01$fit)
tax.trends$year<-year.unscaled
tax.trends<-tax.trends %>% relocate(year)

for (i in 2:length(basins)){
  y.pred<-predict(tax.basin, newdata = data.frame("year_s" = year.pred, "basin_f" = basins[i]), type = "response", se.fit = TRUE)
  tax.trends<-cbind(tax.trends, y.pred$fit)
}
colnames(tax.trends)[2:ncol(tax.trends)]<-basins
tax.trends<-tax.trends %>% gather(key = "basin", value = "fit", 2:ncol(tax.trends))
tax.trends$index<-"taxonomic"


## -- Visualize mean trend and trends per basin

colnames(pred.mean.tax)<-c("year", "fit.mean", "se.mean", "index")

ggplot() +
  ylim (0 ,1) +
  geom_line(data = tax.trends, aes(x = year, y = fit, group = basin), linewidth = 1, color = "gray") +
  geom_line(data = pred.mean.tax, aes(x = year, y = fit.mean), color = "#008837", linewidth = 3, linetype = "longdash")+
  labs(y = "β-diversity", x = "Year") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


## -- Trends in ꞵ-diversity linked to mean EQR

tax.eqr.mean<-glmmTMB(tax_bray ~ year_s * eqr.mean + taxi + (1|basin_f),
                        family = beta_family(link = "logit"),
                        REML = TRUE,
                        data = dat)
summary(tax.eqr.mean)
performance(tax.eqr.mean)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = tax.eqr.mean, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

eqr.mean.pred<-predict(tax.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = mean(dat$eqr.mean), "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean.pred.tax<-data.frame(eqr.mean.pred)
seq(from = min(dat$eqr.mean), to = max(dat$eqr.mean), length.out = 5)
eqr.mean01<-predict(tax.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -1.61, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean02<-predict(tax.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -0.67, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean03<-predict(tax.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 0.27, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean04<-predict(tax.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 1.20, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean05<-predict(tax.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 2.13, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)

eqr.mean.pred.tax<-cbind(eqr.mean.pred.tax, eqr.mean01, eqr.mean02, eqr.mean03, eqr.mean04, eqr.mean05)
colnames(eqr.mean.pred.tax)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
eqr.mean.pred.tax<-eqr.mean.pred.tax %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                    fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # eqr 0.2
                                                    fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # eqr 0.4
                                                    fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # eqr 0.6
                                                    fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # eqr 0.8
                                                    fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # eqr 1.0
eqr.mean.pred.tax$year<-year.unscaled
eqr.mean.pred.tax<-eqr.mean.pred.tax %>% relocate(year)

eqr.mean.pred.tax %>% ggplot() +
  ylim(0.7,0.9) +
  geom_line(aes(y = fit.1, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#B22222") + # more negative slope (degradation)
  geom_line(aes(y = fit.2, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.5, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#00468B") + # more positive slope (improvement)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to mean EQR  (recovering communities)

tax.deg<-c("CR02", "DE06", "DE07", "DE09", "FI01", "FR03", "FR04", "GR01", "GR02", "HU01", "IR01", "SP07", "SW01", "SW02", "SW05")

tax.eqr.rec<-glmmTMB(tax_bray ~ year_s * eqr.mean + taxi + (1|basin_f),
                     family = beta_family(link = "logit"),
                     REML = TRUE,
                     data = dat[!(dat$basin_f %in% tax.deg), ])
summary(tax.eqr.rec)

# Predicted year effect given values of predictor

eqr.rec.pred<-predict(tax.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = mean(dat$eqr.mean), "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.rec.pred.tax<-data.frame(eqr.rec.pred)
seq(from = min(dat$eqr.mean), to = max(dat$eqr.mean), length.out = 5)
eqr.rec01<-predict(tax.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -1.61, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.rec02<-predict(tax.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -0.67, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.rec03<-predict(tax.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 0.27, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.rec04<-predict(tax.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 1.20, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.rec05<-predict(tax.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 2.13, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)

eqr.rec.pred.tax<-cbind(eqr.rec.pred.tax, eqr.rec01, eqr.rec02, eqr.rec03, eqr.rec04, eqr.rec05)
colnames(eqr.rec.pred.tax)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
eqr.rec.pred.tax<-eqr.rec.pred.tax %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                              fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # eqr 0.2
                                              fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # eqr 0.4
                                              fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # eqr 0.6
                                              fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # eqr 0.8
                                              fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # eqr 1.0
eqr.rec.pred.tax$year<-year.unscaled
eqr.rec.pred.tax<-eqr.rec.pred.tax %>% relocate(year)

eqr.rec.pred.tax %>% ggplot() +
  ylim(0.7,0.9) +
  geom_line(aes(y = fit.1, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#B22222") + # more negative slope (degradation)
  geom_line(aes(y = fit.2, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.5, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#00468B") + # more positive slope (improvement)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to temperature trend

tax.temp.change<-glmmTMB(tax_bray ~ year_s * temp.change + taxi + (1|basin_f),
                         family = beta_family(link = "logit"),
                         REML = TRUE,
                         data = dat)
summary(tax.temp.change)
performance(tax.temp.change)

# Test temp.changeoral autocorrelation

Esqr <- simulateResiduals(fittedModel = tax.temp.change, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

temp.change.pred.tax<-predict(tax.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = mean(dat$temp.change), "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change.pred.tax<-data.frame(temp.change.pred.tax)
seq(from = min(dat$temp.change), to = max(dat$temp.change), length.out = 5)
temp.change01<-predict(tax.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = -2.02, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change02<-predict(tax.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = -0.96, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change03<-predict(tax.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 0.09, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change04<-predict(tax.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 1.14, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change05<-predict(tax.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 2.19, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)

temp.change.pred.tax<-cbind(temp.change.pred.tax, temp.change01, temp.change02, temp.change03, temp.change04, temp.change05)
colnames(temp.change.pred.tax)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
temp.change.pred.tax<-temp.change.pred.tax %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                        fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # temp.change 0.2
                                        fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # temp.change 0.4
                                        fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # temp.change 0.6
                                        fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # temp.change 0.8
                                        fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # temp.change 1.0
temp.change.pred.tax$year<-year.unscaled
temp.change.pred.tax<-temp.change.pred.tax %>% relocate(year)

temp.change.pred.tax %>% ggplot() +
  ylim(0.7,0.9) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower slope (less warming or even cooling)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher slope (more warming)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to mean temperature trend (recovering communities)

tax.temp.rec<-glmmTMB(tax_bray ~ year_s * temp.change + taxi + (1|basin_f),
                      family = beta_family(link = "logit"),
                      REML = TRUE,
                      data = dat[!(dat$basin_f %in% tax.deg), ])
summary(tax.temp.rec)

# Predicted year effect given values of predictor

temp.rec.pred.tax<-predict(tax.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = mean(dat$temp.change), "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.rec.pred.tax<-data.frame(temp.rec.pred.tax)
seq(from = min(dat$temp.change), to = max(dat$temp.change), length.out = 5)
temp.change01<-predict(tax.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = -2.02, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change02<-predict(tax.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = -0.96, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change03<-predict(tax.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 0.09, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change04<-predict(tax.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 1.14, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change05<-predict(tax.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 2.19, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)

temp.rec.pred.tax<-cbind(temp.rec.pred.tax, temp.change01, temp.change02, temp.change03, temp.change04, temp.change05)
colnames(temp.rec.pred.tax)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
temp.rec.pred.tax<-temp.rec.pred.tax %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # temp.change 0.2
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # temp.change 0.4
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # temp.change 0.6
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # temp.change 0.8
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # temp.change 1.0
temp.rec.pred.tax$year<-year.unscaled
temp.rec.pred.tax<-temp.rec.pred.tax %>% relocate(year)

temp.rec.pred.tax %>% ggplot() +
  ylim(0.7,0.9) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower slope (less warming or even cooling)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher slope (more warming)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to land cover stress: PC1 is higher forest and lower urban and crops; PC2 is higher urbanization

tax.land<-glmmTMB(tax_bray ~ year_s * PC1 + year_s * PC2 + taxi + (1|basin_f),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(tax.land)
performance(tax.land)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = tax.land, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

land.pred.tax<-predict(tax.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = mean(dat$PC2), "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
land.pred.tax<-data.frame(land.pred.tax)
seq(from = min(dat$PC1), to = max(dat$PC1), length.out = 5)
seq(from = min(dat$PC2), to = max(dat$PC2), length.out = 5)
# Predict for a gradient of urbanization (PC2)
land01<-predict(tax.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -1.44, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
land02<-predict(tax.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -0.28, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
land03<-predict(tax.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 0.88, "taxi" = mean(dat$taxi),  "basin_f" = NA), type = "response", se.fit = TRUE)
land04<-predict(tax.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 2.03, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
land05<-predict(tax.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 3.18, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)

land.pred.tax<-cbind(land.pred.tax, land01, land02, land03, land04, land05)
colnames(land.pred.tax)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
land.pred.tax<-land.pred.tax %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                        fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # higher forest
                                        fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2),
                                        fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3),
                                        fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4),
                                        fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # higher urban
land.pred.tax$year<-year.unscaled
land.pred.tax<-land.pred.tax %>% relocate(year)

land.pred.tax %>% ggplot() +
  ylim(0.7,0.9) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") + 
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower stress (more forest)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") + 
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") + 
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") + 
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") + 
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") + 
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") + 
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher stress (more urban)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to land cover stress (Recovering communities)

tax.land.rec<-glmmTMB(tax_bray ~ year_s * PC1 + year_s * PC2 + taxi + (1|basin_f),
                      family = beta_family(link = "logit"),
                      REML = TRUE,
                      data = dat[!(dat$basin_f %in% tax.deg), ])
summary(tax.land.rec)

# Predicted year effect given values of predictor

land.rec.pred.tax<-predict(tax.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = mean(dat$PC2), "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
land.rec.pred.tax<-data.frame(land.rec.pred.tax)
seq(from = min(dat$PC1), to = max(dat$PC1), length.out = 5)
seq(from = min(dat$PC2), to = max(dat$PC2), length.out = 5)
# Predict for a gradient of urbanization (PC2)
land01<-predict(tax.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -1.44, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
land02<-predict(tax.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -0.28, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
land03<-predict(tax.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 0.88, "taxi" = mean(dat$taxi),  "basin_f" = NA), type = "response", se.fit = TRUE)
land04<-predict(tax.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 2.03, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)
land05<-predict(tax.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 3.18, "taxi" = mean(dat$taxi), "basin_f" = NA), type = "response", se.fit = TRUE)

land.rec.pred.tax<-cbind(land.rec.pred.tax, land01, land02, land03, land04, land05)
colnames(land.rec.pred.tax)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
land.rec.pred.tax<-land.rec.pred.tax %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # higher forest
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2),
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3),
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4),
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # higher urban
land.rec.pred.tax$year<-year.unscaled
land.rec.pred.tax<-land.rec.pred.tax %>% relocate(year)

land.rec.pred.tax %>% ggplot() +
  ylim(0.7,0.9) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") + 
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower stress (more forest)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") + 
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") + 
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") + 
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") + 
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") + 
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") + 
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher stress (more urban)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



### --- Biological trait composition

## -- Mean trend

bio.mean<-glmmTMB(bio_bray ~ year_s + bioi + (1|basin_f),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(bio.mean)
Esqr <- simulateResiduals(fittedModel = bio.mean, plot = FALSE)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Esqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = FALSE)
plotResiduals(Esqr, quantreg = TRUE, smoothScatter = FALSE)
plotResiduals(Esqr, form = dat$year_s, xlab = "Rank-transformed X")

# Predicted values

pred.mean.bio<-predict(bio.mean, newdata = data.frame("year_s" = year.pred, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
pred.mean.bio<-data.frame(pred.mean.bio)
pred.mean.bio$year<-year.unscaled
pred.mean.bio<-pred.mean.bio %>% relocate(year)
pred.mean.bio$index<-"biological"


## -- Trends (over time) in homogenization per basin

bio.basin<-glmmTMB(bio_bray ~ basin_f + basin_f:year_s,
                   family = beta_family(link = "logit"),
                   REML = TRUE,
                   data = dat)
summary(bio.basin)

# Predicted values

y.pred.be01<-predict(bio.basin, newdata = data.frame("year_s" = year.pred, "basin_f" = "BE01"), type = "response", se.fit = TRUE)
bio.trends<-data.frame(y.pred.be01$fit)
bio.trends$year<-year.unscaled
bio.trends<-bio.trends %>% relocate(year)

for (i in 2:length(basins)){
  y.pred<-predict(bio.basin, newdata = data.frame("year_s" = year.pred, "basin_f" = basins[i]), type = "response", se.fit = TRUE)
  bio.trends<-cbind(bio.trends, y.pred$fit)
}
colnames(bio.trends)[2:ncol(bio.trends)]<-basins
bio.trends<-bio.trends %>% gather(key = "basin", value = "fit", 2:ncol(bio.trends))
bio.trends$index<-"biological"


## -- Visualize mean trend and trends per basin

colnames(pred.mean.bio)<-c("year", "fit.mean", "se.mean", "index")

ggplot() +
  ylim (0 ,1) +
  geom_line(data = bio.trends, aes(x = year, y = fit, group = basin), linewidth = 1, color = "gray") +
  geom_line(data = pred.mean.bio, aes(x = year, y = fit.mean), color = "#5e3c99", linewidth = 3, linetype = "longdash")+
  labs(y = "β-diversity", x = "Year") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))



## -- Trends in ꞵ-diversity linked to mean EQR

bio.eqr.mean<-glmmTMB(bio_bray ~ year_s * eqr.mean + bioi + (1|basin_f),
                      family = beta_family(link = "logit"),
                      REML = TRUE,
                      data = dat)
summary(bio.eqr.mean)
performance(bio.eqr.mean)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = bio.eqr.mean, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

eqr.mean.pred<-predict(bio.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = mean(dat$eqr.mean), "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean.pred.bio<-data.frame(eqr.mean.pred)
eqr.mean01<-predict(bio.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -1.61, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean02<-predict(bio.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -0.67, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean03<-predict(bio.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 0.27, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean04<-predict(bio.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 1.20, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean05<-predict(bio.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 2.13, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)

eqr.mean.pred.bio<-cbind(eqr.mean.pred.bio, eqr.mean01, eqr.mean02, eqr.mean03, eqr.mean04, eqr.mean05)
colnames(eqr.mean.pred.bio)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
eqr.mean.pred.bio<-eqr.mean.pred.bio %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # eqr 0.2
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # eqr 0.4
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # eqr 0.6
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # eqr 0.8
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # eqr 1.0
eqr.mean.pred.bio$year<-year.unscaled
eqr.mean.pred.bio<-eqr.mean.pred.bio %>% relocate(year)

eqr.mean.pred.bio %>% ggplot() +
  ylim(0.3,0.6) +
  geom_line(aes(y = fit.1, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#B22222") + # lower EQR
  geom_line(aes(y = fit.2, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.5, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#00468B") + # higher EQR
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to mean EQR (Recovering communities)

bio.deg<-c("DE06", "DE07", "DE08", "FI01", "FR04", "FR06", "GR01", "GR02", "HU01", "IR01", "SP07")

bio.eqr.rec<-glmmTMB(bio_bray ~ year_s * eqr.mean + bioi + (1|basin_f),
                     family = beta_family(link = "logit"),
                     REML = TRUE,
                     data = dat[!(dat$basin_f %in% bio.deg), ])
summary(bio.eqr.rec)

# Predicted year effect given values of predictor

eqr.rec.pred<-predict(bio.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = mean(dat$eqr.mean), "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.rec.pred.bio<-data.frame(eqr.rec.pred)
eqr.mean01<-predict(bio.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -1.61, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean02<-predict(bio.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -0.67, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean03<-predict(bio.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 0.27, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean04<-predict(bio.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 1.20, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean05<-predict(bio.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 2.13, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)

eqr.rec.pred.bio<-cbind(eqr.rec.pred.bio, eqr.mean01, eqr.mean02, eqr.mean03, eqr.mean04, eqr.mean05)
colnames(eqr.rec.pred.bio)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
eqr.rec.pred.bio<-eqr.rec.pred.bio %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                              fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # eqr 0.2
                                              fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # eqr 0.4
                                              fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # eqr 0.6
                                              fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # eqr 0.8
                                              fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # eqr 1.0
eqr.rec.pred.bio$year<-year.unscaled
eqr.rec.pred.bio<-eqr.rec.pred.bio %>% relocate(year)

eqr.rec.pred.bio %>% ggplot() +
  ylim(0.3,0.6) +
  geom_line(aes(y = fit.1, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#B22222") + # lower EQR
  geom_line(aes(y = fit.2, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.5, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#00468B") + # higher EQR
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to temperature trend

bio.temp.change<-glmmTMB(bio_bray ~ year_s * temp.change + bioi + (1|basin_f),
                         family = beta_family(link = "logit"),
                         REML = TRUE,
                         data = dat)
summary(bio.temp.change)
performance(bio.temp.change)

# Test temp.changeoral autocorrelation

Esqr <- simulateResiduals(fittedModel = bio.temp.change, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

temp.change.pred.bio<-predict(bio.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = mean(dat$temp.change), "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change.pred.bio<-data.frame(temp.change.pred.bio)
temp.change01<-predict(bio.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = -2.02, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change02<-predict(bio.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = -0.96, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change03<-predict(bio.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 0.09, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change04<-predict(bio.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 1.14, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change05<-predict(bio.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 2.19, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)

temp.change.pred.bio<-cbind(temp.change.pred.bio, temp.change01, temp.change02, temp.change03, temp.change04, temp.change05)
colnames(temp.change.pred.bio)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
temp.change.pred.bio<-temp.change.pred.bio %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                      fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # temp.change 0.2
                                                      fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # temp.change 0.4
                                                      fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # temp.change 0.6
                                                      fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # temp.change 0.8
                                                      fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # temp.change 1.0
temp.change.pred.bio$year<-year.unscaled
temp.change.pred.bio<-temp.change.pred.bio %>% relocate(year)

temp.change.pred.bio %>% ggplot() +
  ylim(0.3,0.6) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower slope (less warming or even cooling)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher slope (more warming)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to mean temperature trend (Recovering communities)

bio.temp.rec<-glmmTMB(bio_bray ~ year_s * temp.change + bioi + (1|basin_f),
                      family = beta_family(link = "logit"),
                      REML = TRUE,
                      data = dat[!(dat$basin_f %in% bio.deg), ])
summary(bio.temp.rec)

# Predicted year effect given values of predictor

temp.rec.pred.bio<-predict(bio.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = mean(dat$temp.change), "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.rec.pred.bio<-data.frame(temp.rec.pred.bio)
temp.change01<-predict(bio.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = -2.02, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change02<-predict(bio.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = -0.96, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change03<-predict(bio.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 0.09, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change04<-predict(bio.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 1.14, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change05<-predict(bio.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 2.19, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)

temp.rec.pred.bio<-cbind(temp.rec.pred.bio, temp.change01, temp.change02, temp.change03, temp.change04, temp.change05)
colnames(temp.rec.pred.bio)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
temp.rec.pred.bio<-temp.rec.pred.bio %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # temp.change 0.2
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # temp.change 0.4
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # temp.change 0.6
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # temp.change 0.8
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # temp.change 1.0
temp.rec.pred.bio$year<-year.unscaled
temp.rec.pred.bio<-temp.rec.pred.bio %>% relocate(year)

temp.rec.pred.bio %>% ggplot() +
  ylim(0.3,0.6) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower slope (less warming or even cooling)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher slope (more warming)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to land cover stress

bio.land<-glmmTMB(bio_bray ~ year_s * PC1 + year_s * PC2 + bioi + (1|basin_f),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(bio.land)
performance(bio.land)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = bio.land, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

land.pred.bio<-predict(bio.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = mean(dat$PC2), "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land.pred.bio<-data.frame(land.pred.bio)
land01<-predict(bio.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -1.44, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land02<-predict(bio.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -0.28, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land03<-predict(bio.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 0.88, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land04<-predict(bio.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 2.03, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land05<-predict(bio.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 3.18, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)

land.pred.bio<-cbind(land.pred.bio, land01, land02, land03, land04, land05)
colnames(land.pred.bio)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
land.pred.bio<-land.pred.bio %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                        fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # higher forest
                                        fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2),
                                        fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3),
                                        fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4),
                                        fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # higher urban
land.pred.bio$year<-year.unscaled
land.pred.bio<-land.pred.bio %>% relocate(year)

land.pred.bio %>% ggplot() +
  ylim(0.3,0.6) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower stress (more forest)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") + 
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") + 
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") + 
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher stress (more urban)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to land cover stress (Recovering communities)

bio.land.rec<-glmmTMB(bio_bray ~ year_s * PC1 + year_s * PC2 + bioi + (1|basin_f),
                      family = beta_family(link = "logit"),
                      REML = TRUE,
                      data = dat[!(dat$basin_f %in% bio.deg), ])
summary(bio.land.rec)

# Predicted year effect given values of predictor

land.rec.pred.bio<-predict(bio.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = mean(dat$PC2), "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land.rec.pred.bio<-data.frame(land.rec.pred.bio)
land01<-predict(bio.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -1.44, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land02<-predict(bio.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -0.28, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land03<-predict(bio.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 0.88, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land04<-predict(bio.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 2.03, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)
land05<-predict(bio.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 3.18, "bioi" = mean(dat$bioi), "basin_f" = NA), type = "response", se.fit = TRUE)

land.rec.pred.bio<-cbind(land.rec.pred.bio, land01, land02, land03, land04, land05)
colnames(land.rec.pred.bio)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
land.rec.pred.bio<-land.rec.pred.bio %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # higher forest
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2),
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3),
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4),
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # higher urban
land.rec.pred.bio$year<-year.unscaled
land.rec.pred.bio<-land.rec.pred.bio %>% relocate(year)

land.rec.pred.bio %>% ggplot() +
  ylim(0.3,0.6) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower stress (more forest)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") + 
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") + 
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") + 
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher stress (more urban)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



#### --- Ecological composition

## -- Mean trend

eco.mean<-glmmTMB(eco_bray ~ year_s + ecoi + (1|basin_f),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(eco.mean)
Esqr <- simulateResiduals(fittedModel = eco.mean, plot = FALSE)
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Esqr, testUniformity = TRUE, 
           testOutliers = TRUE, testDispersion = FALSE)
plotResiduals(Esqr, quantreg = TRUE, smoothScatter = FALSE)
plotResiduals(Esqr, form = dat$year_s, xlab = "Rank-transformed X")

# Predicted values

pred.mean.eco<-predict(eco.mean, newdata = data.frame("year_s" = year.pred, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
pred.mean.eco<-data.frame(pred.mean.eco)
pred.mean.eco$year<-year.unscaled
pred.mean.eco<-pred.mean.eco %>% relocate(year)
pred.mean.eco$index<-"ecological"


## -- Trends (over time) in homogenization per basin

eco.basin<-glmmTMB(eco_bray ~ basin_f + basin_f:year_s,
                   family = beta_family(link = "logit"),
                   REML = TRUE,
                   data = dat)
summary(eco.basin)

# Predicted values

y.pred.be01<-predict(eco.basin, newdata = data.frame("year_s" = year.pred, "basin_f" = "BE01"), type = "response", se.fit = TRUE)
eco.trends<-data.frame(y.pred.be01$fit)
eco.trends$year<-year.unscaled
eco.trends<-eco.trends %>% relocate(year)

for (i in 2:length(basins)){
  y.pred<-predict(eco.basin, newdata = data.frame("year_s" = year.pred, "basin_f" = basins[i]), type = "response", se.fit = TRUE)
  eco.trends<-cbind(eco.trends, y.pred$fit)
}
colnames(eco.trends)[2:ncol(eco.trends)]<-basins
eco.trends<-eco.trends %>% gather(key = "basin", value = "fit", 2:ncol(eco.trends))
eco.trends$index<-"ecological"


## -- Visualize mean trend and trends per basin

colnames(pred.mean.eco)<-c("year", "fit.mean", "se.mean", "index")

ggplot() +
  ylim(0 ,1) +
  geom_line(data = eco.trends, aes(x = year, y = fit, group = basin), linewidth = 1, color = "gray") +
  geom_line(data = pred.mean.eco, aes(x = year, y = fit.mean), color = "#e66101", linewidth = 3, linetype = "solid") +
  labs(y = "β-diversity", x = "Year") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))



## -- Trends in ꞵ-diversity linked to mean EQR

eco.eqr.mean<-glmmTMB(eco_bray ~ year_s * eqr.mean + ecoi + (1|basin_f),
                      family = beta_family(link = "logit"),
                      REML = TRUE,
                      data = dat)
summary(eco.eqr.mean)
performance(eco.eqr.mean)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = eco.eqr.mean, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

eqr.mean.pred<-predict(eco.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = mean(dat$eqr.mean), "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean.pred.eco<-data.frame(eqr.mean.pred)
eqr.mean01<-predict(eco.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -1.61, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean02<-predict(eco.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -0.67, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean03<-predict(eco.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 0.27, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean04<-predict(eco.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 1.20, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean05<-predict(eco.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 2.13, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)

eqr.mean.pred.eco<-cbind(eqr.mean.pred.eco, eqr.mean01, eqr.mean02, eqr.mean03, eqr.mean04, eqr.mean05)
colnames(eqr.mean.pred.eco)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
eqr.mean.pred.eco<-eqr.mean.pred.eco %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # eqr 0.2
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # eqr 0.4
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # eqr 0.6
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # eqr 0.8
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # eqr 1.0
eqr.mean.pred.eco$year<-year.unscaled
eqr.mean.pred.eco<-eqr.mean.pred.eco %>% relocate(year)

eqr.mean.pred.eco %>% ggplot() +
  ylim(0.18,0.4) +
  geom_line(aes(y = fit.1, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#B22222") + # lower EQR
  geom_line(aes(y = fit.2, x = year), color = "#E66100", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.5, x = year), color = "#00468B", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#00468B") + # higher EQR
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to mean EQR (Recovering communities)

eco.deg<-c("CR02", "DE04", "DE06", "DE07", "DE08", "DE09", "FI01", "FR04", "GR01", "GR02", "GR04", "HU01", "IR01", "SP07", "SW05", "UK01")

eco.eqr.rec<-glmmTMB(eco_bray ~ year_s * eqr.mean + ecoi + (1|basin_f),
                     family = beta_family(link = "logit"),
                     REML = TRUE,
                     data = dat[!(dat$basin_f %in% eco.deg), ])
summary(eco.eqr.rec)

# Predicted year effect given values of predictor

eqr.rec.pred<-predict(eco.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = mean(dat$eqr.mean), "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.rec.pred.eco<-data.frame(eqr.rec.pred)
eqr.mean01<-predict(eco.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -1.61, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean02<-predict(eco.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -0.67, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean03<-predict(eco.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 0.27, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean04<-predict(eco.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 1.20, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean05<-predict(eco.eqr.rec, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 2.13, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)

eqr.rec.pred.eco<-cbind(eqr.rec.pred.eco, eqr.mean01, eqr.mean02, eqr.mean03, eqr.mean04, eqr.mean05)
colnames(eqr.rec.pred.eco)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
eqr.rec.pred.eco<-eqr.rec.pred.eco %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                              fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # eqr 0.2
                                              fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # eqr 0.4
                                              fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # eqr 0.6
                                              fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # eqr 0.8
                                              fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # eqr 1.0
eqr.rec.pred.eco$year<-year.unscaled
eqr.rec.pred.eco<-eqr.rec.pred.eco %>% relocate(year)

eqr.rec.pred.eco %>% ggplot() +
  ylim(0.18,0.4) +
  geom_line(aes(y = fit.1, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#B22222") + # lower EQR
  geom_line(aes(y = fit.2, x = year), color = "#E66100", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.5, x = year), color = "#00468B", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#00468B") + # higher EQR
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to temperature trend

eco.temp.change<-glmmTMB(eco_bray ~ year_s * temp.change + ecoi + (1|basin_f),
                         family = beta_family(link = "logit"),
                         REML = TRUE,
                         data = dat)
summary(eco.temp.change)
performance(eco.temp.change)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = eco.temp.change, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

temp.change.pred.eco<-predict(eco.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = mean(dat$temp.change), "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change.pred.eco<-data.frame(temp.change.pred.eco)
temp.change01<-predict(eco.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = -2.02, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change02<-predict(eco.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = -0.96, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change03<-predict(eco.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 0.09, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change04<-predict(eco.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 1.14, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change05<-predict(eco.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 2.19, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)

temp.change.pred.eco<-cbind(temp.change.pred.eco, temp.change01, temp.change02, temp.change03, temp.change04, temp.change05)
colnames(temp.change.pred.eco)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
temp.change.pred.eco<-temp.change.pred.eco %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                      fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # temp.change 0.2
                                                      fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # temp.change 0.4
                                                      fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # temp.change 0.6
                                                      fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # temp.change 0.8
                                                      fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # temp.change 1.0
temp.change.pred.eco$year<-year.unscaled
temp.change.pred.eco<-temp.change.pred.eco %>% relocate(year)

temp.change.pred.eco %>% ggplot() +
  ylim(0.18,0.4) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower slope (less warming or even cooling)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher slope (more warming)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to mean temperature trend (Recovering communities)

eco.temp.rec<-glmmTMB(eco_bray ~ year_s * temp.change + ecoi + (1|basin_f),
                      family = beta_family(link = "logit"),
                      REML = TRUE,
                      data = dat[!(dat$basin_f %in% eco.deg), ])
summary(eco.temp.rec)

# Predicted year effect given values of predictor

temp.rec.pred.eco<-predict(eco.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = mean(dat$temp.change), "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.rec.pred.eco<-data.frame(temp.rec.pred.eco)
temp.change01<-predict(eco.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = -2.02, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change02<-predict(eco.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = -0.96, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change03<-predict(eco.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 0.09, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change04<-predict(eco.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 1.14, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change05<-predict(eco.temp.rec, newdata = data.frame("year_s" = year.pred, "temp.change" = 2.19, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)

temp.rec.pred.eco<-cbind(temp.rec.pred.eco, temp.change01, temp.change02, temp.change03, temp.change04, temp.change05)
colnames(temp.rec.pred.eco)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
temp.rec.pred.eco<-temp.rec.pred.eco %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # temp.change 0.2
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # temp.change 0.4
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # temp.change 0.6
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # temp.change 0.8
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # temp.change 1.0
temp.rec.pred.eco$year<-year.unscaled
temp.rec.pred.eco<-temp.rec.pred.eco %>% relocate(year)

temp.rec.pred.eco %>% ggplot() +
  ylim(0.18,0.4) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower slope (less warming or even cooling)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher slope (more warming)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to land cover stress

eco.land<-glmmTMB(eco_bray ~ year_s * PC1 + year_s * PC2 + ecoi + (1|basin_f),
                  family = beta_family(link = "logit"),
                  REML = TRUE,
                  data = dat)
summary(eco.land)
performance(eco.land)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = eco.land, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = dat$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(dat$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

land.pred.eco<-predict(eco.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = mean(dat$PC2), "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land.pred.eco<-data.frame(land.pred.eco)
land01<-predict(eco.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -1.44, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land02<-predict(eco.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -0.28, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land03<-predict(eco.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 0.88, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land04<-predict(eco.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 2.03, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land05<-predict(eco.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 3.18, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)

land.pred.eco<-cbind(land.pred.eco, land01, land02, land03, land04, land05)
colnames(land.pred.eco)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
land.pred.eco<-land.pred.eco %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                        fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # higher forest
                                        fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2),
                                        fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3),
                                        fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4),
                                        fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # higher urban
land.pred.eco$year<-year.unscaled
land.pred.eco<-land.pred.eco %>% relocate(year)

land.pred.eco %>% ggplot() +
  ylim(0.18,0.4) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower stress (more forest)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") + 
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") + 
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") + 
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher stress (more urban)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in ꞵ-diversity linked to land cover stress (Recovering communities)

eco.land.rec<-glmmTMB(eco_bray ~ year_s * PC1 + year_s * PC2 + ecoi + (1|basin_f),
                      family = beta_family(link = "logit"),
                      REML = TRUE,
                      data = dat[!(dat$basin_f %in% eco.deg), ])
summary(eco.land.rec)

# Predicted year effect given values of predictor

land.rec.pred.eco<-predict(eco.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = mean(dat$PC2), "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land.rec.pred.eco<-data.frame(land.rec.pred.eco)
land01<-predict(eco.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -1.44, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land02<-predict(eco.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = -0.28, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land03<-predict(eco.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 0.88, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land04<-predict(eco.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 2.03, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)
land05<-predict(eco.land.rec, newdata = data.frame("year_s" = year.pred, "PC1" = mean(dat$PC1), "PC2" = 3.18, "ecoi" = mean(dat$ecoi), "basin_f" = NA), type = "response", se.fit = TRUE)

land.rec.pred.eco<-cbind(land.rec.pred.eco, land01, land02, land03, land04, land05)
colnames(land.rec.pred.eco)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
land.rec.pred.eco<-land.rec.pred.eco %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # higher forest
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2),
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3),
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4),
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # higher urban
land.rec.pred.eco$year<-year.unscaled
land.rec.pred.eco<-land.rec.pred.eco %>% relocate(year)

land.rec.pred.eco %>% ggplot() +
  ylim(0.18,0.4) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower stress (more forest)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") + 
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") + 
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") + 
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher stress (more urban)
  labs(y = "ꞵ-diversity", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



### --- Visualize trends per basin and biodiversity facet

basin.trends<-merge(bio.trends, tax.trends, by = c("basin", "year", "index", "fit"), all.x = TRUE, all.y = TRUE)
basin.trends<-merge(basin.trends, eco.trends, by = c("basin", "year", "index", "fit"), all.x = TRUE, all.y = TRUE)

basin.trends %>% ggplot(aes(x = year, y = fit, group = index, color = index)) +
  ylim(0 ,1) +
  geom_line(linewidth = 2) +
  labs(y = NULL, x = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 30)) +
  scale_color_manual(values = c("#e66101", "#5e3c99", "#008837")) +
  facet_wrap(~ basin, ncol = 9) +
  theme(strip.background = element_blank(), strip.text = element_blank())




### --- Extract estimates of basin trends from models

basin.estimates<-as.data.frame(coef(summary(tax.basin))$cond)
colnames(basin.estimates)<-c("tax_Estimate", "tax_StdError", "tax_zValue", "tax_pValue")
basin.estimates$term<-rownames(basin.estimates)

bio.slopes<-as.data.frame(coef(summary(bio.basin))$cond)
colnames(bio.slopes)<-c("bio_Estimate", "bio_StdError", "bio_zValue", "bio_pValue")
bio.slopes$term<-rownames(bio.slopes)

eco.slopes<-as.data.frame(coef(summary(eco.basin))$cond)
colnames(eco.slopes)<-c("eco_Estimate", "eco_StdError", "eco_zValue", "eco_pValue")
eco.slopes$term<-rownames(eco.slopes)

basin.estimates<-merge(basin.estimates, bio.slopes, by = c("term"))
basin.estimates<-merge(basin.estimates, eco.slopes, by = c("term"))

# Select basin:year

basin.estimates<-basin.estimates %>% filter(grepl("year_s", term))
basin.estimates$basin<-basins

# Calculate CI (alpha = 0.05) from SE

basin.estimates<-basin.estimates %>% mutate(tax.low = basin.estimates$tax_Estimate - (basin.estimates$tax_StdError*1.96),
                                            tax.high = basin.estimates$tax_Estimate + (basin.estimates$tax_StdError*1.96),
                                            bio.low = basin.estimates$bio_Estimate - (basin.estimates$bio_StdError*1.96),
                                            bio.high = basin.estimates$bio_Estimate + (basin.estimates$bio_StdError*1.96),
                                            eco.low = basin.estimates$eco_Estimate - (basin.estimates$eco_StdError*1.96),
                                            eco.high = basin.estimates$eco_Estimate + (basin.estimates$eco_StdError*1.96),)

# Plots

basin.estimates %>% ggplot(aes(x = tax_Estimate, y = reorder(basin, tax_Estimate))) +
  xlim(-2, 2) +
  geom_point(size = 2, color = "#008837") +
  geom_errorbar(aes(xmin = tax.low, xmax = tax.high), color = "#008837") +
  scale_y_discrete(limits = rev) +
  labs(x = "Slope", y = "Basin ID") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  annotate("rect", xmin = c(-2), xmax = c(0), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = c("lightgray"))

basin.estimates %>% ggplot(aes(x = bio_Estimate, reorder(basin, bio_Estimate))) +
  xlim(-1, 1) +
  geom_point(size = 2, color = "#5e3c99") +
  geom_errorbar(aes(xmin = bio.low, xmax = bio.high), color = "#5e3c99") +
  scale_y_discrete(limits = rev) +
  labs(x = "Slope", y = "Basin ID") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  annotate("rect", xmin = c(-1), xmax = c(0), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = c("lightgray"))

basin.estimates %>% ggplot(aes(x = eco_Estimate, reorder(basin, eco_Estimate))) +
  xlim(-1, 1) +
  geom_point(size = 2, color = "#e66101") +
  geom_errorbar(aes(xmin = eco.low, xmax = eco.high), color = "#e66101") +
  scale_y_discrete(limits = rev) +
  labs(x = "Slope", y = "Basin ID") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  annotate("rect", xmin = c(-1), xmax = c(0), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = c("lightgray"))




### --- "Johnson-Neyman point" or the "simple slopes crossover point" : the point where the slope changes direction

# 1. Mean EQR

# Get marginal effects of year_s at different eqr.mean values
me.tax.eqr.mean<-slopes(tax.eqr.mean, variables = "year_s", newdata = datagrid(eqr.mean = seq(min(dat$eqr.mean), max(dat$eqr.mean), length.out = 100), taxi = mean(dat$taxi), year_s = mean(dat$year_s), basin_f = NA))
me.tax.eqr.mean<-as.data.frame(me.tax.eqr.mean)

me.bio.eqr.mean<-slopes(bio.eqr.mean, variables = "year_s", newdata = datagrid(eqr.mean = seq(min(dat$eqr.mean), max(dat$eqr.mean), length.out = 100), bioi = mean(dat$bioi), year_s = mean(dat$year_s), basin_f = NA))
me.bio.eqr.mean<-as.data.frame(me.bio.eqr.mean)

me.eco.eqr.mean<-slopes(eco.eqr.mean, variables = "year_s", newdata = datagrid(eqr.mean = seq(min(dat$eqr.mean), max(dat$eqr.mean), length.out = 100), ecoi = mean(dat$ecoi), year_s = mean(dat$year_s), basin_f = NA))
me.eco.eqr.mean<-as.data.frame(me.eco.eqr.mean)

me.eqr.mean<-cbind(me.tax.eqr.mean[, c("estimate", "conf.low", "conf.high", "eqr.mean")], # eqr.mean only needed once
                   me.bio.eqr.mean[, c("estimate", "conf.low", "conf.high")],
                   me.eco.eqr.mean[, c("estimate", "conf.low", "conf.high")])
colnames(me.eqr.mean)<-c("tax.estimate", "tax.conf.low", "tax.conf.high", "eqr.mean",
                         "bio.estimate", "bio.conf.low", "bio.conf.high",
                         "eco.estimate", "eco.conf.low", "eco.conf.high")
scale(env$eqr)
me.eqr.mean$eqr.unscaled<-me.eqr.mean$eqr.mean * 0.2832634 + 0.6706811

me.eqr.mean %>% ggplot() +
  ylim(-0.018, 0.018) +
  geom_line(aes(x = eqr.unscaled, y = tax.estimate), color = "#008837", linewidth = 3, linetype = "longdash") +
  geom_line(aes(x = eqr.unscaled, y = bio.estimate), color = "#5e3c99", linewidth = 3, linetype = "longdash") +
  geom_line(aes(x = eqr.unscaled, y = eco.estimate), color = "#e66101", linewidth = 3, linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Mean EQR", y = "ꞵ-diversity trend") +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_reverse()



# 2. Temperature change

# Get marginal effects of year_s at different temp.change values
me.tax.temp.change<-slopes(tax.temp.change, variables = "year_s", newdata = datagrid(temp.change = seq(min(dat$temp.change), max(dat$temp.change), length.out = 100), taxi = mean(dat$taxi), year_s = mean(dat$year_s), basin_f = NA))
me.tax.temp.change<-as.data.frame(me.tax.temp.change)

me.bio.temp.change<-slopes(bio.temp.change, variables = "year_s", newdata = datagrid(temp.change = seq(min(dat$temp.change), max(dat$temp.change), length.out = 100), bioi = mean(dat$bioi), year_s = mean(dat$year_s), basin_f = NA))
me.bio.temp.change<-as.data.frame(me.bio.temp.change)

me.eco.temp.change<-slopes(eco.temp.change, variables = "year_s", newdata = datagrid(temp.change = seq(min(dat$temp.change), max(dat$temp.change), length.out = 100), ecoi = mean(dat$ecoi), year_s = mean(dat$year_s), basin_f = NA))
me.eco.temp.change<-as.data.frame(me.eco.temp.change)

me.temp.change<-cbind(me.tax.temp.change[, c("estimate", "conf.low", "conf.high", "temp.change")], # temp.change only needed once
                      me.bio.temp.change[, c("estimate", "conf.low", "conf.high")],
                      me.eco.temp.change[, c("estimate", "conf.low", "conf.high")])
colnames(me.temp.change)<-c("tax.estimate", "tax.conf.low", "tax.conf.high", "temp.change",
                            "bio.estimate", "bio.conf.low", "bio.conf.high",
                            "eco.estimate", "eco.conf.low", "eco.conf.high")
scale(env$temp.change)
me.temp.change$temp.unscaled<-me.temp.change$temp.change * 0.02638092 + 0.04086026

me.temp.change %>% ggplot() +
  ylim(-0.015, 0.015) +
  geom_line(aes(x = temp.unscaled, y = tax.estimate), color = "#008837", linewidth = 3, linetype = "longdash") +
  geom_line(aes(x = temp.unscaled, y = bio.estimate), color = "#5e3c99", linewidth = 3, linetype = "solid") +
  geom_line(aes(x = temp.unscaled, y = eco.estimate), color = "#e66101", linewidth = 3, linetype = "longdash") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Temperature trend", y = "ꞵ-diversity trend") +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# 2. Land cover pressure (urbanization)

# Get marginal effects of year_s at different land values
me.tax.land<-slopes(tax.land, variables = "year_s", newdata = datagrid(PC2 = seq(min(dat$PC2), max(dat$PC2), length.out = 100), PC1 = mean(dat$PC1), taxi = mean(dat$taxi), year_s = mean(dat$year_s), basin_f = NA))
me.tax.land<-as.data.frame(me.tax.land)

me.bio.land<-slopes(bio.land, variables = "year_s", newdata = datagrid(PC2 = seq(min(dat$PC2), max(dat$PC2), length.out = 100), PC1 = mean(dat$PC1), bioi = mean(dat$bioi), year_s = mean(dat$year_s), basin_f = NA))
me.bio.land<-as.data.frame(me.bio.land)

me.eco.land<-slopes(eco.land, variables = "year_s", newdata = datagrid(PC2 = seq(min(dat$PC2), max(dat$PC2), length.out = 100), PC1 = mean(dat$PC1), ecoi = mean(dat$ecoi), year_s = mean(dat$year_s), basin_f = NA))
me.eco.land<-as.data.frame(me.eco.land)

me.land<-cbind(me.tax.land[, c("estimate", "conf.low", "conf.high", "PC2")], # PC2 only needed once
               me.bio.land[, c("estimate", "conf.low", "conf.high")],
               me.eco.land[, c("estimate", "conf.low", "conf.high")])
colnames(me.land)<-c("tax.estimate", "tax.conf.low", "tax.conf.high", "PC2",
                     "bio.estimate", "bio.conf.low", "bio.conf.high",
                     "eco.estimate", "eco.conf.low", "eco.conf.high")

me.land %>% ggplot() +
  ylim(-0.022, 0.022) +
  geom_line(aes(x = PC2, y = tax.estimate), color = "#008837", linewidth = 3, linetype = "longdash") +
  geom_line(aes(x = PC2, y = bio.estimate), color = "#5e3c99", linewidth = 3, linetype = "solid") +
  geom_line(aes(x = PC2, y = eco.estimate), color = "#e66101", linewidth = 3, linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Land cover", y = "ꞵ-diversity trend") +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())



### --- Residual variation

dat.res<-dat[, c("basin", "year")]


# Residuals (mean trend models) vs distance among sites

dat.res<-dat.res %>% mutate(tax.mean.res = residuals(tax.mean, type = c("response")),
                            bio.mean.res = residuals(bio.mean, type = c("response")),
                            eco.mean.res = residuals(eco.mean, type = c("response")))

# 1. Mean and maximum distances between sampling sites

# Function to calculate mean and max distances for each basin

geo.distances<-data.frame()

for (b in basins) {
  df_basin<-subset(env, basin == b)
  
  coords<-as.matrix(df_basin[, c("longitude", "latitude")])
  dist_mat<-as.matrix(dist(coords))
  dists<-dist_mat[upper.tri(dist_mat)]
  
  geo.distances<-rbind(geo.distances, data.frame(
    basin = b,
    mean_distance = mean(dists),
    max_distance = max(dists)
  ))
}

dat.res<-merge(dat.res, geo.distances, by = c("basin"), all.x = TRUE)

## Mean distance

dat.res %>% ggplot(aes(x=log(mean_distance), y = log(tax.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Mean distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

dat.res %>% ggplot(aes(x=log(mean_distance), y = log(bio.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Mean distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

dat.res %>% ggplot(aes(x=log(mean_distance), y = log(eco.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Mean distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 15)) +

  
## Maximum distance

dat.res %>% ggplot(aes(x=log(max_distance), y = log(tax.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Maximum distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +

dat.res %>% ggplot(aes(x=log(max_distance), y = log(bio.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Maximum distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +

dat.res %>% ggplot(aes(x=log(max_distance), y = log(eco.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Maximum distance between sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5))


### --- Summary descriptors

tb1<-alpha %>% group_by(basin) %>% summarise(no.sites = length(unique(site)))
median(tb1$no.sites)

tb2<-dat %>% group_by(basin) %>% summarise(no.years = length(unique(year)))
median(tb2$no.years)

dat.res<-merge(dat.res, tb1, by = c("basin"), all.x = TRUE)
dat.res<-merge(dat.res, tb2, by = c("basin"), all.x = TRUE)

## Number of sites

dat.res %>% ggplot(aes(x = log(no.sites), y = log(tax.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Number of sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +
  
  dat.res %>% ggplot(aes(x = log(no.sites), y = log(bio.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Number of sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +
  
  dat.res %>% ggplot(aes(x = log(no.sites), y = log(eco.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Number of sites", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +
  
  ## Number of years
  
  dat.res %>% ggplot(aes(x = log(no.years), y = log(tax.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Time series length", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +
  
  dat.res %>% ggplot(aes(x = log(no.years), y = log(bio.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Time series length", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5)) +
  
  dat.res %>% ggplot(aes(x = log(no.years), y = log(eco.mean.res))) +
  geom_point() +
  theme_classic() +
  labs(x = "Time series length", y = "Fitted residuals") +
  theme(text = element_text(size = 14.5))




###### ----- Ecopart: mechanisms driving ꞵ-diversity change -----

# Import data set

dat1<-read.csv(file.choose(), header = TRUE) # Ecopart.csv

## How environmental variables affect additions and subtractions

dat1<-merge(dat1, dat[, c("basin", "year", "eqr.i", "eqr.mean", "temp.change", "PC1", "PC2")],
            by = c("basin", "year"), all.x = TRUE)

## -- Taxonomic additions and subtractions

# 1. Mean EQR

tax.add.eqr.mean<-glmmTMB(tax_addition ~ eqr.mean + (1 |basin),
                            family = gaussian,
                            data = dat1)
summary(tax.add.eqr.mean)

tax.sub.eqr.mean<-glmmTMB(tax_subtraction ~ eqr.mean + (1 |basin),
                            family = gaussian,
                            data = dat1)
summary(tax.sub.eqr.mean)

# This gives predicted values across the range of eqr.mean, marginalizing over random effects
tax.add.eqr.preds<-data.frame(ggpredict(tax.add.eqr.mean, terms = "eqr.mean", type = "fixed"))  # type = "re" if you want conditional on random effects
scale(land.mean$eqr.mean)
tax.add.eqr.preds$mean.eqr<-tax.add.eqr.preds$x * 0.2189595 + 0.6935694
tax.sub.eqr.preds<-data.frame(ggpredict(tax.sub.eqr.mean, terms = "eqr.mean", type = "fixed"))
tax.eqr.preds<-cbind(tax.add.eqr.preds, tax.sub.eqr.preds)
colnames(tax.eqr.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add", "mean.eqr",
                           "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

tax.eqr.preds %>% ggplot() +
  geom_line(aes(x = mean.eqr, y = predicted.add), color = "#A1B56C", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.eqr, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = mean.eqr, y = predicted.sub), color = "#7E587E", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.eqr, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  labs(x = "Mean EQR", y = "Change in ꞵ-diversity") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(-0.05, 0.05)) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  scale_x_reverse()


# 2. Temperature change

tax.add.temp.change<-glmmTMB(tax_addition ~ temp.change + (1 |basin),
                            family = gaussian,
                            data = dat1)
summary(tax.add.temp.change)

tax.sub.temp.change<-glmmTMB(tax_subtraction ~ temp.change + (1 |basin),
                            family = gaussian,
                            data = dat1)
summary(tax.sub.temp.change)

tax.add.temp.preds<-data.frame(ggpredict(tax.add.temp.change, terms = "temp.change", type = "fixed"))  # type = "re" if you want conditional on random effects
scale(env$temp.change)
tax.add.temp.preds$mean.temp<-tax.add.temp.preds$x * 0.02638092 + 0.04086026
tax.sub.temp.preds<-data.frame(ggpredict(tax.sub.temp.change, terms = "temp.change", type = "fixed"))
tax.temp.preds<-cbind(tax.add.temp.preds, tax.sub.temp.preds)
colnames(tax.temp.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add", "mean.temp",
                           "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

tax.temp.preds %>% ggplot() +
  geom_line(aes(x = mean.temp, y = predicted.add), color = "#A1B56C", linetype = "solid", linewidth = 3) +
  geom_ribbon(aes(x = mean.temp, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = mean.temp, y = predicted.sub), color = "#7E587E", linetype = "solid", linewidth = 3) +
  geom_ribbon(aes(x = mean.temp, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  labs(x = "Temperature trend", y = "Change in ꞵ-diversity") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(-0.05, 0.05)) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


# 3. Land cover stress

tax.add.land<-glmmTMB(tax_addition ~ PC1 + PC2 + (1 |basin),
                             family = gaussian,
                             data = dat1)
summary(tax.add.land)

tax.sub.land<-glmmTMB(tax_subtraction ~ PC1 + PC2 + (1 |basin),
                             family = gaussian,
                             data = dat1)
summary(tax.sub.land)

tax.add.land.preds<-data.frame(ggpredict(tax.add.land, terms = "PC2", type = "fixed"))  # type = "re" if you want conditional on random effects
tax.sub.land.preds<-data.frame(ggpredict(tax.sub.land, terms = "PC2", type = "fixed"))
tax.land.preds<-cbind(tax.add.land.preds, tax.sub.land.preds)
colnames(tax.land.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add",
                            "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

tax.land.preds %>% ggplot() +
  geom_line(aes(x = x.add, y = predicted.add), color = "#A1B56C", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = x.add, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = x.sub, y = predicted.sub), color = "#7E587E", linetype = "solid", linewidth = 3) +
  geom_ribbon(aes(x = x.sub, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  labs(x = "Land cover", y = "Change in ꞵ-diversity") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(-0.05, 0.05)) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))



## -- Biological traits additions and subtractions

# 1. Mean EQR

bio.add.eqr.mean<-glmmTMB(bio_addition ~ eqr.mean + (1 |basin),
                          family = gaussian,
                          data = dat1)
summary(bio.add.eqr.mean)

bio.sub.eqr.mean<-glmmTMB(bio_subtraction ~ eqr.mean + (1 |basin),
                          family = gaussian,
                          data = dat1)
summary(bio.sub.eqr.mean)

bio.add.eqr.preds<-data.frame(ggpredict(bio.add.eqr.mean, terms = "eqr.mean", type = "fixed"))  # type = "re" if you want conditional on random effects
bio.add.eqr.preds$mean.eqr<-bio.add.eqr.preds$x * 0.2189595 + 0.6935694
bio.sub.eqr.preds<-data.frame(ggpredict(bio.sub.eqr.mean, terms = "eqr.mean", type = "fixed"))
bio.eqr.preds<-cbind(bio.add.eqr.preds, bio.sub.eqr.preds)
colnames(bio.eqr.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add", "mean.eqr",
                           "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

bio.eqr.preds %>% ggplot() +
  geom_line(aes(x = mean.eqr, y = predicted.add), color = "#A1B56C", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.eqr, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = mean.eqr, y = predicted.sub), color = "#7E587E", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.eqr, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  labs(x = "Mean EQR", y = "Change in ꞵ-diversity") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(-0.01, 0.01)) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  scale_x_reverse()


# 2. Temperature change

bio.add.temp.change<-glmmTMB(bio_addition ~ temp.change + (1 |basin),
                             family = gaussian,
                             data = dat1)
summary(bio.add.temp.change)

bio.sub.temp.change<-glmmTMB(bio_subtraction ~ temp.change + (1 |basin),
                             family = gaussian,
                             data = dat1)
summary(bio.sub.temp.change)

bio.add.temp.preds<-data.frame(ggpredict(bio.add.temp.change, terms = "temp.change", type = "fixed"))  # type = "re" if you want conditional on random effects
bio.add.temp.preds$mean.temp<-bio.add.temp.preds$x * 0.02638092 + 0.04086026
bio.sub.temp.preds<-data.frame(ggpredict(bio.sub.temp.change, terms = "temp.change", type = "fixed"))
bio.temp.preds<-cbind(bio.add.temp.preds, bio.sub.temp.preds)
colnames(bio.temp.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add", "mean.temp",
                            "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

bio.temp.preds %>% ggplot() +
  geom_line(aes(x = mean.temp, y = predicted.add), color = "#A1B56C", linetype = "solid", linewidth = 3) +
  geom_ribbon(aes(x = mean.temp, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = mean.temp, y = predicted.sub), color = "#7E587E", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.temp, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  labs(x = "Temperature trend", y = "Change in ꞵ-diversity") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(-0.01, 0.01)) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


# 3. Land cover stress

bio.add.land<-glmmTMB(bio_addition ~ PC1 + PC2 + (1 |basin),
                      family = gaussian,
                      data = dat1)
summary(bio.add.land)

bio.sub.land<-glmmTMB(bio_subtraction ~ PC1 + PC2 + (1 |basin),
                      family = gaussian,
                      data = dat1)
summary(bio.sub.land)

bio.add.land.preds<-data.frame(ggpredict(bio.add.land, terms = "PC2", type = "fixed"))  # type = "re" if you want conditional on random effects
bio.sub.land.preds<-data.frame(ggpredict(bio.sub.land, terms = "PC2", type = "fixed"))
bio.land.preds<-cbind(bio.add.land.preds, bio.sub.land.preds)
colnames(bio.land.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add",
                            "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

bio.land.preds %>% ggplot() +
  geom_line(aes(x = x.add, y = predicted.add), color = "#A1B56C", linetype = "solid", linewidth = 3) +
  geom_ribbon(aes(x = x.add, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = x.sub, y = predicted.sub), color = "#7E587E", linetype = "solid", linewidth = 3) +
  geom_ribbon(aes(x = x.sub, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  labs(x = "Land cover", y = "Change in ꞵ-diversity") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(-0.01, 0.01)) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))



## -- Ecological traits additions and subtractions

# 1. Mean EQR

eco.add.eqr.mean<-glmmTMB(eco_addition ~ eqr.mean + (1 |basin),
                          family = gaussian,
                          data = dat1)
summary(eco.add.eqr.mean)

eco.sub.eqr.mean<-glmmTMB(eco_subtraction ~ eqr.mean + (1 |basin),
                          family = gaussian,
                          data = dat1)
summary(eco.sub.eqr.mean)

eco.add.eqr.preds<-data.frame(ggpredict(eco.add.eqr.mean, terms = "eqr.mean", type = "fixed"))  # type = "re" if you want conditional on random effects
eco.add.eqr.preds$mean.eqr<-eco.add.eqr.preds$x * 0.2189595 + 0.6935694
eco.sub.eqr.preds<-data.frame(ggpredict(eco.sub.eqr.mean, terms = "eqr.mean", type = "fixed"))
eco.eqr.preds<-cbind(eco.add.eqr.preds, eco.sub.eqr.preds)
colnames(eco.eqr.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add", "mean.eqr",
                           "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

eco.eqr.preds %>% ggplot() +
  geom_line(aes(x = mean.eqr, y = predicted.add), color = "#A1B56C", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.eqr, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = mean.eqr, y = predicted.sub), color = "#7E587E", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.eqr, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  labs(x = "Mean EQR", y = "Change in ꞵ-diversity") +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  coord_cartesian(ylim = c(-0.005, 0.005)) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  scale_x_reverse()


# 2. Temperature change

eco.add.temp.change<-glmmTMB(eco_addition ~ temp.change + (1 |basin),
                             family = gaussian,
                             data = dat1)
summary(eco.add.temp.change)

eco.sub.temp.change<-glmmTMB(eco_subtraction ~ temp.change + (1 |basin),
                             family = gaussian,
                             data = dat1)
summary(eco.sub.temp.change)

eco.add.temp.preds<-data.frame(ggpredict(eco.add.temp.change, terms = "temp.change", type = "fixed"))  # type = "re" if you want conditional on random effects
eco.add.temp.preds$mean.temp<-eco.add.temp.preds$x * 0.02638092 + 0.04086026
eco.sub.temp.preds<-data.frame(ggpredict(eco.sub.temp.change, terms = "temp.change", type = "fixed"))
eco.temp.preds<-cbind(eco.add.temp.preds, eco.sub.temp.preds)
colnames(eco.temp.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add", "mean.temp",
                            "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

eco.temp.preds %>% ggplot() +
  geom_line(aes(x = mean.temp, y = predicted.add), color = "#A1B56C", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.temp, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = mean.temp, y = predicted.sub), color = "#7E587E", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = mean.temp, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  coord_cartesian(ylim = c(-0.005, 0.005)) +
  labs(x = "Temperature trend", y = "Change in ꞵ-diversity") +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


# 3. Land cover stress

eco.add.land<-glmmTMB(eco_addition ~ PC1 + PC2 + (1 |basin),
                      family = gaussian,
                      data = dat1)
summary(eco.add.land)

eco.sub.land<-glmmTMB(eco_subtraction ~ PC1 + PC2 + (1 |basin),
                      family = gaussian,
                      data = dat1)
summary(eco.sub.land)

eco.add.land.preds<-data.frame(ggpredict(eco.add.land, terms = "PC2", type = "fixed"))  # type = "re" if you want conditional on random effects
eco.sub.land.preds<-data.frame(ggpredict(eco.sub.land, terms = "PC2", type = "fixed"))
eco.land.preds<-cbind(eco.add.land.preds, eco.sub.land.preds)
colnames(eco.land.preds)<-c("x.add", "predicted.add", "std.error.add", "conf.low.add", "conf.high.add", "group.add",
                            "x.sub", "predicted.sub", "std.error.sub", "conf.low.sub", "conf.high.sub", "group.sub")

eco.land.preds %>% ggplot() +
  geom_line(aes(x = x.add, y = predicted.add), color = "#A1B56C", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = x.add, ymin = conf.low.add, ymax = conf.high.add), fill = "#A1B56C", alpha = 0.5) +
  geom_line(aes(x = x.sub, y = predicted.sub), color = "#7E587E", linetype = "longdash", linewidth = 3) +
  geom_ribbon(aes(x = x.sub, ymin = conf.low.sub, ymax = conf.high.sub), fill = "#7E587E", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  labs(x = "Land cover", y = "Change in ꞵ-diversity") +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  coord_cartesian(ylim = c(-0.005, 0.005)) +
  theme_bw() +
  theme(text = element_text(size = 35), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))




##### ----- Individual traits: additions and subtractions ----- #####

## -- Biological traits

biotraits<-read.csv(file.choose(), header = TRUE) # Ecopart Biological traits.csv
traits<-unique(biotraits$biotrait)
biotraits<-biotraits[!is.na(biotraits$additive), ]

# Additions
model<-glmmTMB(additive ~ 1 + (1|basin),
               family = gaussian,
               data = subset(biotraits, biotrait == traits[1]))
bio.add.coef<-data.frame(coef(summary(model))$cond)

for (i in 2:length(traits)){
  model<-glmmTMB(additive ~ 1 + (1|basin),
                 family = gaussian,
                 data = subset(biotraits, biotrait == traits[i]))
  coef<-data.frame(coef(summary(model))$cond)
  bio.add.coef<-rbind(bio.add.coef, coef)
}
bio.add.coef$biotrait<-traits
bio.add.coef$component<-"additions"

# Subtractions
model<-glmmTMB(subtractive ~ 1 + (1|basin),
               family = gaussian,
               data = subset(biotraits, biotrait == traits[1]))
bio.sub.coef<-data.frame(coef(summary(model))$cond)

for (i in 2:length(traits)){
  model<-glmmTMB(subtractive ~ 1 + (1|basin),
                 family = gaussian,
                 data = subset(biotraits, biotrait == traits[i]))
  coef<-data.frame(coef(summary(model))$cond)
  bio.sub.coef<-rbind(bio.sub.coef, coef)
}
bio.sub.coef$biotrait<-traits
bio.sub.coef$component<-"subtractions"

bio.coef<-rbind(bio.add.coef, bio.sub.coef)
bio.coef<-bio.coef %>% mutate(CI05 = Estimate - (Std..Error*1.96),
                              CI95 = Estimate + (Std..Error*1.96))
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
  xlim(-0.0015, 0.0015) +
  geom_point(aes(color = significance), size = 3) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = significance),
                linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c(">8 cm", "4-8 cm", "2-4 cm", "1-2 cm", "0.5-1 cm", "0.25-0.5 cm", "<0.25 cm",
                                            "None", "Housing", "Dormant eggs", "Diapause", "Cocoons",
                                            "Other", "Parasites", "Xylophagous", "Miners", "Predators",
                                            "Passive filter feeders", "Active filter feeders", "Gatherers", "Shedders", "Grazers",
                                            ">1 per year", "1 per year", "<1 per year",
                                            "Aquatic passive", "Aquatic active", "Aerial passive", "Aerial active",
                                            "Adult", "Nymph", "Larva", "Egg")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        legend.position = "none") +
  scale_color_manual(values = c("#A1B56C", "gray", "#7E587E")) +
  labs(y = NULL, )



## -- Ecological traits

ecotraits<-read.csv(file.choose(), header = TRUE) # Ecopart Ecological traits.csv
traits<-unique(ecotraits$ecotrait)
ecotraits<-ecotraits[!is.na(ecotraits$additive), ]
ecotraits<-ecotraits[!is.na(ecotraits$subtractive), ]

# Additions
model<-glmmTMB(additive ~ 1 + (1|basin),
               family = gaussian,
               data = subset(ecotraits, ecotrait == traits[1]))
eco.add.coef<-data.frame(coef(summary(model))$cond)

for (i in 2:length(traits)){
  model<-glmmTMB(additive ~ 1 + (1|basin),
                 family = gaussian,
                 data = subset(ecotraits, ecotrait == traits[i]))
  coef<-data.frame(coef(summary(model))$cond)
  eco.add.coef<-rbind(eco.add.coef, coef)
}
eco.add.coef$ecotrait<-traits
eco.add.coef$component<-"additions"


# Subtractions
model<-glmmTMB(subtractive ~ 1 + (1|basin),
               family = gaussian,
               data = subset(ecotraits, ecotrait == traits[1]))
eco.sub.coef<-data.frame(coef(summary(model))$cond)

for (i in 2:length(traits)){
  model<-glmmTMB(subtractive ~ 1 + (1|basin),
                 family = gaussian,
                 data = subset(ecotraits, ecotrait == traits[i]))
  coef<-data.frame(coef(summary(model))$cond)
  eco.sub.coef<-rbind(eco.sub.coef, coef)
}
eco.sub.coef$ecotrait<-traits
eco.sub.coef$component<-"subtractions"

eco.coef<-rbind(eco.add.coef, eco.sub.coef)
eco.coef<-eco.coef %>% mutate(CI05 = Estimate - (Std..Error*1.96),
                              CI95 = Estimate + (Std..Error*1.96))
eco.coef$significance<-ifelse(eco.coef$Pr...z..<0.05, eco.coef$component, "not significant")

eco.coef$ecotrait<-factor(eco.coef$ecotrait,
                          levels = c("tachtemp_eut", "tachtemp_psy", "tachtemp_the",
                                     "microhab_arg", "microhab_pel", "microhab_psa", "microhab_aka", "microhab_lit", "microhab_phy",
                                     "microhab_pom", "microhab_oth",
                                     "tachsap_x", "tachsap_o", "tachsap_b", "tachsap_a", "tachsap_p"))

eco.coef %>% ggplot(aes(x = Estimate, y = ecotrait, group = significance)) +
  xlim(-0.0015, 0.0015) +
  geom_point(aes(color = significance), size = 3) +
  geom_errorbar(aes(xmin = CI05, xmax = CI95, color = significance),
                linewidth = 0.8) +
  scale_y_discrete(limits = rev, labels = c("Polysaprobic", "α-Mesosaprobic", "ꞵ-Mesosaprobic", "Oligosaprobic", "Xenosaprobic",
                                            "Other", "POM", "Phytal", "Lithral", "Sense", "Psammal", "Pelal", "Argyllal",
                                            "Thermophilic", "Psychrophilic", "Eurythermal")) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        legend.position = "none") +
  scale_color_manual(values = c("#A1B56C", "gray", "#7E587E")) +
  labs(y = NULL)



#### ---- Taxonomic richness along anthropogenic stress gradients

env.basin<-dat[!duplicated(dat$basin), c("basin", "eqr.mean", "temp.change", "PC1", "PC2")]
alpha<-merge(alpha, env.basin, by = c("basin"), all.x = TRUE)

## -- Trends in taxon richness linked to mean EQR

rich.eqr.mean<-glmmTMB(tax.richness ~ year_s * eqr.mean + (1|basin_f),
                       family = nbinom2,
                       REML = TRUE,
                       data = alpha)
summary(rich.eqr.mean)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = rich.eqr.mean, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = alpha$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(alpha$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

eqr.mean.pred<-predict(rich.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = mean(alpha$eqr.mean), "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean.pred.tax<-data.frame(eqr.mean.pred)
seq(from = min(alpha$eqr.mean), to = max(alpha$eqr.mean), length.out = 5)
eqr.mean01<-predict(rich.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -1.61, "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean02<-predict(rich.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = -0.67, "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean03<-predict(rich.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 0.27, "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean04<-predict(rich.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 1.20, "basin_f" = NA), type = "response", se.fit = TRUE)
eqr.mean05<-predict(rich.eqr.mean, newdata = data.frame("year_s" = year.pred, "eqr.mean" = 2.13, "basin_f" = NA), type = "response", se.fit = TRUE)

eqr.mean.pred.tax<-cbind(eqr.mean.pred.tax, eqr.mean01, eqr.mean02, eqr.mean03, eqr.mean04, eqr.mean05)
colnames(eqr.mean.pred.tax)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
eqr.mean.pred.tax<-eqr.mean.pred.tax %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # eqr 0.2
                                                fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # eqr 0.4
                                                fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # eqr 0.6
                                                fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # eqr 0.8
                                                fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # eqr 1.0
eqr.mean.pred.tax$year<-year.unscaled
eqr.mean.pred.tax<-eqr.mean.pred.tax %>% relocate(year)

eqr.mean.pred.tax %>% ggplot() +
  ylim(5,45) +
  geom_line(aes(y = fit.1, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#B22222") + # more negative slope (degradation)
  geom_line(aes(y = fit.2, x = year), color = "#E66100", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.5, x = year), color = "#00468B", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#00468B") + # more positive slope (improvement)
  labs(y = "Richness", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 30), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in taxon richness linked to mean temperature trend

rich.temp.change<-glmmTMB(tax.richness ~ year_s * temp.change + (1|basin_f),
                          family = nbinom2,
                          REML = TRUE,
                          data = alpha)
summary(rich.temp.change)

# Test temp.changeoral autocorrelation

Esqr <- simulateResiduals(fittedModel = rich.temp.change, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = alpha$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(alpha$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

temp.change.pred.bio<-predict(rich.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = mean(alpha$temp.change), "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change.pred.bio<-data.frame(temp.change.pred.bio)
temp.change01<-predict(rich.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = -2.02, "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change02<-predict(rich.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = -0.96, "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change03<-predict(rich.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 0.09, "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change04<-predict(rich.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 1.14, "basin_f" = NA), type = "response", se.fit = TRUE)
temp.change05<-predict(rich.temp.change, newdata = data.frame("year_s" = year.pred, "temp.change" = 2.19, "basin_f" = NA), type = "response", se.fit = TRUE)

temp.change.pred.bio<-cbind(temp.change.pred.bio, temp.change01, temp.change02, temp.change03, temp.change04, temp.change05)
colnames(temp.change.pred.bio)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
temp.change.pred.bio<-temp.change.pred.bio %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                                      fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # temp.change 0.2
                                                      fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # temp.change 0.4
                                                      fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # temp.change 0.6
                                                      fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # temp.change 0.8
                                                      fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # temp.change 1.0
temp.change.pred.bio$year<-year.unscaled
temp.change.pred.bio<-temp.change.pred.bio %>% relocate(year)

temp.change.pred.bio %>% ggplot() +
  ylim(5,45) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower slope (less warming or even cooling)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "longdash") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher slope (more warming)
  labs(y = "Richness", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 30), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



## -- Trends in taxon richness linked to land cover stress

rich.land<-glmmTMB(tax.richness ~ year_s * PC1 + year_s * PC2 + (1|basin_f),
                   family = nbinom2,
                   REML = TRUE,
                   data = alpha)
summary(rich.land)

# Test temporal autocorrelation

Esqr <- simulateResiduals(fittedModel = rich.land, plot = FALSE)
rsims<-recalculateResiduals(Esqr, group = alpha$year_s, seed = 123, rotation = "estimated") # recalculate for repeated measures
testTemporalAutocorrelation(rsims, time = unique(alpha$year_s),
                            alternative = c("two.sided"), plot = F)

# Predicted year effect given values of predictor

land.pred.bio<-predict(rich.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(alpha$PC1), "PC2" = mean(alpha$PC2), "basin_f" = NA), type = "response", se.fit = TRUE)
land.pred.bio<-data.frame(land.pred.bio)
land01<-predict(rich.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(alpha$PC1), "PC2" = -1.44, "basin_f" = NA), type = "response", se.fit = TRUE)
land02<-predict(rich.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(alpha$PC1), "PC2" = -0.28, "basin_f" = NA), type = "response", se.fit = TRUE)
land03<-predict(rich.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(alpha$PC1), "PC2" = 0.88, "basin_f" = NA), type = "response", se.fit = TRUE)
land04<-predict(rich.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(alpha$PC1), "PC2" = 2.03, "basin_f" = NA), type = "response", se.fit = TRUE)
land05<-predict(rich.land, newdata = data.frame("year_s" = year.pred, "PC1" = mean(alpha$PC1), "PC2" = 3.18, "basin_f" = NA), type = "response", se.fit = TRUE)

land.pred.bio<-cbind(land.pred.bio, land01, land02, land03, land04, land05)
colnames(land.pred.bio)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
land.pred.bio<-land.pred.bio %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                                        fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # higher forest
                                        fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2),
                                        fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3),
                                        fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4),
                                        fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # higher urban
land.pred.bio$year<-year.unscaled
land.pred.bio<-land.pred.bio %>% relocate(year)

land.pred.bio %>% ggplot() +
  ylim(5,45) +
  geom_line(aes(y = fit.1, x = year), color = "#00468B", linewidth = 3, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#00468B") + # lower stress (more forest)
  geom_line(aes(y = fit.2, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#91BFDB") + 
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") + 
  geom_line(aes(y = fit.4, x = year), color = "#E66100", linewidth = 2, linetype = "solid") + 
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#E66100") + 
  geom_line(aes(y = fit.5, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#B22222") + # higher stress (more urban)
  labs(y = "Richness", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 30), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))



##### ----- Model EQR trends over time by initial EQR

eqr.trend<-glmmTMB(eqr ~ year_s*eqr.i.site + (1|site_f), family = gaussian, data = env[!is.na(env$eqr), ])
summary(eqr.trend)

# Predicted values

eqr.pred<-predict(eqr.trend, newdata = data.frame("year_s" = year.pred, "eqr.i.site" = mean(env$eqr.i.site, na.rm = TRUE), "site_f" = NA), type = "response", se.fit = TRUE)
eqr.pred<-data.frame(eqr.pred)
eqr.pred02<-predict(eqr.trend, newdata = data.frame("year_s" = year.pred, "eqr.i.site" = 0.2, "site_f" = NA), type = "response", se.fit = TRUE)
eqr.pred04<-predict(eqr.trend, newdata = data.frame("year_s" = year.pred, "eqr.i.site" = 0.4, "site_f" = NA), type = "response", se.fit = TRUE)
eqr.pred06<-predict(eqr.trend, newdata = data.frame("year_s" = year.pred, "eqr.i.site" = 0.6, "site_f" = NA), type = "response", se.fit = TRUE)
eqr.pred08<-predict(eqr.trend, newdata = data.frame("year_s" = year.pred, "eqr.i.site" = 0.8, "site_f" = NA), type = "response", se.fit = TRUE)
eqr.pred10<-predict(eqr.trend, newdata = data.frame("year_s" = year.pred, "eqr.i.site" = 1.0, "site_f" = NA), type = "response", se.fit = TRUE)


eqr.pred<-cbind(eqr.pred, eqr.pred02, eqr.pred04, eqr.pred06, eqr.pred08, eqr.pred10)
colnames(eqr.pred)<-c("fit", "se.fit", "fit.1", "se.fit.1", "fit.2", "se.fit.2", "fit.3", "se.fit.3", "fit.4", "se.fit.4", "fit.5", "se.fit.5")
eqr.pred<-eqr.pred %>% mutate(fit.ci05 = fit - (1.96 * se.fit), fit.ci95 = fit + (1.96 * se.fit), # mean
                              fit1.ci05 = fit.1 - (1.96 * se.fit.1), fit1.ci95 = fit.1 + (1.96 * se.fit.1), # eqr 0.2
                              fit2.ci05 = fit.2 - (1.96 * se.fit.2), fit2.ci95 = fit.2 + (1.96 * se.fit.2), # eqr 0.4
                              fit3.ci05 = fit.3 - (1.96 * se.fit.3), fit3.ci95 = fit.3 + (1.96 * se.fit.3), # eqr 0.6
                              fit4.ci05 = fit.4 - (1.96 * se.fit.4), fit4.ci95 = fit.4 + (1.96 * se.fit.4), # eqr 0.8
                              fit5.ci05 = fit.5 - (1.96 * se.fit.5), fit5.ci95 = fit.5 + (1.96 * se.fit.5)) # eqr 1.0

eqr.pred$year<-year.unscaled
eqr.pred<-eqr.pred %>% relocate(year)

eqr.pred %>% ggplot() +
  #ylim(0,1) +
  geom_line(aes(y = fit.1, x = year), color = "#B22222", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit1.ci05, ymax = fit1.ci95), alpha = 0.05, fill = "#B22222") +
  geom_line(aes(y = fit.2, x = year), color = "#E66100", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit2.ci05, ymax = fit2.ci95), alpha = 0.05, fill = "#E66100") +
  geom_line(aes(y = fit.3, x = year), color = "#FFD700", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit3.ci05, ymax = fit3.ci95), alpha = 0.05, fill = "#FFD700") +
  geom_line(aes(y = fit.4, x = year), color = "#91BFDB", linewidth = 2, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit4.ci05, ymax = fit4.ci95), alpha = 0.05, fill = "#91BFDB") +
  geom_line(aes(y = fit.5, x = year), color = "#00468B", linewidth = 3, linetype = "solid") +
  geom_ribbon(aes(x = year, ymin = fit5.ci05, ymax = fit5.ci95), alpha = 0.05, fill = "#00468B") +
  labs(y = "EQR", x = "Year", title = NULL) +
  theme_bw() +
  theme(text = element_text(size = 30), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_line(linewidth = 0.8))

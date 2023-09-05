library(lme4)
library(lmerTest)
library(sjPlot)
library(tidyverse)
library(data.table)
library(MuMIn)
library(car)

# read in gams
l <- list.files("data/sppSpecific_gamOutputsCSVs/", full.names = T)

gam_df <- l %>% 
  map_df(~fread(.))

gam_df <- gam_df %>% 
  filter(!is.na(fit)) 

# calculate phenoMetrics
phenoMet <- gam_df %>% 
  group_by(Site, validName) %>% 
  summarise(peak = doy[which.max(fit)],
            se = se[which.max(fit)],
            duration = sum(fit > 0.5))

# read in urbanization metrics
urb <- read.csv('data/urbanStressors/urbanization_gradient_SITE.csv') %>% 
  mutate(Site = str_to_title(Site))

phenoMet <- left_join(phenoMet, urb) %>% 
  ungroup()

traits <- read.csv("data/traits/traits_totalSppList.csv") %>% 
  mutate(minWingspan = if_else(notes == "WS",
                               true = minWingspan/2,
                               false = as.double(minWingspan)),
         maxWingspan = if_else(notes == "WS",
                               true = maxWingspan/2,
                               false = as.double(maxWingspan))) 

traits <- traits %>% 
  mutate(test = minWingspan - maxWingspan)

traits <-
  mutate(rowwise(traits),
         medWingspan = median(c(minWingspan, maxWingspan)))

mdf <- phenoMet %>% 
  left_join(traits, by = "validName")

### transform things before 
mdf <- mdf %>% 
  mutate(peak = if_else(condition = peak > 365, true = peak - 365, false = peak))

# read in tempNiche data
tempNiche <- read.csv("data/traits//geographicdata_fixed.csv")

mdf <- mdf %>% 
  left_join(tempNiche, by = c("validName" = "species"))

# read in temp data
temp <- read.csv("data/urbanStressors/temp_gradient.csv")
mdf <- left_join(mdf, temp)

mdf <- mdf %>% 
  mutate(rel_temp = (mean_temp - 1.02) * -1) %>% 
  mutate(rel_hsi = (mean_hsi - 1.33) * -1)


mdf_raw <- mdf

mdf <- mdf %>% 
  mutate(Dev_10 = scale(Dev_10),
         Dev_1 = scale(Dev_1),
         rel_temp = scale(rel_temp),
         max_lat = scale(max_lat),
         maxWingspan = scale(maxWingspan),
         bio1_mean = scale(bio1_mean)) %>% 
  select(peak, Site, validName, Dev_1, Dev_10, rel_temp, voltinism, maxWingspan, bio1_mean)

mdf <- na.omit(mdf)

# predictors = 
# 1)hostPlantSpecialization, 2)max_lat, 3) rel_temp/Dev1/Dev10, 4)voltinism,
# 5) maxWingspan

# species must be found at at least 2 sites
# sites must have at least 2 species
sites <- mdf %>% 
  group_by(Site) %>% 
  summarise(nSpp = length(unique(validName))) ## there is no baca

mdf <- mdf %>% 
  filter(Site != "Baca")

spp <- mdf %>% 
  group_by(validName) %>% 
  summarise(nSites = length(unique(Site))) %>% 
  filter(nSites > 1)

mdf <- mdf %>% 
  filter(validName %in% spp$validName)

mdf <- na.omit(mdf)

# peak modeling
# Dev1
#checking for collinearity among predictors (so no interactions included)
m_col <- lmer(formula = peak ~ Dev_1 + voltinism + maxWingspan + bio1_mean + 
            (1|validName),
          data = mdf, REML = F)
performance::check_collinearity(m_col)

m <- lmer(formula = peak ~ Dev_1 + voltinism + maxWingspan + bio1_mean + 
            Dev_1:voltinism +
            Dev_1:maxWingspan +
            Dev_1:bio1_mean +
            (1|validName),
          data = mdf, REML = F)

performance::check_collinearity(m)
summary(m)
confint(m)

# Dev10 model
m2_col <- lmer(formula = peak ~ Dev_10 + voltinism +  maxWingspan + bio1_mean + 
             (1|validName),
           data = mdf, REML = F)
performance::check_collinearity(m2_col)

m2 <- lmer(formula = peak ~ Dev_10 + voltinism +  maxWingspan + bio1_mean +
             Dev_10:voltinism +
             Dev_10:maxWingspan +
             Dev_10:bio1_mean +
             (1|validName),
           data = mdf, REML = F)

performance::check_collinearity(m2)
summary(m2)
confint(m2)

#rel temp modeling
m3_col <- lmer(formula = peak ~ rel_temp + voltinism  + maxWingspan + bio1_mean + 
             (1|validName),
           data = mdf, REML = F)
performance::check_collinearity(m3_col)

m3 <- lmer(formula = peak ~ rel_temp + voltinism  + maxWingspan + bio1_mean + 
             rel_temp:voltinism +
             rel_temp:maxWingspan +
             rel_temp:bio1_mean +
             (1|validName),
           data = mdf, REML = F)

performance::check_collinearity(m3)
summary(m3)
confint(m3)

AICc(m, m2, m3)# m3 lowest
model.sel(m, m2, m3)
Weights(AICc(m, m2, m3))# m3 clear top model

r.squaredGLMM(m3)

summary(m3)

# sjPlot::tab_model(m3, file = "tabOutputs/spp_specific_LMM.doc")
peak_temp_plot <- plot_model(m3, terms = "rel_temp", type = "pred")
peak_bio1Mean_plot <- plot_model(m3, terms = "bio1_mean", type = "pred")
plot_model(m3, terms = "voltinism", type = "pred")
plot_model(m3, terms = "maxWingspan", type = "pred")
peak_bstmp_plot <- plot_model(m3, terms = c("rel_temp", "maxWingspan [-0.75,0,0.75]"), type = "pred")
peak_int_plot <- plot_model(m3, terms = c("rel_temp", "bio1_mean [-1,0,1]"), type = "pred")

## make manuscript figures
# first for peak
peak_int_plot_df <- peak_int_plot$data
peak_bio1Mean_plot_df <- peak_bio1Mean_plot$data
peak_temp_plot_df <- peak_temp_plot$data
peak_bstmp_plot_df <- peak_bstmp_plot$data

a <- ggplot(peak_temp_plot_df) +
  geom_point(mdf, mapping = aes(x = rel_temp, y = peak), alpha = 0.3) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.25) +
  geom_line(aes(x = x, y = predicted)) +
  labs(y = "Peak DOY", x = "Relative temperature") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))

b <- ggplot(peak_bio1Mean_plot_df) +
  geom_point(mdf, mapping = aes(x = bio1_mean, y = peak), alpha = 0.3) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.25) +
  geom_line(aes(x = x, y = predicted)) +
  labs(y = "Peak DOY", x = "Temperature niche") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))

mdf_plot <- mdf %>% 
  mutate(group = case_when(
    bio1_mean < -0.5 ~ "-1",
    bio1_mean >= -0.5 & bio1_mean <= 0.5 ~ "0",
    bio1_mean > 0.5 ~ "1"
  ))

library(ggnewscale)

c <- ggplot(peak_int_plot_df) +
  geom_point(mdf_plot, mapping = aes(x = rel_temp, y = peak, color = group), alpha = 0.3) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_line(aes(x = x, y = predicted, color = group)) +
  labs(y = "Peak DOY", x = "Relative temperature", 
       fill = "Temp niche", color = "Temp niche") +
  scale_fill_manual(values = c("#2A9D8F", "#D5A220", "#E76F51"),
                    labels = c("-1 (Cold-adapted)", "0 (Average)", "1 (Warm-adapted)")) +
  scale_color_manual(values = c("#2A9D8F", "#D5A220", "#E76F51"),
                     labels = c("-1 (Cold-adapted)", "0 (Average)", "1 (Warm-adapted)")) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16))

c

mdf_plot2 <- mdf %>% 
  mutate(group = case_when(
    maxWingspan < -0.5 ~ "-0.75",
    maxWingspan >= -0.5 & maxWingspan <= 0.5 ~ "0",
    maxWingspan > 0.5 ~ "0.75"
  ))

d <- ggplot(peak_bstmp_plot_df) +
  geom_point(mdf_plot2, mapping = aes(x = rel_temp, y = peak, color = group), alpha = 0.3) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_line(aes(x = x, y = predicted, color = group)) +
  labs(y = "Peak DOY", x = "Relative temperature", 
       fill = "Body size", color = "Body size") +
  scale_fill_manual(values = c("#E56399","#28262C","#DE6E4B"), 
                    labels = c("-0.75 (Small)", "0 (Average)", "0.75 (Large)")) +
  scale_color_manual(values = c("#E56399","#28262C","#DE6E4B"), 
                     labels = c("-0.75 (Small)", "0 (Average)", "0.75 (Large)")) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  theme(legend.position = "bottom",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16))
d

cp2 <- cowplot::plot_grid(c,d, labels = c("A", "B"), nrow = 2, ncol = 1)
ggsave(plot = cp2, filename = "figOutputs/Fig3_sppSpecificIntercations.png", dpi = 450,
       width = 5.5, height = 7)


## write tables for three competing models
m_df <- as.data.frame(summary(m)$coefficients)
m_ci <- as.data.frame(confint(m)[3:10,])
m_df$Lower.CI <- as.vector(m_ci$`2.5 %`)
m_df$Upper.CI <- as.vector(m_ci$`97.5 %`)
m_df <- rownames_to_column(m_df)

write.csv(m_df, "tabOutputs/Dev1_sppSpecificEffects.csv", row.names = F)

m2_df <- as.data.frame(summary(m2)$coefficients)
m2_ci <- as.data.frame(confint(m2)[3:10,])
m2_df$Lower.CI <- as.vector(m2_ci$`2.5 %`)
m2_df$Upper.CI <- as.vector(m2_ci$`97.5 %`)
m2_df <- rownames_to_column(m2_df)


write.csv(m2_df, "tabOutputs/Dev10_sppSpecificEffects.csv", row.names = F)

m3_df <- as.data.frame(summary(m3)$coefficients)
m3_ci <- as.data.frame(confint(m3)[3:10,])
m3_df$Lower.CI <- as.vector(m3_ci$`2.5 %`)
m3_df$Upper.CI <- as.vector(m3_ci$`97.5 %`)
m3_df <- rownames_to_column(m3_df)

write.csv(m3_df, "tabOutputs/relTemp_sppSpecificEffects.csv", row.names = F)

r.squaredGLMM(m3)

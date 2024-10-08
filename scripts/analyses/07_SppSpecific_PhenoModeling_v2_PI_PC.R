library(lme4)
library(lmerTest)
library(sjPlot)
library(tidyverse)
library(data.table)
library(MuMIn)
library(car)

# read csv with site by species combinations with GAMs with properly fit GAMs that have clear peaks
usableSppSites <- read.csv("data/phenoData/usableSppBySites.csv")

# read in gams
l_nb <- list.files("data/sppSpecific_gamOutputsCSVs_negBinom/", full.names = T)

gam_df_nb <- l_nb %>% 
  map_df(~read_csv(.x)) 

gam_df_nb <- gam_df_nb %>% 
  filter(!is.na(fit)) %>% 
  mutate(SiteSpp = paste(Site, validName, sep = "_"))

## filter to SiteSpp 
gam_df_nb <- gam_df_nb %>% 
  filter(SiteSpp %in% usableSppSites$SiteSpp)

# re-assign gam_df_nb as gam_df
gam_df <- gam_df_nb

# calculate phenoMetrics
phenoMet <- gam_df %>% 
  group_by(Site, validName) %>% 
  summarise(peak = doy[which.max(fit)],
            se = se[which.max(fit)],
            duration = sum(fit > 0.5))

# read in urbanization metrics
urb <- read.csv('gisData/landcover/imperviousSurface.csv') %>% 
  mutate(Site = str_to_title(Site))
can <- read.csv('gisData/landcover/canopySurface.csv') %>% 
  mutate(Site = str_to_title(Site))

phenoMet <- left_join(phenoMet, urb) %>% 
  left_join(can) %>% 
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

### transform things before modeling
mdf <- mdf %>% 
  mutate(peak = if_else(condition = peak > 365, true = peak - 365, false = peak))

# read in tempNiche data
tempNiche <- read.csv("data/traits//geographicdata_fixed.csv")

mdf <- mdf %>% 
  left_join(tempNiche, by = c("validName" = "species"))

mdf_raw <- mdf

mdf <- mdf %>% 
  mutate(propImpervious300m = scale(propImpervious300m),
         propCanopy300m = scale(propCanopy300m),
         maxWingspan = scale(maxWingspan),
         bio1_mean = scale(bio1_mean)) %>% 
  select(peak, Site, validName, propImpervious300m, propCanopy300m, maxWingspan, bio1_mean)

mdf <- na.omit(mdf)

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
# Impervious surface
#checking for collinearity among predictors (so no interactions included)
m_col <- lmer(formula = peak ~ propImpervious300m + maxWingspan + bio1_mean + 
                (1|validName),
              data = mdf, REML = F)
performance::check_collinearity(m_col)

m <- lmer(formula = peak ~ propImpervious300m + maxWingspan + bio1_mean + 
            propImpervious300m:maxWingspan +
            propImpervious300m:bio1_mean +
            (1|validName),
          data = mdf, REML = F)

performance::check_collinearity(m)
summary(m)
confint(m)

#canopy cover
m3_col <- lmer(formula = peak ~ propCanopy300m  + maxWingspan + bio1_mean + 
                 (1|validName),
               data = mdf, REML = F)
performance::check_collinearity(m3_col)

m3 <- lmer(formula = peak ~ propCanopy300m  + maxWingspan + bio1_mean + 
             propCanopy300m:maxWingspan +
             propCanopy300m:bio1_mean +
             (1|validName),
           data = mdf, REML = F)

performance::check_collinearity(m3)
summary(m3)
confint(m3)

r.squaredGLMM(m3)

summary(m3)

# sjPlot::tab_model(m3, file = "tabOutputs/spp_specific_LMM.doc")
# SJ PLOT MODEL PREDICTIONS FOR Impervious surface
PI_BS_plot <- plot_model(m, terms = c("propImpervious300m", "maxWingspan [-1,0,1]"), type = "pred")
PI_TN_plot <- plot_model(m, terms = c("propImpervious300m", "bio1_mean [-1,0,1]"), type = "pred")
PI_BS_plot_df <- PI_BS_plot$data
PI_TN_plot_df <- PI_TN_plot$data

CP_BS_plot <- plot_model(m3, terms = c("propCanopy300m", "maxWingspan [-1,0,1]"), type = "pred")
CP_TN_plot <- plot_model(m3, terms = c("propCanopy300m", "bio1_mean [-1,0,1]"), type = "pred")
CP_BS_plot_df <- CP_BS_plot$data
CP_TN_plot_df <- CP_TN_plot$data
## make manuscript figures
# Impervious surface
mdf_plot <- mdf %>% 
  mutate(group = case_when(
    bio1_mean < -0.5 ~ "-1",
    bio1_mean >= -0.5 & bio1_mean <= 0.5 ~ "0",
    bio1_mean > 0.5 ~ "1"
  ))

mdf_plot2 <- mdf %>% 
  mutate(group = case_when(
    maxWingspan < -0.5 ~ "-0.75",
    maxWingspan >= -0.5 & maxWingspan <= 0.5 ~ "0",
    maxWingspan > 0.5 ~ "0.75"
  ))

library(ggnewscale)

### backtransform variables for easier interpretation 
mean_propImpervious <- mean(mdf_raw$propImpervious300m)
sd_propImpervious <- sd(mdf_raw$propImpervious300m)

mean_propCanopy <- mean(mdf_raw$propCanopy300m)
sd_propCanopy <- sd(mdf_raw$propCanopy300m)

mean(mdf_raw$maxWingspan)
sd(mdf_raw$maxWingspan)


anova(re_only, m)
anova(re_only, m3)

mean(mdf_raw$bio1_mean, na.rm = T)
mean(mdf_raw$bio1_mean, na.rm = T) + sd(mdf_raw$bio1_mean, na.rm = T)
mean(mdf_raw$bio1_mean, na.rm = T) - sd(mdf_raw$bio1_mean, na.rm = T)

a <- ggplot(PI_TN_plot_df) +
  geom_jitter(mdf_plot, mapping = aes(x = mean_propImpervious + propImpervious300m * sd_propImpervious,
                                      y = peak), alpha = 0.3) + 
  geom_ribbon(aes(x = mean_propImpervious + x * sd_propImpervious, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_line(aes(x = mean_propImpervious + x * sd_propImpervious, y = predicted, color = group)) +
  labs(y = "Peak DOY", x = "Proportion impervious surface (300 m)", 
       fill = "Temp niche", color = "Temp niche") +
  scale_fill_manual(values = c("#2A9D8F", "#D5A220", "#E76F51"),
                    labels = c( paste0("10.9", "\u00B0", "C", "\n", "(Cold-adapted)"), 
                                paste0("13.9", "\u00B0", "C", "\n", "(Average)"),
                                paste0("16.8", "\u00B0", "C", "\n", "(Warm-Adapted)"))) +
  scale_color_manual(values = c("#2A9D8F", "#D5A220", "#E76F51"),
                     labels = c( paste0("10.9", "\u00B0", "C", "\n", "(Cold-adapted)"), 
                                 paste0("13.9", "\u00B0", "C", "\n", "(Average)"),
                                 paste0("16.8", "\u00B0", "C", "\n", "(Warm-Adapted)")))  +
  theme_classic() +
  theme(legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16))
a

b <- ggplot(PI_BS_plot_df) +
  geom_jitter(mdf_plot2, mapping = aes(x = mean_propImpervious + propImpervious300m * sd_propImpervious,
                                       y = peak), alpha = 0.3) + 
  geom_ribbon(aes(x = mean_propImpervious + x * sd_propImpervious,
                  ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_line(aes(x = mean_propImpervious + x * sd_propImpervious,
                y = predicted, color = group)) +
  labs(y = "Peak DOY", x = "Proportion impervious surface (300 m)", 
       fill = "Body size", color = "Body size") +
  scale_fill_manual(values = c("#E56399","#28262C","#DE6E4B"), 
                    labels = c("14.7mm\n(Small)", "22.6mm\n(Average)", "30.6mm\n(Large)")) +
  scale_color_manual(values = c("#E56399","#28262C","#DE6E4B"), 
                     labels = c("14.7mm (Small)", "22.6mm\n(Average)", "30.6mm\n(Large)")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16))
b

# Canopy cover

c <- ggplot(CP_TN_plot_df) +
  geom_jitter(mdf_plot, mapping = aes(x = mean_propCanopy + propCanopy300m * sd_propCanopy,
                                      y = peak), alpha = 0.3) + 
  geom_ribbon(aes(x = mean_propCanopy + x * sd_propCanopy, 
                  ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_line(aes(x = mean_propCanopy + x * sd_propCanopy, y = predicted, color = group)) +
  labs(y = "Peak DOY", x = "Proportion canopy cover (300 m)", 
       fill = "Temp niche", color = "Temp niche") +
  scale_fill_manual(values = c("#2A9D8F", "#D5A220", "#E76F51"),
                    labels = c( paste0("10.9", "\u00B0", "C", "\n", "(Cold\nadapted)"), 
                                paste0("13.9", "\u00B0", "C", "\n", "(Average)"),
                                paste0("16.8", "\u00B0", "C", "\n", "(Warm\nAdapted)"))) +
  scale_color_manual(values = c("#2A9D8F", "#D5A220", "#E76F51"),
                     labels = c( paste0("10.9", "\u00B0", "C", "\n", "(Cold\nadapted)"), 
                                 paste0("13.9", "\u00B0", "C", "\n", "(Average)"),
                                 paste0("16.8", "\u00B0", "C", "\n", "(Warm\nAdapted)")))  +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11,
                                    face = "bold"))

c

d <- ggplot(CP_BS_plot_df) +
  geom_jitter(mdf_plot2, mapping = aes(x = mean_propCanopy + propCanopy300m * sd_propCanopy,
                                       y = peak), alpha = 0.3) + 
  geom_ribbon(aes(x = mean_propCanopy + x * sd_propCanopy,
                  ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_line(aes(x = mean_propCanopy + x * sd_propCanopy,
                y = predicted, color = group)) +
  labs(y = "Peak DOY", x = "Proportion canopy cover (300 m)", 
       fill = "Body size", color = "Body size") +
  scale_fill_manual(values = c("#E56399","#28262C","#DE6E4B"), 
                    labels = c("14.7mm\n(Small)", "22.6mm\n(Average)", "30.6mm\n(Large)")) +
  scale_color_manual(values = c("#E56399","#28262C","#DE6E4B"), 
                     labels = c("14.7mm\n(Small)", "22.6mm\n(Average)", "30.6mm\n(Large)")) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, 
                                    face = "bold"))
d

cp2 <- cowplot::plot_grid(a,b,c,d, labels = c("A", "B", "C", "D"), 
                          nrow = 2, ncol = 2)
ggsave(plot = cp2, filename = "figOutputs/resubmission/resubmissionSept/Fig2_sppSpecificIntercations_v2_negBinom.png", 
       dpi = 450,
       width = 8.3, height = 8.3)


## write tables for three competing models
m_df <- as.data.frame(summary(m)$coefficients)
m_ci <- as.data.frame(confint(m)[3:8,])
m_df <- m_df %>% 
  mutate(Lower_CI = m_ci$`2.5 %`,
         Upper_CI = m_ci$`97.5 %`)

write.csv(m_df, "tabOutputs/propImperviousEffects_sppSpecific_negBinom.csv", row.names = F)



m3_df <- as.data.frame(summary(m3)$coefficients)
m3_ci <- as.data.frame(confint(m3)[3:8,])
m3_df$Lower.CI <- as.vector(m3_ci$`2.5 %`)
m3_df$Upper.CI <- as.vector(m3_ci$`97.5 %`)
m3_df <- rownames_to_column(m3_df)

write.csv(m3_df, "tabOutputs/propCanopyEffects_sppSpecific_negBinom.csv", row.names = F)

r.squaredGLMM(m3)

AICc(m, m3)# m3 lowest
model.sel(m, m3)
Weights(AICc(m, m3))# m3 clear top model

r.squaredGLMM(m3)

r.squaredGLMM(m)
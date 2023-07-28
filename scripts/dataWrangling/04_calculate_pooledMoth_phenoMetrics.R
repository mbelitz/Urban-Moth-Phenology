library(dplyr)
library(ggplot2)

# read in phenology gams
frass_gams <- read.csv("data/data_products/frassGams_perDay_randomEffect.csv") %>% 
  mutate(unit = "Frass") %>% 
  mutate(rel_abundance = (fit * -1) / (max(fit) * -1))
macro_gams <- read.csv("data/data_products/macroMothGams.csv") %>% 
  mutate(unit = "Macro-moths") %>% 
  mutate(rel_abundance = fit / max(fit))
micro_gams <- read.csv('data/data_products/microMothGams.csv') %>% 
  mutate(unit = "Micro-moths") %>% 
  mutate(rel_abundance = fit / max(fit))

## plot all three on top of each other
tdf <- rbind(frass_gams, macro_gams, micro_gams)

ggplot(tdf, mapping = aes(x = doy2, y = fit, color = unit)) +
  geom_smooth() +
  facet_wrap(~Site) +
  theme_bw()


# grab the peak dates
# start with macro moths where there is always only one peak and it will be before doy 300
macroPhenoMet <- macro_gams %>% 
  filter(doy2 < 300) %>% 
  group_by(Site) %>% 
  summarise(peak = doy2[which.max(fit)],
            se = se[which.max(fit)])

# now do this for micro moths, where we want peak1 and peak 2
# peak 1 is before <= 200 & peak 2 is after 200
microPhenoMet1 <- filter(micro_gams, doy2 <= 200) %>% 
  group_by(Site) %>% 
  summarise(peak1 = doy2[which.max(fit)],
            se1 = se[which.max(fit)])
microPhenoMet2 <- filter(micro_gams, doy2 > 200 & doy2 <= 300) %>% 
  group_by(Site) %>% 
  summarise(peak2 = doy2[which.max(fit)],
            se2 = se[which.max(fit)])
microPhenoMet <- left_join(microPhenoMet1,microPhenoMet2)


# now onto frass, 2 peaks again
# # peak 1 is before <= 172 & peak 2 is after 172
# RIST and DEMI don't have first peak -- interesting b/c these are the most 
# pine dominated forests
frassPhenoMet1 <- filter(frass_gams, doy2 <= 172) %>% 
  group_by(Site) %>% 
  summarise(peak1 = doy2[which.max(fit)],
            se1 = se[which.max(fit)]) %>% 
  filter(Site != "RIST",
         Site != "DEMI")
frassPhenoMet2 <- filter(frass_gams, doy2 > 172) %>% 
  group_by(Site) %>% 
  summarise(peak2 = doy2[which.max(fit)],
            se2 = se[which.max(fit)])
frassPhenoMet <- left_join(frassPhenoMet2, frassPhenoMet1) 

# write df of phenometrics
write.csv(x = microPhenoMet,
          file = "data/phenoData/microMothsPhenoMetrics.csv", row.names = F)
write.csv(x = macroPhenoMet,
          file = "data/phenoData/macroMothsPhenoMetrics.csv", row.names = F)
write.csv(x = frassPhenoMet,
          file = "data/phenoData/frassPhenoMetrics_perDay_trapRE.csv", row.names = F)
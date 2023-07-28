library(dplyr)
library(ggplot2)
library(lubridate)
library(mgcv)

# read in adult data
adult_df <- read.csv("data/field_observations/surveyDateSheet.csv") %>% 
  mutate(doy = yday(mdy(eventDate)),
         year = year(mdy(eventDate))) %>% 
  mutate(doy2 = if_else(condition = year == 2019,
                        true = doy,
                        false = doy + 365))

adult_gb <- adult_df %>% 
  group_by(doy2, location) %>% 
  summarise(macroMoths = sum(macroMoths),
            microMoths = sum(microMoths))

ggplot(adult_gb, aes(x = doy2, y = microMoths, color = location)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 24))

# read in lunar illumination
lunar.phase <- read.csv("data/urbanStressors/lunarIllumination.csv") %>% 
  mutate(Date = mdy(Date)) %>% 
  mutate(doy = yday(Date),
         year = year(Date)) %>% 
  mutate(doy2 = if_else(condition = year == 2019,
                        true = doy,
                        false = doy + 365))

adult_gb <- left_join(adult_gb, lunar.phase, by = "doy2")

# get adult df for each location
auca_adult_df <- filter(adult_gb, location == "AUCA")
rist_adult_df <- filter(adult_gb, location == "RIST")
prcr_adult_df <- filter(adult_gb, location == "PRCR")

bowa_adult_df <- filter(adult_gb, location == "BOWA")
biva_adult_df <- filter(adult_gb, location == "BIVA")
demi_adult_df <- filter(adult_gb, location == "DEMI")

cofr_adult_df <- filter(adult_gb, location == "COFR")
baca_adult_df <- filter(adult_gb, location == "BACA")
joma_adult_df <- filter(adult_gb, location == "JOMA")


# function for gam for each location
gam_pred <- function(site_df){
  
  g <- gam(microMoths ~ s(doy, k = 6, bs = "cr") + s(lunar.phase, k = 3, bs = "cr"),  
           data = site_df, family = poisson)
  
  g_pred <- predict.gam(g, newdata=data.frame(doy=unique(60:424),
                                              lunar.phase = rep(0.5, length(unique(60:424)))), 
                        type="response", se.fit = T)
  
  g_pred_df <- data.frame(doy2 = unique(60:424),
                          fit = g_pred$fit,
                          se = g_pred$se.fit)
  g_pred_df <- g_pred_df %>% 
    mutate(fit = ifelse(fit < 0, 0, fit))
  
  
}

# AUCA # 
auca_gam_pred <- gam_pred(site_df = auca_adult_df) %>% 
  mutate(Site = "AUCA")

auca_plot <- ggplot() +
  geom_point(data = auca_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = auca_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = auca_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Rural site: Austin Cary Forest (AUCA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))
auca_plot

# RIST #
rist_gam_pred <- gam_pred(site_df = rist_adult_df) %>% 
  mutate(Site = "RIST")

rist_plot <- ggplot() +
  geom_point(data = rist_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = rist_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = rist_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Rural site: Longleaf flatwoods reserve (RIST)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))
rist_plot

# PRCR #
prcr_gam_pred <- gam_pred(site_df = prcr_adult_df) %>% 
  mutate(Site = "PRCR")

prcr_plot <- ggplot() +
  geom_point(data = prcr_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = prcr_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = prcr_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Rural site: Prairie Creek Preserve (PRCR)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))
prcr_plot

# BOWA #
bowa_gam_pred <- gam_pred(site_df = bowa_adult_df) %>% 
  mutate(Site = "BOWA")

bowa_plot <- ggplot() +
  geom_point(data = bowa_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = bowa_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = bowa_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Suburban site: Boulware springs park (BOWA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))
bowa_plot

# biva #
biva_gam_pred <- gam_pred(site_df = biva_adult_df) %>% 
  mutate(Site = "BIVA")

biva_plot <- ggplot() +
  geom_point(data = biva_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = biva_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = biva_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Suburban site: Bivens Arm Nature Park (BIVA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))
biva_plot


# demi #
demi_gam_pred <- gam_pred(site_df = demi_adult_df) %>% 
  mutate(Site = "DEMI")

demi_plot <- ggplot() +
  geom_point(data = demi_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = demi_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = demi_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Suburban site: Devil's Millhopper State Park (DEMI)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

demi_plot

# cofr #
cofr_gam_pred <- gam_pred(site_df = cofr_adult_df) %>% 
  mutate(Site = "COFR")

cofr_plot <- ggplot() +
  geom_point(data = cofr_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = cofr_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = cofr_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Urban site: Cofrin Nature Park (COFR)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))
cofr_plot

# joma #
joma_gam_pred <- gam_pred(site_df = joma_adult_df) %>% 
  mutate(Site = "JOMA")

joma_plot <- ggplot() +
  geom_point(data = joma_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = joma_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = joma_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Urban site: John Mahon Nature Park (JOMA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))
joma_plot

# baca #
baca_gam_pred <- gam_pred(site_df = baca_adult_df) %>% 
  mutate(Site = "BACA")

baca_plot <- ggplot() +
  geom_point(data = baca_adult_df, aes(x = doy2, y = microMoths)) +
  geom_line(data = baca_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = baca_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Micro-moths") +
  ggtitle("Urban site: Bartram/Carr Woods (BACA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))
baca_plot


## plot them all together
cp <- cowplot::plot_grid(baca_plot, cofr_plot, joma_plot,
                         demi_plot, bowa_plot, biva_plot,
                         prcr_plot, rist_plot, auca_plot, 
                         nrow = 3, ncol = 3)
cp

ggsave(plot = cp, filename = "figOutputs/microMothGams.png", dpi = 450,
       width = 12, height = 8)


# grab phenometrics
# first let's make a dataframe of everything
tdf <- rbind(auca_gam_pred, prcr_gam_pred, rist_gam_pred,
             demi_gam_pred, bowa_gam_pred, biva_gam_pred,
             baca_gam_pred, cofr_gam_pred, joma_gam_pred)


write.csv(tdf, file = "data/gamOutputs/microMothGams.csv", row.names = F)
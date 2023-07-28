library(dplyr)
library(ggplot2)
library(lubridate)
library(mgcv)


frass_df <- read.csv("data/data_products/frass_biomass_reweigh_doyDiff.csv")

frass_df <- frass_df %>% 
  mutate(MassPerDay = Mass/diff)

frass_df <- frass_df %>% 
  mutate(MassPerDay = log(MassPerDay + 0.001))

ggplot(frass_df, aes(x = doy2, y = MassPerDay, color = as.factor(Trap))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 36),
              se = F) +
  theme_classic() +
  facet_wrap(~Site)

frass_df$Trap <- as.factor(frass_df$Trap)

# get frass df for each site
auca_frass_df <- filter(frass_df, Site == "AUCA")
rist_frass_df <- filter(frass_df, Site == "RIST")
prcr_frass_df <- filter(frass_df, Site == "PRCR")

bowa_frass_df <- filter(frass_df, Site == "BOWA")
biva_frass_df <- filter(frass_df, Site == "BIVA")
demi_frass_df <- filter(frass_df, Site == "DEMI")

cofr_frass_df <- filter(frass_df, Site == "COFR")
baca_frass_df <- filter(frass_df, Site == "BACA")
joma_frass_df <- filter(frass_df, Site == "JOMA")

# gam function for each location
gam_pred <- function(site_df, k){
  
  g <- gam(MassPerDay ~ s(doy2, bs = "cr", k = k) + s(Trap, bs = "re"),
           data = site_df)
  
  g_pred <- predict.gam(g, newdata=data.frame(doy2=unique(60:424), Trap = 6), 
                        type="response", se.fit = T, re = NA)
  
  g_pred_df <- data.frame(doy2 = unique(60:424),
                          fit = g_pred$fit,
                          se = g_pred$se.fit)
}
# AUCA # 
auca_gam_pred <- gam_pred(site_df = auca_frass_df, k = 8) %>% 
  mutate(Site = "AUCA")

auca_plot <- ggplot() +
  geom_point(data = auca_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = auca_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = auca_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass (mg/day)") +
  ggtitle("Rural site: Austin Cary Forest (AUCA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

auca_plot

# RIST #
rist_gam_pred <- gam_pred(site_df = rist_frass_df, k = 8) %>% 
  mutate(Site = "RIST")

rist_plot <- ggplot() +
  geom_point(data = rist_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = rist_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = rist_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass") +
  ggtitle("Rural site: Longleaf flatwoods reserve (RIST)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

rist_plot

# PRCR #
prcr_gam_pred <- gam_pred(site_df = prcr_frass_df, k = 8) %>% 
  mutate(Site = "PRCR")

prcr_plot <- ggplot() +
  geom_point(data = prcr_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = prcr_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = prcr_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass") +
  ggtitle("Rural site: Prairie Creek Preserve (PRCR)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

prcr_plot


# BOWA #
bowa_gam_pred <- gam_pred(site_df = bowa_frass_df, k = 9) %>% 
  mutate(Site = "BOWA")

bowa_plot <- ggplot() +
  geom_point(data = bowa_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = bowa_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = bowa_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass") +
  ggtitle("Suburban site: Boulware springs park (BOWA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

bowa_plot

# biva #
biva_gam_pred <- gam_pred(site_df = biva_frass_df, k = 8) %>% 
  mutate(Site = "BIVA")

biva_plot <- ggplot() +
  geom_point(data = biva_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = biva_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = biva_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass") +
  ggtitle("Suburban site: Bivens Arm Nature Park (BIVA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

biva_plot

# demi #
demi_gam_pred <- gam_pred(site_df = demi_frass_df, k = 8) %>% 
  mutate(Site = "DEMI")

demi_plot <- ggplot() +
  geom_point(data = demi_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = demi_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = demi_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass") +
  ggtitle("Suburban site: Devil's Millhopper State Park (DEMI)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

demi_plot


# cofr #
cofr_gam_pred <- gam_pred(site_df = cofr_frass_df, k = 8) %>% 
  mutate(Site = "COFR")

cofr_plot <- ggplot() +
  geom_point(data = cofr_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = cofr_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = cofr_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass") +
  ggtitle("Urban site: Cofrin Nature Park (COFR)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

cofr_plot


# joma #
joma_gam_pred <- gam_pred(site_df = joma_frass_df, k = 8) %>% 
  mutate(Site = "JOMA")

joma_plot <- ggplot() +
  geom_point(data = joma_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = joma_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = joma_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass") +
  ggtitle("Urban site: John Mahon Nature Park (JOMA)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"))

joma_plot


# baca #
baca_gam_pred <- gam_pred(site_df = baca_frass_df, k = 8) %>% 
  mutate(Site = "BACA")

baca_plot <- ggplot() +
  geom_point(data = baca_frass_df, aes(x = doy2, y = MassPerDay)) +
  geom_line(data = baca_gam_pred, aes(x = doy2, y = fit)) +
  geom_ribbon(data = baca_gam_pred, aes(x = doy2, ymin = fit - se, ymax = fit + se),
              alpha = 0.2) +
  labs(x = "DOY", y = "Frass") +
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

ggsave(plot = cp, filename = "figOutputs/ch3/frassGams_perDay_randomEffect.png", dpi = 450,
       width = 12, height = 8)


# grab phenometrics
# first let's make a dataframe of everything
tdf <- rbind(auca_gam_pred, prcr_gam_pred, rist_gam_pred,
             demi_gam_pred, bowa_gam_pred, biva_gam_pred,
             baca_gam_pred, cofr_gam_pred, joma_gam_pred)
write.csv(tdf, file = "data/data_products/frassGams_perDay_randomEffect.csv", row.names = F)
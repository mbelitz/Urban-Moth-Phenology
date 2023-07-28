library(dplyr)
library(ggplot2)
library(lubridate)
library(mgcv)
library(stringr)

# read in adult dataset
moth_df <- read.csv("data/field_observations/adultDataSet_validNames.csv") %>% 
  distinct(id, .keep_all = T) %>% # this removes duplicate photos of the same individual
  mutate(year = year(mdy(eventDate))) %>% 
  mutate(doy = if_else(
    condition = year == 2019,
    true = yday(mdy(eventDate)),
    false = 365 + yday(mdy(eventDate))
  ))

# how many species do we have enough data to estimate gams for?
enoughSites <- moth_df %>% 
  group_by(validName) %>% 
  summarise(totalObs = n(), nSites = length(unique(Site))) %>% 
  filter(nSites >= 2) %>% 
  filter(totalObs > nSites)

# find validNames with NAs in them
enoughSites <- enoughSites %>% 
  mutate(hasNA = str_detect(validName, "NA")) 
enoughSites <- enoughSites %>% 
  filter(hasNA == F)

#filter initial dataset to these species
moth_df_filter <- moth_df %>% 
  filter(validName %in% enoughSites$validName) %>% 
  filter(!is.na(doy))

# what days had surveys for each site
surveyDates <- read.csv("data/field_observations/surveyDateSheet.csv") %>% 
  select(eventDate, location)

surveyDates$eventDate <- lubridate::mdy(surveyDates$eventDate)
surveyDates <- surveyDates %>% 
  rename(Site = location)
surveyDates <- surveyDates %>% 
  mutate(Site = case_when(Site == "BACA" ~ "Baca",
                          Site == "JOMA" ~ "Joma",
                          Site == "AUCA" ~ "Auca",
                          Site == "BIVA" ~ "Biva",
                          Site == "DEMI" ~ "Demi",
                          Site == "COFR" ~ "Cofr",
                          Site == "RIST" ~ "Rist",
                          Site == "PRCR" ~ "Prcr",
                          Site == "BOWA" ~ "Bowa"))

# don't forget lunar data
lunar.phase <- read.csv("data/urbanStressors/lunarIllumination.csv") %>% 
  mutate(Date = mdy(Date))
surveyDates <- left_join(surveyDates, lunar.phase, by = c("eventDate" = "Date"))
surveyDates <- surveyDates %>% 
  mutate(year = year(eventDate)) %>% 
  mutate(doy = if_else(
    condition = year == 2019,
    true = yday(eventDate),
    false = 365 + yday(eventDate)
  )) %>% 
  select(doy, lunar.phase) %>% 
  distinct(doy, .keep_all = T)

# gam function
gam_pred <- function(site, spp, k = 6){
  
  spp_df <- moth_df_filter %>% 
    filter(validName == spp,
           Site == site)
  
  mdf <- spp_df %>%
    group_by(doy) %>% 
    summarise(count = n()) 
  
  mdf <- left_join(surveyDates, mdf, by = "doy") %>% 
    mutate(count = if_else(
      is.na(count), 
      true = 0, 
      false = as.double(count)))%>% 
    mutate(validName = spp)
  
  g <- gam(count ~ s(doy, k = k, bs = "cr") + s(lunar.phase, k = 3, bs = "cr"),  
           data = mdf, family = poisson)
  
  g_pred <- predict.gam(g, newdata=data.frame(doy=unique(69:424),
                                              lunar.phase = rep(0.2, length(unique(69:424)))), 
                        type="response", se.fit = T)
  
  g_pred_df <- data.frame(doy = unique(69:424),
                          fit = g_pred$fit,
                          se = g_pred$se.fit)
  g_pred_df <- g_pred_df %>% 
    mutate(fit = ifelse(fit < 0, 0, fit)) %>% 
    mutate(Site = site,
           validName = spp)
  
  
}

save_gam_csv <- function(g){
  
  if(nrow(filter(g, fit > 0))){
    
    fp_csvs <- file.path('data/sppSpecific_gamOutputsCSVs/')
    fp_figures <- file.path('figOutputs/gams/')
    
    write.csv(x = g, file = paste0(fp_csvs, '/',"GAM_", 
                                   g$Site[1], "_", g$validName[1],
                                   ".csv"), row.names = F)
    
  } else(print("no data for this species"))
  
}

plot_gam <- function(spp, site, g) {
  
  if(nrow(filter(g, fit > 0))) {
    
    rdf <- moth_df_filter %>% 
      filter(validName == spp,
             Site == site) %>% 
      group_by(doy) %>% 
      summarise(count = n()) %>% 
      mutate(validName = spp,
             Site = site)
    
    p <- ggplot() +
      geom_point(rdf, mapping = aes(x = doy, y = count)) +
      geom_line(g, mapping = aes(x = doy, y = fit)) +
      labs(x = "DOY", y = "Abundance") +
      theme_bw() +
      ggtitle(site)
    
  } else {
    p <- ggplot() +
      ggtitle(site) +
      theme_void() }
}

## lapply function
gam_pipeline <- function(binomial){
  
  baca_gam <- gam_pred(site = "Baca", spp = binomial, k = 6)
  save_gam_csv(g = baca_gam)
  baca_plot <- plot_gam(spp = binomial, site = "Baca", g = baca_gam)
  
  joma_gam <- gam_pred(site = "Joma", spp = binomial, k = 6)
  save_gam_csv(g = joma_gam)
  joma_plot <- plot_gam(spp = binomial, site = "Joma", g = joma_gam)
  
  cofr_gam <- gam_pred(site = "Cofr", spp = binomial, k = 6)
  save_gam_csv(g = cofr_gam)
  cofr_plot <- plot_gam(spp = binomial, site = "Cofr", g = cofr_gam)
  
  biva_gam <- gam_pred(site = "Biva", spp = binomial, k = 6)
  save_gam_csv(g = biva_gam)
  biva_plot <- plot_gam(spp = binomial, site = "Biva", g = biva_gam)
  
  bowa_gam <- gam_pred(site = "Bowa", spp = binomial, k = 6)
  save_gam_csv(g = bowa_gam)
  bowa_plot <- plot_gam(spp = binomial, site = "Bowa", g = bowa_gam)
  
  demi_gam <- gam_pred(site = "Demi", spp = binomial, k = 6)
  save_gam_csv(g = demi_gam)
  demi_plot <- plot_gam(spp = binomial, site = "Demi", g = demi_gam)
  
  auca_gam <- gam_pred(site = "Auca", spp = binomial, k = 6)
  save_gam_csv(g = auca_gam) 
  auca_plot <- plot_gam(spp = binomial, site = "Auca", g = auca_gam)
  
  prcr_gam <- gam_pred(site = "Prcr", spp = binomial, k = 6)
  save_gam_csv(g = prcr_gam)
  prcr_plot <- plot_gam(spp = binomial, site = "Prcr", g = prcr_gam)
  
  rist_gam <- gam_pred(site = "Rist", spp = binomial, k = 6)
  save_gam_csv(g = rist_gam)
  rist_plot <- plot_gam(spp = binomial, site = "Rist", g = rist_gam)
  
  cp <- cowplot::plot_grid(baca_plot, joma_plot, cofr_plot,
                           biva_plot, bowa_plot, demi_plot,
                           rist_plot, auca_plot, prcr_plot)
  
  bw <- stringr::str_replace(binomial, " ", "_")
  ggsave(filename = paste0("figOutputs/gamsPoisson/", bw, ".png"), plot = cp, 
         width = 7, height = 7)
  
}

## get list of species with enough data and start working through these
spp_list <- unique(enoughSites$validName)
spp_list

# This part of the workflow was not 100% automated 
# What we did was 1) fit a GAM for every species in the list [1:141]
# 2) examined the predicted outputs and removed sites from output where either there
# were not enough data to fit the GAM or the GAM did not produce an unambigious peak

# Here's the example code for the first species on the list
# gam_pipeline(binomial = spp_list[1]) # remove Auca, Prcr, Baca

# on the GitHub we will release the species X site GAMs that passed our test

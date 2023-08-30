library(dplyr)
library(ggplot2)
library(sjPlot)
library(MuMIn)
library(lme4)

# read in frass pheno data
frass <- read.csv("data/phenoData/frassPhenoMetrics_perDay_trapRE.csv")
micro <- read.csv("data/phenoData/microMothsPhenoMetrics.csv")
macro <- read.csv("data/phenoData/macroMothsPhenoMetrics.csv")

#join with urbanization data
urb <- read.csv("data/urbanStressors/urbanization_gradient_SITE.csv")
temp <- read.csv('data/urbanStressors/temp_gradient_SITE.csv')
light <- read.csv('data/urbanStressors/lightData_SITE.csv')

frass <- frass %>% 
  left_join(urb) %>% 
  left_join(temp) %>% 
  left_join(light) %>% 
  mutate(rel_temp = (mean_temp - 1.02)*-1)

micro <- micro %>% 
  left_join(urb) %>% 
  left_join(temp) %>% 
  left_join(light) %>% 
  mutate(rel_temp = (mean_temp - 1.02)*-1)

macro <- macro %>% 
  left_join(urb) %>% 
  left_join(temp) %>% 
  left_join(light) %>% 
  mutate(rel_temp = (mean_temp - 1.02)*-1)%>% 
  mutate(meanLight = log(meanLight + 1))

# make model for frass, 4 predictors are dev 10, dev 1, light and temp
frassPeak1_dev1 <- lm(formula = peak1 ~ Dev_1,
                      data = frass)
frassPeak1_dev10 <- lm(formula = peak1 ~ Dev_10,
                       data = frass)
frassPeak1_temp <- lm(formula = peak1 ~ rel_temp,
                      data = frass)

Weights(AICc(frassPeak1_dev1, frassPeak1_dev10, frassPeak1_temp))
model.sel(frassPeak1_dev1, frassPeak1_dev10, frassPeak1_temp)
# get BH correction for three outputs
summary(frassPeak1_dev1) #7.138
confint(frassPeak1_dev1, 'Dev_1', level=0.95) #-9.5 - 23.78
summary(frassPeak1_dev10) 
summary(frassPeak1_temp)

frass_plot1 <- ggplot() +
  geom_point(frass, mapping = aes(x = Dev_1, y = peak1)) +
  geom_line(frass_plot_df1, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(frass_plot_df1, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Proportion developed (1-km)", y = "DOY: First peak") +
  ggtitle("Frass") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))
frass_plot1

# now for frass plot 2
frassPeak2_dev1 <- lm(formula = peak2 ~ Dev_1,
                      data = frass)
frassPeak2_dev10 <- lm(formula = peak2 ~ Dev_10,
                       data = frass)
frassPeak2_temp <- lm(formula = peak2 ~ rel_temp,
                      data = frass)


Weights(AICc(frassPeak2_dev1, frassPeak2_dev10, frassPeak2_temp))
model.sel(frassPeak2_dev1, frassPeak2_dev10, frassPeak2_temp)
summary(frassPeak2_dev10) # -34.1
confint(frassPeak1_dev10) # -19.4 - 35.3

summary(frassPeak2_dev1) # 0.33
confint(frassPeak2_dev1) # -50.4 - 19.5

summary(frassPeak2_temp) # 0.2
confint(frassPeak2_temp) #-56.9 - 15.6

frass_plot_df2 <- plot_model(frassPeak2_dev10, type = "pred", terms = "Dev_10")$data

frass_plot2 <- ggplot() +
  geom_jitter(frass, mapping = aes(x = Dev_10, y = peak2)) +
  geom_line(frass_plot_df2, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(frass_plot_df2, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Proportion developed (10-km)", y = "DOY: Second peak") +
  theme_classic()
frass_plot2

summary(frassPeak2_dev10)

## micro moth time
microPeak1_dev1 <- lm(formula = peak1 ~ Dev_1,
                      data = micro)
microPeak1_dev10 <- lm(formula = peak1 ~ Dev_10,
                       data = micro)
microPeak1_temp <- lm(formula = peak1 ~ rel_temp,
                      data = micro)

Weights(AICc(microPeak1_dev1, microPeak1_dev10, microPeak1_temp))
model.sel(microPeak1_dev1, microPeak1_dev10, microPeak1_temp)

summary(microPeak1_temp) # -38.1
confint(microPeak1_temp) # -70.2 - -5.9

summary(microPeak1_dev1) # -5.7
confint(microPeak1_dev1) # -48.6 - 37.2

summary(microPeak1_dev10) #-0.04
confint(microPeak1_dev10) # -63.1 - 63.1


micro_plot_df1 <- plot_model(microPeak1_temp, type = "pred", terms = "rel_temp")$data

micro_plot1 <- ggplot() +
  geom_jitter(micro, mapping = aes(x = rel_temp, y = peak1)) +
  geom_line(micro_plot_df1, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(micro_plot_df1, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Relative temperature of site", y = "DOY: First peak") +
  ggtitle("Micro-moths") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) 
micro_plot1

# now for micro plot 2
microPeak2_dev1 <- lm(formula = peak2 ~ Dev_1,
                      data = micro)
microPeak2_dev10 <- lm(formula = peak2 ~ Dev_10,
                       data = micro)
microPeak2_temp <- lm(formula = peak2 ~ rel_temp,
                      data = micro)

Weights(AICc(microPeak2_dev1, microPeak2_dev10, microPeak2_temp))
model.sel(microPeak2_dev1, microPeak2_dev10, microPeak2_temp)

micro_plot_df2 <- plot_model(microPeak2_dev10, type = "pred", terms = "Dev_10")$data

summary(microPeak2_dev10) # 25.23
confint(microPeak2_dev10) # -68.6 - 118.9

summary(microPeak2_dev1) # 1.6
confint(microPeak2_dev1) #-64.4 - 67.7

summary(microPeak2_temp) # 1.26
confint(microPeak2_temp) #-70.2 - 72.7

micro_plot2 <- ggplot() +
  geom_jitter(micro, mapping = aes(x = Dev_10, y = peak2)) +
  geom_line(micro_plot_df2, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(micro_plot_df2, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Proportion developed (10-km)", y = "DOY: Second peak") +
  theme_classic()
micro_plot2

### macro time
macroPeak1_dev1 <- lm(formula = peak ~ Dev_1,
                      data = macro)
macroPeak1_dev10 <- lm(formula = peak ~ Dev_10,
                       data = macro)
macroPeak1_temp <- lm(formula = peak ~ rel_temp,
                      data = macro)

Weights(AICc(macroPeak1_dev1, macroPeak1_dev10, macroPeak1_temp))

model.sel(macroPeak1_dev1, macroPeak1_dev10, macroPeak1_temp)

summary(macroPeak1_temp) # -55.8
confint(macroPeak1_temp) # -101.8 - -9.8

summary(macroPeak1_dev1) #17.6
confint(macroPeak1_dev1) #-78.3 - 43.0

summary(macroPeak1_dev10) # -17.6
confint(macroPeak1_dev10) # -107.7 - 72.4

macro_plot_df1 <- plot_model(macroPeak1_temp, type = "pred", terms = "rel_temp")$data

macro_plot1 <- ggplot() +
  geom_point(macro, mapping = aes(x = rel_temp, y = peak)) +
  geom_line(macro_plot_df1, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(macro_plot_df1, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Relative tempearture", y = "Peak DOY") +
  ggtitle("Macro-moths") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) 
macro_plot1

## correlation plot of the urbanization variables
cor_df <- macro %>% 
  select(rel_temp, Dev_10, Dev_1) %>% 
  rename(Rel.Temp = rel_temp,
         Dev.1km = Dev_1,
         Dev.10km = Dev_10)

r <- cor(cor_df)
library(ggcorrplot)
cor_plot <- ggcorrplot(r, hc.order = T, type = "lower", lab = T)

c1 <- ggplot() +
  #  geom_line(aes(cor_df$Rel.Temp, cor_df$Dev.10km)) +
  geom_point(aes(cor_df$Rel.Temp, cor_df$Dev.10km), size = 0.75, shape = 1) +
  labs(x = "Rel.Temp", y = "Dev.10") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  annotate("text", label = "r=0.55", x = 0.15, y = 0.75, size = 4) +
  theme_classic() +
  theme(axis.title = element_text(size = 10))

c2 <- ggplot() +
  #  geom_line(aes(cor_df$Rel.Temp, cor_df$Dev.1km)) +
  geom_point(aes(cor_df$Rel.Temp, cor_df$Dev.1km), size = 0.75, shape = 1) +
  labs(x = "Rel.Temp", y = "Dev.1") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  annotate("text", label = "r=0.55", x = 0.15, y = 0.75, size = 4) +
  theme_classic() +
  theme(axis.title = element_text(size = 10))

c3 <- ggplot() +
  #  geom_line(aes(cor_df$Dev.10, cor_df$Dev.1km)) +
  geom_point(aes(cor_df$Dev.10, cor_df$Dev.1km), size = 0.75, shape = 1) +
  labs(x = "Dev.10", y = "Dev.1") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  annotate("text", label = "r=0.92", x = 0.15, y = 0.75, size = 4) +
  theme_classic() +
  theme(axis.title = element_text(size = 10))

cp_cor <- cowplot::plot_grid(c2,c1, c3)
cp_cor

library(ggpubr)
ga <- ggpubr::ggarrange(frass_plot1, micro_plot1, macro_plot1,
                        frass_plot2, micro_plot2, cor_plot, 
                        labels = c("A", "B", "C", "D", "E", "F"))

ga

ggsave(ga, filename = "figOutputs/Fig2_PooledPhenoResponses.png", width = 9, height = 5)

### write table
lm_to_table <- function(x){
  
  df <- data.frame(Predictor = colnames(x$model)[2],
                   Estimate = summary(x$coefficients)[[1]],
                   Lower.CI = confint(x)[2,1],
                   Upper.CI = confint(x)[2,2],
                   R2 = MuMIn::r.squaredGLMM(x)[[1]])
  
  return(df)
  
}


tbl <- rbind(
  lm_to_table(frassPeak1_dev1),
  lm_to_table(frassPeak1_dev10),
  lm_to_table(frassPeak1_temp),                         
                         
  lm_to_table(frassPeak2_dev1),
  lm_to_table(frassPeak2_dev10),
  lm_to_table(frassPeak2_temp),                        

  lm_to_table(microPeak1_dev1),                         
  lm_to_table(microPeak1_dev10),                         
  lm_to_table(microPeak1_temp),  
  
  lm_to_table(microPeak2_dev1),                         
  lm_to_table(microPeak2_dev10),                         
  lm_to_table(microPeak2_temp),  
  
  lm_to_table(macroPeak1_dev1),                         
  lm_to_table(macroPeak1_dev10),                         
  lm_to_table(macroPeak1_temp))

tbl$Model <- c(rep("Frass", 6),
               rep("Micro-moth", 6),
               rep("Macro-moth", 3))

tbl$Peak <- c(1,1,1,2,2,2,
              1,1,1,2,2,2,
              1,1,1)

write.csv(tbl, file = "tabOutputs/univariateModels.csv", row.names = F)

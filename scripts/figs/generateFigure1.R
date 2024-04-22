library(dplyr)
library(ggplot2)
library(sjPlot)
library(MuMIn)
library(lme4)

# read in frass pheno data
frass <- read.csv("data/frassPhenoMetrics_perDay_trapRE.csv")
micro <- read.csv("data/microMothsPhenoMetrics_negBinom.csv")
macro <- read.csv("data/macroMothsPhenoMetrics_negBinom.csv")

#join with urbanization data
urb <- read.csv("data/imperviousSurface.csv")
canopy <- read.csv("data/canopySurface.csv")

frass <- frass %>% 
  left_join(urb) %>% 
  left_join(canopy)

micro <- micro %>% 
  left_join(urb) %>% 
  left_join(canopy)

macro <- macro %>% 
  left_join(urb) %>% 
  left_join(canopy)

# make model for frass, 4 predictors are propImpervious 300m,  propCanopyCover30m
frassPeak1_PI300m <- lm(formula = peak1 ~ propImpervious300m,
                        data = frass)
frassPeak1_CP300m <- lm(formula = peak1 ~ propCanopy300m,
                        data = frass)

Weights(AICc(frassPeak1_PI300m, frassPeak1_CP300m))
model.sel(frassPeak1_PI300m, frassPeak1_CP300m)
summary(frassPeak1_PI300m) 

# Impervious surface
frass_peak1_plotPI_300m_df <- plot_model(frassPeak1_PI300m, type = "pred", terms = "propImpervious300m")$data

frass_peak1_plotDev <- ggplot() +
  geom_point(frass, mapping = aes(x = propImpervious300m, y = peak1)) +
  geom_line(frass_peak1_plotPI_300m_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(frass_peak1_plotPI_300m_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Impervious surface", y = "DOY: First peak") +
  ggtitle("Frass") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
frass_peak1_plotDev

# Canopy cover
frass_peak1_plotCP_300m_df <- plot_model(frassPeak1_CP300m, type = "pred", terms = "propCanopy300m")$data

frass_peak1_plotCP <- ggplot() +
  geom_point(frass, mapping = aes(x = propCanopy300m, y = peak1)) +
  geom_line(frass_peak1_plotCP_300m_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(frass_peak1_plotCP_300m_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Canopy cover", y = "DOY: First peak") +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
frass_peak1_plotCP

# now for frass peak 2
# make model for frass, 4 predictors are propImpervious 300m,  propCanopyCover30m
frassPeak2_PI300m <- lm(formula = peak2 ~ propImpervious300m,
                        data = frass)
frassPeak2_CP300m <- lm(formula = peak2 ~ propCanopy300m,
                        data = frass)

Weights(AICc(frassPeak2_PI300m, frassPeak2_CP300m))
model.sel(frassPeak2_PI300m, frassPeak2_CP300m)
summary(frassPeak2_PI300m) 

# Impervious surface
frass_peak2_plotPI_300m_df <- plot_model(frassPeak2_PI300m, type = "pred", terms = "propImpervious300m")$data

frass_peak2_plotDev <- ggplot() +
  geom_point(frass, mapping = aes(x = propImpervious300m, y = peak2)) +
  geom_line(frass_peak2_plotPI_300m_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(frass_peak2_plotPI_300m_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Impervious surface", y = "DOY: Second peak") +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
frass_peak2_plotDev

# Canopy cover
frass_peak2_plotCP_300m_df <- plot_model(frassPeak2_CP300m, type = "pred", terms = "propCanopy300m")$data

frass_peak2_plotCP <- ggplot() +
  geom_point(frass, mapping = aes(x = propCanopy300m, y = peak2)) +
  geom_line(frass_peak2_plotCP_300m_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(frass_peak2_plotCP_300m_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Canopy cover", y = "DOY: Second peak") +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
frass_peak2_plotCP
################################################################################
################################ micro moth time ###############################
################################################################################
microPeak1_PI300m <- lm(formula = peak1 ~ propImpervious300m,
                        data = micro)
microPeak1_CP300m <- lm(formula = peak1 ~ propCanopy300m,
                        data = micro)

Weights(AICc(microPeak1_PI300m, microPeak1_CP300m))
model.sel(microPeak1_PI300m, microPeak1_CP300m)
summary(microPeak1_PI300m) 

# Impervious surface
micro_peak1_plotPI_300m_df <- plot_model(microPeak1_PI300m, type = "pred", terms = "propImpervious300m")$data

micro_peak1_plotDev <- ggplot() +
  geom_point(micro, mapping = aes(x = propImpervious300m, y = peak1)) +
  geom_line(micro_peak1_plotPI_300m_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(micro_peak1_plotPI_300m_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Impervious surface", y = "DOY: First peak") +
  ggtitle("Micro-moths") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
micro_peak1_plotDev

# Canopy cover
micro_peak1_plotCP_300m_df <- plot_model(microPeak1_CP300m, type = "pred", terms = "propCanopy300m")$data

micro_peak1_plotCP <- ggplot() +
  geom_point(micro, mapping = aes(x = propCanopy300m, y = peak1)) +
  geom_line(micro_peak1_plotCP_300m_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(micro_peak1_plotCP_300m_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Canopy cover", y = "DOY: First peak") +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
micro_peak1_plotCP

# now for micro peak 2
# make model for micro, 4 predictors are propImpervious 300m,  propCanopyCover30m
microPeak2_PI300m <- lm(formula = peak2 ~ propImpervious300m,
                        data = micro)
microPeak2_CP300m <- lm(formula = peak2 ~ propCanopy300m,
                        data = micro)

Weights(AICc(microPeak2_PI300m, microPeak2_CP300m))
model.sel(microPeak2_PI300m, microPeak2_CP300m)
summary(microPeak2_PI300m) 

# Impervious surface
micro_peak2_plotPI_300m_df <- plot_model(microPeak2_PI300m, type = "pred", terms = "propImpervious300m")$data

micro_peak2_plotDev <- ggplot() +
  geom_point(micro, mapping = aes(x = propImpervious300m, y = peak2)) +
  geom_line(micro_peak2_plotPI_300m_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(micro_peak2_plotPI_300m_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Impervious surface", y = "DOY: Second peak") +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
micro_peak2_plotDev

# Canopy cover
micro_peak2_plotCP_300m_df <- plot_model(microPeak2_CP300m, type = "pred", terms = "propCanopy300m")$data

micro_peak2_plotCP <- ggplot() +
  geom_point(micro, mapping = aes(x = propCanopy300m, y = peak2)) +
  geom_line(micro_peak2_plotCP_300m_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(micro_peak2_plotCP_300m_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Canopy cover", y = "DOY: Second peak") +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
micro_peak2_plotCP

################################################################################
################################# macro time####################################
################################################################################
macroPeak1_PI3 <- lm(formula = peak ~ propImpervious300m,
                     data = macro)
macroPeak1_CP3 <- lm(formula = peak ~ propCanopy300m,
                     data = macro)

Weights(AICc(macroPeak1_PI3, macroPeak1_CP3))

model.sel(macroPeak1_PI3, macroPeak1_CP3)


macro_plotTemp_df <- plot_model(macroPeak1_PI3, type = "pred", terms = "propImpervious300m")$data

macro_peak_plotPI <- ggplot() +
  geom_point(macro, mapping = aes(x = propImpervious300m, y = peak)) +
  geom_line(macro_plotTemp_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(macro_plotTemp_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Impervious surface", y = "Peak DOY") +
  ggtitle("Macro-moths") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
macro_peak_plotPI

# canopy cover
macro_plotCP_df <- plot_model(macroPeak1_CP3, type = "pred", terms = "propCanopy300m")$data

macro_peak_plotCP <- ggplot() +
  geom_point(macro, mapping = aes(x = propCanopy300m, y = peak)) +
  geom_line(macro_plotCP_df, mapping = aes(x = x, y = predicted)) +
  geom_ribbon(macro_plotCP_df, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  labs(x = "Canopy cover", y = "Peak DOY") +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 18),
        axis.text=element_text(size=13.5))
macro_peak_plotCP


# make multi paneled plot
source('scripts/01_plotStudyArea.R')

library(ggpubr)

# Arrange together
p_left <- ggarrange(imperviousSurfacePlot, canopyCoverPlot, ncol = 1, labels = c("A", "B"),
                    label.y = 0.8, font.label = list(size = 18))
ga1 <- ggpubr::ggarrange(frass_peak1_plotDev, micro_peak1_plotDev, macro_peak_plotPI,
                         frass_peak2_plotDev, micro_peak2_plotDev, ggplot() + theme_void(), 
                         labels = c("C", "G", "K", "D", "H"), font.label = list(size = 18))

ga2 <- ggpubr::ggarrange(frass_peak1_plotCP, micro_peak1_plotCP, macro_peak_plotCP,
                         frass_peak2_plotCP, micro_peak2_plotCP, ggplot() + theme_void(),
                         labels = c("E", "I", "L", "F", "J"), font.label = list(size = 18))

p_right <- ggpubr::ggarrange(ga1, ga2, nrow = 2, ncol = 1)

p_all <- ggarrange(p_left, p_right, widths = c(1,1.7))

ggsave(p_all, filename = "figs/Fig2_PooledPhenoResponsesV3_negBinom.png", 
       width = 20, height = 15, dpi = 450)


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
  lm_to_table(frassPeak1_PI300m),
  lm_to_table(frassPeak1_CP300m),                  
  
  lm_to_table(frassPeak2_PI300m),
  lm_to_table(frassPeak2_CP300m),                        
  
  lm_to_table(microPeak1_PI300m),
  lm_to_table(microPeak1_CP300m),  
  
  lm_to_table(microPeak2_PI300m),
  lm_to_table(microPeak2_CP300m),  
  
  lm_to_table(macroPeak1_PI3),                         
  lm_to_table(macroPeak1_CP3))

tbl$Model <- c(rep("Frass", 4),
               rep("Micro-moth", 4),
               rep("Macro-moth", 2))

tbl$Peak <- c(1,1,2,2,
              1,1,2,2,
              1,1)

write.csv(tbl, file = "tabOutputs/univariateModels.csv", row.names = F)
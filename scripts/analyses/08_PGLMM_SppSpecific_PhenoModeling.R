library(phyr)
library(ape)

source("scripts/analyses/07_SppSpecific_PhenoModeling.R")

# function to capitalize spp names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

tt <- read.tree("data/traits/insect_tree_wBranches_allSpecies.tre")
tt$tip.label <- stringr::str_replace(tt$tip.label, pattern = "_", " ")
tt$tip.label <- firstup(tt$tip.label)

mdf_phylo <- left_join(mdf, data.frame(validName = tt$tip.label, Phylo = "Yes"))

#imput missing phylo with Heterocampa obliqua
mdf_phylo <- mdf_phylo %>% 
  mutate(pseudoValidName = if_else(is.na(Phylo), true = "Heterocampa obliqua",
                                   false = validName))

sppNotInAnlysis <- data.frame(validName = tt$tip.label) %>% 
  filter(!validName %in% mdf_phylo$validName)

tt <- ape::drop.tip(tt, tip = sppNotInAnlysis$validName)

# peak model
summary(m3)
#Formula: peak ~ rel_temp + bio1_mean + rel_temp:bio1_mean + (1 | validName)

peak_pglmm <- pglmm(peak ~ rel_temp + voltinism  + maxWingspan + 
                      hostPlantSepcialization + bio1_mean + bio1_sd +
                      rel_temp:voltinism +
                      rel_temp:maxWingspan +
                      rel_temp:hostPlantSepcialization +
                      rel_temp:bio1_mean +
                      rel_temp:bio1_sd +
                      (1|pseudoValidName), 
                    data = mdf_phylo, 
                    cov_ranef = list(pseudoValidName = tt), 
                    bayes = TRUE)

summary(peak_pglmm)
plot_bayes(peak_pglmm)


#PGLMM has same effects as LMM so will report LMM
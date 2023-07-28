library(phyr)
library(ape)

source("scripts/Modeling/speciesSpecificPhenoModels_Poisson.R")

# function to capitalize spp names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

tt <- read.tree("data/phylogeny/insect_tree_wBranches_allSpecies.tre")
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
summary(m3_top)
#Formula: peak ~ rel_temp + bio1_mean + rel_temp:bio1_mean + (1 | validName)

peak_pglmm <- pglmm(formula = peak ~ rel_temp + voltinism + maxWingspan + bio1_mean + 
                      rel_temp:bio1_mean + rel_temp:maxWingspan + 
                      (1 | pseudoValidName__), 
                    data = mdf_phylo, 
                    cov_ranef = list(pseudoValidName = tt), 
                    bayes = TRUE)

summary(peak_pglmm)
plot_bayes(peak_pglmm)


#PGLMM has same effects as LMM so will report LMM
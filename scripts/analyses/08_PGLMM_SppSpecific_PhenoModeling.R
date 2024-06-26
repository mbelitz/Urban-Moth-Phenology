library(phyr)
library(ape)

source("scripts/analyses/07_SppSpecific_PhenoModeling_v2_PI_PC.R")

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

# impervious surface model
summary(m)
#Formula: peak ~ rel_temp + bio1_mean + rel_temp:bio1_mean + (1 | validName)

is_pglmm <- pglmm(formula = peak ~ propImpervious300m + maxWingspan + bio1_mean + 
                      propImpervious300m:maxWingspan +
                      propImpervious300m:bio1_mean +
                      (1|pseudoValidName__), 
                    data = mdf_phylo, 
                    cov_ranef = list(pseudoValidName = tt), 
                    bayes = TRUE)

summary(is_pglmm)
plot_bayes(is_pglmm)

# canopy cover model
summary(m)
#Formula: peak ~ rel_temp + bio1_mean + rel_temp:bio1_mean + (1 | validName)

cc_pglmm <- pglmm(formula = peak ~ propCanopy300m + maxWingspan + bio1_mean + 
                    propCanopy300m:maxWingspan +
                    propCanopy300m:bio1_mean +
                    (1|pseudoValidName__), 
                  data = mdf_phylo, 
                  cov_ranef = list(pseudoValidName = tt), 
                  bayes = TRUE)

summary(cc_pglmm)
plot_bayes(cc_pglmm)


#PGLMM has same effects as LMM so will report LMM
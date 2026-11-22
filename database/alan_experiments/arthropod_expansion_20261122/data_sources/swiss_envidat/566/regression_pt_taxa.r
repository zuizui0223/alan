# load packages ----
library(dplyr)
library(conflicted)
library(emmeans)
library(lmerTest)
library(DHARMa)
library(regclass)
library(car)
library(caTools)
library(performance)
library(see)
library(multcompView)
library(lsmeans)
library(multcomp)
library(GLMMadaptive)
library(patchwork)
library(sjmisc)
library(writexl)
library(rio)
library(MuMIn)
library(AER)
library(optimx)
library(tidyverse)
library(data.table)
library(readxl)
library(ggplot2)
library(conflicted)
library(writexl)
require(grid)
library(ggpubr)
library(lme4)
library(ggsignif)
library(gmodels)
library(modelbased)
library(modelsummary)

conflict_prefer("filter", "dplyr")
conflict_prefer("margin", "ggplot2")
conflict_prefer("summarise", "dplyr")
conflict_prefer("lmer", "lmerTest")
conflict_prefer("filter", "dplyr")
conflict_prefer("step", "lmerTest")
conflict_prefer("residuals", "GLMMadaptive")

# load data
C <- read.csv("1_Data/PT/ALANeX_PT_all_clean_control.csv", sep = ",")

# filter data (including control term)
CC <- C %>% 
  filter(temp %in% c("4000", "3000", "2000", "control")) %>%
  filter(level %in% c("50", "100", "control")) %>%
  filter(issues %in% c("0"))

# filter data (excluding control term) 
RR <- C %>% 
  filter(temp %in% c("2000", "3000", "4000")) %>%
  filter(level %in% c("50", "100")) %>% 
  filter(issues %in% c("0"))

## PT - Taxa
#### Nematocera (incl. control) ----
lmNEc <- lmer(log(Nematocera + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=CC)

## contrasts
nem_kc <- estimate_means(lmNEc, at = "control")
estimate_contrasts(lmNEc, contrast = "control")

#### Lepidoptera (incl. control)----
lmLEc <- lmer(log(Lepidoptera + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=CC)

## contrasts 
lep_kc <- estimate_means(lmLEc, at = "control")
estimate_contrasts(lmLEc, contrast = "control")

#### Brachycera (incl. control)----
lmBRc <- lmer(log(Brachycera + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=CC)

## contrasts
bra_kc <- estimate_means(lmBRc, at = "control")
estimate_contrasts(lmBRc, contrast = "control")

#### Coleoptera (incl. control)----
lmCOc <- lmer(log(Coleoptera + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=CC)

## contrasts
cor_kc <- estimate_means(lmCOc, at = "control")
estimate_contrasts(lmCOc, contrast = "control")

#### Hymenoptera (incl. control)----
lmHYc <- lmer(log(Hymenoptera + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=CC)

## contrasts
hym_kc <- estimate_means(lmHYc, at = "control")
estimate_contrasts(lmHYc, contrast = "control")

#### Hemiptera (incl. control)----
lmHEc <- lmer(log(Hemiptera + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=CC)

## contrasts
hem_kc <- estimate_means(lmHEc, at = "control")
estimate_contrasts(lmHEc, contrast = "control")

#### Neuropteroidea (incl. control)----
lmNEUc <- lmer(log(Neuropteroidea + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=CC)

## contrasts
neu_kc <- estimate_means(lmNEUc, at = "control")
estimate_contrasts(lmNEUc, contrast = "control")

#### Trichoptera (incl. control)----
lmTRc <- lmer(log(Trichoptera + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=CC)

## contrasts
tri_kc <- estimate_means(lmTRc, at = "control")
estimate_contrasts(lmTRc, contrast = "control")


#### Nematocera (excl. control)----
lmNE <- lmer(log(Nematocera + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=RR)

## contrasts
### individual treatments
nem_k <- estimate_means(lmNE, at = "temp")
estimate_contrasts(lmNE, contrast = "temp")

nem_l <- estimate_means(lmNE, at = "level")
estimate_contrasts(lmNE, contrast = "level")

nem_d <- estimate_means(lmNE, at = "shape")
estimate_contrasts(lmNE, contrast = "shape")

### interactions
nem_kl <- estimate_means(lmNE, at = "temp:level")
estimate_contrasts(lmNE, at = "temp", contrast ="level")

nem_dl <- estimate_means(lmNE, at = "shape:level")
estimate_contrasts(lmNE, at = "shape", contrast ="level")

nem_kd <- estimate_means(lmNE, at = "temp:shape")
estimate_contrasts(lmNE, at = "temp", contrast ="shape")

#### Lepidoptera (excl. control)----
lmLE <- lmer(log(Lepidoptera + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=RR)


## contrasts
### individual treatments
lep_k <- estimate_means(lmLE, at = "temp")
estimate_contrasts(lmLE, contrast = "temp")

lep_l <- estimate_means(lmLE, at = "level")
estimate_contrasts(lmLE, contrast = "level")

lep_d <- estimate_means(lmLE, at = "shape")
estimate_contrasts(lmLE, contrast = "shape")

### interactions
lep_kl <- estimate_means(lmLE, at = "temp:level")
estimate_contrasts(lmLE, at = "temp", contrast ="level")

lep_dl <- estimate_means(lmLE, at = "shape:level")
estimate_contrasts(lmLE, at = "shape", contrast ="level")

lep_kd <- estimate_means(lmLE, at = "temp:shape")
estimate_contrasts(lmLE, at = "temp", contrast ="shape")


#### Brachycera (excl. control)----
lmBR <- lmer(log(Brachycera + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=RR)

## contrasts
#### individual treatments
bra_k <- estimate_means(lmBR, at = "temp")
estimate_contrasts(lmBR, contrast = "temp")

bra_l <- estimate_means(lmBR, at = "level")
estimate_contrasts(lmBR, contrast = "level")

bra_d <- estimate_means(lmBR, at = "shape")
estimate_contrasts(lmBR, contrast = "shape")

### interactions
bra_kl <- estimate_means(lmBR, at = "temp:level")
estimate_contrasts(lmBR, at = "temp", contrast ="level")

bra_dl <- estimate_means(lmBR, at = "shape:level")
estimate_contrasts(lmBR, at = "shape", contrast ="level")

bra_kd <- estimate_means(lmBR, at = "temp:shape")
estimate_contrasts(lmBR, at = "temp", contrast ="shape")



#### Coleoptera (excl. control)----
lmCO <- lmer(log(Coleoptera + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot), data=RR)

## contrasts
### individual treatments
cor_k <- estimate_means(lmCO, at = "temp")
estimate_contrasts(lmCO, contrast = "temp")

cor_l <- estimate_means(lmCO, at = "level")
estimate_contrasts(lmCO, contrast = "level")

cor_d <- estimate_means(lmCO, at = "shape")
estimate_contrasts(lmCO, contrast = "shape")

### interactions
cor_kl <- estimate_means(lmCO, at = "temp:level")
estimate_contrasts(lmCO, at = "temp", contrast ="level")

cor_dl <- estimate_means(lmCO, at = "shape:level")
estimate_contrasts(lmCO, at = "shape", contrast ="level")

cor_kd <- estimate_means(lmCO, at = "temp:shape")
estimate_contrasts(lmCO, at = "temp", contrast ="shape")


#### Hymenoptera (excl. control)----
lmHY <- lmer(log(Hymenoptera + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=RR)


## contrasts
### individual treatments
hym_k <- estimate_means(lmHY, at = "temp")
estimate_contrasts(lmHY, contrast = "temp")

hym_l <- estimate_means(lmHY, at = "level")
estimate_contrasts(lmHY, contrast = "level")

hym_d <- estimate_means(lmHY, at = "shape")
estimate_contrasts(lmHY, contrast = "shape")

### interactions
hym_kl <- estimate_means(lmHY, at = "temp:level")
estimate_contrasts(lmHY, at = "temp", contrast ="level")

hym_dl <- estimate_means(lmHY, at = "shape:level")
estimate_contrasts(lmHY, at = "shape", contrast ="level")

hym_kd <- estimate_means(lmHY, at = "temp:shape")
estimate_contrasts(lmHY, at = "temp", contrast ="shape")


#### Hemiptera (excl. control)----
lmHE <- lmer(log(Hemiptera + 1) ~ temp * level + temp * shape + shape * level + (1|week_ID) + (1|site/plot), data=RR)

## contrasts
### individual treatments
hem_k <- estimate_means(lmHE, at = "temp")
estimate_contrasts(lmHE, contrast = "temp")

hem_l <- estimate_means(lmHE, at = "level")
estimate_contrasts(lmHE, contrast = "level")

hem_d <- estimate_means(lmHE, at = "shape")
estimate_contrasts(lmHE, contrast = "shape")

### interactions
hem_kl <- estimate_means(lmHE, at = "temp:level")
estimate_contrasts(lmHE, at = "temp", contrast ="level")

hem_dl <- estimate_means(lmHE, at = "shape:level")
estimate_contrasts(lmHE, at = "shape", contrast ="level")

hem_kd <- estimate_means(lmHE, at = "temp:shape")
estimate_contrasts(lmHE, at = "temp", contrast ="shape")


#### Neuropteroidea (excl. control)----
lmNEU <- lmer(log(Neuropteroidea + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=RR)

## contrasts
### individual treatments
neu_k <- estimate_means(lmNEU, at = "temp")
estimate_contrasts(lmNEU, contrast = "temp")

neu_l <- estimate_means(lmNEU, at = "level")
estimate_contrasts(lmNEU, contrast = "level")

neu_d <- estimate_means(lmNEU, at = "shape")
estimate_contrasts(lmNEU, contrast = "shape")

### interactions
neu_kl <- estimate_means(lmNEU, at = "temp:level")
estimate_contrasts(lmNEU, at = "temp", contrast ="level")

neu_dl <- estimate_means(lmNEU, at = "shape:level")
estimate_contrasts(lmNEU, at = "shape", contrast ="level")

neu_kd <- estimate_means(lmNEU, at = "temp:shape")
estimate_contrasts(lmNEU, at = "temp", contrast ="shape")




#### Trichoptera (excl. control)----
lmTR <- lmer(log(Trichoptera + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site / plot) , data=RR)

## contrasts
### individual treatments
tri_k <- estimate_means(lmTR, at = "temp")
estimate_contrasts(lmTR, contrast = "temp")

tri_l <- estimate_means(lmTR, at = "level")
estimate_contrasts(lmTR, contrast = "level")

tri_d <- estimate_means(lmTR, at = "shape")
estimate_contrasts(lmTR, contrast = "shape")

### interactions
tri_kl <- estimate_means(lmTR, at = "temp:level")
estimate_contrasts(lmTR, at = "temp", contrast ="level")

tri_dl <- estimate_means(lmTR, at = "shape:level")
estimate_contrasts(lmTR, at = "shape", contrast ="level")

tri_kd <- estimate_means(lmTR, at = "temp:shape")
estimate_contrasts(lmTR, at = "temp", contrast ="shape")


##### ANOVA Tables & R2 & AIC ----
## incl. Control----
# Nematocera
x <- Anova(lmNEc, type=3)
r.squaredGLMM(lmNEc)
AIC(lmNEc)

anova_nem <- as.data.frame(x)
write_xlsx(anova_nem, "NematoceraC.xlsx")

# Lepidoptera
x <- Anova(lmLEc, type=3)
r.squaredGLMM(lmLEc)
AIC(lmLEc)

anova_lep <- as.data.frame(x)
write_xlsx(anova_lep, "LepidopteraC.xlsx")

# Brachycera
x <- Anova(lmBRc, type=3)
r.squaredGLMM(lmBRc)
AIC(lmBRc)

anova_bra <- as.data.frame(x)
write_xlsx(anova_bra, "BrachyceraC.xlsx")

# Coleoptera
x <- Anova(lmCOc, type=3)
r.squaredGLMM(lmCOc)
AIC(lmCOc)

anova_col <- as.data.frame(x)
write_xlsx(anova_col, "ColeopteraC.xlsx")

# Hymenoptera
x <- Anova(lmHYc, type=3)
r.squaredGLMM(lmHYc)
AIC(lmHYc)

anova_hym <- as.data.frame(x)
write_xlsx(anova_hym, "HymenopteraC.xlsx")

# Hemiptera
x <- Anova(lmHEc, type=3)
r.squaredGLMM(lmHEc)
AIC(lmHEc)

anova_hen <- as.data.frame(x)
write_xlsx(anova_hen, "HemipteraC.xlsx")

# Neuropteroidea
x <- Anova(lmNEUc, type=3)
r.squaredGLMM(lmNEUc)
AIC(lmNEUc)

anova_neu <- as.data.frame(x)
write_xlsx(anova_neu, "NeuropteroideaC.xlsx")

# Trichoptera
x <- Anova(lmTRc, type=3)
r.squaredGLMM(lmTRc)
AIC(lmTRc)

anova_tri <- as.data.frame(x)
write_xlsx(anova_tri, "TrichopteraC.xlsx")



## excl. Control----
# Nematocera
x <- Anova(lmNE, type=3)
r.squaredGLMM(lmNE)
AIC(lmNE)

anova_nem <- as.data.frame(x)
write_xlsx(anova_nem, "Nematocera.xlsx")

# Lepidoptera
x <- Anova(lmLE, type=3)
r.squaredGLMM(lmLE)
AIC(lmLE)

anova_lep <- as.data.frame(x)
write_xlsx(anova_lep, "Lepidoptera.xlsx")

# Brachycera
x <- Anova(lmBR, type=3)
r.squaredGLMM(lmBR)
AIC(lmBR)

anova_bra <- as.data.frame(x)
write_xlsx(anova_bra, "Brachycera.xlsx")

# Coleoptera
x <- Anova(lmCO, type=3)
r.squaredGLMM(lmCO)
AIC(lmCO)

anova_col <- as.data.frame(x)
write_xlsx(anova_col, "Coleoptera.xlsx")

# Hymenoptera
x <- Anova(lmHY, type=3)
r.squaredGLMM(lmHY)
AIC(lmHY)

anova_hym <- as.data.frame(x)
write_xlsx(anova_hym, "Hymenoptera.xlsx")

# Hemiptera
x <- Anova(lmHE, type=3)
r.squaredGLMM(lmHE)
AIC(lmHE)

anova_hen <- as.data.frame(x)
write_xlsx(anova_hen, "Hemiptera.xlsx")

# Neuropteroidea
x <- Anova(lmNEU, type=3)
r.squaredGLMM(lmNEU)
AIC(lmNEU)

anova_neu <- as.data.frame(x)
write_xlsx(anova_neu, "Neuropteroidea.xlsx")

# Trichoptera
x <- Anova(lmTR, type=3)
r.squaredGLMM(lmTR)
AIC(lmTR)

anova_tri <- as.data.frame(x)
write_xlsx(anova_tri, "Trichoptera.xlsx")


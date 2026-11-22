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
K <- read.csv("1_Data/TRF/ALANeX_TRF_all_clean_control.csv", sep = ",")

# filter data (including control term) 
KK <- K %>% 
  filter(temp %in% c("2000", "3000", "4000", "control")) %>%
  filter(level %in% c("50", "100", "control")) %>%
  filter(issues %in% c("0"))

# filter data (excluding control term) 
DD <- K %>% 
  filter(temp %in% c("2000", "3000", "4000")) %>%
  filter(level %in% c("50", "100")) %>% 
  filter(issues %in% c("0"))

## TRF - Taxa
#### Collembola (incl. control) ----
lmCOLc <- lmer(log(Collembola + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts
coll_kc <- estimate_means(lmCOLc, at = "control")
estimate_contrasts(lmCOLc, contrast = "control")

#### Acarina (incl. control)----
lmACc <- lmer(log(Acarina + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts 
aca_kc <- estimate_means(lmACc, at = "control")
estimate_contrasts(lmACc, contrast = "control")

#### Formicidae (incl. control)----
lmFOc <- lmer(log(Formicidae + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts
for_kc <- estimate_means(lmFOc, at = "control")
estimate_contrasts(lmFOc, contrast = "control")

#### Arachnida (incl. control)----
lmARc <- lmer(log(Arachnida + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts
ara_kc <- estimate_means(lmARc, at = "control")
estimate_contrasts(lmARc, contrast = "control")

#### Coleo_Rest (incl. control)----
lmCORc <- lmer(log(Coleo_Rest + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts
col_kc <- estimate_means(lmCORc, at = "control")
estimate_contrasts(lmCORc, contrast = "control")

#### Coleo_Pred (incl. control)----
lmCOPc <- lmer(log(Coleo_Pred + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts
cop_kc <- estimate_means(lmCOPc, at = "control")
estimate_contrasts(lmCOPc, contrast = "control")

#### Myriapoda (incl. control)----
lmMYc <- lmer(log(Myriapoda + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts
myr_kc <- estimate_means(lmMYc, at = "control")
estimate_contrasts(lmMYc, contrast = "control")

#### Hemiptera (incl. control)----
lmHEMc <- lmer(log(Hemiptera + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts
hen_kc <- estimate_means(lmHEMc, at = "control")
estimate_contrasts(lmHEMc, contrast = "control")


#### Collembola (excl. control)----
lmCOL <- lmer(log(Collembola + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=DD)

## contrasts
### individual treatments
coll_k <- estimate_means(lmCOL, at = "temp")
estimate_contrasts(lmCOL, contrast = "temp")

coll_l <- estimate_means(lmCOL, at = "level")
estimate_contrasts(lmCOL, contrast = "level")

coll_d <- estimate_means(lmCOL, at = "shape")
estimate_contrasts(lmCOL, contrast = "shape")

### interactions
coll_kl <- estimate_means(lmCOL, at = "temp:level")
estimate_contrasts(lmCOL, at = "temp", contrast ="level")

coll_dl <- estimate_means(lmCOL, at = "shape:level")
estimate_contrasts(lmCOL, at = "shape", contrast ="level")

coll_kd <- estimate_means(lmCOL, at = "temp:shape")
estimate_contrasts(lmCOL, at = "temp", contrast ="shape")

#### Acarina (excl. control)----
lmAC <- lmer(log(Acarina + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=DD)


## contrasts
### individual treatments
aca_k <- estimate_means(lmAC, at = "temp")
estimate_contrasts(lmAC, contrast = "temp")

aca_l <- estimate_means(lmAC, at = "level")
estimate_contrasts(lmAC, contrast = "level")

aca_d <- estimate_means(lmAC, at = "shape")
estimate_contrasts(lmAC, contrast = "shape")

### interactions
aca_kl <- estimate_means(lmAC, at = "temp:level")
estimate_contrasts(lmAC, at = "temp", contrast ="level")

aca_dl <- estimate_means(lmAC, at = "shape:level")
estimate_contrasts(lmAC, at = "shape", contrast ="level")

aca_kd <- estimate_means(lmAC, at = "temp:shape")
estimate_contrasts(lmAC, at = "temp", contrast ="shape")


#### Formicidae (excl. control)----
lmFO <- lmer(log(Formicidae + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=DD)

## contrasts
#### individual treatments
for_k <- estimate_means(lmFO, at = "temp")
estimate_contrasts(lmFO, contrast = "temp")

for_l <- estimate_means(lmFO, at = "level")
estimate_contrasts(lmFO, contrast = "level")

for_d <- estimate_means(lmFO, at = "shape")
estimate_contrasts(lmFO, contrast = "shape")

### interactions
for_kl <- estimate_means(lmFO, at = "temp:level")
estimate_contrasts(lmFO, at = "temp", contrast ="level")

for_dl <- estimate_means(lmFO, at = "shape:level")
estimate_contrasts(lmFO, at = "shape", contrast ="level")

for_kd <- estimate_means(lmFO, at = "temp:shape")
estimate_contrasts(lmFO, at = "temp", contrast ="shape")



#### Arachnida (excl. control)----
lmAR <- lmer(log(Arachnida + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot), data=DD)

## contrasts
### individual treatments
ara_k <- estimate_means(lmAR, at = "temp")
estimate_contrasts(lmAR, contrast = "temp")

ara_l <- estimate_means(lmAR, at = "level")
estimate_contrasts(lmAR, contrast = "level")

ara_d <- estimate_means(lmAR, at = "shape")
estimate_contrasts(lmAR, contrast = "shape")

### interactions
ara_kl <- estimate_means(lmAR, at = "temp:level")
estimate_contrasts(lmAR, at = "temp", contrast ="level")

ara_dl <- estimate_means(lmAR, at = "shape:level")
estimate_contrasts(lmAR, at = "shape", contrast ="level")

ara_kd <- estimate_means(lmAR, at = "temp:shape")
estimate_contrasts(lmAR, at = "temp", contrast ="shape")




#### Coleo_Rest (excl. control)----
lmCOR <- lmer(log(Coleo_Rest + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=DD)


## contrasts
### individual treatments
col_k <- estimate_means(lmCOR, at = "temp")
estimate_contrasts(lmCOR, contrast = "temp")

col_l <- estimate_means(lmCOR, at = "level")
estimate_contrasts(lmCOR, contrast = "level")

col_d <- estimate_means(lmCOR, at = "shape")
estimate_contrasts(lmCOR, contrast = "shape")

### interactions
col_kl <- estimate_means(lmCOR, at = "temp:level")
estimate_contrasts(lmCOR, at = "temp", contrast ="level")

col_dl <- estimate_means(lmCOR, at = "shape:level")
estimate_contrasts(lmCOR, at = "shape", contrast ="level")

col_kd <- estimate_means(lmCOR, at = "temp:shape")
estimate_contrasts(lmCOR, at = "temp", contrast ="shape")



#### Coleo_Pred (excl. control)----
lmCOP <- lmer(log(Coleo_Pred + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot), data=DD)

## contrasts
### individual treatments
cop_k <- estimate_means(lmCOP, at = "temp")
estimate_contrasts(lmCOP, contrast = "temp")

cop_l <- estimate_means(lmCOP, at = "level")
estimate_contrasts(lmCOP, contrast = "level")

cop_d <- estimate_means(lmCOP, at = "shape")
estimate_contrasts(lmCOP, contrast = "shape")

### interactions
cop_kl <- estimate_means(lmCOP, at = "temp:level")
estimate_contrasts(lmCOP, at = "temp", contrast ="level")

cop_dl <- estimate_means(lmCOP, at = "shape:level")
estimate_contrasts(lmCOP, at = "shape", contrast ="level")

cop_kd <- estimate_means(lmCOP, at = "temp:shape")
estimate_contrasts(lmCOP, at = "temp", contrast ="shape")


#### Myriapoda (excl. control)----
lmMY <- lmer(log(Myriapoda + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=DD)

## contrasts
### individual treatments
myr_k <- estimate_means(lmMY, at = "temp")
estimate_contrasts(lmMY, contrast = "temp")

myr_l <- estimate_means(lmMY, at = "level")
estimate_contrasts(lmMY, contrast = "level")

myr_d <- estimate_means(lmMY, at = "shape")
estimate_contrasts(lmMY, contrast = "shape")

### interactions
myr_kl <- estimate_means(lmMY, at = "temp:level")
estimate_contrasts(lmMY, at = "temp", contrast ="level")

myr_dl <- estimate_means(lmMY, at = "shape:level")
estimate_contrasts(lmMY, at = "shape", contrast ="level")

myr_kd <- estimate_means(lmMY, at = "temp:shape")
estimate_contrasts(lmMY, at = "temp", contrast ="shape")




#### Hemiptera (excl. control)----
lmHEM <- lmer(log(Hemiptera + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site / plot) , data=DD)

## contrasts
### individual treatments
hen_k <- estimate_means(lmHEM, at = "temp")
estimate_contrasts(lmHEM, contrast = "temp")

hen_l <- estimate_means(lmHEM, at = "level")
estimate_contrasts(lmHEM, contrast = "level")

hen_d <- estimate_means(lmHEM, at = "shape")
estimate_contrasts(lmHEM, contrast = "shape")

### interactions
hen_kl <- estimate_means(lmHEM, at = "temp:level")
estimate_contrasts(lmHEM, at = "temp", contrast ="level")

hen_dl <- estimate_means(lmHEM, at = "shape:level")
estimate_contrasts(lmHEM, at = "shape", contrast ="level")

hen_kd <- estimate_means(lmHEM, at = "temp:shape")
estimate_contrasts(lmHEM, at = "temp", contrast ="shape")







##### ANOVA Tables & R2 & AIC ----
## incl. Control----
# Collembola
x <- Anova(lmCOLc, type=3)
r.squaredGLMM(lmCOLc)
AIC(lmCOLc)

anova_coll <- as.data.frame(x)
write_xlsx(anova_coll, "CollembolaC.xlsx")

# Acarina
x <- Anova(lmACc, type=3)
r.squaredGLMM(lmACc)
AIC(lmACc)

anova_aca <- as.data.frame(x)
write_xlsx(anova_aca, "AcarinaC.xlsx")

# Formicidae
x <- Anova(lmFOc, type=3)
r.squaredGLMM(lmFOc)
AIC(lmFOc)

anova_for <- as.data.frame(x)
write_xlsx(anova_for, "FormicidaeC.xlsx")

# Arachnida
x <- Anova(lmARc, type=3)
r.squaredGLMM(lmARc)
AIC(lmARc)

anova_ara <- as.data.frame(x)
write_xlsx(anova_ara, "ArachnidaC.xlsx")

# Coleoptera pred.
x <- Anova(lmCOPc, type=3)
r.squaredGLMM(lmCOPc)
AIC(lmCOPc)

anova_cop <- as.data.frame(x)
write_xlsx(anova_cop, "Coleo_PredC.xlsx")

# Coleoptera rest
x <- Anova(lmCORc, type=3)
r.squaredGLMM(lmCORc)
AIC(lmCORc)

anova_cor <- as.data.frame(x)
write_xlsx(anova_cor, "Coleo_RestC.xlsx")

# Myriapoda
x <- Anova(lmMYc, type=3)
r.squaredGLMM(lmMYc)
AIC(lmMYc)

anova_myr <- as.data.frame(x)
write_xlsx(anova_myr, "MyriapodaC.xlsx")

# Hemiptera
x <- Anova(lmHEMc, type=3)
r.squaredGLMM(lmHEMc)
AIC(lmHEMc)

anova_hem <- as.data.frame(x)
write_xlsx(anova_hem, "HemipteraC.xlsx")



## excl. Control----
# Collembola
x <- Anova(lmCOL, type=3)
r.squaredGLMM(lmCOL)
AIC(lmCOL)

anova_coll <- as.data.frame(x)
write_xlsx(anova_coll, "Collembola.xlsx")

# Acarina
x <- Anova(lmAC, type=3)
r.squaredGLMM(lmAC)
AIC(lmAC)

anova_aca <- as.data.frame(x)
write_xlsx(anova_aca, "Acarina.xlsx")

# Formicidae
x <- Anova(lmFO, type=3)
r.squaredGLMM(lmFO)
AIC(lmFO)

anova_for <- as.data.frame(x)
write_xlsx(anova_for, "Formicidae.xlsx")

# Arachnida
x <- Anova(lmAR, type=3)
r.squaredGLMM(lmAR)
AIC(lmAR)

anova_ara <- as.data.frame(x)
write_xlsx(anova_ara, "Arachnida.xlsx")

# Coleoptera pred.
x <- Anova(lmCOP, type=3)
r.squaredGLMM(lmCOP)
AIC(lmCOP)

anova_cop <- as.data.frame(x)
write_xlsx(anova_cop, "Coleo_Pred.xlsx")

# Coleoptera rest
x <- Anova(lmCOR, type=3)
r.squaredGLMM(lmCOR)
AIC(lmCOR)

anova_cor <- as.data.frame(x)
write_xlsx(anova_cor, "Coleo_Rest.xlsx")

# Myriapoda
x <- Anova(lmMY, type=3)
r.squaredGLMM(lmMY)
AIC(lmMY)

anova_myr <- as.data.frame(x)
write_xlsx(anova_myr, "Myriapoda.xlsx")

# Hemiptera
x <- Anova(lmHEM, type=3)
r.squaredGLMM(lmHEM)
AIC(lmHEM)

anova_hem <- as.data.frame(x)
write_xlsx(anova_hem, "Hemiptera.xlsx")


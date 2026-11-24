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

#### TRF & PT - All Groups (incl. control) ----
# load data
K <- read.csv("1_Data/PT&TRF/PT_TRF_all.csv", sep = ",")

KK <- K %>% 
  filter(temp %in% c("2000", "3000", "4000", "control")) %>%     ##### OHNE KONTROLLE; KONTROLLE RAUSGEFILTERT
  filter(level %in% c("50", "100", "control")) %>%
  filter(issues %in% c("0"))


lmAllc <- lmer(log(SumAll + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)


## contrasts
All_kc <- estimate_means(lmAllc, at = "control")
estimate_contrasts(lmAllc, contrast = "control")

#### PT - All Groups (incl. control) ----
### load data
C <- read.csv("1_Data/PT/ALANeX_PT_all_clean_control.csv", sep = ",")

CC <- C %>% 
  filter(temp %in% c("4000", "3000", "2000", "control")) %>%     ##### OHNE KONTROLLE; KONTROLLE RAUSGEFILTERT
  filter(level %in% c("50", "100", "control")) %>%
  filter(issues %in% c("0"))

lmPTc <- lmer(log(SumAll + 1) ~ control + temp + level + shape + (1|week_ID) + (1|site/plot) , data=CC)


## contrasts
PT_kc <- estimate_means(lmPTc, at = "control")
estimate_contrasts(lmPTc, contrast = "control")

#### TRF - All Groups (incl. control) ----
### load data
K <- read.csv("1_Data/TRF/ALANeX_TRF_all_clean_control.csv", sep = ",")

KK <- K %>% 
  filter(temp %in% c("2000", "3000", "4000", "control")) %>%     ##### OHNE KONTROLLE; KONTROLLE RAUSGEFILTERT
  filter(level %in% c("50", "100", "control")) %>%
  filter(issues %in% c("0"))

lmTRFc <- lmer(log(SumAll + 1) ~ control + temp + level + shape + (1|year/week_ID) + (1|site/plot) , data=KK)

## contrasts
TRF_kc <- estimate_means(lmTRFc, at = "control")
estimate_contrasts(lmTRFc, contrast = "control")




#### TRF & PT - All Groups (excl. control) ----
# load data
R <- read.csv("1_Data/PT&TRF/PT_TRF_all.csv", sep = ",")

RR <- R %>% 
  filter(temp %in% c("4000", "3000", "2000")) %>%     ##### OHNE KONTROLLE; KONTROLLE RAUSGEFILTERT
  filter(level %in% c("50", "100")) %>%
  filter(issues %in% c("0"))


lmALL <- lmer(log(SumAll + 1) ~ temp*level + temp*shape + shape*level + (1|week_ID) + (1|site/plot) , data=RR)


## contrasts
### individual treatments
All_k <- estimate_means(lmALL, at = "temp")
estimate_contrasts(lmALL, contrast = "temp")

All_l <- estimate_means(lmALL, at = "level")
estimate_contrasts(lmALL, contrast = "level")

All_d <- estimate_means(lmALL, at = "shape")
estimate_contrasts(lmALL, contrast = "shape")

### interactions
All_kl <- estimate_means(lmALL, at = "temp:level")
estimate_contrasts(lmALL, at = "temp", contrast ="level")

All_dl <- estimate_means(lmALL, at = "shape:level")
estimate_contrasts(lmALL, at = "shape", contrast ="level")

All_kd <- estimate_means(lmALL, at = "temp:shape")
estimate_contrasts(lmALL, at = "temp", contrast ="shape")


#### PT - All Groups (excl. control) ----
# load data 
R <- read.csv("1_Data/PT/ALANeX_PT_all_clean_control.csv", sep = ",")

RR <- R %>% 
  filter(temp %in% c("4000", "3000", "2000")) %>% 
  filter(level %in% c("50", "100")) %>%
  filter(issues %in% c("0"))

lmPT <- lmer(log(SumAll + 1) ~ temp*level + temp*shape + shape*level + (1|week_ID) + (1|site/plot) , data=RR)

## contrasts
##### individual treatments
PT_k <- estimate_means(lmPT, at = "temp")
estimate_contrasts(lmPT, contrast = "temp")

PT_l <- estimate_means(lmPT, at = "level")
estimate_contrasts(lmPT, contrast = "level")

PT_d <- estimate_means(lmPT, at = "shape")
estimate_contrasts(lmPT, contrast = "shape")

##### interactions
PT_kl <- estimate_means(lmPT, at = "temp:level")
estimate_contrasts(lmPT, at = "temp", contrast ="level")

PT_dl <- estimate_means(lmPT, at = "shape:level")
estimate_contrasts(lmPT, at = "shape", contrast ="level")

PT_kd <- estimate_means(lmPT, at = "temp:shape")
estimate_contrasts(lmPT, at = "temp", contrast ="shape")






#### TRF - All Groups (excl. control) ----
# load data 
D <- read.csv("1_Data/TRF/ALANeX_TRF_all_clean_control.csv", sep = ",")

DD <- D %>% 
  filter(temp %in% c("2000", "3000", "4000")) %>%     ##### OHNE KONTROLLE; KONTROLLE RAUSGEFILTERT
  filter(level %in% c("50", "100")) %>%
  filter(issues %in% c("0"))

lmTRF <- lmer(log(SumAll + 1) ~ temp * level + temp * shape + shape * level + (1|year/week_ID) + (1|site/plot) , data=DD)


## contrasts
#### individual treatments 
TRF_k <- estimate_means(lmTRF, at = "temp")
estimate_contrasts(lmTRF, contrast = "temp")

TRF_l <- estimate_means(lmTRF, at = "level")
estimate_contrasts(lmTRF, contrast = "level")

TRF_d <- estimate_means(lmTRF, at = "shape")
estimate_contrasts(lmTRF, contrast = "shape")

#### interactions
TRF_kl <- estimate_means(lmTRF, at = "temp:level")
estimate_contrasts(lmTRF, at = "temp", contrast ="level")

TRF_dl <- estimate_means(lmTRF, at = "shape:level")
estimate_contrasts(lmTRF, at = "shape", contrast ="level")

TRF_kd <- estimate_means(lmTRF, at = "temp:shape")
estimate_contrasts(lmTRF, at = "temp", contrast ="shape")


##### ANOVA Tables & R2 & AIC ----
## incl. Control----
# PT & TRF
x <- Anova(lmAllc, type=3)
r.squaredGLMM(lmAllc)
AIC(lmAllc)

anova_Allc <- as.data.frame(x)
write_xlsx(anova_Allc, "PT&TRFc.xlsx")

# PT
x <- Anova(lmPTc, type=3)
r.squaredGLMM(lmPTc)
AIC(lmPTc)

anova_PTc <- as.data.frame(x)
write_xlsx(anova_PTc, "PTc.xlsx")

# TRF
x <- Anova(lmTRFc, type=3)
r.squaredGLMM(lmTRFc)
AIC(lmTRFc)

anova_TRFc <- as.data.frame(x)
write_xlsx(anova_TRFc, "TRFc.xlsx")


## excl. Control----
# PT & TRF
x <- Anova(lmALL, type=3)
r.squaredGLMM(lmALL)
AIC(lmALL)

anova_ALL <- as.data.frame(x)
write_xlsx(anova_ALL, "PT&TRF.xlsx")

# PT
x <- Anova(lmPT, type=3)
r.squaredGLMM(lmPT)
AIC(lmPT)

anova_PT <- as.data.frame(x)
write_xlsx(anova_PT, "PT.xlsx")

# TRF
x <- Anova(lmTRF, type=3)
r.squaredGLMM(lmTRF)
AIC(lmTRF)

anova_TRF <- as.data.frame(x)
write_xlsx(anova_TRF, "TRF.xlsx")


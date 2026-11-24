my_colors <- c("4000" = "#FFF59D",   
               "3000" = "#FDD835",   
               "2000" = "#FF6F00",   
               "control" = "#212121",
               "50" = "#3F51B5",
               "100" = "#9FA8DA",
               "diff" = "#2E7D32",
               "std" = "#A5D6A7")  


# load packages ----
library(permute)
library(lattice)
library(vegan)
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

# filter dataset ----
NMDS_TRF <- read.csv("1_Data/NMDS/AlANeX_TRF_NMDS.csv")
NMDS_TRF <- NMDS_TRF %>% filter(temp %in% c("2000", "3000", "4000")) %>% 
  filter(level %in% c("100", "50")) %>% 
  filter(issues %in% c("0"))

## year and week_ID as factors ----
NMDS_TRF$year <- as.factor(NMDS_TRF$year)
NMDS_TRF$week_ID <- as.factor(NMDS_TRF$week_ID)

## remove all rows with 0 ----
AL_TRF <- NMDS_TRF[,10:17]
tot <- rowSums(AL_TRF)
AL_TRF <- AL_TRF[tot > 0, ]

## create a matrix for the nmds analysis ----
AL_TRFm <- data.matrix(AL_TRF)
#AL_TRFm <- log(AL_TRFm + 1)
nmds_TRF <- metaMDS(AL_TRFm, distance = "bray", trymax = 200)

stressplot(nmds_TRF)
plot(nmds_TRF)

#### remove all rows with 0 from NMDS_TRF dataset to then merge the original dataset with the nmds values calculated ----
NMDS_TRF <- NMDS_TRF[tot > 0, ]
nmds_TRF_scores <- as.data.frame(scores(nmds_TRF)$sites)
NMDS_TRF <- data.frame(NMDS_TRF, nmds_TRF_scores)


NMDS_TRF$temp <- factor(NMDS_TRF$temp, levels = c("4000", "3000", "2000","control","dark"))
NMDS_TRF$shape <- factor(NMDS_TRF$shape, levels = c("std", "diff", "control", "dark"))

##### temp plot ----
TRFtemp <- ggplot(NMDS_TRF, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = temp)) + 
  stat_ellipse(aes(colour = temp), type = "norm", level = 0.95, linewidth = 1) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = c(0.99, 0.99), 
        legend.justification = c(0.99, 0.99),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 13, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  #theme(legend.background = element_rect(colour = "black", fill = "white")) + 
  labs(x = "NMDS1", colour = "DL", y = "NMDS2", shape = "Type")  + 
  scale_colour_manual(values = my_colors,  name = "Color Temp.", labels = c("4000 K", "3000 K", "2200 K"))+
  xlim(-0.6, 0.6) + ylim(-0.6, 0.6)

##### level plot ----
TRFlevel <- ggplot(NMDS_TRF, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = level)) + 
  stat_ellipse(aes(colour = level), type = "norm", level = 0.95, linewidth = 1) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = c(0.99, 0.99), 
        legend.justification = c(0.99, 0.99),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 13, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  #theme(legend.background = element_rect(colour = "black", fill = "white")) + 
  labs(x = "NMDS1", colour = "DL", y = "NMDS2", shape = "Type")  + 
  scale_colour_manual(values = my_colors,  name = "Light level", labels = c("100 %", "50 %"))+
  xlim(-0.6, 0.6) + ylim(-0.6, 0.6)

##### shape plot ----
TRFshape <- ggplot(NMDS_TRF, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = shape)) + 
  stat_ellipse(aes(colour = shape), type = "norm", level = 0.95, linewidth = 1) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = c(0.99, 0.99), 
        legend.justification = c(0.99, 0.99),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 13, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  #theme(legend.background = element_rect(colour = "black", fill = "white")) + 
  labs(x = "NMDS1", colour = "DL", y = "NMDS2", shape = "Type")  + 
  scale_colour_manual(values = my_colors,  name = "Lum. shape", labels = c("standard", "diffused")) +
  xlim(-0.6, 0.6) + ylim(-0.6, 0.6)

figure1 <- ggarrange(TRFtemp + TRFlevel + TRFshape) +
  theme(plot.margin = unit(c(1, 1, 1.5, 1), "cm"))

ggarrange( figure1, figure2,  ncol = 1,  nrow = 2) +
  theme(plot.margin = unit(c(4, 2, 2, 3), "cm"))

ggsave(path = "3_Plots", "TRF_Comboplot_NMDS_StackedBar.pdf", width = 55, height = 45, units = c("cm"), dpi = 600)


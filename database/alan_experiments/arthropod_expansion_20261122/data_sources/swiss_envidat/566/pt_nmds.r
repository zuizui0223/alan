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

conflict_prefer("filter", "dplyr")
conflict_prefer("margin", "ggplot2")
conflict_prefer("summarise", "dplyr")

# filter dataset ----
NMDS_PT <- read.csv("1_Data/NMDS/ALANeX_PT_NMDS.csv")
NMDS_PT <- NMDS_PT %>% filter(temp %in% c("2000", "3000", "4000")) %>% 
  filter(level %in% c("100", "50")) %>% 
  filter(issues %in% c("0"))

## year and week_ID as factors ----
NMDS_PT$year <- as.factor(NMDS_PT$year)
NMDS_PT$week_ID <- as.factor(NMDS_PT$week_ID)

## remove all rows with 0 ----
AL_PT <- NMDS_PT[,9:16]
tot <- rowSums(AL_PT)
AL_PT <- AL_PT[tot > 0, ]

## create a matrix for the nmds analysis ----
AL_PTm <- data.matrix(AL_PT)
AL_PTm <- na.omit(AL_PTm)
nmds_PT <- metaMDS(AL_PTm, distance = "bray")

stressplot(nmds_PT)
plot(nmds_PT)

#### remove all rows with 0 from ALAN dataset to then merge the original dataset with the nmds values calculated ----
NMDS_PT <- NMDS_PT[tot > 0, ]
nmds_PT_scores <- as.data.frame(scores(nmds_PT)$sites)
NMDS_PT <- data.frame(NMDS_PT, nmds_PT_scores)

NMDS_PT$temp <- factor(NMDS_PT$temp, levels = c("4000", "3000", "2000","control","dark"))
NMDS_PT$shape <- factor(NMDS_PT$shape, levels = c("std", "diff", "control", "dark"))

##### temp plot ----
PTtemp <- ggplot(NMDS_PT, aes(x = NMDS1, y = NMDS2)) + 
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
  xlim(-0.8, 0.8) + ylim(-0.8, 0.8)

##### level plot ----
PTlevel <- ggplot(NMDS_PT, aes(x = NMDS1, y = NMDS2)) + 
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
  xlim(-0.8, 0.8) + ylim(-0.8, 0.8)

##### shape plot ----
PTshape <- ggplot(NMDS_PT, aes(x = NMDS1, y = NMDS2)) + 
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
  xlim(-0.8, 0.8) + ylim(-0.8, 0.8)


figure1 <- ggarrange(PTtemp + PTlevel + PTshape)

ggarrange( figure1, figure2,  ncol = 1,  nrow = 2) +
  theme(plot.margin = unit(c(4, 2, 2, 3), "cm"))

ggsave(path = "3_Plots", "XX.pdf", width = 55, height = 40, units = c("cm"), dpi = 600)

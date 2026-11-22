library(vegan)

NMDS_PT <- read.csv("1_Data/NMDS/ALANeX_PT_NMDS.csv")
NMDS_PT <- NMDS_PT %>% filter(temp %in% c("2000", "3000", "4000")) %>% 
  filter(level %in% c("100", "50")) %>% 
  filter(issues %in% c("0"))

## year and week_ID as factors ----
NMDS_PT$year <- as.factor(NMDS_PT$year)
NMDS_PT$week_ID <- as.factor(NMDS_PT$week_ID)

## remove all rows with 0 ----
AL_PT <- NMDS_PT[,9:15]
tot <- rowSums(AL_PT)
AL_PT <- AL_PT[tot > 0, ]

## create a matrix for the nmds analysis ----
AL_PTm <- data.matrix(AL_PT)
AL_PTm <- na.omit(AL_PTm)

# Environmental data: select the relevant columns
environmental_data <- NMDS_PT %>% select(temp, level, shape)


# Perform PERMANOVA
adonis2(AL_PTm ~ temp + level + shape, 
                            data = environmental_data, 
                            method = "bray", 
                            permutations = 999)


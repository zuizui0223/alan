library(vegan)

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

# Environmental data: select the relevant columns
environmental_data <- NMDS_TRF %>% select(temp, level, shape)


# Perform PERMANOVA
adonis2(AL_TRFm ~ temp + level + shape, 
                            data = environmental_data, 
                            method = "bray", 
                            permutations = 999)


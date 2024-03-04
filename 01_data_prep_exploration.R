### Data cleaning + exploration script ###

pacman::p_load( 'rstan', 'brms','tidybayes', 'dplyr','bayesplot', 'rstanarm', 
                'kableExtra', 'patchwork', 'loo', 'emmeans', 'broom', 'broom.mixed', 'cmdstanr',
                'GGally', 'MetBrewer', 'rcartocolor', 'ggrepel', 'devEMF', 'flextable', 'officer',
                'grid', 'ggplot2', 'forcats', 'tidyverse')

setwd("C:/Users/msmi0005/Documents/FLX_mosquitofish_mesocosm") #FOR PC
#setwd("/Users/marcusmichelangeli/Desktop/FLX_mosquitofish_mesocosm/")

#set file path
data_path = "./data/"

#-------------------------------------------------------------------------#
#> 1. Fish size data ######################################################
#-------------------------------------------------------------------------#

#load data
fish_size_data <- read.csv(paste0(data_path, "fish_size_data.csv"))    
glimpse(fish_size_data)

#some missing values in Treatment column
fish_size_data <- fish_size_data %>%
  mutate(Treatment = case_when(
    Mesocosm == 'T1' ~ 'High',
    Mesocosm == 'T3' ~ 'Control',
    Mesocosm %in% c('T7', 'T11') ~ 'Low',
    TRUE ~ Treatment
  ))

#change ordering of factor levels
fish_size_data$Treatment <- forcats::fct_relevel(fish_size_data$Treatment, "Control", "Low", "High")
levels(fish_size_data$Treatment)
fish_size_data$Stage <- forcats::fct_relevel(fish_size_data$Stage, "Pre-treatment", "Post-treatment")
levels(fish_size_data$Stage)


#Get average size of fish pre-treatment for each mesocosm
fish_meso_size <- 
  fish_size_data %>% 
  group_by(Mesocosm) %>% 
  mutate(Avg_length = mean(Length, na.rm = T),
         Avg_weight = mean(Weight, na.rm = T)) %>% 
  ungroup() %>% 
  filter(Stage == 'Pre-treatment') %>% 
  group_by(Mesocosm) %>%
  reframe(Avg_length = Avg_length,
          Avg_weight = Avg_weight,
          Avg_length_pre = mean(Length, na.rm = T),
          Avg_weight_pre = mean(Weight, na.rm = T)) %>% 
  distinct(Mesocosm, .keep_all = TRUE)

#save file as RDS
saveRDS(fish_size_data, file = paste0(data_path, "fish_size_data.rds"))


#--------------------------------------------------------------------------#
#> 2. Algae data ##########################################################
#--------------------------------------------------------------------------#

#load data
algae_data <- read.csv(paste0(data_path, "algae_data.csv"), sep = ";")
glimpse(algae_data)

#change ordering of factor levels
algae_data$Stage <- forcats::fct_relevel(algae_data$Stage, "Establishment", "Fish only", "Fish_Exposure", "Exposure only")
levels(algae_data$Stage)

algae_data$Treatment <- forcats::fct_relevel(algae_data$Treatment, "Control", "Low", "High")
levels(algae_data$Treatment)

#convert Date
algae_data$Date <- as.Date(algae_data$Date, format = "%d/%m/%Y")

#create a new column called sample_side to indicate on which side a sample was taken from in a mesocosm
algae_data <- algae_data %>%
  mutate(Sample_side = case_when(
    stringr::str_detect(stringr::str_sub(Sample_ID, 3, 3), "L") | stringr::str_detect(stringr::str_sub(Sample_ID, 4, 4), "L") ~ "Left",
    stringr::str_detect(stringr::str_sub(Sample_ID, 3, 3), "R") | stringr::str_detect(stringr::str_sub(Sample_ID, 4, 4), "R") ~ "Right",
    TRUE ~ NA_character_
  ))

#Create a new column to define the sample number, where the number refers to the date order
algae_data <- algae_data %>% 
  mutate(Sample_date = dense_rank(Date)) %>% 
  mutate(Sample_ID = paste0(Sample_ID, "_", Sample_date))

algae_data %>% distinct(length(Sample_ID)) #we have no overlapping IDs

#Now I'm going to get the average values of each algae parameter for each mesocosm per sampling date (i.e. get the average between the left and right tile)
algae_data <- algae_data %>%
  group_by(Mesocosm, Date, Sample_type) %>%
  mutate(GPP_avg = mean(GPP),
         Chl_a_avg = mean(Chl_a),
         Chl_b_avg = mean(Chl_b),
         Chl_c_avg = mean(Chl_c),
         CR_avg = mean(CR))

#Add average fish length and weight data for each mesocosm
algae_data <- merge(algae_data, fish_meso_size, by = c('Mesocosm'))

#standardize weight and length
algae_data$Avg_length_std = as.numeric(scale(algae_data$Avg_length))
algae_data$Avg_weight_std = as.numeric(scale(algae_data$Avg_weight))
algae_data$Avg_length_pre_std = as.numeric(scale(algae_data$Avg_length_pre))
algae_data$Avg_weight_pre_std = as.numeric(scale(algae_data$Avg_weight_pre))

#sampling_date should be continuous and standardize
algae_data$Sample_date <- as.numeric(algae_data$Sample_date)
algae_data$Sample_date_std <- as.numeric(scale(algae_data$Sample_date))

#I'm going to create 2 seperate datasets. 1 for samples taken from tiles (biofilm) and one taken from the water column (water). For tile dataset, I will only analyse the average of the right and left tile, so i can remove one side

biofilm_algae_data = algae_data %>%
  filter(Sample_type == 'Tile') %>% 
  filter(Sample_side == 'Right')

water_algae_data = algae_data %>% 
  filter(Sample_type == 'Water column') %>% 
  select(-Sample_side, -NPP)

#save file as RDS
saveRDS(biofilm_algae_data, file = paste0(data_path, "biofilm_algae_data.rds"))
saveRDS(water_algae_data, file = paste0(data_path, "water_algae_data.rds"))

#--------------------------------------------------------------------------#
#> 3. Insect chemical analysis #############################################
#--------------------------------------------------------------------------#

#load data
insect_chem_data <- read.csv(paste0(data_path, "insect_chem_analysis.csv"))
glimpse(insect_chem_data)

#make flx_ng numeric and change <LOQ to NA for coding purposes
insect_chem_data <- insect_chem_data %>% 
  mutate(FLX_ng = as.numeric(ifelse(FLX_ng == "<LOQ", 'NA', FLX_ng)))

#rename treatment variables and re-order them. 
insect_chem_data <- insect_chem_data %>% 
  filter(Treatment != 'Blank' ) %>%  #filter out the blank from the dataset
  mutate(Treatment = recode_factor(Treatment, 
                                   "CTRL" = "Control",
                                   "Low" = "Low",
                                   "High" = "High"))
levels(insect_chem_data$Treatment)

#save file as RDS
saveRDS(insect_chem_data, file = paste0(data_path, "insect_chem_data.rds"))

#--------------------------------------------------------------------------#
#> 4. Water nutrients data #############################################
#--------------------------------------------------------------------------#

#load data
nutrient_data <- read.csv(paste0(data_path, "nutrient_analysis.csv"))
glimpse(nutrient_data)

#remove notes column
nutrient_data <- nutrient_data %>% 
  select(-Notes)

#Fix the Sample_no to remove the hypen from the T-1
nutrient_data$Sample_no <- gsub("T-", "T", nutrient_data$Sample_no)

#Convert date into a date format
nutrient_data$Date <- as.Date(nutrient_data$Date, format = "%d/%m/%Y")

#Add a treatment column
nutrient_data <- nutrient_data %>%
  mutate(Treatment = case_when(
    Sample_no %in% c("T1", "T4", "T9", "T12") ~ "High",
    Sample_no %in% c("T5", "T2", "T11", "T7") ~ "Low",
    Sample_no %in% c("T8", "T6", "T3", "T10") ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "Low", "High")))

levels(nutrient_data$Treatment)

#Add a stage column, indicating at what stage in the experiment the nutrients were sampled
nutrient_data <- nutrient_data %>%
  mutate(Stage = case_when(
    Date >= as.Date("2020-09-25") & Date <= as.Date("2020-10-15") ~ "Establishment",
    Date >= as.Date("2020-10-16") & Date <= as.Date("2020-11-05") ~ "Fish only",
    Date >= as.Date("2020-11-06") & Date <= as.Date("2020-11-26") ~ "Fish_Exposure",
    Date >= as.Date("2020-11-27") & Date <= as.Date("2020-12-17") ~ "Exposure only",
    TRUE ~ NA_character_
  )) %>% 
  mutate(Stage = factor(Stage, 
                        levels = c("Establishment", 
                                   "Fish only", 
                                   "Fish_Exposure",
                                   "Exposure only")))

#Need to change measurements from character to numeric values, Im also adding a column which tells us whether there detectable concentration or not
nutrient_data <- 
  nutrient_data %>% 
  mutate(NH3_mgL_det = ifelse(NH3_mgL == '<0.001', 0, 1),
         FRP_mgL_det = ifelse(FRP_mgL == '<0.001', 0, 1),
         Nox_mgL_det = ifelse(Nox_mgL == '<0.001', 0, 1)) %>% 
  mutate_at(vars(c("NH3_mgL", "FRP_mgL", "Nox_mgL")),
            ~ ifelse(.x == "<0.001", 0.001/sqrt(2), as.numeric(.x))) %>% 
  mutate_at(vars(c("NH3_mgL", "FRP_mgL", "Nox_mgL")), round, digits = 4)

#Need to consider how to deal with non-detections in the analysis. No straight answer. 

nutrient_data <- nutrient_data %>% 
  mutate(Sample_date = dense_rank(Date)) 

#filter out sample data so that we only have the last sample date for each stage of the experiment. This is to maintain consistency between stages  as it precedes dates for the other measured parameters

#note that we have two nutrient sampling points for the fish_only stage. The first sampling point is around the time that we added nutrients so it should be quite high...I want to check this. 
nutrient_data <- nutrient_data %>% 
  filter(!Sample_date == '1') %>% 
  filter(!Sample_date == '2') %>% 
  filter(!Sample_date == '4')

#filter out blanks and JMR
nutrient_data <- nutrient_data %>% 
  filter(!Sample_no == 'Blank') %>% 
  filter(!Sample_no == 'JMR')
#note that sample no in this case also refers to mesocosm number
#I'm going to rename this column for consistency 

nutrient_data <- nutrient_data %>% 
  rename('Mesocosm' = 'Sample_no')

#Add average fish length and weight data for each mesocosm
nutrient_data <- merge(nutrient_data, fish_meso_size, by = c('Mesocosm'))

#standardize weight and length
nutrient_data$Avg_length_std = as.numeric(scale(nutrient_data$Avg_length))
nutrient_data$Avg_weight_std = as.numeric(scale(nutrient_data$Avg_weight))
nutrient_data$Avg_length_pre_std = as.numeric(scale(nutrient_data$Avg_length_pre))
nutrient_data$Avg_weight_pre_std = as.numeric(scale(nutrient_data$Avg_weight_pre))

#sampling_date should be continuous and standardize
nutrient_data$Sample_date <- as.numeric(nutrient_data$Sample_date)
nutrient_data$Sample_date_std <- as.numeric(scale(nutrient_data$Sample_date))

#save file as RDS
saveRDS(nutrient_data, file = paste0(data_path, "nutrient_data.rds"))

#--------------------------------------------------------------------------#
#> 5. Plankton data #############################################
#--------------------------------------------------------------------------#

#load data
plankton_data <- read.csv(paste0(data_path, "plankton_data.csv"))
glimpse(plankton_data)

#convert date to date format
plankton_data$Date <- as.Date(plankton_data$Date, format = "%d/%m/%Y")

#Add sample_date
#Create a new column to define the sample number, where the number refers to the date order
plankton_data <- plankton_data %>% 
  mutate(Sample_date = dense_rank(Date))

#Add a treatment column
plankton_data <- plankton_data %>%
  mutate(Treatment = case_when(
    Mesocosm %in% c("T1", "T4", "T9", "T12") ~ "High",
    Mesocosm %in% c("T5", "T2", "T11", "T7") ~ "Low",
    Mesocosm %in% c("T8", "T6", "T3", "T10") ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "Low", "High")))

levels(plankton_data$Treatment)

#change ordering of Stage factor levels
plankton_data$Stage <- forcats::fct_relevel(plankton_data$Stage, "Establishment", "Fish only", "Fish_Exposure", "Exposure only")
levels(plankton_data$Stage)

#change NAs to 0
plankton_data <- plankton_data %>%
  mutate(across(c(Hyrdidae:Arachinida_Acarina), ~ifelse(is.na(.), 0, .)))

#create a column counting the number of different species and abundance in a sample
plankton_data$num_species_found <- rowSums(plankton_data[, 4:25] > 0)
plankton_data$total_abundance <- rowSums(plankton_data[, 4:25])

#combine columns within the same order into a single 'Order' column
plankton_data_wide <- 
  plankton_data %>%
  pivot_longer(cols = -c(Date, Mesocosm, Stage, Treatment, num_species_found, total_abundance, Sample_date), 
               names_to = "Class", 
               values_to = "Count") %>%
  separate(Class, c("Order"), sep = "_") %>%
  pivot_wider(names_from = Order, values_from = Count, values_fn = sum, values_fill = 0)

head(plankton_data_wide)

#create a column for number of orders found
plankton_data_wide$num_order_found <- rowSums(plankton_data_wide[, 8:21] > 0)

#Add average fish length and weight data for each mesocosm
plankton_data_wide <- merge(plankton_data_wide, fish_meso_size, by = c('Mesocosm'))

#Create a sample id within stage
plankton_data_wide <- plankton_data_wide %>%
  group_by(Stage) %>%
  mutate(Sample_num = 
           ifelse(Stage == "Establishment", ifelse(Sample_date == min(Sample_date), 1, 2),
                  ifelse(Stage == "Fish only", ifelse(Sample_date == min(Sample_date), 1, 2),
                         ifelse(Stage == "Fish_Exposure", ifelse(Sample_date == min(Sample_date), 1, 2),
                                ifelse(Stage == "Exposure only", ifelse(Sample_date == min(Sample_date), 1, 2), NA)))))


#standardize weight and length
plankton_data_wide$Avg_length_std = as.numeric(scale(plankton_data_wide$Avg_length))
plankton_data_wide$Avg_weight_std = as.numeric(scale(plankton_data_wide$Avg_weight))
plankton_data_wide$Avg_length_pre_std = as.numeric(scale(plankton_data_wide$Avg_length_pre))
plankton_data_wide$Avg_weight_pre_std = as.numeric(scale(plankton_data_wide$Avg_weight_pre))

#sampling_date should be continuous and standardize
plankton_data_wide$Sample_date <- as.numeric(plankton_data_wide$Sample_date)
plankton_data_wide$Sample_date_std <- as.numeric(scale(plankton_data_wide$Sample_date))

#save file as RDS
saveRDS(plankton_data_wide, file = paste0(data_path, "plankton_data_wide.rds"))

#--------------------------------------------------------------------------#
#> 6. Water chemistry data  #############################################
#--------------------------------------------------------------------------#

#load data
water_chem_data <- read.csv(paste0(data_path, "water_chem_data.csv"))
glimpse(water_chem_data)

#Fix the Sample to remove the the /XX 
water_chem_data <- 
  water_chem_data %>% 
  separate(Sample, into = c("Mesocosm", "X"), sep = "/") %>% 
  select(-X)

#convert date to date format
water_chem_data$Date <- as.Date(water_chem_data$Date, format = "%d/%m/%Y")

#change ordering of treatment factor level to make more sense
water_chem_data$Treatment <- forcats::fct_relevel(water_chem_data$Treatment, "Control", "Low", "High")
levels(water_chem_data$Treatment)

#Add a stage of Experiment column into the dataset
water_chem_data <- water_chem_data %>%
  mutate(Stage = case_when(
    Date >= as.Date("2020-11-06") & Date <= as.Date("2020-11-26") ~ "Fish_Exposure",
    Date >= as.Date("2020-11-27") & Date <= as.Date("2020-12-17") ~ "Exposure only",
    TRUE ~ NA_character_
  )) %>% 
  mutate(Stage = factor(Stage, 
                        levels = c("Fish_Exposure",
                                   "Exposure only"), 
                        ordered = TRUE))

#save file as RDS
saveRDS(water_chem_data, file = paste0(data_path, "water_chem_data.rds"))

#--------------------------------------------------------------------------#
#> 7. Water parameter data  #############################################
#--------------------------------------------------------------------------#

#load data
water_param_data <- read.csv(paste0(data_path, "water_parameter_data.csv"))
glimpse(water_param_data)  
  
#need to rename the mesocosm column (there is a spelling mistake)
water_param_data <- 
  water_param_data %>% 
  rename(Mesocosm = Mesoccosm)

#convert date to date format
water_param_data$Date <- as.Date(water_param_data$Date, format = "%d/%m/%Y")

#convert time to time format (AEDT)
#first I need to create a datetime column
water_param_data$datetime <-  paste(water_param_data$Date, water_param_data$Time)
water_param_data$Time <- as.POSIXct(water_param_data$datetime, format="%Y-%m-%d %H:%M", tz = "Australia/Sydney")

#Add a treatment column
water_param_data <- water_param_data %>%
  mutate(Treatment = case_when(
    Mesocosm %in% c("T1", "T4", "T9", "T12") ~ "High",
    Mesocosm %in% c("T5", "T2", "T11", "T7") ~ "Low",
    Mesocosm %in% c("T8", "T6", "T3", "T10") ~ "Control",
    TRUE ~ NA_character_
  )) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "Low", "High")))

levels(water_param_data$Treatment)

#change ordering of Stage factor level to make more sense
water_param_data$Stage <- forcats::fct_relevel(water_param_data$Stage, 
                                               "Establishment", "Fish only", "Fish_Exposure", "Exposure only")
levels(water_param_data$Stage)

#I want to create a sample ID to more easily distinguish samples from each other
water_param_data$Week <- paste("W", water_param_data$Week, sep = "")

water_param_data$Sample_ID <- paste(water_param_data$Mesocosm, 
                                    water_param_data$Week, 
                                    water_param_data$Time_of_day,
                                    water_param_data$Probe,
                                    water_param_data$Sample,
                                    sep = "_")

#I might just make a new column called sample_data for consistency between datasets
water_param_data <- 
  water_param_data %>% 
  mutate(Sample_date = dense_rank(as.integer(substring(Week, 2))))

#Some weird DO estimate (way too large)...guessing they are mistakes in data entry. Im going to make them NAs
water_param_data <- 
  water_param_data %>% 
  mutate(DO_mgL = ifelse(DO_mgL > 50, NA, DO_mgL))

#Im going to create 2 datasets
#1. will contain Temp and pH data
#2. the other will contain the DO data which will make wide format to calculate CR

#Get the averages etc.s here for the data analysis
#Need to still calculate respiration rate.

#1. Isolate temp and pH data

temp_ph_data <- water_param_data %>% 
  select(-DO_perc, -DO_mgL, -Turb) %>% 
  filter(!Time_of_day == 'Sunrise') %>%   
  filter(!Time_of_day == 'Sunset') %>% 
  group_by(Treatment, Stage, Mesocosm, Sample_date, Time_of_day) %>% 
  reframe(
    Temp_avg = mean(Temp, na.rm = T),
    pH_avg = mean(pH, na.rm = T)
  )
  

#2. We also need to calculate community respiration using dissolved oxygen measurements at sunset vs sunrise.

DO_CR_data <- 
  water_param_data %>% 
  select(-Temp, -pH, -Turb) %>% 
  #get the average
  pivot_longer(
    cols = c('DO_mgL', 'DO_perc'),
    names_to = 'DO_measure',
    values_to = 'Value') %>% 
  filter(!is.na(Value)) %>% 
  group_by(Treatment, Stage, Mesocosm, Sample_date, Time_of_day, DO_measure) %>% 
  #Get average of the DO measurement for each mesocosm and sampling date
  reframe(
    Value = mean(Value)
  ) %>% 
  pivot_wider(names_from = c(Time_of_day, DO_measure),
              values_from = Value) %>% 
  #Calculate respiration
  mutate(CR_DO_mgL = Sunset_DO_mgL - Sunrise_DO_mgL,
         CR_DO_perc = Sunset_DO_perc - Sunrise_DO_perc)
  
  
#save file as RDS
saveRDS(temp_ph_data, file = paste0(data_path, "temp_ph_data.rds"))
saveRDS(DO_CR_data, file = paste0(data_path, "DO_CR_data.rds"))
  
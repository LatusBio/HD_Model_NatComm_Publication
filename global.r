library(readr)

########################################
# Below are the variables that will be used in the model
########################################

######## HD ISS #############
# Create some logic to calculate and display HDISS stages based
# Read in phenotypic data from EnrollHD
TotalFunctionalCapacity_Score <- read_csv('UHDRS_Total_Functional_Capacity_AIscraped_5yearRange_individualYears.csv')
#TotalFunctionalCapacity_Score[is.na(TotalFunctionalCapacity_Score)] <- 1
Motor_Score <- read_csv('meanMotorScore_AIscraped_5yearRange_individualYears.csv')
#Motor_Score[is.na(Motor_Score)] <- 50
# Read in HD-ISS age associated thresholds 
HD_ISS_Thresholds_ByAge <- read_csv('HD_ISS_Thresholds_ByAge.csv')

# Values for testing that the function works
#CurrentAge <- 49
#cag <- 40
#brainVolume <- .8
#motorScore <- 17
#functionalCapacity <- 13

# Create function to calculate HD-ISS stage
HD_ISS_Stage_Calc_FUN <- function(CurrentAge, cag, brainVolume, motorScore, functionalCapacity) {
  #Pull age related threshold from key
  # Get appropriate row
  HD_ISS_Thresholds_ByAge_trim1 <- HD_ISS_Thresholds_ByAge[CurrentAge <= HD_ISS_Thresholds_ByAge$Age,]
  HD_ISS_Thresholds_ByAge_trim2 <- HD_ISS_Thresholds_ByAge_trim1[CurrentAge + 4 >= HD_ISS_Thresholds_ByAge_trim1$Age,]
  
  # Define stage 0
  if (cag <= 35) {
    return(0)
    #print(0)
  } else if (cag >= 36 & brainVolume >= .95 & motorScore <= HD_ISS_Thresholds_ByAge_trim2$TMS[1] & functionalCapacity >= HD_ISS_Thresholds_ByAge_trim2$TFC[1]) {
    return(0)
    #print(0)
    # Define stage 1
  } else if (cag >= 36 & brainVolume <= .95 & motorScore <= HD_ISS_Thresholds_ByAge_trim2$TMS[1] & functionalCapacity >= HD_ISS_Thresholds_ByAge_trim2$TFC[1]) {
    return(1)
    #print(1)
    # Define stage 2
  } else if (cag >= 36 & brainVolume <= .95 & motorScore > HD_ISS_Thresholds_ByAge_trim2$TMS[1] & functionalCapacity >= HD_ISS_Thresholds_ByAge_trim2$TFC[1]) {
    return(2)
    #print(2)
    # Define stage 3
  } else if (cag >= 36 & brainVolume <= .95 & motorScore > HD_ISS_Thresholds_ByAge_trim2$TMS[1] & functionalCapacity < HD_ISS_Thresholds_ByAge_trim2$TFC[1]) {
    return(3)
    #print(3)
  } else {
    # Fallback: brainVolume >= .95 but motor/functional abnormal → Stage 0
    return(0)
  }
}


##### Therapeutic variables #########
# Estimated reduction in annual rate of repeat expansion for cells transduced with therapeutic.  Example: "4" will yield a four fold lower rate of repeat expansion.
#highlightTransducedCells <- TRUE
#therapy_on <- TRUE
# therapeuticModifier <- 0.5
# transductionRate <- 0.5 # a value of .5 will yield a transduction rate of 50%
# age_at_intervention <- 60
# germlineHDRepeatLength <- 40
# 
# ##### Time variables ##### 
# maxLifetime <- 75
# 
# ##### Brain and Cell variables ##### 
# nNeurons_CaudatePutamen <- 3000
# #nNeurons_CaudatePutamen <- 110000000
# nMSNs_CaudatePutamen <- nNeurons_CaudatePutamen * .925 # it is estimated that between 90 and 95% of all neurons in caudate and putamen are MSNs
# Caudate_Volume_HealthControl_Human <- c(8:9) # Hobbs et al AJNR 2010
# Caudate_Volume_EarlyHD_Human <- c(3:7)

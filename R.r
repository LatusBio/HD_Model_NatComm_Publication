########################################
# Below are the variables that will be used in the model
########################################

##### Therapeutic variables #########
# Estimated reduction in annual rate of repeat expansion for cells transduced with therapeutic.  Example: "4" will yield a four fold lower rate of repeat expansion.
highlightTransducedCells <- TRUE
therapy_on <- TRUE
therapeuticModifier <- 4
transductionRate <- .5 # a value of .5 will yield a transduction rate of 50%
age_at_intervention <- 60
germlineHDRepeatLength <- 40

##### Time variables ##### 
maxLifetime <- 75

##### Brain and Cell variables ##### 
nNeurons_CaudatePutamen <- 3000
#nNeurons_CaudatePutamen <- 110000000
nMSNs_CaudatePutamen <- nNeurons_CaudatePutamen * .925 # it is estimated that between 90 and 95% of all neurons in caudate and putamen are MSNs
Caudate_Volume_HealthControl_Human <- c(8:9) # Hobbs et al AJNR 2010
Caudate_Volume_EarlyHD_Human <- c(3:7)

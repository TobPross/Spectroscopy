
#  Spectroscopy script version 6.8
#
#  The script covers the complete work flow
#  from raw spectral files to predicted values in 12 steps:
#
#  Step 1  Define operation folders
#  Step 2  Splice and convert .asd files to .txt files
#  Step 3  Aggregate repeated scans to single .txt files
#  Step 4  Import spectral data and set calibration parameters
#  Step 5  Creating a testset
#  Step 6  Run calibration
#  Step 7  Select model
#  Step 8  Run validation
#  Step 9  Inspect model parameters
#  Step 10 Outlier detection
#  Step 11 Export Model
#  Step 12 Apply Model(s)

# Load libraries
library(plantspec)
library(Hmisc)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(ggplot2)

##### Step 1: Define operation folders -----------------------------------------
# This script accesses the defined folders and writes on the disk.
# Making a backup of all files before starting is highly recommended.

PATH_SPECTRA_RAW          <- ("C:/Spectroscopy/Spectra_for_calibration_raw")      ## Path that contains raw unprocessed spectral .asd files
PATH_SPECTRA_SPLICED      <- ("C:/Spectroscopy/Spectra_for_calibration_spliced")  ## spliced spectra will be written to this path
PATH_SPECTRA_CALIBRATION  <- ("C:/Spectroscopy/Spectra_for_calibration")          ## averaged spectra will be written to this path
PATH_SPECTRA_PREDICT_RAW  <- ("C:/Spectroscopy/Spectra_to_be_predicted_raw")      ## Path that contains spectral files that 
PATH_SPECTRA_PREDICT      <- ("C:/Spectroscopy/Spectra_to_be_predicted")          ## Path that contains spectral files that 
PATH_MODELS               <- ("C:/Spectroscopy/Prediction_models")                ## Prediction models will be written to this path
PATH_RESULTS              <- ("C:/Spectroscopy/Results")                          ## Predicted values will be written to this path

##### Step 2: Splice and convert .asd files to .txt files ----------------------
# Only necessary if you use ASD data files, otherwise skip this step
library(asdreader)  
library(spectacles)

SPLICING = TRUE   ## Spicing is recommended for instruments with more than one sensor
PLOTTING = TRUE   ## Plotting is recommended to visually quick-check the splicing results
FILES_LIST <- list.files(path = PATH_SPECTRA_RAW)
run <- length(FILES_LIST)

for (i in 1:run) {
  
  my_spectrum <- get_spectra(paste0(PATH_SPECTRA_RAW, "/", FILES_LIST[i]))
  
  if (SPLICING == FALSE) {
    spectral_numbers <- dimnames(my_spectrum)[[2]]
    spectral_data    <- as.numeric(my_spectrum)
    
    output_frame <- data.frame(spectral_numbers, spectral_data)
    FILENAME <- paste0(FILES_LIST[i], ".txt")
    
    if (PLOTTING ==TRUE) {plot(output_frame, main=FILENAME)}
    
     write.table(output_frame, paste0(PATH_SPECTRA_SPLICED, "/", FILENAME) , row.names = FALSE, col.names=FALSE)
  }
  
  if (SPLICING == TRUE) {
    my_wl <- as.numeric(dimnames(my_spectrum)[[2]])
    my_id <- FILES_LIST[i]
    my_nir <- as.numeric(my_spectrum)
    
    my_s <- Spectra(wl = my_wl, nir = my_nir, id = my_id)
    
    # Splicing points for ASD FieldSpec4 = list(c(750, 1000), c(1950, 1800))
    my_s_spliced <- splice(my_s, locations = list(c(750, 1000), c(1950, 1800))) 
    
    spectral_numbers <- my_s_spliced@wl
    spectral_data    <- as.numeric(my_s_spliced@nir[1,])
    output_frame <- data.frame(spectral_numbers, spectral_data)
    FILENAME <- paste0(FILES_LIST[i], ".txt")
    
    # Splicing points are marked in blue 
    if (PLOTTING ==TRUE) {plot(output_frame, main=FILENAME, axes=FALSE)
      axis(side=1, at=seq(350, 2500, by=100))
      axis(side=2, at=seq(0, 1, by=0.1))
      box()
      abline(v=1000, col="blue")
      abline(v=1800, col="blue")
    }
      write.table(output_frame, paste0(PATH_SPECTRA_SPLICED, "/", FILENAME) , row.names = FALSE, col.names=FALSE)
    }
  }

# Check for completion
OUTPUT_LIST <- list.files(path = PATH_SPECTRA_SPLICED)
OUTPUT_LIST <- gsub(".txt", "", OUTPUT_LIST)
if (identical(FILES_LIST, OUTPUT_LIST)==TRUE) {
  print("transcoding complete")  
} else {
  # Generate a list of missing files
  MATCH_LIST <- match(FILES_LIST, OUTPUT_LIST)
  LOST_INDEX <- which(is.na(MATCH_LIST))
  MISSING_FILES <- FILES_LIST[LOST_INDEX]
  write.csv(MISSING_FILES, file=paste0(PATH_RESULTS, "/","missing.csv"))
  print(paste0("not all files could be trascoded, manually check ",PATH_RESULTS, "/","missing.csv"))
}

##### Step 3: Aggregate repeated scans to single .txt files --------------------
# only required if you have more than one spectral file per sample, otherwise skip this step
# requirement: - One folder that contains all spectral files to convert "/PATH_SPECTRA_SPLICED/"
#              - a csv file that contains one entry for each spectral file an information to which sample the file belongs   

  SPECTRA <- readSpectra(filelist = list.files(PATH_SPECTRA_SPLICED, full.names = T), wave_unit = "nm", measurement_unit = "reflectance", sep="", header=F)
 #plot.spectra.matrix(SPECTRA)

    # The reference file must include a File and Sample variable
  DATA_REFERENCE <- read.csv("C:/Spectroscopy/Reference/Calibration_List.csv", header = T, sep="")
    head(DATA_REFERENCE) # DATA_REFERENCE[,1] should contain the file names; DATA_REFERENCE[,2] should contain the sample to which the file belongs

# Check if the Reference matches the Spectra    
  identical(DATA_REFERENCE[,1], dimnames(SPECTRA)[[1]]) # TRUE

# Averaging spectra on $Sample level  
  SPECTRA <- as.spectra.matrix(SPECTRA)
  SPECTRA_AV  <- averageSpectra(spec=SPECTRA, by=DATA_REFERENCE[,2])

  writeSpectra(SPECTRA_AV, path=PATH_SPECTRA_CALIBRATION)

##### Step 4:  Import spectral data and set calibration parameters -------------
# requirement: - One folder that contains all spectral files (1 spectral file = 1 sample, spliced and aggregated if necessary) "/PATH_SPECTRA_CALIBRATION/"
#              - One csv file that contains one entry for each spectral file with the according reference value
#              - Samples and reference values must be in the same order 

SPECTRA <- readSpectra(filelist = list.files(PATH_SPECTRA_CALIBRATION, full.names = T), wave_unit = "nm", measurement_unit = "reflectance")
 
# manually check spectral data
plot.spectra.matrix(SPECTRA)

# check for missing spectra
any(is.na(SPECTRA)) #FALSE

# Import Reference Dataset
DATA_REFERENCE <- read.csv("C:/Spectroscopy/Reference/Reference_Traits.csv", sep=";", header=T)
# check which reference data is available
head(DATA_REFERENCE)

# Create new Variable for clarity
REFERENCE <- DATA_REFERENCE[,11] ## Select the reference variable here
REFERENCE                        ## Double check if the correct reference variable is selected

REFERENCE[1:100]

# check variable distribution
qqnorm(REFERENCE, pch = 1, frame = FALSE)
qqline(REFERENCE, col = "red", lwd = 2)
hist(REFERENCE)


# Set Parameter for region list
segments_number    <- 6     ## Meaningful values are 3 - 10 (6 is a good solution for most applications) Increasing this value above 8 is not recommended. Doing so will result in high system saturation and system instability.
smallest_segment   <- 72    ## Meaningful values are 42 - 250 (72 is a good solution for most applications)
minimum_frequency  <- 350   ## Starting frequency of the spectrometer (350 for ASD FieldSpec4 #  3800 for Bruker MPA)
maximum_frequency  <- 2500  ## Ending frequency of the spectrometer (2500 for ASD FieldSpec4 # 12400 for Bruker MPA)
maximum_rank       <- 12    ## Meaningful values are 5-15, 12 is a good solution for most applications
wavesequence <- seq(from =minimum_frequency+smallest_segment, to=maximum_frequency-smallest_segment, by=1)

# Define function to sample and validate frequencies
sample_valid_frequencies <- function() {
  wavesequence <- seq(from = minimum_frequency + smallest_segment, 
                        to = maximum_frequency - smallest_segment, 
                        by = 1)
  repeat {
    custom_sequence <- sample(wavesequence, segments_number - 1)
    custom_sequence <- sort(custom_sequence)
    if (all(diff(custom_sequence) > smallest_segment)) {
      return(c(minimum_frequency, custom_sequence, maximum_frequency))
    }
  }
}

##### Step 5: Create a test-set ------------------------------------------------
# Select a training set for model fitting. Will give a TRUE/FALSE vector.
# specifying TRUE for training/calibration data and FALSE for test/validation set data.
# p defines the percentage of VALIDATION samples: 0.3 will give 30% validation and 70% calibration samples.
trainingset <- !(subdivideDataset(spectra = SPECTRA, 
                                  method = "PCAKS",            ## Validation Method. KS, PCAKS or random are good solutions for most applications.
                                  component = REFERENCE,       ## Component is required for several methods.
                                  p = 0.5,                     ## Size of validation Set. Meaningful values are 0.1 - 0.9 (0.5 is a good solution for most applications).
                                  type = "validation",
                                  output="logical"))

#Manually check calibration set distribution
DATA_REFERENCE$Level_ID[trainingset]

##### Step 6: Run calibration --------------------------------------------------
# Create cluster for parallel computing and select number of parallel tested iteration (number of parallel cores)
cores=detectCores(logical = FALSE)   ## Automatically detect the number of physical cores
cl <- makeCluster(cores-1)           ## The use of (cores-1) is recommended. Alternatively insert the number of  cores here.
registerDoParallel(cl)               ## Register Cluster
on.exit(stopCluster(cl))             ## Stop Cluster on exit

# Select number of total tested iterations
n <- (cores-1) * 5


Sys.time()

# start calibration 
finalMatrix <- foreach(i=1:n, .combine=rbind, .packages=c('foreach', 'doParallel', 'dplyr', 'tidyr', 'plantspec', 'Hmisc')) %dopar% {
  
  # Create random sequence of NIR regions within the preset parameter and assemble region list
   full_sequence <- sample_valid_frequencies()
   region_list <- mapply(function(x, y) c(x, y), 
                        head(full_sequence, -1), 
                        tail(full_sequence, -1), 
                        SIMPLIFY = FALSE)
  
  # Optimize model parameter 
  Optimized_Par <-  optimizePLS(component = REFERENCE,
                                spectra = SPECTRA,
                                max_comps=maximum_rank,                   
                                region_list = region_list,
                                training_set = trainingset,
                                parallel = FALSE)
  
  # Output to finalMatrix
  return(Optimized_Par)
}

#stop cluster
stopCluster(cl)
Sys.time()

##### Step 7: Select model -----------------------------------------------------
# Best model is selected by lowest RMSE
run <- dim(finalMatrix)

# Use this section for test-set validation
best_RMSE <- finalMatrix[1,1][[1]]$RMSEP[1]   
for (i in 1:run[1])
{
  if (finalMatrix[i,1][[1]]$RMSEP[1] <= best_RMSE){
    best_RMSE <- finalMatrix[i,1][[1]]$RMSEP[1]
    
    a <- finalMatrix[i,1][[1]]
    b <- finalMatrix[i,2][[1]]
    c <- finalMatrix[i,3][[1]]
    Selected_Model <- list(a,b,c)
    names(Selected_Model) <- c("optimization_results", "param_subsets", "param_preproc")
  }
}

# Use this section for cross validation
'best_RMSE <- finalMatrix[1,1][[1]]$RMSECV[1]
for (i in 1:run[1])
{
  if (finalMatrix[i,1][[1]]$RMSECV[1] <= best_RMSE){
    best_RMSE <- finalMatrix[i,1][[1]]$RMSECV[1]
    
    a <- finalMatrix[i,1][[1]]
    b <- finalMatrix[i,2][[1]]
    c <- finalMatrix[i,3][[1]]
    Selected_Model <- list(a,b,c)
    names(Selected_Model) <- c("optimization_results", "param_subsets", "param_preproc")
  }
}'

##### Step 8:  Model validation ------------------------------------------------
# Select parameter for Model Calibration 
Best_Model <- calibrate(component = REFERENCE, 
                        spectra = SPECTRA, 
                        optimal_params = Selected_Model, 
                        optimal_model = 1,                 # use optimal_model = 1 for model with lowest RMSE. You might want to manually check the top 10 models and decide to chose a different one
                        validation = "testset",           
                        training_set = trainingset,        
                        max_comps = maximum_rank)          # it is possible to manually change the model rank for test purposes at this point 

##### Step 9:  Inspect model parameters ----------------------------------------
# Inspect model parameters
Best_Model$R2_Cal
Best_Model$R2_Val
Best_Model$RMSEP
Best_Model$rank               
Best_Model$regions
Best_Model$preproc
#Best_Model$data$component

# View Model
plot.PLScalibration(Best_Model, ncomp = Best_Model$rank)    # check calibration accuracy select value from Best_Model$rank here, or modify carefully

PERDICTIONS <- predictPLS(object= Best_Model, newdata = SPECTRA)  # Combined calibration and validation plot
plot (as.numeric(PERDICTIONS) ~ REFERENCE)                        #
abline(lm((as.numeric(PERDICTIONS) ~ REFERENCE)))                 # 
abline(a=0, b=1, col="blue")                                      # adding x=y line for comparison

summary(lm(as.numeric(PERDICTIONS) ~ REFERENCE))$r.squared        # combined R2 of calibration and validation. For quick assessment only. Typically R2_Cal and R2_Val are better quality indicators.
coef(lm(as.numeric(PERDICTIONS) ~ REFERENCE))[1]


# Combined calibration and validation plot, seperated by species
SPECTRA_NEW <- readSpectra(filelist = list.files(PATH_SPECTRA_CALIBRATION, full.names = T), wave_unit = "nm", measurement_unit = "reflectance")
PERDICTIONS <- predictPLS(object= Best_Model, newdata = SPECTRA_NEW)  
Species <- (DATA_REFERENCE$Species)
plot_frame <-data.frame(PERDICTIONS, Best_Model$data$component, Best_Model$training_set, Species)
names(plot_frame) <- c("Predicted", "Measured", "Calibration", "Species")
x1 <- min(plot_frame$Measured)
x2 <- max(plot_frame$Measured)

library(ggplot2)
p <- ggplot(plot_frame, aes(Measured, Predicted))
p + geom_point(aes(colour = Species, shape=Calibration),  size=0.7)  +
  geom_segment(aes(x = x1, y = x1, xend = x2, yend = x2), colour="blue", linetype = 2, size=0.2) +
  geom_smooth(aes(colour = Calibration),  method = "lm", se = FALSE, size=0.4) +
  scale_shape_manual(values = c(17, 16)) +
  theme_classic() +
  scale_x_continuous(name="Measured") +
  scale_y_continuous(name="Predicted")

# Check optimal model Rank
SPECTRA_CAL <- readSpectra(filelist = list.files(PATH_SPECTRA_CALIBRATION, full.names = T)[trainingset], wave_unit = "nm", measurement_unit = "reflectance")
REFERENCE_CAL <- REFERENCE[trainingset]
SPECTRA_VAL <- readSpectra(filelist = list.files(PATH_SPECTRA_CALIBRATION, full.names = T)[!trainingset], wave_unit = "nm", measurement_unit = "reflectance")
REFERENCE_VAL <- REFERENCE[!trainingset]

Ranks_Overview <- matrix(nrow=maximum_rank, ncol = 5)
colnames(Ranks_Overview) <-c("Rank","R_SQ_Cal","RMSE_CAL", "R_SQ_VAL", "RMSE_VAL")

for (i in 1:maximum_rank){
  
  k <- i+1
  Ranks_Overview[i,1] <- k    # cannot create model with rank = 1
  Best_Model$rank   <- k
  
  PERDICTIONS_CAL       <- predictPLS(object= Best_Model, newdata = SPECTRA_CAL)                           # Predict Calibration data with Rank = i
  Ranks_Overview[i,2] <- (summary(lm(as.numeric(PERDICTIONS_CAL) ~ REFERENCE_CAL))$r.squared)              # get R2 for Calibration data
  residuals_cal         <- as.vector(resid(summary.lm(lm(as.numeric(PERDICTIONS_CAL) ~ REFERENCE_CAL))))   # extract residuals from lm
  Ranks_Overview[i,3] <- sqrt(mean(residuals_cal^2))                                                       # get RMSE for Calibration data
  
  PERDICTIONS_VAL       <- predictPLS(object= Best_Model, newdata = SPECTRA_VAL)                           # Predict Validation data with Rank = i
  Ranks_Overview[i,4]    <- (summary(lm(as.numeric(PERDICTIONS_VAL) ~ REFERENCE_VAL))$r.squared)           # get R2 for Validation data
  residuals_val         <- as.vector(resid(summary.lm(lm(as.numeric(PERDICTIONS_VAL) ~ REFERENCE_VAL))))   # extract residuals from lm
  Ranks_Overview[i,5] <- sqrt(mean(residuals_val^2))                                                       # get RMSE for Validation data
}

library("reshape2")
Ranks_Overview <- as.data.frame(Ranks_Overview)
Ranks_Overview2 <- melt(Ranks_Overview, id="Rank")
Ranks_Overview3 <- cbind(Ranks_Overview2, substr(Ranks_Overview2[,2], 1, 4))
colnames(Ranks_Overview3)[4] <- "category"

ggplot(data=Ranks_Overview3,
       aes(x=Rank, y=value, colour=variable)) +
  geom_line()+ 
  facet_wrap(vars(category), scales = "free")

##### Step 10: Outlier detection -----------------------------------------------
# Removing outliers is not recommended without further indication. 
# This has to be done manually by removing them from the INPUT_PATH and the Reference.csv file.

DIFF <- abs(as.numeric(PERDICTIONS) - REFERENCE)
FValue_DIFF <- DIFF
DIFF_2 <- DIFF^2 

Limit999 <- qf(.999, df1=1, df2=length(DIFF))
Limit99  <- qf(.99, df1=1, df2=length(DIFF))
Limit95  <- qf(.95, df1=1, df2=length(DIFF)) 
Limit90  <- qf(.90, df1=1, df2=length(DIFF))

for(i in 1:length(DIFF)){
  Value1 <- (length(DIFF) -1) * DIFF[i]^2
  Value2 <- (sum(DIFF_2) - DIFF[i])
  
  FValue_DIFF[i] <- Value1 / Value2}

plot(FValue_DIFF)

abline(h=Limit999,col="red")   # alpha = 0.999
abline(h=Limit99, col="blue")  # alpha = 0.99
abline(h=Limit95, col="green") # alpha = 0.95
abline(h=Limit90, col="black") # alpha = 0.90

outlier_index_999 <- (which(FValue_DIFF>Limit999))
outlier_index_99  <- (which(FValue_DIFF>Limit99))
outlier_index_95  <- (which(FValue_DIFF>Limit95))
outlier_index_90  <- (which(FValue_DIFF>Limit90))

FValue_DIFF[outlier_index_999]
FValue_DIFF[outlier_index_99]
FValue_DIFF[outlier_index_95]
FValue_DIFF[outlier_index_90]

# Names of outlier Spectra
dimnames(SPECTRA)[[1]][outlier_index_999]
dimnames(SPECTRA)[[1]][outlier_index_99]
dimnames(SPECTRA)[[1]][outlier_index_95]
dimnames(SPECTRA)[[1]][outlier_index_90]


##### Step 11: Export Model ----------------------------------------------------
save(Best_Model, file = paste0(PATH_MODELS, "/", "Prediction_v1.RData"))

##### Step 12 Apply prediction models to spectral files ------------------------
# requirement: - One folder that contains all spectral files on which the prediction models will be applied "/PATH_SPECTRA_PREDICT/"
#              - One folder that contains all prediction models  "/PATH_MODELS/"
#              - One folder where the results will be written "/PATH_RESULTS/"
# It is possible to apply only one or multiple prediction models at once  

library(devtools)
library(plantspec)
library(Hmisc)

# read spectral files
SPECTRA <- readSpectra(filelist = list.files(PATH_SPECTRA_PREDICT, full.names = T), wave_unit = "nm", measurement_unit = "reflectance", sep="")

plot.spectra.matrix(SPECTRA)

MODELS_LIST <- list.files(path = PATH_MODELS)
FILES_LIST <- list.files(path = PATH_SPECTRA_PREDICT) 

Predictions_Matrix <- matrix(nrow = length(FILES_LIST), ncol = length(MODELS_LIST))
colnames(Predictions_Matrix) <- as.list(MODELS_LIST)
rownames(Predictions_Matrix) <- as.list(FILES_LIST)
Run <- length(MODELS_LIST)

for (i in 1:Run){
  load(paste0(PATH_MODELS,"/", MODELS_LIST[i]))
  PERDICTIONS <- predictPLS(object= Best_Model, newdata = SPECTRA)
  Predictions_Matrix[,i] <- as.numeric(PERDICTIONS)
}

write.csv(Predictions_Matrix, file=paste0(PATH_RESULTS,"/Results.csv"), row.names = TRUE, na = "NA")

##### End of script ------------------------------------------------------------
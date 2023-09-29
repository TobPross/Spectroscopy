<h2 align="center"> Performance optimization</h2>

## General remarks

This script will run multiple parallel iterations of partial least squares regression (plsr) with randomly chosen variables to find the best calibration model. The runtime of the script depends greatly on the complexity of the dataset, in particular on the number of calibration and validation samples, the structure of the wave sequence and the number of ranks of the calibration model. It is important to note that allocating more time and calculation power does not necessarily increase the model quality but it might help to find the optimal model. Allocating more resources to the calibration will converge the found model slowly to the optimal model. 
  
### Calibration type (Test-set validation vs cross-validation) and training set 
Test-set validation is significantly faster in calibration and validation than cross-validation (approximately factor 5-10, depending on the dataset). Additionally, test-set validation has a lower risk of accidentally overfitting the model. However, cross-validation is a very good solution for small sample sizes (n<80). 

 ```md
 trainingset <- !(subdivideDataset(spectra = SPECTRA, 
                                   method = "random", 
                                   component = REFERENCE, 
                                   p = 0.5, 
                                   type = "validation", 
                                   output="logical")) 
```
The command "subdivideDataset" divides the dataset in training set (calibration samples) and test set (validation samples). The calibration set should account for the complete variation in your dataset, as well as for the variation that is expected in future predictions. There is no consensus on which is the best selection method for a test set and what is the optimal ratio of calibration to validation samples are. Kennard Stones algorithm and PCA are commonly used (here 'KS' and 'PCAKS'). Using 'random' is also legit, in particular it would be a good idea to use multiple 'random' sets once you found the optimal model parameters (see below) in order to test model stability. If you analyse your dataset for multiple traits, it is advised to use the same calibration set for all traits. If your dataset contains multiple species, check if both calibration and validation set contain all species. If a species occurs only once in your dataset, it is recommended to include it in the calibration dataset. For medium sized datasets (n=150) it is advised to have more calibration samples (for example 70% calibration, 30% validation). For very large datasets (n=500) it is also okay to have a balanced ratio (50/50), in particular if you are worried about model stability. Having more validation then calibration samples is uncommon, but it might be a good idea if you would like to check if your model is overfitted. 

 If you wish to use cross validation (recommended for small datasets) you will need to change the following code sections:
 - in step 5 you need to define a testset even if you wish to use cross-validation.
 (for example use 'random' 50/50) 
 - in step 7 replace all "$RMSEP" by "$RMSECV" 
 - in step 8 delete the line "training_set = trainingset," 
 - in step 8 change "validation = "testset"," to "validation = "LOO"," 

## Wave sequence and complexity 

 ```md
segments_number    <- 6
smallest_segment   <- 72
minimum_frequency  <- 350
maximum_frequency  <- 2500
maximum_rank       <- 12
```
The 'segment_ number’ is the number of spectral regions that are analysed in every run (model size). This variable has the highest impact on the resource demand. Increasing this value will exponentially increase the memory demand of the analysis. Keep in mind that running the script on multiple CPU cores will also multiply the memory demand. If the computer runs out of memory it will use a swap file on the hard disk which is orders of magnitude slower. If it runs out of swap space the process might crash. If a system has limited memory but a strong CPU, it is advised to reduce this value and run multiple models instead. If you wish to run very large models (9+) it is advised to run it on only a few or single CPU core. When in doubt, it is advised to carefully monitor the memory utilisation over the complete calibration process (memory demand will increase in the process). Reducing the wave sequence complexity by increasing the 'minimum_frequency' or decreasing the 'maximum_frequency' will in most cases not decrease resource demand but might in some cases affect model quality. For example, if there are noisy regions at the edge of the spectrum, manually exclude them could increase model quality. However it is not advised to remove spectral regions without indication. The ‘maximum_rank’ is the highest number of factors of the plsr regression that will be tested (model complexity). This variable also has high impact on the resource demand. Typically, prediction models for leaf traits have a rank of approximately 7-12. Increasing this value too much (20+) will significantly increase the time that is needed for the calculation. 

##  Parallel analysis 
  ```md
 cores=detectCores(logical = FALSE) 
 cl <- makeCluster(cores) 
 registerDoParallel(cl) 
 on.exit(stopCluster(cl)) 
 finalMatrix <- foreach(i=1:300, .combine=rbind) %dopar% { . } 
 ```
The core of this script is the repeated parallel analysis of different randomly drawn configurations of wave sequences. By manipulating these factors, you define how many of these configurations are tested. The more configurations you test, the higher is the likelihood that you will come close to an "optimal" configuration. The number of cores can be automatically detected "cores=detectCores(logical = FALSE)" or manually set "cores<-5". The number in the loop should be devisable by the number of cores that you use. For example if you want to test 300 configurations on 5 cores, put "foreach(i=1:300." and the script will run 60 iterations of 5 parallel runs. Before you start your main analysis it would be a good idea to make a test run with one iteration "cl <- makeCluster(5) .foreach(i=1:5." and check the time via System.time(). It is advised that you run on 50 percent of the cores or up to cores-1 (if your system has 6 cores, use 3-5).

## Model complexity
 ```md
 Optimized_Par <- optimizePLS(component = REFERENCE,
                              spectra = SPECTRA, 
                              max_comps = maximum_rank,
                              training_set = trainingset,
                              region_list = Region_List,
                              parallel = FALSE)
 ```
The 'max_comps' variable defines the highest possible model rank that is tested in the calibration process. It is okay to slightly overestimate this value, as the calibration process will most likely result in a model that has a lower rank then the 'max_comps'. However, if the resulting model rank is equal to 'max_comps' it might be a good idea to rerun the analyses with a higher 'max_comps' value. In this case it is advised to manually check all model ranks which should reveal a peak (meaning that increasing the model rank to a certain point should increase model quality but beyond that, model quality should decrease), otherwise the model might be overfitted. The rank of a model is typically between 5 and 15 depending on the trait. This variable has high impact on the resource demand. Decreasing this value will increase calculation speed, but might negatively affect model quality.

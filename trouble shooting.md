
<h2 align="center"> Trouble shooting </h2>

## General remarks

For trouble shooting it is best to reduce model complexity parameters to values that are slightly above the minimal meaningful values. This will speed up the diagnosis and will still allow to test for the full functionality of the script. 
For example:
```md
segments_number    <- 3
smallest_segment   <- 72
minimum_frequency  <- 350
maximum_frequency  <- 2500
maximum_rank       <- 3
```



## Reading Spectral Files error
 ```md
 Error in names(x) <- value : 'names' attribute [3] must be the same length as the vector [2]
 ```
This error might occur if the wrong separator is used. Use sep="," to clarify.

```md
Error in (function (file, header = FALSE, sep = "", quote = "\"'", dec = ".",  : 
  no lines available in input
```
This means that some of the files in the /PATH_SPECTRA_CALIBRATION/ cannot be read. This occurs for example if there are other than spectral files in the folder (for example R script or reference file). Manually check the folder.

```md
Warning message:
In (function (file, header = FALSE, sep = "", quote = "\"'", dec = ".",  :
  incomplete final line found by readTableHeader on /PATH/FILENAME
```
This means that some of the files in the /PATH_SPECTRA_CALIBRATION/ have incomplete spectral information. Manually check the given file.
```md
Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
  line 2151 did not have 2 elements
```
This means that a spectral file has an empty field (Please note that ‘empty field’ does not mean NA in this context). Manually check the file. If some of the spectral data is missing, you can specify it as ‘NA’. While this might cause problems when creating a testset (see below), the calibration appears to work fine.

## Spectral Files not recognised  
The script uses a proprietary data type called 'spectra.matrix'. In some rare cases, the spectral information is not recognized.
To fix this, try to run:
 ```md
SPECTRA <- as.spectra.matrix(SPECTRA)
 ```
## Creating a Test-Set error
```md
Error in svd(x, nu = 0, nv = k) : infinite or missing values in 'x'
```
This means that there is an invalid value for in one of the reference files. Try loading a smaller dataset to isolate the problem. Alternatively change the method to ‘random’

## Creating a Test-Set crash or stuck
This process sometimes crashes or gets stuck in an endless loop. For example if you try to create two even numbered dataset using an uneven number of samples or in similar circumstances. Try using method = "random", or change from "p = 0.5", to "p = 0.49" for troubleshooting.

## Calibration error
The calibration process runs as part as a 'dplyr' parallel pipeline, which makes tracing errors difficult. In case of a crash it is advised to run the code without the '%dopar% {}' to track the problem (set "i<-1"). 

 ```md
 Error "`$<-.data.frame`(`*tmp*`, "spectra", value = c(0, 0, 0, 0, 0, : replacement has x rows, data has y"
  ```
This error occurs when the number of calibration spectra does not match the number of reference values. This occurs for example if you load files with header while they don’t contain one.


```md
 Error in `$<-.data.frame`(`*tmp*`, "F", value = numeric(0)) : replacement has x rows, data has y
 (while x = 0 and y is a numerical value)
```
This error occurs when the wave sequence does not match the range of the spectral files. Check if you selected the correct minimum and maximum frequency. Also manually check the spectral files for their spectral range. Also check if your files contain header information and if this is considered while loading the file.

## Validation error
An error occurs when the wrong variable type is addressed in validation (i.e. when test set variables are used for cross validation and vice versa). Check the calibration type section above for correct variable names.

## Bad Calibration Results
Run the optimization without the "%dopar% {.}" to manually observe the calibration events. A typical symptom of this is, that the vast majority of the calibration events gives rank 1 or 2 as the best model rank, even if you selected a higher "max_rank". This means that the algorithm cannot find a good correlation between the spectra and the reference files. Doublecheck the reference data, are these physiologically meaningful data or is there a high amount of outliers? Is it possible that there was a line skip in the reference file? It might be a good idea to test a sub-dataset (for example only the first 20 files + reference) to isolate the problem. If this works, try to enlarge this dataset until you find the outliers. Finally keep in mind that not everything can be calibrated, there might just not be a signal of what you are searching for in the spectra. Therefore it might be a good idea to start your analysis with a trait that is known to be predictable, for example LDMC for FieldSpec or nitrogen for Bruker MPA - this way you can make sure that the process itself works.

## Bad Validation Results
Typically, the validation R2 is slightly worse than the calibration R2. A difference of 0.05 to 0.1 is normal, a 0.2 difference is concerning (for example R2_cal=0.85; R2_val=0.65). If there is a big discrepancy between these values, this is a sign that the model is overfitted. In this case you should have a look in the 'finalMatrix' if there is a model on position 2, 3, 4 etc. that has a better R2_Cal to R2_Val ratio, this might a better choice even if the RMSE is worse. Furthermore you should analyse the test-set. Are all species in calibration and validation set? Are the values evenly distributed? Consider increasing the number of validation samples. Consider using a different algorithm for sample selection.


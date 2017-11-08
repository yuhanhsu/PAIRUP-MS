##########################################################################################
### Config file for running QC_submission_script.r (for C8-pos data)
##########################################################################################

### INPUT SETTINGS

dir.src <- "../src" # path to src directory

file.raw <- "example_raw_data_C8-pos.txt" # raw metabolite data file
file.sampleRunWindow <- "example_sample_run_window.txt" # What file to read for defining windows

met.type <- "C8-pos" # profilng method name
met.column <- "Metabolite" # column containing metabolite names
IS.name <- "Known_C8-pos_IS" # internal standard name (used for IS normalization)
PP.label <- "HPP" # prefix of sample IDs of 1st set of control samples (e.g. pooled plasma, used for normalization)
biorec.label <- "BPP" # prefix of sample IDs of 2nd set of control samples (e.g. commercial, used for fine-tuning)


### OUTPUT SETTINGS

# output directory name
dir.out <- paste("QC_output_", met.type, sep="") 

# final output file name
file.final_output <- paste(dir.out, "/X_QCed_", met.type,".txt", sep="")


### PARAMETER SETTINGS

# specify the required parameters:
sp.scale <- 0.5 # df of spline will be minimun of this number multiplied by number of unique PP, and the number of unique PP itself. So make it less than 1 to be effective. The smaller sp.scale, the less overfitting of the spline.
do.extrap <- FALSE # set to false if you do not want to extrapolate beyond enpoint PP samples within a fitted spline
sp.PPnum <- 5 # number of PP samples to fit each spline to
sp.window <- 40 # number of samples in each spline
tol.spline <- 15 # parameter in smooth.spline. A tolerance for same-ness or uniqueness of the x values. The values are binned into bins of size tol and values which fall into the same bin are regarded as the same

plot.it <- TRUE # do you want to plot any spline fits
plot.allsplines <- FALSE # whether plot spline fits for all metabolites
n.splineplot <- 5 # number of plots to produce of spline fit, randomly sampled (if plot.allsplines==F)

cv.PP <- 0.1 # Threshold to filter out peaks after normalization based on the CV of their PP values
sd.point <- 4 # Threshold to filter out extreme outlier points for each metabolite
sd.window <- 3 # Threshold to filter out outlier windows for each metabolite

plot.final <- TRUE # do you want to plot final QC'ed emtabolite data


##########################################################################################
### Function for determining if a pair of m/z values could be a match
### input: x and y = m/z values 
###	   x_ion and y_ion = "POS" or "NEG" ionization mode
###	   mode = "noAdduct" (match x and y directly without checking for adduct ions) 
###	          or "adduct" (allow matches between all adduct m/z values of x and y)
###	   w = window to allow for matching to a value
### pre-defined variables: posAdducts and negAdducts = list of adduct masses to consider
### (called by: "MatchSignalsByImp.r" and "MatchSignalsByPredRt.r")
###
### author: Yu-Han Hsu
##########################################################################################

matchMZ <- function(x, y, x_ion, y_ion, mode, w) {
	"function for determining if two m/z values could correspond to the same compound (with same mass)"

	# DEFINE M/Z DIFFERENCES TO CHECK
	if (x_ion == "POS" & y_ion == "POS") { # POS-POS matching
		if (mode == "adduct") {
			delAdduct <- c(0,c(dist(posAdducts))) # 0 + all pairwise differences between adduct masses
		} else {
			delAdduct <- c(0)
		}
	} else if (x_ion == "NEG" & y_ion == "NEG") { # NEG-NEG matching
		if (mode == "adduct") {
			delAdduct <- c(0,c(dist(negAdducts)))
		} else {
			delAdduct <- c(0)
		}
	} else if (x_ion == "POS" & y_ion == "NEG") { # POS-NEG matching
		if (mode == "adduct") {
			delAdduct <- c(outer(posAdducts,negAdducts,"-")) # pairwise differences b/w POS and NEG adduct masses
		} else {
			delAdduct <- c(1.007276 * 2) # difference b/w M+H and M-H
		}
	} else if (x_ion == "NEG" & y_ion == "POS") { # NEG-POS matching
		if (mode == "adduct" ) {
			delAdduct <- c(outer(negAdducts,posAdducts,"-")) # pairwise differences b/w NEG and POS adduct masses
		} else {
			delAdduct <- c(-1.007276 * 2) # difference b/w M-H and M+H
		}
	}

	# CHECK IF M/Z VALUES ARE WITHIN DEFINED DIFFERENCES +- ERROR WINDOW
	delMZ <- x - y
	isMatch <- FALSE	
	if (x_ion == y_ion) { # POS-POS or NEG-NEG: sign of difference does not matter, check abs diff b/w values
		delMZ <- abs(delMZ)
	}

	for (d in delAdduct) {
		if (delMZ > (d - w) & delMZ < (d + w)) {
			isMatch <- TRUE
			break
		}
	}	
	
	return(isMatch)
}


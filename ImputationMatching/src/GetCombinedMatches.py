##########################################################################################
### Create "combined" matches by merging "noAdduct" and "adduct" matching results
### (always give priority to "noAdduct" entry if it exists)
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser

# read in input/output file paths
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
noAdductFile = config.get("Files","noAdductFile")
adductFile = config.get("Files","adductFile")
outFilePath = config.get("Files","outFilePath")

# read in and store all "noAdduct" matches
noAddDict = {} # Metabolite1 as key, entire entry/line as value
with open(noAdductFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		noAddDict[line.split('\t')[0]] = line

print len(noAddDict)

# read in "adduct" matches, replace with "noAdduct" match if it exists, output
with open(adductFile,'r') as inFile, open(outFilePath,'w') as outFile:
	outFile.write(next(inFile))
	for line in inFile:
		met1 = line.split('\t')[0]
		if noAddDict.get(met1) is not None:
			outFile.write(noAddDict[met1])
		else:
			outFile.write(line)

print 'SCRIPT COMPLETED'

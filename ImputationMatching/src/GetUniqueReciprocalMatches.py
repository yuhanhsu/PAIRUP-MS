##########################################################################################
### Get "UNIQUE" and "RECIPROCAL" Dataset1-Datset2 matches
### (1) UNIQUE: for Met1 -> Met2  matches that match to the same Met2,
###	pick the one with best Match_Corr (for Imp matches) or DeltaRT (for PredRt)
### (2) RECIPROCAL: require Met1 -> Met2 and Met2 -> Met1 matches to be the same
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
matchApp = config.get("Parameters","matchApp")
match12File = config.get("Files","match12File")
match21File = config.get("Files","match21File")
uniqFile = config.get("Files","uniqFile")
recipFile = config.get("Files","recipFile")

# store match21 info
match21_dict = {} # Met2 as key, Met1 as value
with open(match21File,'r') as inFile:
	header = next(inFile)

	for line in inFile:
		fields = line.split()
		met2 = fields[0]
		met1 = fields[5]

		match21_dict[met2] = met1


# store match12 info, output reciprocal matches
match12_dup2_dict = {} # Met2 as key, list of lines as value
with open(match12File,'r') as inFile, open(recipFile,'w') as outFile:
	header = next(inFile)
	outFile.write(header)

	for line in inFile:
		fields = line.split()
		met1 = fields[0]
		met2 = fields[5]

		match12_dup2_dict.setdefault(met2,[]).append(line)

		if match21_dict.get(met2) == met1:
			outFile.write(line)


# output UNIQUE Met1-Met2 matches
with open(uniqFile,'w') as outFile:
	outFile.write(header)
	for m in match12_dup2_dict:
		
		if matchApp == 'PredRt':
			bestRt = float('inf')
			bestLine = ''
			for line in match12_dup2_dict[m]:
				rt = abs(float(line.split()[12]))
				if rt < bestRt:
					bestRt = rt
					bestLine = line
			
			outFile.write(bestLine)

		else:
			bestCorr = -1
			bestLine = ''
			for line in match12_dup2_dict[m]:
				corr = float(line.split()[12])
				if corr > bestCorr:
					bestCorr = corr
					bestLine = line
		
			outFile.write(bestLine)

print 'SCRIPT COMPLETED'

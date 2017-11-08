##########################################################################################
### Convert Dataset1 signal lists to Dataset2 (shared or matched) signal lists
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser

# read in input/output file paths
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
sharedFile = config.get("Input Files","sharedFile")
matchedFile = config.get("Input Files","matchedFile")
listFiles = config.get("Input Files","listFiles").split(',')
outMetFile = config.get("Output Files","outMetFile")
outListFiles = config.get("Output Files","outListFiles").split(',')


# ----------------------------------------------------------------------------
# read in shared knowns
matchDict = {} # Dataset1 signal as key, Dataset2 signal as value
with open(sharedFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		if fields[0] != 'NA' and fields[1] != 'NA':
			matchDict[fields[0]] = fields[1]

nShared = len(matchDict)
print '# of shared knowns: %d' % nShared

# read in matched signals (optional)
if matchedFile != 'NA':
	with open(matchedFile,'r') as inFile:
		next(inFile)
		for line in inFile:
			fields = line.strip().split('\t')
			if (fields[0] != 'NA' and fields[1] != 'NA'
				and fields[0] not in matchDict): # only store if Met1 not a shared known
				matchDict[fields[0]] = fields[1]

print '# of matched signals: %d' % (len(matchDict) - nShared)
print 'total # of mapped (shared + matched) signals: %d' % len(matchDict)

# output list of all mapped (shared or matched) Dataset2 signals (optional)
if outMetFile != 'NA':
	with open(outMetFile,'w') as outFile:
		for x in matchDict.keys():
			outFile.write(matchDict[x] + '\n')

# read in each input file with Dataset1 signal list(s), convert to Dataset2 signals and output convereted list(s)
for i in range(len(listFiles)):
	with open(listFiles[i],'r') as inFile, open(outListFiles[i],'w') as outFile:

		for line in inFile:
			outList = [matchDict[x] for x in line.strip().split(';') 
					if x in matchDict]

			outFile.write(';'.join(outList)+'\n')
		
print 'SCRIPT COMPLETED'

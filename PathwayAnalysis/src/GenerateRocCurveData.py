##########################################################################################
### Generate data file (from an annotation matrix) for plotting ROC curve in R
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser

# read in input/out file paths
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
setFile = config.get("Input Files","setFile")
annFile = config.get("Input Files","annFile")
filterFile = config.get("Input Files","filterFile")
outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in list of met sets to include (if not NA)
filterList = set() # unordered list of met sets
if filterFile != 'NA':
	with open(filterFile,'r') as inFile:
		next(inFile)
		for line in inFile:
			fields = line.strip().split('\t')
			if float(fields[1]) < 0.05:
				filterList.add(fields[0])

# read in met sets -> known mets mapping
metList = set() # unordered list of all mets included in met sets
setMetDict = {} # SET_ID as key, unordered list of mets as value
with open(setFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		mets = set(fields[1].split(';'))
		
		if ((filterFile == 'NA') or
		    (fields[0] in filterList)):
			setMetDict[fields[0]] = mets
			metList.update(mets)

print '# of sets: %d' % len(setMetDict)
print '# of mets: %d' % len(metList)

# read in met set scores for mets in metList (from annotation matrix)
metSetScoreDict = {} # met as key, list of set score (SpearmanR) as value
setList = [] # list of sets in order they appear in matrix
setInds = [] # indices of sets in setHeaders
with open(annFile,'r') as inFile: 
	setHeaders = next(inFile).strip().split('\t')[1:]
	if filterFile == 'NA':
		setInds = range(len(setHeaders))
	else:
		setInds = [i for i in range(len(setHeaders))
				if setHeaders[i] in filterList]

	setList = [setHeaders[i] for i in setInds]

	for line in inFile:
		fields = line.strip().split('\t')
		if fields[0] in metList:
			scores = [fields[i+1] for i in setInds]
			metSetScoreDict[fields[0]] = scores

print '# of mets with set scores: %d' % len(metSetScoreDict)
print '# of sets with scores: %d' % len(metSetScoreDict[list(metList)[0]])

# output: score + label (0=met in set, 1=met not in set)
with open(outPath,'w') as outFile:
	outFile.write('Score\tLabel\n')

	for i in range(len(setList)):
		mset = setList[i]
		setMets = setMetDict[mset]

		for met in list(metList):
			score = metSetScoreDict[met][i]
			label = 1
			if met in setMets:
				label = 0

			outFile.write('%s\t%d\n' % (score,label))

print 'SCRIPT COMPLETED'

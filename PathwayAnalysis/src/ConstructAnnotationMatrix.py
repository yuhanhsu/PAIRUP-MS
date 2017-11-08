##########################################################################################
### Construct metabolite signal x metaboltie set annotation matrix
### 
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser, re, math
from scipy.stats import ranksums, spearmanr

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
setFile = config.get("Input Files","setFile")
matrixFile = config.get("Input Files","matrixFile")
mcFile = config.get("Input Files","mcFile")
outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in met set annotations
setDict = {} # SET_ID as key, unordered list of met names as value
setList = [] # list of SET_IDs
setMetList = set() # unordered list of all mets found in met sets
with open(setFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		mets = set(fields[1].split(';'))
		setDict[fields[0]] = mets
		setList.append(fields[0])	
		setMetList.update(mets)
		
print '# of met sets to use: %d' % len(setList)
print '# of mets in met sets: %d' % len(setMetList)

# read in list of MCs
mcList = []
with open(mcFile,'r') as inFile:
	for line in inFile:
		mcList.append(int(line.strip().split('MC')[1]) - 1) # -1 to get 0-based MC index

print '# of MCs to use: %d' % len(mcList)

# read in signal x MC score matrix
mcDict = {} # met name as key, list of mcScores as value
with open(matrixFile,'r') as mcFile:
	next(mcFile)
	for line in mcFile:
		fields = line.strip().split('\t')
		metName = fields[0]
		tempScores = fields[1:]
		mcScores = map(float,[tempScores[i] for i in mcList])
		mcDict[metName] = mcScores	

print '# of metabolites with abundance data: %d' % len(mcDict)
print '# of MC scores for each metabolite: %d' % len(mcDict[metName])


# ----------------------------------------------------------------------------
# function to calculate met set x MC matrix
def calcSetMcScores(leftOutMet=None):	
	setMcScores = {}  # SET_ID as key, list of MC z-scores as value
	
	for setID in setList:
		# known metabolites in set
		metsInSet = set(setDict[setID])

		# known metabolites not in set
		metsOthers = set(setMetList - metsInSet)

		# leave out leftOutMet
		if leftOutMet in metsInSet:
			metsInSet.remove(leftOutMet)	
		elif leftOutMet in metsOthers:
			metsOthers.remove(leftOutMet)

		# MC scores for mets in set vs. not in set
		scoresInSet = [] # nested lists [met][mcScore]
		scoresOthers = []
		[scoresInSet.append(mcDict[met]) for met in list(metsInSet)]
		[scoresOthers.append(mcDict[met]) for met in list(metsOthers)]

		# rank sums test for each MC
		mcZscores = []
		for mc in range(len(mcList)):
			z, p = ranksums([x[mc] for x in scoresInSet],
					[y[mc] for y in scoresOthers])

			mcZscores.append(z)

		setMcScores[setID] = mcZscores

	return setMcScores


# -----------------------------------------------------------------------------
# Calculate signal x met set annotation matrix

print 'Calculating Matrix...'

with open(outPath,'w') as outFile:

	outFile.write('Signal\t' + '\t'.join(setList) + '\n')
	i = 0
	
	# for each met with met set annotations,
	# do leave-one-out met set x MC matrix calculation,
	# followed by Spearman correlation between met and sets
	# for all other mets/signals, compute complete met set x MC matrix,
	# then Spearman corr

	setScores = calcSetMcScores() # complete met set x MC matrix
	
	for met in mcDict:
		outFile.write(met)
		
		if met in setMetList:
			
			# met set x MC matrix without using data for current met
			tempSetScores = calcSetMcScores(met)
			for setID in setList:
				r,p = spearmanr(mcDict[met],tempSetScores[setID])
				outFile.write('\t%.4g' % r)
		else:	
			for setID in setList:
				r,p = spearmanr(mcDict[met], setScores[setID])	
				outFile.write('\t%.4g' % r)
		
		outFile.write('\n')

		i += 1
		if i % 100 == 0:
			print i

print 'SCRIPT COMPLETED'


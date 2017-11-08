##########################################################################################
### Calculate enrichment of original met set metabolites in reconstituted met sets
### (enrichment rank-sum p-value used as "label confidence score")
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser
from scipy.stats import ranksums

# read in input/out file paths
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
setFile = config.get("Input Files","setFile")
annFile = config.get("Input Files","annFile")
outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in set -> known mets mapping
metList = set() # unordered list of all mets included in sets
setMetDict = {} # SET_ID as key, unordered list of mets as value
with open(setFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		mets = set(fields[1].split(';'))
		setMetDict[fields[0]] = mets
		metList.update(mets)

print '# of sets: %d' % len(setMetDict)
print '# of mets: %d' % len(metList)

# read in reconstituted met set scores for mets in metList
# (from signal x met set annotation matrix)
metSetScoreDict = {} # met as key, list of set scores (absolute SpearmanR) as value
setList = [] # list of sets in order they appear in matrix
with open(annFile,'r') as inFile:
	setList = next(inFile).strip().split('\t')[1:]
	for line in inFile:
		fields = line.strip().split('\t')
		if fields[0] in metList:
			metSetScoreDict[fields[0]] = map(abs,map(float,fields[1:]))

print '# of mets with set scores: %d' % len(metSetScoreDict)
print '# of sets with scores : %d' % len(metSetScoreDict[list(metList)[0]])


# ----------------------------------------------------------------------------
# perform rank-sum test to compare met set scores of known metabolites
# originally in the set vs. not in the set; output p-value
with open(outPath,'w') as outFile:
	outFile.write('SET_ID\tRankSumP\n')

	for i in range(len(setList)):
		mset = setList[i]
		setMets = setMetDict[mset] # mets in set
		notSetMets = metList - setMets # mets not in set

		goodScores = [metSetScoreDict[m][i] for m in list(setMets)]
		badScores = [metSetScoreDict[m][i] for m in list(notSetMets)]

		z,p = ranksums(goodScores,badScores)
		outFile.write('%s\t%.3g\n' % (mset,p))

print 'SCRIPT COMPLETED'

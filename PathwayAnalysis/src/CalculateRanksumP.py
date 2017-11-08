##########################################################################################
### Calcualte pathway enrichment for list(s) of metabolite signals using rank-sum test
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser
from scipy.stats import ranksums

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
nStart = int(config.get("Parameters","nStart"))
nEnd = int(config.get("Parameters","nEnd"))
nSignals = int(config.get("Parameters","nSignals"))
listPath = config.get("Input Files","listPath")
annotationPath = config.get("Input Files","annotationPath")
metPath = config.get("Input Files","metPath")
outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in list of all signals
metabList = [] # all signals to include in analysis
with open(metPath,'r') as inFile:
	for line in inFile:
		metabList.append(line.strip())

# read in ranked signal lists to analyze
goodMetabsDict = {} # assoc signals as key, list of iter# (line# in listPath file) as value
with open(listPath,'r') as inFile:
	i = 0
	for line in inFile:

		if i in range(nStart,nEnd):
			metabs = sorted(line.strip().split(';')[:nSignals])
			metabKey = ';'.join(metabs)
			goodMetabsDict.setdefault(metabKey,[]).append(i)

		i += 1
		if i == nEnd:
			break

print '# of metabolite lists: %d' % (nEnd-nStart)
print '# of unique lists: %d' % len(goodMetabsDict)
print '# of good signals: %d' % nSignals

# read in signal x metabolite set annotation matrix
metSetDict = {} # metID as key, list of set scores (abs rho) as value
setList = [] # list of sets in order they appear in matrix
with open(annotationPath,'r') as inFile:
	setList = next(inFile).strip().split('\t')[1:]

	for line in inFile:
		fields = line.strip().split('\t')
		metID = fields[0]

		if metID in set(metabList):
			absCors = map(abs,map(float,fields[1:]))
			metSetDict[metID] = absCors

print 'total # of metabolite sets: %d' % len(setList)
print 'total # of signals: %d' % len(metSetDict)


# ----------------------------------------------------------------------------
# calculate met set enrichment rank-sum p-values
print 'Calculating enrichment for metabolite sets:'
iterSetPvalues = {} # iter# as key, list of set p-values as value
for entry in goodMetabsDict:
	goodMetabs = entry.split(';')
	iters = goodMetabsDict[entry]

	badMetabs = list(metabList) # this should treat duplicates seperately
	for met in goodMetabs:
		 badMetabs.remove(met)

	for i in range(len(setList)):

		# rank-sum test
		goodScores = [metSetDict[met][i] for met in goodMetabs]
		badScores = [metSetDict[met][i] for met in badMetabs]
		z, p = ranksums(goodScores,badScores)

		for j in iters:
			iterSetPvalues.setdefault(j,[]).append('%.3g' % p)

	print 'entry done'

# output p-values
with open(outPath, 'w') as outFile:
	for i in range(nStart,nEnd):
		outFile.write('\t'.join(iterSetPvalues[i]) + '\n')

print 'SCRIPT COMPLETED'

##########################################################################################
### Identify enriched metabolite sets for each MC
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser
from scipy.stats import ranksums
from random import sample

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
nIters = int(config.get("Parameters","nIters"))
setFile = config.get("Input Files","setFile")
matrixFile = config.get("Input Files","matrixFile")
outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in met set -> mets mapping
metList = set() # unordered list of known mets with met set annotations
setList = [] # list of met sets (SET_ID)
setMetDict = {} # SET_ID as key, list of mets as value
with open(setFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		mets = set(fields[1].split(';'))
		setMetDict[fields[0]] = mets
		setList.append(fields[0])
		metList.update(mets)

print '# of met sets: %d' % len(setList)
print '# of unique metabolites: %d' % len(metList)

# read in met signal x MC score matrix
metMcDict = {} # met as key, list of MC scores as value
with open(matrixFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		if fields[0] in metList:
			metMcDict[fields[0]] = map(float,fields[1:])

print '# of mets with MC scores: %d ' % len(metMcDict)
print '# of MCs: %d' % len(metMcDict[list(metList)[0]])


# ----------------------------------------------------------------------------
# calculate MC-met set enrichment rank-sum p-values
mcPvalues = [] # list of list of p-values: mcPvalues[mc][setID]
obsPvalues = [] # list of all observed p-values
nullPvalues = [] # list of all null p-values

for mc in range(len(metMcDict[list(metList)[0]])):
	tempPvalues = []
	
	# rank-sum test to compare MC scores of mets in vs. not in met set
	for ms in setList:
		setMets = setMetDict[ms] # known mets in met set
		notSetMets = metList - setMets # known mets not in met set

		goodScores = [metMcDict[met][mc] for met in list(setMets)]
		badScores = [metMcDict[met][mc] for met in list(notSetMets)]

		z,p = ranksums(goodScores,badScores)
		tempPvalues.append(p)	

		# null iterations: compare set of random mets vs. others
		for i in range(nIters):
			nullMets = set(sample(metList,len(setMets)))
			notNullMets = metList - nullMets

			goodScores = [metMcDict[met][mc] for met in list(nullMets)]
			badScores = [metMcDict[met][mc] for met in list(notNullMets)]
			
			z,p = ranksums(goodScores,badScores)
			nullPvalues.append(p)

	mcPvalues.append(tempPvalues)
	obsPvalues.extend(tempPvalues)

	# log progress
	print mc+1

# calcualte FDR
thresholds = sorted([i*10**-j for i in range(1,10) for j in range(0,11)])[:-22]
thresholds.append(1) # don't care about entries above nominal p-value sig.
fdrDict = {} # p-value threshold as key, FDR as value
for t in thresholds:
	obsCount = sum(x <= t for x in obsPvalues)
	nullCount = sum(x <= t for x in nullPvalues)/float(nIters)

	if obsCount == 0:
		fdrDict[t] = 1
	else:
		fdrDict[t] = min(nullCount/obsCount,1)

print 'FDR calculation done'


# ----------------------------------------------------------------------------
# for each MC, output: (1) best met set(s) and (2) all sets with FDR < 0.05
with open(outPath,'w') as outFile:
	outFile.write('MC\tBestP\tBestFDR\tBestSets\tFDR0.05_Sets\n')
	for mc in range(len(mcPvalues)):
		mcFdrs = []
		for i in range(len(setList)):
			for t in thresholds:
				if mcPvalues[mc][i] <= t:
					mcFdrs.append(fdrDict[t])
					break

		bestP = min(mcPvalues[mc])
		bestInds = [i for i in range(len(setList)) 
				if mcPvalues[mc][i] == bestP]
		bestFdr = mcFdrs[bestInds[0]]

		bestP = mcPvalues[mc][bestInds[0]]
		bestSets = '|'.join([setList[i] for i in bestInds])

		sigInds = [i for i in range(len(setList))
				if mcFdrs[i] < 0.05]
		sigSets = ''
		if len(sigInds) > 0: 
			sigSets = '|'.join([setList[i] for i in sigInds])

		outFile.write('MC%d\t%.3g\t%.3g\t%s\t%s\n' 
			% (mc+1,bestP,bestFdr,bestSets,sigSets))

print 'SCRIPT COMPLETED'

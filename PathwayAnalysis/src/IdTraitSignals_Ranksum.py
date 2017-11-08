##########################################################################################
### Identify signals associated with a dichotomous/binary trait using rank-sum test
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser
from scipy.stats import ranksums

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
metFile = config.get("Input Files","metFile")
traitFile = config.get("Input Files","traitFile")
trait = config.get("Parameters","trait")
pThreshold = float(config.get("Parameters","pThreshold"))
outTableFile = config.get("Output Files","outTableFile")
outListFile = config.get("Output Files","outListFile")


# ----------------------------------------------------------------------------
# read in metabolite data
sampleOrder = [] # order of samples (rows)
metOrder = [] # order of signals (columns)
metDict = {} # signal as key, list of sample values as value (order saved in sampleOrder)

with open(metFile,'r') as mFile:
	header = next(mFile).strip().split('\t')
	metOrder = header[1:len(header)]
	
	for line in mFile:
		fields = line.strip().split('\t')
		sampleOrder.append(fields[0])
		for i in range(len(metOrder)):
			metDict.setdefault(metOrder[i],[]).append(float(fields[i+1]))

print '# of samples: %d' % len(sampleOrder)
print '# of signals: %d' % len(metOrder)

# read in trait data
traitSamplesDict = {} # trait value as key, unordered list of sample IDs as value 
with open(traitFile,'r') as pFile:
	
	# get column index of trait
	header = next(pFile).strip().split('\t')
	pInd = header.index(trait)

	# store trait
	for line in pFile:
		fields = line.strip().split('\t')
		pValue = fields[pInd]
		if pValue != 'NA': # do not store sample if trait value is missing
			traitSamplesDict.setdefault(pValue,set()).add(fields[0])
		
# store sample indices for each trait group in sampleOrder
traitLabels = traitSamplesDict.keys()
xInds = []
yInds = []
for i in range(len(sampleOrder)):
	if sampleOrder[i] in traitSamplesDict[traitLabels[0]]:
		xInds.append(i)
	elif sampleOrder[i] in traitSamplesDict[traitLabels[1]]:
		yInds.append(i)

print '# of %s samples: %d' % (traitLabels[0],len(xInds))
print '# of %s samples: %d' % (traitLabels[1],len(yInds))


# ----------------------------------------------------------------------------
# perform rank-sum test for each signal
results = [] # nested list of [signal, z, p]
for met in metOrder:
	metValues = metDict[met]
	xValues = [metValues[i] for i in xInds]
	yValues = [metValues[i] for i in yInds]

	z, p = ranksums(xValues, yValues)

	results.append([met, z, p])

# sort results by p-value
results.sort(key=lambda x: x[2])

# output: (1) rank-sum stats for all signals and (2) significant signals
with open(outTableFile,'w') as outFile1, open(outListFile,'w') as outFile2:
	outFile1.write('Signal\tRankSumZ\tP\n')

	sigList = []
	for result in results:
		outFile1.write('%s\t%.4g\t%.4g\n' % (result[0],result[1],result[2]))

		if result[2] < pThreshold:
			sigList.append(result[0])
	
	outFile2.write(';'.join(sigList) + '\n')

print 'SCRIPT COMPLETED'

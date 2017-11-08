##########################################################################################
### Generate null lists of signals ranked by rank-sum p-values with permuted binary traits
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser
from random import shuffle
from scipy.stats import ranksums

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
metFile = config.get("Input Files","metFile")
traitFile = config.get("Input Files","traitFile")
trait = config.get("Parameters","trait")
N = int(config.get("Parameters","N"))
outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in metabolite data
sampleOrder = [] # order of samples (rows)
metOrder = [] # order of signals (columns)
metDict = {} # signal as key, list of sample values as value (order saved in sampleOrder)

with open(metFile,'r') as mFile:
	header = next(mFile).strip().split('\t')
	metOrder = header[1:]

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
		
# store number of samples with each trait label
traitLabels = traitSamplesDict.keys()
xNum = len(traitSamplesDict[traitLabels[0]])
yNum = len(traitSamplesDict[traitLabels[1]])

print '# of %s samples: %d' % (traitLabels[0],xNum)
print '# of %s samples: %d' % (traitLabels[1],yNum)

# store sample indices for samples with trait data available
sampInds = []
for i in range(len(sampleOrder)):
	if (sampleOrder[i] in traitSamplesDict[traitLabels[0]] 
		or sampleOrder[i] in traitSamplesDict[traitLabels[1]]):
		sampInds.append(i)

print 'total # of samples: %d' % len(sampInds)


# ----------------------------------------------------------------------------
# generate N lists of signals ranked by association with null traits
with open(outPath,'w') as outFile:
	for i in range(N):
		results = [] # nested list of [signal, z, p]
		
		# generate "null" sample labels
		shuffle(sampInds)
		xInds = sampInds[:xNum]
		yInds = sampInds[xNum:]

		for met in metOrder:
			metValues = metDict[met]
			xValues = [metValues[i] for i in xInds]
			yValues = [metValues[i] for i in yInds]

			z, p = ranksums(xValues, yValues)
			results.append([met, z, p])

		results.sort(key=lambda x: x[2])

		outFile.write(';'.join([x[0] for x in results]) + '\n')

		if (i % 100) == 0:
			print i

print 'SCRIPT COMPLETED'

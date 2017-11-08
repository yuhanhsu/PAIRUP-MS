##########################################################################################
### Generate null lists of signals ranked by Pearson p-values with permuted trait values
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser
from random import shuffle
from scipy.stats import pearsonr

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
metFile = config.get("Input Files","metFile")
traitFile = config.get("Input Files","traitFile")
trait = config.get("Parameters","trait")
N = int(config.get("Parameters","N"))
outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in trait data
traitDict = {} # sampleID as key, trait value as value
with open(traitFile,'r') as pFile:

	# get column index of trait
	header = next(pFile).strip().split('\t')
	pInd = header.index(trait)

	# store trait
	for line in pFile:
		fields = line.strip().split('\t')
		pValue = fields[pInd]
		if pValue != 'NA': # do not store sample if trait value is missing
			traitDict[fields[0]] = float(pValue)

# read in metabolite data
sampleOrder = [] # order of samples (rows)
metOrder = [] # order of signals (columns)
metDict = {} # signal as key, list of sample values as value (order saved in sampleOrder)

with open(metFile,'r') as mFile:
	header = next(mFile).strip().split('\t')
	metOrder = header[1:]

	for line in mFile:
		fields = line.strip().split('\t')
		if fields[0] in traitDict: # only want samples with trait data
			sampleOrder.append(fields[0])
			for i in range(len(metOrder)):
				metDict.setdefault(metOrder[i],[]).append(float(fields[i+1]))

print '# of samples: %d' % len(sampleOrder)
print '# of signals: %d' % len(metOrder)


traitValues = [] # vector of trait values for samples stored in metDict
for s in sampleOrder:
	traitValues.append(traitDict[s])

print '# of traitValues: %d' % len(traitValues)


# ----------------------------------------------------------------------------
# generate N lists of signals ranked by association with null traits
with open(outPath,'w') as outFile:
	for i in range(N):
		results = [] # nested list of [signal, r, p]
		
		# generate "null" trait
		shuffle(traitValues)

		for met in metOrder:
			metValues = metDict[met]
			r, p = pearsonr(metValues, traitValues)
			results.append([met,r,p])

		results.sort(key=lambda x: x[2])

		outFile.write(';'.join([x[0] for x in results]) + '\n')

		if (i % 100) == 0:
			print i

print 'SCRIPT COMPLETED'


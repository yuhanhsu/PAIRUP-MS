##########################################################################################
### Identify signals associated with a continuous trait using Pearson correlation test 
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser
from scipy.stats import pearsonr

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
	metOrder = header[1:len(header)]
	
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

print '# of stored trait values: %d' % len(traitValues) 


# ----------------------------------------------------------------------------
# calculate correlation between signals and trait
results = [] # nested list of [signal, r, p]
for met in metOrder:
	metValues = metDict[met]
	r, p = pearsonr(metValues, traitValues)

	results.append([met, r, p])

# sort results by p-value
results.sort(key=lambda x: x[2])

# output: (1) correlation stats for all signals and (2) significant signals
with open(outTableFile,'w') as outFile1, open(outListFile,'w') as outFile2:
	outFile1.write('Signal\tCorrelation\tP\n')

	sigList = []
	for result in results:
		outFile1.write('%s\t%.4g\t%.4g\n' % (result[0],result[1],result[2]))

		if result[2] < pThreshold:
			sigList.append(result[0])

	outFile2.write(';'.join(sigList) + '\n')

print 'SCRIPT COMPLETED'

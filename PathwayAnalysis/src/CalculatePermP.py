##########################################################################################
### Calculate permutation p-values for list(s) of observed rank-sum (or Fisher) p-values
### by comparing against null rank-sum (or Fisher) p-values
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])

obsPath = config.get("Input Files","obsPath")
obsStart = int(config.get("Parameters","obsStart"))
obsEnd = int(config.get("Parameters","obsEnd"))

nullPath = config.get("Input Files","nullPath")
nullStart = int(config.get("Parameters","nullStart"))
nullEnd = int(config.get("Parameters","nullEnd"))

outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in observed p-values
obsPvalues = [] # list of lists of observed met set enrichment p-values
with open(obsPath,'r') as inFile:
	i = 0
	for line in inFile:
		if i in range(obsStart,obsEnd):
			obsPvalues.append(map(float,line.strip().split('\t')))
		
		i += 1
		if i == obsEnd:
			break

nObs = len(obsPvalues)
print '# of observed p-value lists to analyze: %d' % nObs

# read in null p-values
nullPvalues = []
with open(nullPath,'r') as inFile:
	i = 0
	for line in inFile:
		if i in range(nullStart,nullEnd):
			nullPvalues.append(map(float,line.strip().split('\t')))

		i += 1
		if i == nullEnd:
			break

nNull = len(nullPvalues)
print '# of null p-value lists to use: %d' % nNull

# number of metabolite sets (number of p-values in each list)
nSets = len(obsPvalues[0])
print '# of metabolite sets: %d' % nSets


# ----------------------------------------------------------------------------
# calculate and output permutation p-values
print 'Calculating permutation p-values...'
with open(outPath,'w') as outFile:
	for i in range(nObs):
		permPs = []
		for j in range(nSets):
			permPs.append(
				len([x for x in nullPvalues 
				if x[j] <= obsPvalues[i][j]])/float(nNull))

		outFile.write('\t'.join(['%.3g' % x for x in permPs]) + '\n')

		print 'set %d' % i

print 'SCRIPT COMPLETED'


##########################################################################################
### Calculate pathway over-representation for list(s) of known metabolites
### using Fisher's exact test
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser
from scipy.stats import fisher_exact

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
nStart = int(config.get("Parameters","nStart"))
nEnd = int(config.get("Parameters","nEnd"))
nMets = int(config.get("Parameters","nMets"))
listPath = config.get("Input Files","listPath")
annotationPath = config.get("Input Files","annotationPath")
metPath = config.get("Input Files","metPath") 
outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in list of known metabolites to include in analysis
metabList = []
with open(metPath,'r') as inFile:
	for line in inFile:
		metabList.append(line.strip())

print '# of metabolites to include for analysis: %d' % len(metabList)

# read in list of  metabolites to analyze
goodMetabsDict = {} # assoc metabs as key, list of iter# (line# in listPath file) as value 
with open(listPath,'r') as inFile:
	i = 0
	for line in inFile:
		if i in range(nStart,nEnd):
			metabs = sorted(line.strip().split(';')[:nMets])
			metabKey = ';'.join(metabs)
			goodMetabsDict.setdefault(metabKey,[]).append(i)

		i += 1
		if i == nEnd:
			break

print '# of metabolite lists: %d' % (nEnd-nStart)
print '# of unique lists: %d' % len(goodMetabsDict)
print '# of good metabs: %d' % nMets

# read in binary metabolite set annotaions
setMetDict = {} # met set ID as key, SET of met names as value (restrict to metabolites included for analysis)
setList = [] # list of sets in order they appear in file (should be in numeric order, same as in reconstituted matrix) 
with open(annotationPath,'r') as inFile:
	next(inFile)	
	
	for line in inFile:
		fields = line.strip().split('\t')
		setID = fields[0]
		metNames = fields[1].split(';')

		setMetDict[setID] = set(metNames) & set(metabList)
		setList.append(setID)

print 'total # of metabolite sets: %d' % len(setList)


# ----------------------------------------------------------------------------
# calculate pathway over-representation Fisher's exact p-values
print 'Calculating over-representation for metabolite sets:'
iterSetPvalues = {} # iter# as key, list of set p-values as value
for entry in goodMetabsDict:
	goodMetabs = entry.split(';')
	iters = goodMetabsDict[entry]
	
	badMetabs = list(metabList) # this should treat duplicates seperately
	for met in goodMetabs:
		badMetabs.remove(met)

	for i in range(len(setList)):
		setName = setList[i]
		setMets = setMetDict[setName]

		# 2 x 2 contigency table (Trait-associated yes/no x Annotated to be in set yes/no)
		# count duplicates separately
		yy = len([met for met in goodMetabs if met in setMets]) # trait-associated, in set
		yn = len([met for met in goodMetabs if met not in setMets]) # trait-associated, not in set
		ny = len([met for met in badMetabs if met in setMets]) # not trait-asscoiated, in set
		nn = len([met for met in badMetabs if met not in setMets]) # not trait-associated, not in set

		stat, p = fisher_exact([[yy,yn],[ny,nn]],alternative='greater')
	
		for j in iters:
			iterSetPvalues.setdefault(j,[]).append('%.3g' % p)

	print 'entry done'

# output p-values
with open(outPath,'w') as outFile:
	for i in range(nStart,nEnd):
		outFile.write('\t'.join(iterSetPvalues[i]) + '\n')

print 'SCRIPT COMPLETED'

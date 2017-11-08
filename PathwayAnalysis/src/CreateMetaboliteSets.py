##########################################################################################
### Generate metabolite set annotations to be used for pathway reconstitution 
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser, re

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
nMets = int(config.get("Parameters","nMets"))
pathFile = config.get("Input Files","pathFile")
metFile = config.get("Input Files","metFile")
idFile = config.get("Input Files","idFile")
outFilePath = config.get("Output Files","outFilePath")


# -----------------------------------------------------------------------------
# save CPDB pathway -> metabolite IDs mapping
pathDict = {} # source:path as key, [[kegg IDs],[pubchemIDs],[chebiIDs]] as value
keggStr = re.compile('kegg:(C[0-9]+)') # currently ignoring GXXXXX IDs (glycans)
pubchemStr = re.compile('pubchem-compound:([0-9]+)')
chebiStr = re.compile('chebi:([0-9]+)')

with open(pathFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		pathName = fields[1] + ':' + fields[0]
		metIDs = fields[2].split(',')

		keggIDs = []
		pubchemIDs = []
		chebiIDs = []
		for id in metIDs:
			if keggStr.match(id) is not None:
				keggIDs.append(keggStr.match(id).group(1))
			elif pubchemStr.match(id) is not None:
				pubchemIDs.append(pubchemStr.match(id).group(1))
			elif chebiStr.match(id) is not None:
				chebiIDs.append(chebiStr.match(id).group(1))

		if pathDict.get(pathName) is None:
			pathDict[pathName] = [keggIDs,pubchemIDs,chebiIDs]
		else:
			pathDict[pathName][0].extend(keggIDs)
			pathDict[pathName][1].extend(pubchemIDs)
			pathDict[pathName][2].extend(chebiIDs)

print '# of CPDB pathways: %d' % len(pathDict)

# -----------------------------------------------------------------------------
# save list of known metabolites found in metabolite data
metList = []

with open(metFile,'r') as inFile:
	for line in inFile:
		metList.append(line.strip())

print '# of known metabolites in data: %d' % len(metList)

# -----------------------------------------------------------------------------
# generate CPDB pathway -> known met names mapping
goodPathDict = {} # pathway as key, found met names as value

with open(idFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		metName = fields[0]
		kegg = fields[2]
		pubchem = fields[3]
		chebi = fields[4]

		for pathway in pathDict:
			if metName in metList and (kegg in pathDict[pathway][0]
				or pubchem in pathDict[pathway][1]
				or chebi in pathDict[pathway][2]):

				goodPathDict.setdefault(pathway,[]).append(metName)

print '# of CPDB pathways with any known mets: %d' % len(goodPathDict)

# -----------------------------------------------------------------------------
# get unique met sets (w/ corresponding pathway names) with >= nMets metabolites
metsetDict = {} # metabolite set (met names separated by ';') as key, pathway names as value
tempPathCount = 0 # total number of pathways in met sets
goodMetList = [] # list of metabolites included in met sets

for pathway in goodPathDict:
	if len(goodPathDict[pathway]) >= nMets:
		metStr = ';'.join(goodPathDict[pathway])
		metsetDict.setdefault(metStr,[]).append(pathway)
	
		goodMetList.extend(goodPathDict[pathway])

		tempPathCount += 1

print '# of CPDB pathways with >= %d mets: %d' % (nMets,tempPathCount)
print '# of unique met sets with >= %d mets: %d' % (nMets,len(metsetDict))
print '# of unique metabolites in met sets: %d' % len(set(goodMetList))

# -----------------------------------------------------------------------------
# assign SET_IDs to met sets, output ID, mets, pathways
setCount = 1
with open(outFilePath,'w') as outFile:
	outFile.write('SET_ID\tMetabolites\tPathways\n')

	metSets = sorted(metsetDict.keys())
	for ms in metSets:
		setID = 'SET_%d' % setCount
		pathways = '|'.join(metsetDict[ms])

		outList = [setID,ms,pathways]
		outFile.write('\t'.join(outList) + '\n')

		setCount += 1

print 'SCRIPT COMPLETED'

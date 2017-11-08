##########################################################################################
### Calculate FDR for a range of enrichment permutation p-value thresholds
### (by comparing observed vs. null permutation p-values)
### Output final met set enrichment results
###
### author: Yu-Han Hsu
##########################################################################################

import sys, ConfigParser

# read in input/output file paths and parameters
config = ConfigParser.ConfigParser()
config.read(sys.argv[1])

# command-line arguments
obsPath = config.get("Input Files","obsPath")

nullPath = config.get("Input Files","nullPath")
nullStart = int(config.get("Parameters","nullStart"))
nullEnd = int(config.get("Parameters","nullEnd"))

permN = int(config.get("Parameters","permN"))

nomPath = config.get("Input Files","nomPath")

setDictPath = config.get("Input Files","setDictPath")

setValidPath = config.get("Input Files","setValidPath")

outPath = config.get("Output Files","outPath")


# ----------------------------------------------------------------------------
# read in metabolite sets, store set -> CPDB pathways mapping
resultsList = [] # list of lists [SET_ID, nominal p, perm p] as value
setPathwaysDict = {} # dict of SET_ID -> CPDB pathways
with open(setDictPath,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		resultsList.append([fields[0]])
		setPathwaysDict[fields[0]] = fields[2]

# read in metabolite set post-reconstitution label confidence scores
setValidScoresDict = {} # dict of SET_ID -> rank-sum p
with open(setValidPath,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		setValidScoresDict[fields[0]] = fields[1]

# read in nominal p-values
with open(nomPath,'r') as inFile:
	nomPvalues = map(float,next(inFile).strip().split('\t'))

[resultsList[i].append(nomPvalues[i]) for i in range(len(resultsList))]

# read in observed perm p-values
obsPvalues = [] # list of observed perm p-values
with open(obsPath,'r') as inFile:
	obsPvalues = map(float,next(inFile).strip().split('\t'))

[resultsList[i].append(obsPvalues[i]) for i in range(len(resultsList))]

N = len(resultsList)
print '# of p-values/metabolite sets in one analysis: %d' % N

# read in null perm p-values
nullPvalues = [] # list of null perm p-values
with open(nullPath,'r') as inFile:
	i = 0
	for line in inFile:
		if i in range(nullStart,nullEnd):	
			[nullPvalues.append(x) for x in
				map(float,line.strip().split('\t'))]
	
		i += 1
		if i == nullEnd:
			break

M = len(nullPvalues)/N # number of iterations (of null analysis)
print '# of null iterations: %d' % M

# range of p-value thresholds for FDR calculation
thresholds = [float(x)/permN for x in range(permN+1)]
print 'range of p-value thresholds: %d - %d' % (thresholds[0], thresholds[-1])

# calculate FDR for each threshold
fdrDict = {} # p-value threshold as key, FDR as value
for t in thresholds:
	obsCount = sum(x <= t for x in obsPvalues) # = rank?
	nullCount = sum(x <= t for x in nullPvalues)/float(M)

	if obsCount == 0:
		fdrDict[t] = 'NA' #1
	else:
		fdrDict[t]  = nullCount/obsCount
		if fdrDict[t] > 1: # cap FDR at 1
			fdrDict[t] = 1

# adjust FDR so it is monotonically increasing as P increases
minFdr = 1
for t in reversed(thresholds):
	if fdrDict[t] > minFdr:
		fdrDict[t] = minFdr
	else:
		minFdr = fdrDict[t]

print 'FDR calculation done'


# ----------------------------------------------------------------------------
# output metabolite set enrichment results with FDR
# sort sets by enrichment significance
resultsList.sort(key=lambda x: (x[2],x[1]))

with open(outPath,'w') as outFile:
	outFile.write('SET_ID\tNominal_P\tPerm_P\tFDR\tLabel_P\tCPDB_Pathways\n')

	for r in resultsList:
		setid = r[0]
		nomp = '%.3g' % r[1]
		permp = '%.3g' % r[2]
		fdr = '%.3g' % fdrDict[r[2]]
		label = setValidScoresDict[setid]
		pathways = setPathwaysDict[setid]

		outFile.write(setid + '\t' + nomp + '\t' + permp + '\t'
			+ fdr + '\t' + label + '\t' + pathways + '\n')

print 'SCRIPT COMPLETED'

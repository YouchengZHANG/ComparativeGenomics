import sys, re

# acquire first needed data

geneOrderList = []

aHandle = open (sys.argv [1])

lines = aHandle.readlines ()

for aLine in lines:

	aLine = aLine.replace ("\n", "")

	if aLine.startswith (">"):

		#print aLine [0:len (aLine)],
		geneOrderList.append (aLine [0:len (aLine)])
#print geneOrderList

# acquire second needed data

partOfCluster = {}

bHandle = open (sys.argv [2])

lines = bHandle.readlines ()

id = 0

for aLine in lines:

	aLine = aLine.replace ("\n", "")
	#print aLine
	#words = aLine.rstrip('\r').split ("\t")
	words = aLine.split (" ")
	#print words
	for aWord in words:

		if not partOfCluster.has_key (aWord):

			partOfCluster [aWord] = id

	id = id + 1

#print partOfCluster
# put together

for aGene in geneOrderList:

	if partOfCluster.has_key (aGene):

		print partOfCluster [aGene],


# put together and handle the strandness
#for aGene in geneOrderList:
#
#	if partOfCluster.has_key(aGene):
#
#		if aGene.endswith(‘_rev’):
#
#			print(‘-’+ str(partOfCluster [aGene]),)
#
#		else:
#			print(‘+’+ str(partOfCluster [aGene]),)



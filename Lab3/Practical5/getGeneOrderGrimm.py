import sys

# this script reads in all the gene order files
# limits the analysis to those clusters that are present in
# all the supplied genomes
# and renumbers the result

# first, it is necessary to get the list of what is present

maxGenes = int (sys.argv [1])

noFiles = len (sys.argv) - 2

presentGenes = {}
geneLists = []

# read the input files

for i in range (2, noFiles + 2):
	
	# open this file

	fileHandle = open (sys.argv [i])
	lines = fileHandle.readlines ()
	
	# assume the input is a single line
	# of space-separated numbers

	genes = lines [0].split (" ")

	geneLists.append (lines [0]) # save this list

	# keep track on how many times each gene appears in this particular genome

	countsInGenome = {}

	for gene in genes:

		if countsInGenome.has_key (gene):
		
			countsInGenome [gene] = countsInGenome [gene] + 1

		else:

			countsInGenome [gene] = 1
	
	# with that done, for each gene that occurred exactly once
	# see if previously checked genomes included it

	for gene in countsInGenome.keys ():

		if countsInGenome [gene] == 1 and i == 2: # first file must be handled differently

			presentGenes [gene] = 1

		elif countsInGenome [gene] == 1 and presentGenes.has_key (gene):

			presentGenes [gene] = presentGenes [gene] + 1

# remove any genes from the analysis that are not present exactly once in each genome

id = 0

for gene in presentGenes.keys ():

	if not (presentGenes [gene] == noFiles and id < maxGenes):

		del presentGenes [gene]

	else:

		id = id + 1

# assign new indexing from 1 to end

id = 1

newMap = {}

for gene in presentGenes.keys ():

	newMap [gene] = id

	id = id + 1

# print out the revised gene lists

id = 1

for geneList in geneLists:

	print "Genome ", id, ":"

	genes = geneList.split (" ")

	for gene in genes:

		if presentGenes.has_key (gene):

			print newMap [gene],

	id = id + 1

	print
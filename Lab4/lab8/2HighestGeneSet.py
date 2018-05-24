# This script is to get highest counts of overlap of gene sets
import sys

#First one is that my blast result
f1 = open(sys.argv[1])
mygene = f1.readlines()
gene_list =[]
for i in mygene:
	i = i.rstrip("\n")
	gene_list.append(i)

#print (len(gene_list))
#f2 is experimental gene set
f2 = open(sys.argv[2])
gene_set = f2.readlines()

count = []
#Each line in experiment gene set
for group in gene_set:
	line_num = 1
	gene_count = 0
	geneset = set(group.split())
	#print (geneset)
	for gene in geneset:
		#print (gene)
		if gene in gene_list:
			gene_count +=1
	count.append(gene_count)
print (count)




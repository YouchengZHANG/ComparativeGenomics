# This script is to get gene symbols
import sys


f = open (sys.argv[1])
results = f.readlines()

for lines in results:
	#print (lines)
	symbol = lines.split()
	symbol = symbol[1].split('|')
	symbol = symbol[2].split('_')
	print (symbol[0])

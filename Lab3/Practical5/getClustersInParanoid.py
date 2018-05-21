#!/usr/bin/python
# This script can be modified to extract ortholog groups
# from multiple InParanoid runs of one proteome (i.e. 10) vs other proteomes (i.e. 15, 17, 19)
from collections import OrderedDict
import re


def readSQLTable(path):
    groupmap = OrderedDict()
    with open(path , 'r') as pfa15:
        for line in pfa15:
            group, bitscore, species, inparalog, seqname = line.rstrip().split()
            if group not in groupmap:
                groupmap[group] = set([seqname])
            else:
                groupmap[group].add(seqname)
    return(groupmap)


# Sqltable files from InParanoid run
pfa15 = readSQLTable('./sqltable.10.fa.txt.pfa-15.fa.txt.pfa')
pfa17 = readSQLTable('./sqltable.10.fa.txt.pfa-17.fa.txt.pfa')
pfa19 = readSQLTable('./sqltable.10.fa.txt.pfa-19.fa.txt.pfa')

output = []
for key, value in pfa15.iteritems():
    # Initialize with current group
    groups = value
    # Search for other clusters with reference in other sql
    for item in value:
        # Regex matching 10.fa.txt lines (reference sequence)
        # This needs to be adjusted to your reference proteome
        m = re.match('.*10.fa.txt_(orf[0-9]+)', item)
        if m:
            # We have ref, find it in other files
            ref = m.group(1)
            # Use that reference to find other groups
            for sqlmap in [pfa17, pfa19]:
                for val in sqlmap.itervalues():
                    # If reference in another group, add it
                    if item in val:
                        groups = groups.union(val)
    output.append(groups)

# Print ortholog groups
for item in output:
    print '\t'.join(item)


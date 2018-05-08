#!/usr/bin/env python2
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='Fasta file with your genome')
parser.add_argument('glimmer', help='Glimmer prediction file')
parser.add_argument('--translate', help='Translate to protein',
                    action='store_true')
args = parser.parse_args()

# Read the sequence file
seq_record = SeqIO.parse(open(args.fasta), "fasta").next()

outseqrecords = []
# Read the glimmer file, record by record
for inline in file(args.glimmer):
    if '>' in inline:
        seqname = inline.split()[0][1:]
        outfilename = "%s.nfa" % (seqname)
        if args.translate:
            outfilename = "%s.pfa" % (seqname)
        continue
    if "orf" not in inline:
        continue
    orfname, sbegin, send, rf, score = inline.strip().split()
    sbegin = int(sbegin)
    send = int(send)
    rf = int(rf)
    # reverse complement
    if rf < 0:
        sbegin, send = send, sbegin
    sbegin -= 1     # Python indexes start a 0
    score = float(score)
    # split the sequence record
    newseq = seq_record.seq[sbegin:send]
    if rf < 0:
        newseq = newseq.reverse_complement()
    # translate if requested
    if args.translate:
        newseq = newseq.translate(stop_symbol='')
    # add a sequence record to the output
    name = seqname+"_"+orfname
    if rf < 0:
        name = seqname+"_"+orfname+'_rev'
    outseqrecords.append(SeqRecord(newseq,
                                   id=name,
                                   description=''))

SeqIO.write(outseqrecords, open(outfilename, "w"), "fasta")


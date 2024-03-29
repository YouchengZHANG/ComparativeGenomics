


#  ORFsFinder   (by Yuanyuan XI, Youcheng ZHANG)
#  This scripts when runs searches the start and stop codon in a fasta sequence
#  and outputs a multiple fasta file with ORFs sequences and the corresponding ID
#  as well as alternatively outputs the tranlated protein multiple fasta file


#  STEP 1: Import library and add argument
#  Usage: ORFsFinder.py [-h] filename [-t] 
#  filename refers to the input genome, [-t] states the genome types pro:prokaryotes or euk:eukaryotes
#  [--translate] identifies whether outputs the translated protein fasta file
#  Example: $ ./ORFsFinder.py FILENAME -t pro --translate
import re
import argparse
from Bio import SeqIO
ORFsFinder = argparse.ArgumentParser(description='### ORFsFinder for eukaryotes and prokaryotes ###')
ORFsFinder.add_argument('filename',type=str,metavar='FILENAME',default='.',help='input fasta file')
ORFsFinder.add_argument('-t',dest='type',type=str,default='pro',help='pro:prokaryotes or euk:eukaryotes')
ORFsFinder.add_argument('--translate', help='translate to protein', action='store_true')
args = ORFsFinder.parse_args()
print('ORFsFinder program is running...')


#  STEP 2: Predict and extract the ORFs
#  Both forward and reverse ORFs are considered
#  For prokaryotes, minimum length is assumed to be 150bp and Pribnow box content at TSS region is used as a filter
#  For eukaryotes, minimum length is assumed to be 300bp and TATA box content at TSS region is used as a filter
fo = open(args.filename, 'r+')
f = fo.read().splitlines()
dna = f[1]

if args.type == 'pro':
    print('Types confirmed: prokaryotes')
    print('Start finding open reading frames...')
    pattern = re.compile(r'ATG(?:[A-Z]{3}(?<!TAG|TAA|TGA)){50,}(?:TAG|TAA|TGA)')
if args.type == 'euk':
    print('Types confirmed: eukaryotes')
    print('Start finding open reading frames...')
    pattern = re.compile(r'ATG(?:[A-Z]{3}(?<!TAG|TAA|TGA)){100,}(?:TAG|TAA|TGA)')

dna_rev = dna[::-1].translate(str.maketrans("ATGC","TACG"))
orfs_forward_1 = pattern.findall(dna)
orfs_forward_2 = pattern.findall(dna[1:])
orfs_forward_3 = pattern.findall(dna[2:])
orfs_reverse_1 = pattern.findall(dna_rev)
orfs_reverse_2 = pattern.findall(dna_rev[1:])
orfs_reverse_3 = pattern.findall(dna_rev[2:])
orfs_list = list(set(orfs_forward_1 + orfs_forward_2 + orfs_forward_3))
orfs_list_rev = list(set(orfs_reverse_1 + orfs_reverse_2 + orfs_reverse_3))
print(len(orfs_list)) 
print(len(orfs_list + orfs_list_rev))

if args.type == 'pro':
    print('Calculating the promoter content (Pribnow box) at TSS region in predicted orfs...')
    tmp0,tmp1 = [], []
    for orfs in orfs_list:
        s = dna[dna.index(orfs)-15:dna.index(orfs)]
        if ( s.count('A') + s.count('T') ) / 15 >= 0.6:
            tmp0.append(orfs)
    for orfs in orfs_list_rev:
        s = dna_rev[dna_rev.index(orfs)-15:dna_rev.index(orfs)]
        if ( s.count('A') + s.count('T') ) / 15 >= 0.6:
            tmp1.append(orfs)
    orfs_list = tmp0
    orfs_list_rev = tmp1
    print(len(orfs_list)) 
    print(len(orfs_list + orfs_list_rev))

if args.type == 'euk':
    print('Calculating the promoter content (TATA box) at TSS region in predicted orfs...')
    tmp0,tmp1 = [], []
    TATABox = ['TATA']
    for orfs in orfs_list:
        s = dna[dna.index(orfs)-100:dna.index(orfs)-30]
        for i in TATABox:
            if i in s:
                tmp0.append(orfs)
                break
    for orfs in orfs_list_rev:
        s = dna_rev[dna_rev.index(orfs)-100:dna_rev.index(orfs)-30]
        for i in TATABox:
            if i in s:
                tmp1.append(orfs)
                break
    orfs_list = tmp0
    orfs_list_rev = tmp1
    print(len(orfs_list)) 
    print(len(orfs_list + orfs_list_rev))


#####################################################################################
#  STEP 3 (optional): Evaluate prediction performance 
#  ORFsFinder Outcome Evalution program is written in separate script: ORFsEvaluator.py
#  glimmer.predict is used as the true prediction file with specific format
#  Orfs coordinates are compared to measure the prediction outcome

#print('Start evaluating the predicted outcome...')
#print('Obtaining the coordinates of each orf in PREDICTED outcome...')
#tmp2, tmp3 = [], []
#for orfs in orfs_list:
#    tmp2.append([int(dna.index(orfs))+1,int(dna.index(orfs)+len(orfs))])
#for orfs in orfs_list_rev:
#    orfs_rr = orfs[::-1].translate(str.maketrans("ATGC","TACG"))
#    tmp3.append([int(dna.index(orfs_rr)+len(orfs_rr)),int(dna.index(orfs_rr))+1])
#orfs_pred = tmp2 + tmp3

#print('Obtaining the coordinates of each orf in TRUE outcome...')
#truefilename = '34.glimmer.predict' 
#fo_true = open(truefilename,'r+')
#f_true = fo_true.read().splitlines()
#tmp4 = []
#for s in f_true[1:]:
#    tmp4.append([int(s.split()[1]),int(s.split()[2])])
#orfs_true = tmp4

#print('Calculating the confusion matrix... ')
#count_TP, count_FP, count_FN = 0, 0, 0
#for c in orfs_pred: 
#    if c in orfs_true:
#        count_TP += 1 
#for c in orfs_pred:
#    if not c in orfs_true: 
#        count_FP += 1
#for c in orfs_true:
#    if not c in orfs_pred: 
#        count_FN += 1
#dna_length = len(dna)
#avr_length = list4[i]
#count_TN = (2*(dna_length/avr_length)) - (count_TP + count_FP + count_FN)
#ACC = (count_TP + count_TN) / ( count_TP + count_TN + count_FP + count_FN )
#Precision = count_TP / (count_TP + count_FP)
#Recall = count_TP / (count_TP + count_FN)
#F1 = 2*count_TP / (2*count_TP + count_FP + count_FN)
#print('True positive (TP): ' + str(count_TP))
#print('False positive (FP): ' + str(count_FP))
#print('False negative (FN): ' + str(count_FN))
#print('True negative (TN): ' + str(count_TN))
#print('Accuracy (ACC): ' + str(ACC))
#print('Precision (PPV): ' + str(Precision))
#print('Recall (TPR): ' + str(Recall))
#print('F1 (F1 score): ' + str(F1))
#####################################################################################


#  STEP 4: Output and save nucleotide sequence multi-fasta file
#  if --translate is added, also outputs protein multi-fasta file
#  .nfa refers to the nulceotide multi-fasta
#  .pfa refers to the protein sequence multi-fasta ()
print('Creating nulceotide multi-fasta file...')
fn = open(args.filename + '.nfa', 'w')
num = 1
for s in orfs_list:
    fn.write('>'+ str(args.filename)+ '_orf' + str(num).zfill(8) + '\n')
    fn.write( s + '\n')
    num += 1
num = int(len(orfs_list)) + 1  
for s in orfs_list_rev:
    fn.write('>'+ str(args.filename) + '_orf' + str(num).zfill(8) + '_rev' + '\n')
    fn.write( s + '\n')
    num += 1
fn.close()
print(str(args.filename) + '.nfa' + 'has been saved.')

if args.translate:
    fn = open(args.filename + '.nfa', 'r+')
    print('Creating protein sequence multi-fasta file...')
    prot_id, prot_seq = [], []
    for record in SeqIO.parse(fn, 'fasta') :
        ids = record.id
        prot_id.append('>' + str(ids))
        prot_seq.append(str(record.seq.translate(to_stop=True)))

    fp = open(args.filename + '.pfa', 'w')
    for i in range(len(prot_id)):
        fp.write( prot_id[i] + '\n')
        fp.write( prot_seq[i] + '\n')
    fn.close()
    fp.close()
    print(str(args.filename) + '.pfa' + 'has been saved.')
    print('ORFsFinder program is completed.')



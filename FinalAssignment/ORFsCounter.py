


#  ORFsCounter  (by Yuanyuan XI, Youcheng ZHANG)
#  This script when runs calculates the properties of ORFs including
#  the numnber of ORFs, average length, maximum length, minimum length


#  STEP 1: Calculation
filename = ['04.fa.nfa','05.fa.nfa','08.fa.nfa','16.fa.nfa','34.fa.nfa',]
for i in range(5):
    fo = open(filename[i],'r+')
    f = fo.read().splitlines()
    num_ORFs = len(f) / 2
    length = []
    for n in range(1,len(f),2):
        length.append(int(len(f[n])))
    avr_length = sum(length) / num_ORFs
    max_length = max(length)
    min_length = min(length)
    print(str(filename[i]))
    print('The numnber of ORFS: ' + str(num_ORFs))
    print('The average length: ' + str(avr_length))
    print('The maximum length: ' + str(max_length))
    print('The minimum length: ' + str(min_length))


###########################################
#filename = '34_nucleotide.fa'
#fo = open(filename,'r+')
#f = fo.read().splitlines()
#count = 0
#for s in f:
#    if '>' in s:
#        count += 1
#print(count)
    
# 04: 10397
# 05: 957
# 08: 1867
# 16: 6414
# 34: 484 /287    
###########################################



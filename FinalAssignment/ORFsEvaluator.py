


#  ORFsFinder Outcome Evalution   (by Yuanyuan XI, Youcheng ZHANG)
#  This script when runs evalutes the predicted outcome performance 
#  by calculating the confusion matrix including 
#  True positive (TP), False positive (FP), False negative (FN), True negative (TN)
#  Accuracy (ACC), Precision (PPV), Recall (TPR), F1 (F1 score)


#  STEP 1 : Open files and exact the orfs list from predicted multi fasta files
print('Evalution program is running...')
list1 = ['04.fa','05.fa','08.fa','16.fa','34.fa']
list2 = ['04.fa.nfa','05.fa.nfa','08.fa.nfa','16.fa.nfa','34.fa.nfa']
list3 = ['04.glimmer.predict','05.glimmer.predict','08.glimmer.predict','16.glimmer.predict','34.glimmer.predict' ]
list4 = [555,678,633,615,1368]
for i in range(5):
    orfs_list, orfs_list_rev = [], []
    filename = list1[i]
    fo = open(filename, 'r+')
    f = fo.read().splitlines()
    dna = f[1]
    filename2 = list2[i]
    fo2 = open(filename2, 'r+')
    f2 = fo2.read().splitlines()
    for m in range(len(f2)):
        if '>' in f2[m]:
            if not '_rev' in f2[m]:
                orfs_list.append(f2[m+1])
            if '_rev' in f2[m]:
                orfs_list_rev.append(f2[m+1])


#  STEP 2 : Evaluate prediction performance 
#  glimmer.predict is used as the true prediction file with specific format
#  Orfs coordinates are compared to measure the prediction outcome
    print('Start evaluating the predicted outcome...')
    print('Obtaining the coordinates of each orf in PREDICTED outcome...')
    tmp2, tmp3 = [], []
    for orfs in orfs_list:
        tmp2.append([int(dna.index(orfs))+1,int(dna.index(orfs)+len(orfs))])
    for orfs in orfs_list_rev:
        orfs_rr = orfs[::-1].translate(str.maketrans("ATGC","TACG"))
        tmp3.append([int(dna.index(orfs_rr)+len(orfs_rr)),int(dna.index(orfs_rr))+1])
    orfs_pred = tmp2 + tmp3
    #print(len(orfs_pred))

    print('Obtaining the coordinates of each orf in TRUE outcome...')
    truefilename = str(list3[i])
    fo_true = open(truefilename,'r+')
    f_true = fo_true.read().splitlines()
    tmp4 = []
    for s in f_true[1:]:
        tmp4.append([int(s.split()[1]),int(s.split()[2])])
    orfs_true = tmp4
    #print(len(orfs_true))

    print('Calculating the confusion matrix... ')
    print(str(list1[i]))
    count_TP, count_FP, count_FN = 0, 0, 0
    for c in orfs_pred: 
        if c in orfs_true:
            count_TP += 1 
    for c in orfs_pred:
        if not c in orfs_true: 
            count_FP += 1
    for c in orfs_true:
        if not c in orfs_pred: 
            count_FN += 1
    dna_length = len(dna)
    avr_length = list4[i]
    count_TN = (2*(dna_length/avr_length)) - (count_TP + count_FP + count_FN)
    ACC = (count_TP + count_TN) / ( count_TP + count_TN + count_FP + count_FN )
    Precision = count_TP / (count_TP + count_FP)
    Recall = count_TP / (count_TP + count_FN)
    F1 = 2*count_TP / (2*count_TP + count_FP + count_FN)
    print('True positive (TP): ' + str(count_TP))
    print('False positive (FP): ' + str(count_FP))
    print('False negative (FN): ' + str(count_FN))
    print('True negative (TN): ' + str(count_TN))
    print('Accuracy (ACC): ' + str(ACC))
    print('Precision (PPV): ' + str(Precision))
    print('Recall (TPR): ' + str(Recall))
    print('F1 (F1 score): ' + str(F1))

print('Evalution program is completed.')   




##########################################################################
#  TP + FP = num predict 
#  
#     TP  FP  FM
# 04: 399 1652 9997
# 05: 351 590 605
# 08: 726 1423 1136
# 16: 544 1894 5869
# 34: 114 74 369
#
#     TP  FP  FN  TN
# 04: 399 1652 9997
# 05: 351 588 605 1532    2*len(genome)/ave_length
# 08: 726 1423 1136
# 16: 544 1894 5869
# 34: 114 74 369


# 2051  939  2149  2438  188
# 


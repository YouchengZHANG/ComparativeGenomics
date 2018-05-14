


#  Ortholog Extractor and Construction
#  This script when runs firstly creates 'orth_cluster.txt' that contains all the best hits in one cluster file
#  secondly selects 10 different ortholog clusters and extracts the corresponding protein sequences, 
#  then saves each ortholog cluster into 10 separated fasta files, named 'orth_01.fa', 'orth_02.fa', ... , 'orth_10.fa'

Parse01 = open('16_04_parse.txt','r+')
Parse02 = open('16_05_parse.txt','r+')
Parse03 = open('16_08_parse.txt','r+')

f1 = Parse01.read().splitlines()
f2 = Parse02.read().splitlines()
f3 = Parse03.read().splitlines()
#print(f1[0:30])

l1 = [s1.split(' ') for s1 in f1]
l2 = [s2.split(' ') for s2 in f2]
l3 = [s3.split(' ') for s3 in f3]
#print(l1[0:30])

ref1 = [l1[i][0] for i in range(len(f1))]
db1 = [l1[i][2] for i in range(len(f1))]
ref2 = [l2[i][0] for i in range(len(f2))]
db2 = [l2[i][2] for i in range(len(f2))]
ref3 = [l3[i][0] for i in range(len(f2))]
db3 = [l3[i][2] for i in range(len(f2))]
#print(len(ref1),len(ref2),len(ref2))

orth_cluster = [x for x in ref1 if x in ref2 and x in ref3]
#print(len(orth_cluster))

#  Create 'orth_cluster.txt' file
cluster = open('orth_cluster.txt','w')
for i in range(len(orth_cluster)):
	cluster.write(orth_cluster[i] + ' ' + db1[ref1.index(orth_cluster[i])] + ' ' + db2[ref2.index(orth_cluster[i])] + ' ' + db3[ref3.index(orth_cluster[i])] + '\n')
cluster.close()



#  Cluster Genome Extractor
#  This script when runs creates 10 cluster fasta files with ortholog full sequences in each species
#  Select 10 different ortholog clusters and gather corresponding sequences to 10 separated files

Prot01 = open('04_p.fa','r+')
Prot02 = open('05_p.fa','r+')
Prot03 = open('08_p.fa','r+')
Prot04 = open('16_p.fa','r+')

p1 = Prot01.read().splitlines()
p2 = Prot02.read().splitlines()
p3 = Prot03.read().splitlines()
p4 = Prot04.read().splitlines()

orth_select = ['>./16.fa.txt_orf00023_rev', '>./16.fa.txt_orf00025_rev','>./16.fa.txt_orf00209','>./16.fa.txt_orf00325','>./16.fa.txt_orf00404','>./16.fa.txt_orf00446_rev','>./16.fa.txt_orf00774_rev','>./16.fa.txt_orf00896_rev','>./16.fa.txt_orf00977_rev','>./16.fa.txt_orf01008_rev']
for k,i in enumerate(orth_select):
	vars = 'o' + str(int(k)+1)
	names = 'orth_' + str(int(k)+1)
	vars = open(names,'w')
	vars.write(i + '\n' + p4[int(p4.index(i))+1] + '\n' + db1[ref1.index(i)] + '\n' + p1[p1.index(db1[ref1.index(i)])+1] + '\n' + db2[ref2.index(i)] + '\n' + p2[p2.index(db2[ref2.index(i)])+1] + '\n' + db3[ref3.index(i)] + '\n' + p3[p3.index(db3[ref3.index(i)])+1] + '\n')	



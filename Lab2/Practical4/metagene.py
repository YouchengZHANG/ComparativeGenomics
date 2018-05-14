


#  Metagene Extractor and Construction
#  This script when runs creates 'orth_meatagene.fa' that contains all the aligned gene in all cluster files

s_16, s_04, s_05, s_08 = '', '', '', ''

for i in range(1,11):
	names = 'orth_' + str(i) + '_k' + '.fa'
	fo = open(names,'r+')
	f = fo.read().splitlines()
	s_16 += str(f[1])
	s_04 += str(f[3])
	s_05 += str(f[5])
	s_08 += str(f[7])
	fo.close()

fm = open('orth_metagene.fa','w')
fm.write('>16' + '\n' + s_16 + '\n' + '>04' + '\n' + s_04 + '\n' + '>05' + '\n' + s_05 + '\n' + '>08' + '\n' + s_08 + '\n')
fm.close()



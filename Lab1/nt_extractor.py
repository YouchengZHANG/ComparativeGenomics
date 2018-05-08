

def func(filename):
    output = open('34_nucleotide.fa','w')
    fo = open(filename,'r+')
    f = fo.readlines()
    #print(f[0:1252])
    f = f[1252:]
    for i in range(len(f)):
        if '_bp' in f[i]:
            output.write(f[i])
        if 'atcg' in f[i]:
            output.write(f[i])
    output.close


if __name__ == '__main__':
    func('34_nucleotide.out')

#ARNDCQEGHILKMFPSTWYV
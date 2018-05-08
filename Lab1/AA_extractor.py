

def func(filename):
    output = open('34_aminoacid.fa','w')
    fo = open(filename,'r+')
    f = fo.read().splitlines()
    #print(f[0:1252])
    f = f[1252:]
    for i in range(len(f)):
        if not '' == f[i]:
            output.write(f[i] + '\n')
    output.close


if __name__ == '__main__':
    func('34_aminoacid.out')

#ARNDCQEGHILKMFPSTWYV

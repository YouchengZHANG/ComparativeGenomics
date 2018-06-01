# This script is to do statistics about genome.
import sys

#This script is to calculate uniq aa/nt in a fasta file
def uniref_aa(file_name):
    filehandle = open (file_name,'r')
    b = set()
    for line in filehandle:
        if not line.startswith('>'):
            line = line.rstrip('\n')
            b.update(line)
    b = list(b)
    b.sort()
    print (b)


# This script is to calculate nt/aa frequency in a fasta file
def nucl_fre(filename):
    f1 = open(filename,'r')
    freq = dict()
    for line in f1:
        if not line.startswith('>'):
            line = line.strip('\n')
            nucl = list(line)
            #print (len(nucl))
            for nt in nucl:
                if not nt in freq:
                    freq[nt] = 1
                    #print (freq[nt])
                else: 
                    freq[nt] +=1
    print_freq(freq)
    return freq

def dinucl_fre(filename):
    f1 = open(filename, 'r')
    freq = dict()
    for line in f1:
        if not line.startswith('>'):
            line = line.strip('\n')
            #print (len(nucl))
            for i in range(len(line)):
                nt = line[i:i+2]
                #print (nt)
                if not nt in freq:
                    if len(nt) == 1:
                        break
                    else:
                        freq[nt] = 1
                        #print (freq[nt])
                else:
                    freq[nt] +=1
    print_freq(freq)
    return freq

def print_freq(freq_dict):
    sorted_key = sorted(freq_dict.keys())
    for key in sorted_key:
        print("{}\t{}".format(key,freq_dict[key]))
    
if __name__ == '__main__':
    #uniref_aa(sys.argv[1])
    freq1 = nucl_fre(sys.argv[1])
    freq2 = dinucl_fre(sys.argv[1])

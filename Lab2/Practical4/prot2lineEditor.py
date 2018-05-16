


# Edit the raw dataset into real 2 lines dataset

def separate(filename):
    f = open(filename,'r+')
    fr = f.read().splitlines()
    tmp_pid, tmp_seq = '', ''
    lpid, lseq = [],[]

    for n in range(len(fr)):
        if '>' in fr[n]:
            if n == len(fr)-1:
                print('This is the end of the file.')
            else:
                lpid.append(fr[n]+'\n')

        elif not '>' in fr[n]:
            if '>' in fr[n+1]:
                tmp_seq = tmp_seq + fr[n] + '\n'
                lseq.append(tmp_seq)
                tmp_seq = ''
            else:
                tmp_seq = tmp_seq + fr[n]  
    
    o = open('16_p_2line.fa','w')
    for i in range(0,len(lpid)):
       o.write(lpid[i]+lseq[i])

if __name__ == '__main__':
    from datetime import datetime
    start_time = datetime.now()
    print('Program is running...')
    separate('16_p.fa')
    end_time = datetime.now()
    print('Starting from',start_time,'to',end_time)
    print('Running Time: {}'.format(end_time - start_time))

    

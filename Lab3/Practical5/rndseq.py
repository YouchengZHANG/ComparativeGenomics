#! /usr/bin/env python
"prints first argument aa sequences of length second argument to screen"
import random
import sys
number=int(sys.argv[1])
length=int(sys.argv[2])

#number=int(6326)
#length=int(20)

aas="ARNDCEQGHILKMFPSTWYV"
def r20():
    return random.randrange(20)
for i in range(number):
    s=""
    for j in range(length):
        s+=aas[r20()]
    print s

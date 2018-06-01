# This script is to calculate distance among five genomes
import sys
import numpy as np
import math

data = np.genfromtxt(sys.argv[1], dtype=float,skip_header = 1, usecols = range (1,6) )
#data = np.asarray(data)

print (data.shape)
#print (data.shape())
def calculation_mul(data):
    data = np.multiply(data,100)
    result = np.zeros((5,5),dtype = float)
    #print (data)
    for i in range(5):
        for j in range(i+1,5):
            summ = 0
            for col in range(len(list(data))):
                diff = (data[col,j]-data[col,i])**2
                summ +=diff
                #print (data[col,j], data[col,i])
            distance = math.sqrt(summ)
            result[i,j] = result[j,i] = distance
    return result

def calculation_sin(data):
    data = np.multiply(data,100)
    result = np.zeros((5,5),dtype = float)
    #print (data)
    for i in range(5):
        for j in range(i+1,5):
            diff = (data[j]-data[i])**2
            distance = math.sqrt(diff)
            result[i,j] = result[j,i] = distance
    return result
        

if __name__ == "__main__":
    #result = calculation_mul(data)
    result = calculation_sin(data)
    #You could change your file name here
    #np.savetxt("gc_matrix.dis",result, fmt='%.2f', delimiter='\t', newline='\n', header="04fa\t05fa\t08fa\t16fa\t34fa\t")

## this code is used for read a list of number and output the distribution
## input file , number of bin,
## output max, min, the distribution
import scipy
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import StringIO
import re

##################################################
#read the file from the command line
if(len(sys.argv)<2):
    sys.stderr.write('Usage : i. filename\n')
    sys.stderr.write('        ii. number of bins\n')
    exit()
##sys.argv first argv is the python file, if you want to
##input the filename, filename is the 2rd argv
carg = 1
file  = str(sys.argv[carg])
carg += 1
nbins = int(sys.argv[carg])
########################################################

def readfile(file):

    f = open(file,'r')
    # read the first colunms and the second columns
    col1 =[]
    col2 =[]
    for line in f:
        num = line.split()
        col1.append(int(num[0]))  # id of the particles
        col2.append(float(num[1]))  # local density of particles
    col_matrix = np.array([col1,col2])
    #print(col_matrix)
    return col_matrix,nbins
#######################################################
col_matrix, nbins = readfile(file)
#val,bin_edges = np.histogram(col_matrix[1],nbins)
val,bin_edges,patches = plt.hist(col_matrix[1],nbins)
plt.show()

#######################################################
def savefile(val,bin_edges):
    save = 'bin-'+file
    print (save)
    f2 = open(save,'w')
    for i in range(len(val)):
        f2.write('%f %f \n'%(bin_edges[i],val[i]))
    f2.write('%f\n'%bin_edges[-1])
    
savefile(val,bin_edges)

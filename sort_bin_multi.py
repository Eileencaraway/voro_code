## this code is used for read a list of number and output the distribution
## input file , number of bin,
## output max, min, the distribution
import scipy
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
import sys
import StringIO
import re

##################################################
#read the file from the command line
if(len(sys.argv)<2):
    sys.stderr.write('Usage : i. number of bins\n')
    sys.stderr.write('        ii. density\n')
    exit()
##sys.argv first argv is the python file, if you want to
##input the filename, filename is the 2rd argv
carg = 1
nbins = int(sys.argv[carg])
carg +=1
cen_den  = float(sys.argv[carg])
########################################################

def readfile(file):
    f = open(file,'r')
    # read the first colunms and the second columns
    col =[]
    for line in f:
        num = line.split()
        #col1.append(int(num[0]))  # id of the particles, seems id is useless here
        col.append(float(num[1]))  # local density of particles
    col_density = np.array(col)
    #print(col_matrix)
    return col_density  ## output an array
#######################################################
## get a list of file name
FileList = []
for infile in sorted(glob.glob('local_density-*.dat')):
    print infile
    FileList.append(infile)
## define an array to store all the data in all the files
density = np.zeros((len(FileList),1000)) # 1000 particles
for i in range(len(FileList)):
    density[i] = readfile(FileList[i])
#######################################################
# produce a bin to store all the numbers
maxall = np.amax(density)
minall = np.amin(density)
#HalfLength = max(maxall - cen_den, cen_den - minall)
#def_min = cen_den - HalfLength
#def_max = cen_den + HalfLength
#######################################################
val = np.zeros((len(FileList),nbins))
for i in range(len(FileList)):
    val[i],bin_edges = np.histogram(density[i],nbins,range=(minall,maxall))
bin_size = bin_edges[1]-bin_edges[0]
avg_val = np.mean(val,axis=0)
nor_val = avg_val[:]/sum(avg_val)
#######################################################
f2 = open('bin-local-density.txt','w')
for i in range(len(nor_val)):
    f2.write('%f %f \n'%(bin_edges[i],nor_val[i]))
f2.write('%f\n'%bin_edges[-1])

#val,bin_edges,patches = plt.hist(col_matrix[1],nbins)
plt.plot(bin_edges[0:-1]+bin_size/2,nor_val)
plt.xlabel('local density')
plt.ylabel('number of particles')
plt.title('phi%f'%cen_den)
plt.savefig("distribution.png")
#normalize for val

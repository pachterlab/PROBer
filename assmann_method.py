#!/usr/bin/env python

from sys import argv, exit
from math import log

if len(argv) != 5:
    print("Usage: python assmann_method.py length minus_channel_data.dat plus_channel_data.dat output.theta")
    exit(-1)

length = int(argv[1])

def analyzeCounts(input_file):
    """ Load counts and analyze them according to Assmann's paper """

    counts = [0] * (length + 1)
        
    fin = open(input_file, 'r')
    for line in fin:
        line = line.rstrip('\n')
        arr = line.split('\t')
        counts[int(arr[0])] += 1
    fin.close()

    for i in xrange(length + 1):
        if counts[i] == 0:
            counts[i] = 0.01
    
    avg = 0.0
    for c in counts:
        avg += log(c)
    avg /= (length + 1)
    
    result = [0.0] * length
    for i in xrange(length):
        result[i] = counts[i + 1] / avg
        
    return result

M = analyzeCounts(argv[2])
P = analyzeCounts(argv[3])

theta = [max(P[i] - M[i], 0) for i in xrange(length)]
denom = sum(theta)
theta = [str(i / denom) for i in theta]

fout = open(argv[4], 'w')
fout.write(str(denom) + '\t' + '\t'.join(theta) + '\n')
fout.close()

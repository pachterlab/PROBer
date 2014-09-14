#!/usr/bin/env python

from sys import argv, exit
from math import log

if len(argv) != 3:
    print("Usage: python beta2theta.py input.beta output.theta")
    exit(-1)

fin = open(argv[1], 'r')
line = fin.readline().rstrip('\n')
fin.close()

arr = line.split('\t')
length = int(arr[0])

theta = [-log(1.0 - float(i)) for i in arr[1:]]
c = sum(theta)
theta = [str(i / c) for i in theta]

fout = open(argv[2], 'w')
fout.write(str(c) + '\t' + '\t'.join(theta) + '\n')
fout.close()









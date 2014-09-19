#!/usr/bin/env python

from sys import argv, exit
import pysam

if len(argv) != 3:
    print("Usage: python extractHeader.py input.bam output_header.txt")
    exit(-1)

sam_in = pysam.Samfile(argv[1], "rb")
out = open(argv[2], "w")
ref = sam_in.references
lens = sam_in.lengths
size = len(sam_in.references)
for i in xrange(size):
    out.write(ref[i] + '\t' + str(int(lens[i])) + '\n')
out.close()

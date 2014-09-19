#!/usr/bin/env python

from sys import argv, exit
import numpy
import numpy.linalg
import scipy.stats

primer_len = 6

class Transcript:
    def __init__(self, name, len):
        self.name = name
        self.TPM = len - primer_len + 1
        self.cor = 0.0
        self.l1 = 1e20

def calcTPM(header_file, theta_file):
    """Calculate TPM values"""
    trans = []
    fin = open(header_file, 'r')
    for line in fin:
        arr = line[:-1].split('\t')
        trans.append(Transcript(arr[0], float(arr[1])))
    fin.close()

    fin = open(theta_file, 'r')
    size = int(fin.readline()[:-1])
    thetas = fin.readline()[:-1].split('\t')[1:]
    denom = 0.0
    for i in xrange(size):
        trans[i].TPM = float(thetas[i]) / trans[i].TPM
        denom += trans[i].TPM
    for i in xrange(size):
        trans[i].TPM = trans[i].TPM / denom * 1e6
    fin.close()

    return trans

def calcCor(gtfile, estfile):
    """Calculate correlations"""
    fgt = open(gtfile, "r")
    fest = open(estfile, "r")

    fgt.readline()
    fest.readline()
    size = len(trans)

    for i in xrange(size):
        gt = [float(x) for x in fgt.readline()[:-1].split('\t')[2:]]
        est = [float(x) for x in fest.readline()[:-1].split('\t')[2:]]
        if sum(gt) > 0 and sum(est) > 0:
            trans[i].cor = scipy.stats.pearsonr(gt, est)[0]
            trans[i].l1 = numpy.linalg.norm(numpy.array(est) - numpy.array(gt), ord=1) / len(gt)
        if (i + 1) % 10000 == 0:
            print("FIN " + str(i + 1))

    fgt.close()
    fest.close()

def printResults(output_name):
    """Print out results"""

    fout = open(output_name + ".tpm", "w")
    for tran in trans:
        fout.write(tran.name + '\t' + str(tran.TPM) + '\t' + str(tran.cor) + '\t' + str(tran.l1) + '\n')
    fout.close()
    
    fout = open(output_name + ".stat", "w")
    intervals = [(1.0, 10.0), (10.0, 100.0), (100.0, 1000.0), (1000.0, 10000.0), (10000.0, 1000000.0)]
    for tp in intervals:
        cor_arr = [x.cor for x in trans if x.TPM >= tp[0] and x.TPM < tp[1]]
        l1_arr = [x.l1 for x in trans if x.TPM >= tp[0] and x.TPM < tp[1]]
        print("[" + str(tp[0]) + " , " + str(tp[1]) + "]: size = " + str(len(cor_arr)) + ", cor.mean = " + str(numpy.mean(cor_arr)) + ", cor.median = " + str(numpy.median(cor_arr)) + ", avgL1.mean = " + str(numpy.mean(l1_arr)) + ", avgL1.median = " + str(numpy.median(l1_arr)))
        fout.write("[" + str(tp[0]) + " , " + str(tp[1]) + "]: size = " + str(len(cor_arr)) + ", cor.mean = " + str(numpy.mean(cor_arr)) + ", cor.median = " + str(numpy.median(cor_arr)) + ", avgL1.mean = " + str(numpy.mean(l1_arr)) + ", avgL1.median = " + str(numpy.median(l1_arr)) + '\n')

    fout.close()

if len(argv) != 6:
    print("Usage: analyzeRes.py header.txt ground_truth.theta ground_truth.[gamma/beta] estimated.[gamma/beta] output_name")
    exit(-1)

trans = calcTPM(argv[1], argv[2])
print("calcTPM finished!")
calcCor(argv[3], argv[4])
print("calcCor finished!")
printResults(argv[5])




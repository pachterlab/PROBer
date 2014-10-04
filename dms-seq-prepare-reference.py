#!/usr/bin/env python

import os
import sys
import argparse
#import subprocess
import glob

import utils
from utils import MyFormatter, MyParser, runProg

#Parse options and arguments
              
parser = MyParser(formatter_class=MyFormatter)

parser.add_argument("ref_fastas", help = "Either a comma-separated list of Multi-FASTA formatted files OR a directory name. If a directory name is specified, This program will read all files with suffix \".fa\" or \".fasta\" in this directory. The files should contain either the sequences of transcripts or an entire genome, depending on whether the --gtf or --gff3  option is used.", metavar = "reference_fasta_file(s)")
parser.add_argument("ref_name", help = "The name of the reference used. This program will generate several reference-related files that are prefixed by this name. This name can contain path information (e.g. /ref/mm9).", metavar = "reference_name")

group = parser.add_mutually_exclusive_group()
group.add_argument("--gtf", help = "<file> is in GTF format. This program will assume reference_fasta_file(s) contains genome sequences and extract transcript sequences using the gene annotation specified in <file>.", metavar = "<file>")
group.add_argument("--gff3", help = "<file> is in GFF3 format. This program will assume reference_fasta_file(s) contains genome sequences and extract transcript sequences using the gene annotation specified in <file>.", metavar = "<file>")

parser.add_argument("--transcript-to-gene-map", help = "L:Use information from <file> to map from transcript (isoform) ids to gene ids.\nEach line of <file> should be of the form:\n\ngene_id\ttranscript_id\n\nwith the two fields separated by a tab character.", metavar = "<file>", dest = "tran2gene")

parser.add_argument("--bowtie", help = "Build Bowtie indices.", action = "store_true")
parser.add_argument("--bowtie-path", help = "The path to the Bowtie executables.", metavar = "<path>")
parser.add_argument("--bowtie2", help = "Build Bowtie2 indices.", action = "store_true")
parser.add_argument("--bowtie2-path", help = "The path to the Bowtie2 executables.", metavar = "<path>")
parser.add_argument("--hobbes2", help = "Build Hobbes2 indices.", action = "store_true")
parser.add_argument("--hobbes2-path", help = "The path to the Hobbes2 executables.", metavar = "<path>")
parser.add_argument("--hobbes2-p", help = "Number of threads to build Hobbes2 indices", type = int, default = 1, metavar = "<int>")
parser.add_argument("--hobbes2-gram-length", help = "Hobbes2 indices' gram length.", type = int, default = 11, metavar = "<int>")

parser.add_argument("-q", "--quiet", help = "Suppress the output of logging information.", action = "store_true")

args = parser.parse_args()

#Set executable directory

mydir = os.path.realpath(os.path.dirname(sys.argv[0]))
os.environ["PATH"] = mydir + os.pathsep + os.getenv("PATH", ".")
os.environ["PYTHONPATH"] = mydir + os.pathsep + os.getenv("PYTHONPATH", ".")

#Run programs

utils.demo = True
       
ref_name = os.path.expanduser(args.ref_name)

# Prepare GTF file
gtf_file = None
if args.gff3 != None:
    gtf_file = ref_name + os.extsep + "gtf"
    runProg("gff3_to_gtf {} {}".format(args.gff3, gtf_file))
elif args.gtf != None:
    gtf_file = args.gtf

# Obtain all reference files
tmp_list = args.ref_fastas.split(',')
fasta_files = []
for afile in tmp_list:
    afile = os.path.expanduser(afile)
    if os.path.isfile(afile):
        fasta_files.append(afile)
    elif os.path.isdir(afile):
        fasta_files.extend(glob.glob("{}/*.fa".format(afile)))
        fasta_files.extend(glob.glob("{}/*.fasta".format(afile)))
    else:
        print("{} does not exist!".format(afile))
        sys.exit(-1)
if len(fasta_files) <= 0:
    print("reference_fasta_file(s) is empty! Please check if the directory name is correct or sequence files are sufficed with either '.fa' or '.fasta'.")
    sys.exit(-1)
        
command = ""

# Extract references
if gtf_file != None:
    command = "dms-seq-extract-reference-transcripts " + ref_name + " " + ("1" if args.quiet else "0") + " " + gtf_file
    command += (" 1 {}".format(args.tran2gene) if args.tran2gene != None else " 0")
    command += " {}".format(" ".join(fasta_files))
    runProg(command)
else:
    command = "dms-seq-synthesis-reference-transcripts " + ref_name + " " + ("1" if args.quiet else "0")
    command += (" 1 {}".format(args.tran2gene) if args.tran2gene != None else " 0")
    command += " {}".format(" ".join(fasta_files))
    runProg(command)

command = "dms-seq-preref {0}.transcripts.fa 1 {0}".format(ref_name)
if args.quiet:
    command += " -q"
runProg(command)

if args.bowtie:
    command = ((args.bowtie_path + os.sep) if args.bowtie_path != None else "") + "bowtie-build -f"
    if args.quiet:
        command += " -q"
    command += " {0}.n2g.idx.fa {0}".format(ref_name)
    runProg(command)
elif args.bowtie2:
    command = ((args.bowtie2_path + os.sep) if args.bowtie2_path != None else "") + "bowtie2-build -f"
    if args.quiet:
        command += " -q"
    command += " {0}.idx.fa {0}".format(ref_name)
    runProg(command)
elif args.hobbes2:
    command = ((args.hobbes2_path + os.sep) if args.hobbes2_path != None else "") + "hobbes-index"
    command += " --sref {0}.idx.fa -i {0}.hix -g {1!s} -p {2!s}".format(ref_name, args.hobbes2_gram_length, args.hobbes2_p)
    runProg(command)  
     
      

        
    
    

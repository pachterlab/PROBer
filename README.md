README
======

[Bo Li](http://math.berkeley.edu/~bli) \(bli25 at berkeley dot edu\)

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Compilation](#compilation)
* [Usage](#usage)
* [Example](#example)
* [Authors](#authors)
* [Acknowledgements](#acknowledgements)

* * *

## <a name="introduction"></a> Introduction

## <a name="compilation"></a> Compilation

To compile, run

    make

## <a name="usage"></a> Usage

### Prepare Reference Sequences

To prepare reference sequence, you should run
'dms-seq-prepare-reference'. Run

    dms-seq-prepare-reference --help

to get usage information.

### Estimate RNA Structure Parameters

To estimate RNA structure parameters, you should run
'dms-seq-estimate-parameters' twice. The first run processes '-'
channel data and the second run processes '+' channel data. Run

    dms-seq-estimate-parameters --help

to get usage information.

### Simulation

To simulate reads, you should run 'dms-seq-simulate-reads'. Run

    dms-seq-simulate-reads

to get usage information.

## <a name="example"></a> Example

Suppose we have arabidopsis genome and gene annotation in two files:
'TAIR10_chr_all.fa' and 'TAIR10_GFF3_genes.gff'. We choose the
reference name as 'arabidopsis' and are only interested in mRNA and
rRNA. The data we have are single-end reads, with minus channel reads
in 'minus.fq' and plus channel reads in 'plus.fq'. We use Bowtie
aligner to align reads and assume Bowtie executables are under
'/sw/bowtie'. We choose sample name as 'dms_sample'. We use 40
cores. In the end, we simulate 10M single-end reads with output name
'dms_sim'.

The commands are listed below:

dms-seq-prepare-reference --gff3 TAIR10_GFF3_genes.gff --gff3-RNA-pattern mRNA,rRNA --bowtie --bowtie-path /sw/bowtie TAIR10_chr_all.fa arabidosis/arabidosis
dms-seq-estimate-parameters -p 40 --primer-length 6 --size-selection-min 150 --size-selection-max 650 --gamma-init 0.01 --beta-init 0.01 --bowtie-path /sw/bowtie arabidosis/arabidosis dms_sample minus --reads minus.fq
dms-seq-estimate-parameters -p 40 --primer-length 6 --size-selection-min 150 --size-selection-max 650 --gamma-init 0.01 --beta-init 0.01 --bowtie-path /sw/bowtie arabidosis/arabidosis dms_sample plus --reads plus.fq
dms-seq-simulate-reads arabidosis/arabidosis dms_sample.temp/dms_sample_minus.config dms_sample dms_sample.stat/dms_sample_minus.read_model minus 10000000 dms_sim
dms-seq-simulate-reads arabidosis/arabidosis dms_sample.temp/dms_sample_plus.config dms_sample dms_sample.stat/dms_sample_plus.read_model plus 10000000 dms_sim
 
## <a name="authors"></a> Authors

This package is implemented by Bo Li. 

## <a name="acknowledgements"></a> Acknowledgements

Thanks Akshay Tambe, Sharon Aviran, and Lior Pachter for their super
helpful suggestions and feedbacks on implementing this package.

Part of this package's codes are adopted from
[RSEM](http://deweylab.biostat.wisc.edu/rsem).

This package uses the
[Boost C++](http://www.boost.org) and
[samtools](http://samtools.sourceforge.net) libraries.

# PROBer

Quantitative modeling of transcriptome-wide RNA structure-probing experiments 

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
* [License](#license)

* * *

## <a name="introduction"></a> Introduction

## <a name="compilation"></a> Compilation

To compile, run

    ```
    make
    ```
## <a name="usage"></a> Usage

### Prepare Reference Sequences

To prepare reference sequence, you should run
`PROBer prepare`. Run

    ```
    PROBer prepare --help
    ```

to get usage information.

### Estimate RNA Structure Parameters

To estimate RNA structure parameters, you should run
`PROBer estimate`. Run

    ```
    PROBer estimate --help
    ```

to get usage information.

### Simulation

To simulate reads, you should run `PROBer simulate`. Run

    ```
    PROBer simulate --help
    ```

to get usage information.

## <a name="example"></a> Example

Suppose we have arabidopsis genome and gene annotation in two files:
'TAIR10_chr_all.fa' and 'TAIR10_GFF3_genes.gff'. We choose the
reference name as 'arabidopsis' and are only interested in mRNA and
rRNA. The data we have are single-end reads with read length 37bp,
with minus channel reads in 'minus.fq' and plus channel reads in
'plus.fq'. The primer length is 6bp, the size selection range is from
21bp to 526bp. We use Bowtie aligner to align reads and assume Bowtie
executables are under '/sw/bowtie'. We choose sample name as
'test_sample'. We use 40 cores. In the end, we simulate 10M single-end
reads with output name 'test_sim'.

The commands are listed below:

PROBer prepare --gff3 TAIR10_GFF3_genes.gff --gff3-RNA-pattern mRNA,rRNA --bowtie --bowtie-path /sw/bowtie TAIR10_chr_all.fa arabidosis/arabidosis
PROBer estimate -p 40 --primer-length 6 --size-selection-min 21 --size-selection-max 526 --read-length 37 --bowtie-path /sw/bowtie arabidosis/arabidosis test_sample --reads minus.fq plus.fq
PROBer simulate arabidosis/arabidosis test_sample.temp/test_sample_minus.config test_sample minus 10000000 test_sim
PROBer simulate arabidosis/arabidosis test_sample.temp/test_sample_plus.config test_sample plus 10000000 test_sim
 
## <a name="authors"></a> Authors

This package is implemented by Bo Li. 

## <a name="acknowledgements"></a> Acknowledgements

Thanks Akshay Tambe, Sharon Aviran, and Lior Pachter for their super
helpful suggestions and feedbacks on implementing this package. Thanks
Harold Pimentel and Pall Melsted for their help on CMake, website and
markdown documents.

Part of this package's codes are adopted from
[RSEM](http://deweylab.biostat.wisc.edu/rsem).

This package uses the
[Boost C++](http://www.boost.org) and
[samtools](http://samtools.sourceforge.net) libraries.

## <a name="license"></a> License

PROBer is licensed under the [GNU General Public License
v3](http://www.gnu.org/licenses/gpl-3.0.html).

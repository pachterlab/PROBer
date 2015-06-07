---
layout: page
show_meta: false
title: "Getting Started"
header:
   image_fullwidth: "campus.jpg"
permalink: "/getting_started/"
---

The short tutorial below explains how to run __PROBer__ using an imaginary example. Please go over the following steps before you continue:

1. [Download][1] PROBer.
2. [Install][2] PROBer. 
3. [Go over the manual][3].

[1]: {{ site.url }}/download/
[2]: {{ site.url }}/resource/
[3]: {{ site.url }}/manual/

If you want __PROBer__ to align raw reads for you, you also need to install [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) and/or [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).  

#### Example

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

~~~
PROBer prepare --gff3 TAIR10_GFF3_genes.gff --gff3-RNA-pattern mRNA,rRNA --bowtie --bowtie-path /sw/bowtie TAIR10_chr_all.fa arabidosis/arabidosis
PROBer estimate -p 40 --primer-length 6 --size-selection-min 21 --size-selection-max 526 --read-length 37 --bowtie-path /sw/bowtie arabidosis/arabidosis test_sample --reads minus.fq plus.fq
PROBer simulate arabidosis/arabidosis test_sample.temp/test_sample_minus.config test_sample minus 10000000 test_sim
PROBer simulate arabidosis/arabidosis test_sample.temp/test_sample_plus.config test_sample plus 10000000 test_sim
~~~


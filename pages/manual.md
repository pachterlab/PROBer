---
layout: page
show_meta: false
title: "Manual"
header:
   image_fullwidth: "RNA.jpg"
permalink: "/manual/"
---

Type `PROBer` to see the all usage commands:

~~~
usage: PROBer [-h] {prepare,estimate,simulate,iCLIP,version} ...

PROBer is a program to quantify chemical modification profiles for a general set of 'toeprinting' assays.

optional arguments:
  -h, --help            show this help message and exit

commands:
    prepare             This command is used to prepare PROBer references.
    estimate            This command is used to estimate chemical modification
                        profiles.
    simulate            This command is used to simulate reads.
    iCLIP               This command is used to allocate multi-mapping reads
                        for iCLIP data.
    version             Show version information.
~~~

PROBer has four usage commands. They are:

#### prepare

`PROBer prepare` extracts transcript sequences, prepares its reference files, and optionally builds Bowtie/Bowtie2 indices. For iCLIP data, you should turn on `--genome` option and this command will only build genome indices for aligners. Type `PROBer prepare --help` to get its usage information:

~~~
usage: PROBer prepare [-h] [--gtf <file> | --gff3 <file>]
                      [--gff3-RNA-pattern <pattern>]
                      [--transcript-to-gene-map <file>] [--genome] [--bowtie]
                      [--bowtie-path <path>] [--bowtie2]
                      [--bowtie2-path <path>] [-q]
                      reference_fasta_files reference_name

This program lets PROBer to build its references and optionally build
Bowtie/Bowtie2 indices. For iCLIP data, users can use this program to build
Bowtie/Bowtie2 indices for their genomes

positional arguments:
  reference_fasta_file(s)
                        Either a comma-separated list of Multi-FASTA formatted
                        files OR a directory name. If a directory name is
                        specified, This program will read all files with
                        suffix ".fa" or ".fasta" in this directory. The files
                        should contain either the sequences of transcripts or
                        an entire genome, depending on whether the --gtf or
                        --gff3 option is used.
  reference_name        The name of the reference used. This program will
                        generate several reference-related files that are
                        prefixed by this name. This name can contain path
                        information (e.g. /ref/mm9).

optional arguments:
  -h, --help            show this help message and exit
  --gtf <file>          <file> is in GTF format. This program will assume
                        reference_fasta_file(s) contains genome sequences and
                        extract transcript sequences using the gene annotation
                        specified in <file>. (default: None)
  --gff3 <file>         <file> is in GFF3 format. This program will assume
                        reference_fasta_file(s) contains genome sequences and
                        extract transcript sequences using the gene annotation
                        specified in <file>. (default: None)
  --gff3-RNA-pattern <pattern>
                        <pattern> is a comma-separated list of transcript
                        categories, e.g. 'mRNA,rRNA'. Only transcripts that
                        match the <pattern> will be extracted. (default: mRNA)
  --transcript-to-gene-map <file>
                        Use information from <file> to map from transcript (isoform) ids to gene ids.
                        Each line of <file> should be of the form:
                        
                        transcript_id gene_id
                        
                        with the two fields separated by a tab character. (default: None)
  --genome              This option is required and only used for iCLIP data;
                        it allows PROBer to call Bowtie/Bowtie2 to build their
                        indices. (default: False)
  --bowtie              Build Bowtie indices. (default: False)
  --bowtie-path <path>  The path to the Bowtie executables. (default: None)
  --bowtie2             Build Bowtie2 indices. (default: False)
  --bowtie2-path <path>
                        The path to the Bowtie2 executables. (default: None)
  -q, --quiet           Suppress the output of logging information. (default:
                        False)

OUTPUT:
  PROBer reference files prefixed by 'reference_name'.
~~~

#### estimate

`PROBer estimate` quantifies chemical modification profiles using sequenced toeprinting data. Its accepts both single-end and paired-end reads either as unprocessed FASTA/FASTQ files or aligned SAM/BAM/CRAM files. Type `PROBer estimate --help` to get its usage information:

~~~
usage: PROBer estimate [options] reference_name sample_name (--alignments input_plus.(sam|bam|cram) [input_minus.(sam|bam|cram)] | --reads plus_channel_mate1_read_file(s) [plus_channel_mate2_read_file(s)] [minus_channel_mate1_read_file(s) [minus_channel_mate2_read_file(s)]])

DESCRIPTION: This program helps users to align reads and estimate RNA
structure parameters. By default, it requires data from both control and
treatment groups. But it works with only treatment data as well.

positional arguments:
  reference_name        The name of the reference used. Users should have run
                        'PROBer prepare' with this name before running this
                        program.
  sample_name           The output name of this run. All outputs use this name
                        as their prefixes.

optional arguments:
  -h, --help            show this help message and exit
  --time                Output time consumed by each step. (default: False)
  --memory              Output memory used by each step. (default: False)
  -q, --quiet           Suppress the output of logging information. (default:
                        False)

Input:
  Input alignments or reads, options are mutually exclusive. If input are
  alignments, all alignments of a same read should group together and each
  paired-end alignment's two mates should be adjacent.

  --alignments input_plus.(sam/bam/cram) [input_minus.(sam/bam/cram)] [input_plus.(sam/bam/cram) [input_minus.(sam/bam/cram)] ...]
                        Input are alignments in SAM/BAM/CRAM formats. If only
                        one alignment file is provided, PROBer assumes control
                        data are not available. (default: None)
  --reads mate_read_file(s) [mate_read_file(s) ...]
                        Input are read files.
                        plus_channel_mate1_read_file(s) and minus_channel_mate1_read_file(s) are comma-separated lists of files containing single-end reads or first mates of paired-end reads
                        plus_mate2_read_file(s) and minus_mate2_read_file(s), present only if '--paired-end' is enabled, are comma-separated lists of files containing second mates of paired-end reads
                        By default, these files should be in FASTQ format. If '--no-quality-scores' is specified, multi-FASTA format files are expected instead.
                        Minus channel reads may be omitted if no control data are available.
                         (default: None)

Basic options:
  --no-quality-scores   Input reads do not contain quality scores. (default:
                        False)
  --paired-end          Input reads are paired-end reads. (default: False)
  -p <int>, --number-of-threads <int>
                        Number of threads this program can use. (default: 1)
  --output-bam          Output transcript BAM file. (default: False)
  --output-logMAP       Output the log MAP probability, which can be used to
                        select priors. (default: False)
  --keep-intermediate-files
                        If PROBer should keep intermediate files. (default:
                        False)

Structure-seq related:
  Set necessary parameters for generating a config file.

  --primer-length <int>
                        Random primer length. (default: 6)
  --size-selection-min <int>
                        The minimum fragment length that can pass the size
                        selection step. (default: None)
  --size-selection-max <int>
                        The maximum fragment length that can pass the size
                        selection step. (default: None)
  --gamma-init <float>  Initial value for all gammas. (default: 0.0001)
  --beta-init <float>   Initial value for all betas. (default: 0.0001)
  --read-length <int>   Read length before trimming adaptors. (default: None)
  --maximum-likelihood  Use maximum likelihood estimates. (default: False)

Alignment options:
  User can choose from Bowtie and Bowtie2. All reads with more than 200
  alignments will be filtered by this script.

  --bowtie              Use bowtie aligner to align reads, with Bowtie
                        parameters "--norc -p number_of_threads -a -m 200 -S".
                        If "--paired-end" is set, additionaly enable Bowtie
                        parameters "-I 1 -X 1000 --chunkmbs 1024". (default:
                        True)
  --bowtie-path <path>  The path to Bowtie executables. (default: None)
  --bowtie2             Use bowtie2 aligner to align reads, indel alignments
                        enabled, with Bowtie2 parameters "--norc -p
                        number_of_threads -k 201". If "--paired-end" is set,
                        additionaly enable Bowtie2 parameters "-I 1 -X 1000
                        --no-mixed --no-discordant". (default: False)
  --bowtie2-path <path>
                        The path to Bowtie2 executables. (default: None)

OUTPUTS:
  sample_name.expr
    Isoform level expression estimates. The first line contains column names separated by a tab character:
    
    transcript_id length effective_length expected_count_minus expected_count_plus TPM FPKM
    
    transcript_id gives the transcript's name. length is the transcript length. effective_length represents the number of positions that can generate a fragment. It is equal to length - primer_length + 1. expected_count_minus is the sum of posterior probabilities of reads coming from this transcript in the (-) channel. expected_count_plus is the counts from (+) channel. TPM is transcript per million. FPKM is fragment per kilobase per millon reads.

    In the rest lines of the file, each line describes a transcript according to the defined columns.

  sample_name.beta
    Estimated beta parameters for each transcript. The first line contains the total number of transcripts. Then each line describes estimated parameters for a different transcript. Within each line, the first field gives the transcript name, the second field provides the number of estimated beta parameters, which is equal to transcript length - primer length. In the end, estimated beta values at each position were given (from 5' end to 3' end).

  sample_name.gamma
    Estimated gamma parameters for each transcript. The first line contains the total number of transcripts. Then each line describes estimated parameters for a different transcript. Within each line, the first field gives the transcript name, the second field provides the number of estimated gamma parameters, which is equal to transcript length - primer length. In the end, estimated gamma values at each position were given (from 5' end to 3' end).

  sample_name_plus.bam
    Only generated when '--output-bam' option is set.

    It is a BAM-formatted file that contains annotated '+' channel read alignments in transcript coordinates. For each alignable BAM line, The MAPQ field is set to min(100, floor(-10 * log10(1.0 - w) + 0.5)), where w is the posterior probability of that alignment being the true mapping of a read. In addition, a new tag ZW:f:value is added, where the value is a single precision floating number representing the posterior probability. All filtered alignment lines has a ZF:A:! tag to identify that it is filtered. Please note that 'ZW' and 'ZF' tags are reserved for PROBer and users need to make sure the aligner output or input BAM/SAM file does not contain these two tags unless the input BAM file is produced by PROBer and alignment/filtering criteria are not changed. Because this file contains all alignment lines produced by the aligner, it can also be used as a replacement of the aligner generated BAM/SAM file.

  sample_name_minus.bam
    Only generated when '--output-bam' option is set.

    It is a BAM-formatted file that contains annotated '-' channel read alignments in transcript coordinates. The annotation format is exactly the same as the one used in 'sample_name_plus.bam'.

  sample_name.logMAP
    Only generated when '--output-logMAP' option is set.

    This file contains the log MAP probability of the observed data given current parameter settings, which can be used to select appropriate priors.

  sample_name.stat
    This folder contains learned model parameters from data. In the folder, 'sample_name_minus.theta' contains the estimated read generating probabilities from '-' channel. 'sample_name_minus.read_model' contains the estimated sequencing error model from '-' channel. 'sample_name_plus.theta' contains the estimated read generating probabilities from '+' channel. 'sample_name_plus.read_model' contains the estimated sequencing error model from '+' channel. The files contained in this folder can be used for simulation.

  sample_name.temp
    This is a temporary folder contains intermediate files. It will be deleted automatically after the program finishes unless '--keep-intermediate-files' option is on.
~~~

#### simulate

`PROBer simulate` is used to generate simulation data based on model parameters learned from real data. If the real data are paired-end reads, it will simulate paired-end reads. Otherwise, it will simulate single-end reads. Type `PROBer simulate --help` to get its usage information:

~~~
usage: PROBer simulate [-h] [--seed <uint32>] [--no-control]
                       reference_name config_file sample_name channel
                       number_of_reads output_name

This program simulates reads using parameters learned from real data by
program 'estimate'.

positional arguments:
  reference_name   The reference's name, should be same as the one used in
                   programs 'prepare' and 'estimate'.
  config_file      A configuration file containting primer length, size
                   selection min and max fragment size etc.
                   'sample_name.temp/sample_name_minus.config' and
                   'sample_name.temp/sample_name_plus.config' can be used
                   here.
  sample_name      This should be the 'sample_name' used in 'PROBer-estimate-
                   parameters'. No slash should be in the end of this string.
  channel          Which channel to simulate. 'minus' stands for the mock-
                   treated channel and 'plus' stands for the modification-
                   treated channel.
  number_of_reads  Number of reads to simulate.
  output_name      Output files' prefix

optional arguments:
  -h, --help       show this help message and exit
  --seed <uint32>  The seed initializing the random number generator used in
                   the simulation. (default: None)
  --no-control     Indicate if the data used to learn simulation parameters do
                   not have a control. (default: False)

OUTPUT:
  If single-end reads are simulated, this program produces 'output_name_(minus|plus).(fa|fq)'. If paired-end reads are simulated, this program produces 'output_name_(minus|plus)_1.(fa|fq)' and 'output_name_(minus|plus)_2.(fa|fq).
~~~

#### iCLIP

`PROBer iCLIP` allocates multi-mapping reads for iCLIP data. Type `PROBer iCLIP --help` to get its usage information:

~~~
usage: PROBer iCLIP [options] sample_name {--alignments input_alignments.[sam/bam/cram] | --reads mate1_read_file(s) [mate2_read_file(s)]}

This program allocates multi-mapping reads for iCLIP data.

positional arguments:
  sample_name           The output name of this run. All outputs use this name
                        as their prefixes.

optional arguments:
  -h, --help            show this help message and exit
  -q, --quiet           Suppress the output of logging information. (default:
                        False)

Input:
  Input alignments or reads, options are mutually exclusive. If input are
  alignments, all alignments of a same read should group together and each
  paired-end alignment's two mates should be adjacent.

  --alignments alignment_file.[sam/bam/cram]
                        Input are alignments in SAM/BAM/CRAM format. (default:
                        None)
  --reads mate_read_file(s) [mate_read_file(s) ...]
                        Input are comma-separated lists of files containing
                        single-end or paired-end reads. If input are single-
                        end reads, only one list is required. If input are
                        paired-end reads, i.e. '--paired-end' is set, PROBer
                        needs two lists --- one for the first mates and the
                        other for the second mates. By default, these files
                        should be in FASTQ format. If '--no-quality-scores' is
                        specified, multi-FASTA format files are expected
                        instead. (default: None)

Basic options:
  --no-quality-scores   Input reads do not contain quality scores. (default:
                        False)
  --paired-end          Input reads are paired-end reads. (default: False)
  -p <int>, --number-of-threads <int>
                        Number of threads this program can use. (default: 1)
  --keep-intermediate-files
                        If PROBer should keep intermediate files. (default:
                        False)

iCLIP options:
  --half-window-size <int>
                        PROBer will borrow information from adjacent crosslink
                        sites within plus/minus half window size to help
                        allocating multi-mapping reads. (default: 25)
  --rounds <int>        Number of EMS iterations to run. (default: 100)
  --maximum-read-length <int>
                        The maximum possible read length. You may set this
                        option only if '--no-qualities' is set. (default:
                        1000)

Alignment options:
  User can choose from Bowtie and Bowtie2. All reads with more than 100
  alignments will be filtered by this script.

  --bowtie              Use bowtie aligner to align reads, with Bowtie
                        parameters "-p number_of_threads -a -m 100 -S
                        --chunkmbs 1024". If "--paired-end" is set,
                        additionaly enable Bowtie parameters "-I 1 -X 1000".
                        (default: True)
  --bowtie-path <path>  The path to Bowtie executables. (default: None)
  --bowtie2             Use bowtie2 aligner to align reads, indel alignments
                        enabled, with Bowtie2 parameters "-p number_of_threads
                        -k 101". If "--paired-end" is set, additionaly enable
                        Bowtie2 parameters "-I 1 -X 1000 --no-mixed --no-
                        discordant". (default: False)
  --bowtie2-path <path>
                        The path to Bowtie2 executables. (default: None)
  --index-name <name>   The base name for Bowtie/Bowtie2 indices (default:
                        None)
  --keep-alignments     Turn on this option will enable PROBer to keep a copy
                        of aligner-produced alignments in
                        'sample_name.alignments.bam'. (default: False)

OUTPUT:
  sample_name.site_info
    This file contains the expected read counts at each unique crosslink site. Each line describes one site and has the following format:
    
    chr ori pos    n_unique    n_multi
    
    chr is the chromosome name, ori is the orientation (+/-) and pos gives the 0-based genomic coordinate in the '+' strand of chromosome chr. chr, ori, and pos together define the genomic location of the crosslink site and they are separated by single spaces. Then separated by single tabs, n_unique gives the number of uniquely mapped reads, and n_multi provides the expected number of multi-mapping reads at this site.

  sample_name.alignments.bam
    Only generated when '--keep-alignments' option is set.

    This file stores the aligner-produced alignments in BAM format.

  sample_name.stat
    This folder contains learned model parameters from data. In the folder, 'sample_name.model' contains the estimated sequencing model parameters.

  sample_name.temp
    This is a temporary folder contains intermediate files. It will be deleted automatically after the program finishes unless '--keep-intermediate-files' option is on.
~~~

#### version

`PROBer version` displays the current version of the software.

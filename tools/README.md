## Introduction

The executables in this fold are used to interpret PROBer's iCLIP
outputs. You should have one annotation file in GTF format and the
PROBer produced *.site_info in hand.

## Compilation

Type `make`.

## Workflow

Suppose you have the yeast GTF as
`Saccharomyces_cerevisiae.R64-1-1.84.gtf` and PROBer output as
`yeast_PROBer.site_info`.

###### First, you need to parse the GTF using `parseGTF`:

```
Usage: parseGTF input.gtf output_name
```

In the above command, `input.gtf` is your input annotation file and
`output_name` is the file name for parsed results. This program will
create to files: `output_name.anno`, which records the parsed
gene/transcript/exon/intron information, and `output_name.idx`, which
is an index used to fast locate entries in `output_name.anno`.

For our example, type

```
parseGTF Saccharomyces_cerevisiae.R64-1-1.84.gtf yeast
```

and two files `yeast.anno` and `yeast.idx` will be produced.

###### Second, annotate each crosslink site using `annotateSite`:

```
Usage: annotateSite annotation_name input.site_info threshold output.annotated.site_info
```

The first argument is the annotation name you fed into `parseGTF`,
i.e. `yeast`. The second argument is the PROBer produced site_info
file, i.e. `yeast_PROBer.site_info`. The third argument is a
threshold: only crosslink sites with >= threshold total reads (unique
+ multi) will be considered. This argument should be set according to
your sequencing depth. The last argument is the output file.

For example, we can type

```
annotateSite yeast yeast_PROBer.site_info 5 yeast_PROBer.anno.site_info
```

###### Lastly, get gene names associated with each selected crosslink site using `showGeneNames`:

```
Usage: showGeneNames annotation_name input.anno.site_info output.gene.site_info
```

To get the desired output, type

```
showGeneNames yeast yeast_PROBer.anno.site_info yeast_PROBer.gene.site_info
```

In the output file `yeast_PROBer.gene.site_info`, the gene IDs and
names (if available) are shown in the last columns. For one gene, gene
ID is shown first. If GTF file contains gene name, gene name is shown
in the parentheses. If more than one gene associate with the crosslink
site, they are separated by single tabs.






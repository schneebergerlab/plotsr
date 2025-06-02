![python package](https://github.com/schneebergerlab/plotsr/actions/workflows/python-package.yml/badge.svg?event=push)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/plotsr.svg?label=Conda%20downloads)](
https://anaconda.org/bioconda/plotsr)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/plotsr/badges/latest_release_date.svg)](https://anaconda.org/bioconda/plotsr)

## Introduction
Plotsr generates high-quality visualisation of synteny and structural rearrangements between multiple genomes. For this, it uses the genomic structural annotations between multiple chromosome-level assemblies.

![Example](./example/ampril_col0_chr3_6600000_10000000.png)

## Installation
The easiest way to install plotsr is using conda:
```
conda install -c bioconda plotsr 
```
For manual installation, the pre-requisites are:
1. Python >= 3.8
2. Python libraries

These can be installed in a conda environment using:
```
conda install numpy=1.21.2 pandas=1.2.4 matplotlib=3.3.4 setuptools
```
Then download plotsr from github and install:
```
git clone https://github.com/schneebergerlab/plotsr.git
cd plotsr
python -m pip install .
```
After this plotsr should be installed in your conda environment. Test it by printing the help message:
```
plotsr -h
```

## Inputs
### Minimum requirements
The three files below are necessary for plotsr to run:
1. <b>Pairwise structural annotations between genomes (syri.out or syri.bedpe).</b>

For example, if genomes A, B, and C are to be visualised in this order, then structural annotations of A vs B and B vs C genome comparisons would be required (A_B.out and B_C.out).

2. <b>A manually created, tab separated TXT file (genomes.txt)</b>: containing paths to chromosome length files and genome names (e.g. A, B, and C). A third column can be added to customise the visualisation of genomes.

Example:
```
$genomes.txt
#file	name	tags
A.chrlen	A	lw:1.5
B.chrlen	B	lw:1.5
C.chrlen	C	lw:1.5
```

Currently, the following tags are available for genomes:
```
ft = File type (fa/cl for fasta/chromosome_length, default = fa); cl files must be in tsv format with chromosome name in column 1 and chromosome length in column 2; using cl files is much faster than using fasta files
lw = line width
lc = line colour
```
Check the [genomes.txt](./example/genomes.txt) for a working example.

<b><i>NOTE</b>: It is required that the order of the genomes is the same as the order in which genomes are compared. For example, the first genome annotation file uses A as a reference and B as query, and the second genome annotation file uses B as a reference and C as query. Then the genomes.txt file should list the genomes in the order A, B, C.</i>

3. <b>Manually created chromosome length files (A.chrlen, B.chrlen, and C.chrlen)</b>: although these files do not appear in the command line to run plotsr, the paths to these files are entered in the genomes.txt file above.

Example: 
```
$A.chrlen
Chr1	30222096
Chr2	20044946
Chr3	23195252
Chr4	18371560
Chr5	26470424
```

Example files can be found in the [example](./example/) folder. 

### Additional inputs
* <b>Tracks</b>: GFF/BED/bedGraph files to indicate chromosomal regions such as Genes and SNPs, as in the example plot above.
* <b>Markers</b>: BED file containing genomic coordinates to indicate markers such as Inversion 3 and Not aligned 1 in the example plot above.
* <b>Configuration file</b>: Additional parameters (colors, spacing, legends) of the plot can be adjusted by passing a config file to the `--cfg` parameter. Description and default values present in the example [base.cfg](./example/base.cfg) file.   

Example files can be found in the [example](./example/) folder. 

## Command line example
For demonstration, we visualise the structural rearrangements between four accessions of <i>Arabidopsis thaliana</i> (Col-0, L<i>er</i>, Cvi, and Eri). All required files are in the [example](./example/) folder. Below is the list of input files:
| File name|  File Description   |
|----|--------|
| `*.chrlen` | Tables containing chromosome lengths for all genomes |
| `*syri.filtered.out` | Pairwise structural annotation files |
| `genomes.txt` | [Genomes information file](#genomes) |
| `tracks.txt` | [Tracks information file](#visualising-tracks) |
| `markers.bed` | [Markers information file](#visualising-markers) |
| `base.cfg` | Configuration file for adjusting visual properties of the plot |

Below are the commands:
```
# Go to the directory containing all these files.
cd example

# Unzip gene annotation and SNPs files. These would be plotted as tracks.
gzip -d TAIR10_GFF3_genes.gff.gz
gzip -d 1001genomes.snps.sorted.bed.gz

# Running plotsr
plotsr --sr col_lersyri.filtered.out \
       --sr ler_cvisyri.filtered.out \
       --sr cvi_erisyri.filtered.out \
       --genomes genomes.txt \
       --tracks tracks.txt \
       --markers markers.bed \
       --cfg base.cfg \
       -o output_plot.png \
       -S 0.5 -W 7 -H 10 -f 8 
```
These steps create the output_plot.png in the `example` directory.

[plotsr.sh](./example/plotsr.sh) file contains ten different commands corresponding to different modes of visualisation (stacked vs itx mode), different selection of genomic regions (all chromosomes, some chromosomes, or specific region), and different orientation of chromosomes (horizontal vs vertical).

## Full pipeline example: from assemblies to plots
Let's say that we want to visualise genomic differences between four genome assemblies: A.fa, B.fa, C.fa, and D.fa. Further, we want to visualsize the genomes in the order A > B > C > D. Then, we can follow the steps below to visualise the structural rearrangements between these genomes using plotsr:

#### Step 1: Align genomes
* First we align the genomes using a whole-genome alignment tool. Here, we use [minimap2](https://github.com/lh3/minimap2) for aligning and [samtools](https://www.htslib.org/download/) for indexing the output alignment BAM file: 
```
# Align genomes
minimap2 -ax asm5 -t 4 --eqx A.fa B.fa \
 | samtools sort -O BAM - > A_B.bam
samtools index A_B.bam
minimap2 -ax asm5 -t 4 --eqx B.fa C.fa \
 | samtools sort -O BAM - > B_C.bam
samtools index B_C.bam
minimap2 -ax asm5 -t 4 --eqx C.fa D.fa \
 | samtools sort -O BAM - > C_D.bam
samtools index C_D.bam
```

#### Step 2: Find pairwise structural annotations between genomes
* Next we find synteny and structural rearrangements between the genomes using [SyRI](https://github.com/schneebergerlab/syri):
```
# Running syri for finding structural rearrangements between A and B
syri -c A_B.bam -r A.fa -q B.fa -F B --prefix A_B &
# Running syri for finding structural rearrangements between B and C
syri -c B_C.bam -r B.fa -q C.fa -F B --prefix B_C &
# Running syri for finding structural rearrangements between C and D
syri -c C_D.bam -r C.fa -q D.fa -F B --prefix C_D &
```
This will generate A_Bsyri.out, B_Csyri.out, and C_Dsyri.out files that contain the structural annotations between genomes and will be used as input to plotsr.

If other methods are used for finding structural annotations, then their output can be parsed to plotsr using the BEDPE format which should have the following columns:
```
Reference chromosome name
Reference start position
Reference end position
Query chromosome name
Query start position
Query end position
Annotation type
```
Valid values for annotation type: SYN, INV, TRA, INVTR, DUP, INVDP. Here:

| <!-- --> |  <!-- -->   |
|----------|--------|
| SYN      | Syntenic |
| INV      | Inversion |
| TRANS    | Translocation |
| INVTR    | Inverted translocation |
| DUP      | Duplication |
| INVDP    | Inverted duplication |

<b><i>NOTE</b>: The BEDPE file must have syntenic region annotations. These are required to group homologous chromosomes from different genomes. Syntenic regions can only be between homologous chromosomes. In case, syntenic regions between homologous chromosomes are not available, then entire homologous chromosomes can be added as syntenic in the BEDPE file manually to allow clustering of homologous chromosomes by plotsr. While plotting, use the `--nosyn` option to skip plotting of these manually added syntenic regions.  </i>


#### Step 3: Running plotsr
Plotsr can be run using the following command: 
```
plotsr \
    --sr A_Bsyri.out \
    --sr B_Csyri.out \
    --sr C_Dsyri.out \
    --genomes genomes.txt \
    -o output_plot.png
```

## Customising alignment visualisations
Additional column can be added in the input structural annotation files to customise specific alignments. Currently, following tags are available:  
```
cl = colour
lw = line width
z  = vertical location (higher value would plot the alignment over other plot elements)
```

Examples:
```
# Example modified syri.out. Inversions on Chr3 would be black. Inversions on Chr4 would be red and have thick line width 
Chr3	18112802	18114029	-	-	Chr3	18084583	18085805	INV662	-	INV	-	cl:black
Chr3	20464781	20466696	-	-	Chr3	20458463	20460390	INV663	-	INV	-	cl:black
Chr4	1347612	1353808	-	-	Chr4	1437445	1445482	INV664	-	INV	-	cl:red;lw:5;z:4
Chr4	1612606	2782621	-	-	Chr4	1746533	2898561	INV665	-	INV	-	cl:red;lw:5;z:4

# Example BEDPE file
Chr1	1771291	1771585	Chr1	1774045	1774339	INV	cl:black;lw:2;z:4
Chr1	2294260	2296795	Chr1	2297217	2299752	INV	cl:black;lw:2;z:4
Chr1	2455543	2464808	Chr1	2458652	2467917	INV	cl:black;lw:2;z:4
```

<b><i>NOTE</b>: If using alignment customisation, then each row should either have one (or more) of the available tags or have an ```-``` </i>

## Tracks and markers
In addition to structural annotations, plotsr can also be used for visualising tracks for genomics features as well as for marking specific positions in the genomes.

#### Visualising tracks

Feature track information should be in BED or bedGraph format and should correspond to the first genome in visualisation. For example, the [tracks.txt](./example/tracks.txt) contains tracks corresponding to the col-0 genome. Plotsr would then calculate and plot the relative frequency of these features in bins along the chromosomes.
Feature tracks are parsed to plotsr as a tab-separated file containing the path and names for the tracks. The visualisation properties of the tracks can be adjusted by providing a third column containing different tags and corresponding values.

```
$tracks.txt
# file	name	tags
TAIR10_GFF3_genes.gff   Genes   ft:gff;bw:10000;nc:black;ns:8;nf:Arial;lc:blue;lw:4;bc:lightblue;ba:0.5
1001genomes.snps.sorted.bed     SNPs    bw:10000;nc:black;ns:8;nf:Arial;lc:sienna;lw:1;bc:peachpuff;ba:0.5
Giraut2011_centromeres.bed     Centromeres     bw:10000;nc:black;ns:8;nf:Arial;lc:olive;lw:1;bc:palegreen;ba:0.5
```
Currently, the following tags are available for tracks.
```
ft = File type (bed/bedgraph/gff, default = bed)
bw = bin width (default=100000)
nc = name colour
ns = name size
nf = name font
nm = name margin      # Additional margin between name and track. Fraction between [0,1]
lc = line colour
lw = line width
bc = background colour
ba = background alpha
ti = track index      # Numbers starting from 1. Tracks with same index are plotted on top of each other. Tracks with index will be plotted above tracks without index
tt = track type       # f: for plotting a filled plot, l: for plotting a line plot
ta = track alpha      # track transparency. Fraction between [0,1]
```

#### Visualising Markers
Plotsr can mark positions of interest in the genomes. Markers are provided as an extended BED file with five columns: chromosome name, start position, end position, genome name, tags (optional).

```
$markers.bed
#chr	start	end genome_id	tags
Chr3	4035330	4035331	eri	mt:v;mc:black;ms:3;tt:Inversion 1;tp:0.02;ts:8;tf:Arial;tc:black
Chr4	2322547	2322548	ler	mt:^;mc:black;ms:3;tt:Inversion 2;tp:-0.07;ts:8;tf:Arial;tc:black
Chr3	8792851	8792852	col-0	mt:.;mc:red;ms:10;tt:Notal aligned;tp:0.02;ts:8;tf:Arial;tc:black
```
The visualisation properties of the markers can be adjusted by adjusting tag values. Currently, the following tags are available for tracks.
```
mt = marker type
mc = marker colour
ms = marker size
tt = text
tc = text colour
ts = text size
tf = text font
tp = text position
```
Check [markers.txt](./config/marker_point_type.txt) for the list of available markers.

## Citation
If you find plotsr helpful, please [cite](https://doi.org/10.1093/bioinformatics/btac196):

Manish Goel, Korbinian Schneeberger, plotsr: visualizing structural similarities and rearrangements between multiple genomes, Bioinformatics, 2022; btac196, https://doi.org/10.1093/bioinformatics/btac196

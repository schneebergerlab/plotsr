[![Anaconda-Server Badge](https://anaconda.org/bioconda/plotsr/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/plotsr/badges/downloads.svg)](https://anaconda.org/bioconda/plotsr)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/plotsr/badges/latest_release_date.svg)](https://anaconda.org/bioconda/plotsr)

## Introduction
Plotsr generates high-quality visualisation of synteny and structural rearrangements between multiple genomes. For this, it uses the genomic structural annotations between multiple chromosome-level assemblies.


![Example](./example/ampril_col0_chr3_6600000_10000000.png)

## Installation:
The easiest method to install plotsr is using anaconda:
```
conda install -c bioconda plotsr 
```
For manual installation the pre-requisites are:
1. Python >= 3.8
2. Python libraries. These can be installed in a conda environment using:
```
conda install numpy=1.21.2 pandas=1.2.4 matplotlib=3.3.4 setuptools
```
Then download plotsr and install:
```
git clone https://github.com/schneebergerlab/plotsr.git
cd plotsr
python setup.py install
```

After this plotsr should be installed and in your environment. Test it by printing the help message:
```
plotsr -h
```

## Inputs requirements
#### Minimal requirements
1. Chromosome-level assemblies for the genomes to be compared 
2. Pairwise structural annotations between genomes

For example, if genomes A, B, and C are to be visualised in this order, then structural annotations of A vs B and B vs C genome comparisons would be required.

#### Additional inputs
* GFF/BED/bedGraph files for adding tracks to the visualisation, like the tracks for genes and SNPs in the [example](Example) plot above.
* Bed file containing genomic coordinates to add markers, like the markers for Inversion 3, Not aligned 1 in the [example](Example) plot above.

## Example visualisation

As example, we would visualise structural rearrangements between four accessions of <i>Arabidopsis thaliana</i> (Col-0, L<i>er</i>, Cvi, and Eri). All required files are in the [example](./example/) folder. Following is the list of the important input files:
| File name|  File Description   |
|----|--------|
| `*.chrlen` | Table containing chromosome lengths |
| `*syri.filtered.out` | Pairwise structural annotation information between genomes |
| `genomes.txt` | [Genomes information file](#genomes) |
| `tracks.txt` | [Tracks information file](#visualising-tracks) |
| `markers.bed` | [Markers information file](#visualising-markers) |
| `base.cfg` | Configuration file for adjusting visual properties of the plot |

The structural rearrangements between the genomes can be visualised using the following commands:
```
cd example
# Unzip gene annotation and SNPs file. These would be plotted as tracks.
gzip -d TAIR10_GFF3_genes.gff.gz
gzip -d 1001genomes.snps.sorted.bed.gz
# Plot using plotsr
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
This would create the output_plot.png.

[plotsr.sh](./example/plotsr.sh) file contains ten different commands corresponding to different modes of visualisation (stacked vs itx mode), different selection of genomic regions (all chromosomes, some chromosomes, or specific region), and different orientation of chromosomes (horizontal vs stacked).

## Pipeline for visualising genomic differences

Let's say that we want to visualise genomic differences between four genome assemblies: A.fa, B.fa, C.fa, and D.fa. Further, we want to visualsize the genomes in the order A > B > C > D. Then, following steps are involved in visualising structural rearrangements between these genomes using plotsr:

#### Step 1: Align the genomes
* Genomes need to aligned using a whole-genome alignment tool. Here, we align the genomes using [minimap2](https://github.com/lh3/minimap2) and index the alignment BAM file using [samtools](https://www.htslib.org/download/): 
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

#### Step 2: Finding structural annotations between genomes
* Next we need to find synteny and structural rearrangements between the genomes. For this, we use [SyRI](https://github.com/schneebergerlab/syri):
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
|----|--------|
| SYN | Syntenic |
| INV | Inversion |
| TRA | Translocation |
| INVTR | Inverted translocation |
| DUP | Duplication |
| INVDP | Inverted duplication |

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

<a name="genomes">
Here, genomes.txt is a tab-separated file containing the path and names for the genomes. A third column can also be added to customise the visualisation of genomes.
</a>

```
$genomes.txt
#file	name	tags
A.fa	A	lw:1.5
B.fa	B	lw:1.5
C.fa	C	lw:1.5
D.fa	D	lw:1.5
```

Currently, the following tags are available for genomes.

```
ft = File type (fa/cl for fasta/chromosome_length, default = fa); cl files must be in tsv format with chromosome name in column 1 and chromosome length in column 2; using cl files is much faster than using fasta files
lw = line width
lc = line colour
```
Check the [genomes.txt](./example/genomes.txt) for a working example.


<b><i>NOTE</b>: It is required that the order of the genomes is the same as the order in which genomes are compared. For example, if the first genome annotation file uses A as a reference and B as query, and the second genome annotation file uses B as a reference and C as query, then the genomes.txt file should list the genomes in the order A, B, C.</i>

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
lc = line colour
lw = line width
bc = background colour
ba = background alpha
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

## Adjusting other parameters
Additional parameters (colors, spacing, legends) of the plot can be adjusted by parsing a config file to the `--cfg` parameter. Description and default values present in the example [base.cfg](./example/base.cfg) file.   

## Citation:
If you find plotsr helpful, please [cite](https://doi.org/10.1093/bioinformatics/btac196):

`Manish Goel, Korbinian Schneeberger, plotsr: visualizing structural similarities and rearrangements between multiple genomes, Bioinformatics, 2022; btac196, https://doi.org/10.1093/bioinformatics/btac196`
### The package is still in development phase. Please report any issue that you may find. Features requests are also welcome.
## Introduction
Plotsr generates high-quality visualisation of synteny and structural rearrangements between multiple genomes. For this it uses the genomic structural annotations between multiple chromosome-level assemblies.
![Example](./example/ampril_col0_chr3_6600000_10000000.png)

## Installation:
Pre-requisites:
1. Python >= 3.8
2. Python libraries: can be installed in a conda environment using:
```
conda install numpy=1.20.2 pandas=1.2.4 matplotlib=3.3.4
```

Download plotsr:
```
git clone https://github.com/schneebergerlab/plotsr.git
```


## Inputs requirements
#### Minimal requirements
1. Chromosome-level assemblies of the genomes to be compared 
2. Pairwise structural annotations between genomes

For example, if genomes A,B, and C are to be compared, then structural annotations between A and B as well as between B and C would be required.

#### Additional inputs
* BED/bedGraph files for adding tracks to the visualisation, like the tracks for genes and SNPs in the [Example](Example) plot above.
* Bed file containing genomic coordinates to add markers, like the markers for Inv3, Notal1, and Notal2 in the [Example](Example) plot above.

## Running example data
Following are the steps for a typical pipeline to visualise structural annotations between genomes. For this, we would use the data available in the ```example ``` folder.

#### Step 1: Aligning genomes
* In the ```example``` folder, assemblies for three strains of _Arabidopsis thaliana_ are available (Col-0: TAIR10.filtered.fa.gz, L<i>er</i>: ler.filtered.fa.gz, C<i>vi</i>: cvi.filtered.fa.gz).
* We would unzip and then align these genomes using minimap2 (v2.17).
```
# Unzip genomes
gzip -d TAIR10.filtered.fa.gz
gzip -d an1.filtered.fa.gz
gzip -d c24.filtered.fa.gz

# Align genomes
minimap2 -ax asm5 -t 4 --eqx TAIR10.filtered.fa ler.filtered.fa \
 | samtools sort -O BAM - > col_ler.bam
samtools index col_ler.bam
minimap2 -ax asm5 -t 4 --eqx ler.filtered.fa cvi.filtered.fa \
 | samtools sort -O BAM - > ler_cvi.bam
samtools index ler_cvi.bam
```

#### Step 2: Finding structural annotations between genomes
* We use SyRI to get structural annotations between the genomes.
```
# Running syri for finding structural rearrangements between Col-0 and Ler
syri -c col_ler.bam -r TAIR10.filtered.fa -q ler.filtered.fa -F B --prefix col_ler &
# Running syri for finding structural rearrangements between Ler and Cvi
syri -c ler_cvi.bam -r ler.filtered.fa -q cvi.filtered.fa -F B --prefix ler_cvi &
```
This will generate col_lersyri.out and ler_cvisyri.out files that contain the structural annotations between genomes and will be the input to plotsr.

If other methods are used for finding structural annotations, then their output can be parsed to plotsr using the BEDPE format which should have the following coloumns:
```
Reference chromosome name
Reference start position
Reference end position
Query chromosome name
Query start position
Query end position
Annotation type
```
Values for annotation type should one of the following: SYN, INV, TRA, INVTR, DUP, INVDP. Here:

| <!-- --> |  <!-- -->   |
|----|--------|
| SYN | Syntenic |
| INV | Inversion |
| TRA | Translocation |
| INVTR | Inverted translocation |
| DUP | Duplication |
| INVDP | Inverted duplication |



```
SYN = Syntenic

``` 
   
#### Step 3: Running plotsr
Plotsr can be run using the following command: 
```
plotsr.py \
    --sr col_lersyri.out \
    --sr ler_cvisyri.out \
    --genomes genomes.txt \
    -o ampril_horizon.png
```
Here, genomes.txt is a tab-separated file containing path and names for the genomes. A third column can alse be added to customisation visualtion properties of genomes. 
```
$genomes.txt
#file	name	tags
TAIR10.filetered.fa	col-0	lw:1.5
ler.filtered.fa	ler	lw:1.5
cvi.filtered.fa	cvi	lw:1.5
```
Currently, following tags are available for tracks.
```
lw = line width
lc = line colour
```

<b><i>NOTE</b>: It is required that the order of the genomes is same as the order in which genomes are compared. For example, if the first genome annotation file uses GenomeA as reference and GenomeB as query, and the second genome annotation file uses GenomeB as reference and GenomeC as query, then the genomes file should list the genomes in the order GenomeA, GenomeB, GenomeC.</i>

## Tracks and markers
In addition to structural annotations, plotsr can also be used for visualising tracks for genomics features as well as for marking specific positions in the genomes.

#### Visualising tracks
Feature track information should in BED or bedGraph format and should correspond to the first genome in visualisation (here for example: col-0). Plotsr would then calculate and plot the relative frequency of these features in bins along the chromosomes.
Feature tracks are parsed to plotsr as a tab-separated file containing the path and names for the tracks. The visualisation properties of the tracks can be adjusted by providing a third-column containing different tags and corresponding values.
```
$tracks.txt
# file	name	tags
TAIR10_GFF3_genes.bed	Genes	ft:bed;bw:50000;nc:black;ns:8;nf:Arial;lc:black;lw:1;bc:lightgrey;ba:0.5
#TAIR10_GFF3_transposons.bed	TEs	bw:50000;nc:black;ns:8;nf:Arial;lc:black;lw:1;bc:lightgrey;ba:0.5
1001genomes.snps.sorted.bed	SNPs	bw:50000;nc:black;ns:8;nf:Arial;lc:navy;lw:1;bc:aqua;ba:0.5
```
Currently, following tags are available for tracks.
```
bw = bin width (default=100000)
ft = File type (bed/bedgraph, default = bed)
nc = name colour
ns = name size
nf = name font
lc = line colour
lw = line width
bc = background colour
ba = background alpha
```

#### Visualising Markers
Plotsr can mark positions of interest in the genomes. Markers are provided as an extended BED file with five columns: chromosome name, start position, end position, genome name, tags (optional)
```
$markers.bed
#chr	start	end genome_id	tags
Chr3	7354325	7354326	cvi	mt:v;mc:black;ms:3;tt:Inv3;tp:0.02;ts:8;tf:Arial;tc:black
Chr4	4571491	4571492	cvi	mt:v;mc:black;ms:3;tt:Inv1;tp:0.02;ts:8;tf:Arial;tc:black
Chr5	5991438	5991439	c24	mt:^;mc:black;ms:3;tt:Inv2;tp:-0.07;ts:8;tf:Arial;tc:black
Chr3	8792851	8792852	col-0	mt:.;mc:red;ms:10;tt:Notal1;tp:0.02;ts:8;tf:Arial;tc:black
Chr3	8682034	8682035	sha	mt:.;mc:red;ms:10;tt:Notal2;tp:0.02;ts:8;tf:Arial;tc:black
```
The visualisation properties of the markers can be adjusted by adjusting tag values. Currently, following tags are available for tracks.
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
Check [markers.txt](./config/markers.txt) for the list of available markers.


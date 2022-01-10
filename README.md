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
####Minimal requirements
1. Chromosome-level assemblies of the genomes to be compared 
2. Pairwise structural annotations between genomes

For example, if genomes A,B, and C are to be compared, then structural annotations between A and B as well as between B and C would be required.

####Additional inputs
* BED/bedGraph files for adding tracks to the visualisation, like the tracks for genes and SNPs in the [Example](Example) plot above.
* Bed file containing genomic coordinates to add markers, like the markers for Inv3, Notal1, and Notal2 in the [Example](Example) plot above.

## Running example data
Following are the steps for a typical pipeline to visualise structural annotations between genomes. For this, we would use the data available in the ```example ``` folder.

Step 1: Aligning genomes
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

Step 2: Finding structural annotations between genomes
* We use SyRI to get structural annotations between the genomes.
```
# Running syri for finding structural rearrangements between Col-0 and Ler
syri -c col_ler.bam -r TAIR10.filtered.fa -q ler.filtered.fa -F B --prefix col_ler &
# Running syri for finding structural rearrangements between Ler and Cvi
syri -c ler_cvi.bam -r ler.filtered.fa -q cvi.filtered.fa -F B --prefix ler_cvi &
```
This will generate col_lersyri.out and ler_cvisyri.out containing the structural annotations between genomes.

Step 3: Running plotsr
Plotsr can be run using the following command: 
```
plotsr.py \
    --sr col_lersyri.out \
    --sr ler_cvisyri.out \
    --genomes genomes.txt \
    -o ampril_horizon.png
```
Here, genomes.txt is a tab-separated file containig path and names for the genomes.
```
# genomes.txt
TAIR10_Filtered.fasta	col-0
ler.filtered.fa	ler
cvi.filtered.fa	cvi
```


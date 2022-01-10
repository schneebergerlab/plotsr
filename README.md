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

## Inputs requirements
1. Chromosome-level assemblies of the genomes to be compared 
2. Pairwise structural annotations between genomes

For example, if genomes A,B, and C are to be compared, then structural annotations between A and B as well as between B and C would be required.
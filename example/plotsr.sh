#!/bin/bash
## Unzip the annotation (TAIR10_GFF3_genes.gff.gz) and SNPs (1001genomes.snps.sorted.bed.gz) files
gzip -d TAIR10_GFF3_genes.gff.gz
gzip -d 1001genomes.snps.sorted.bed.gz

## Load environment where plotsr is installed
conda activate ENV_WITH_PLOTSR  ## Edit to your environment

## Below are examples of different configuration in which the chromosomes can be seected and visualised
# STACKED
# R1
plotsr --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.5 -o R1 -W 7 -H 10 -f 8 --cfg base.cfg --markers markers.bed
# R2
plotsr --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.5 -o R2 -W 7 -H 10 -f 8 --cfg base.cfg -v --markers markers.bed
# R3
plotsr --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.5 -o R3 -W 7 -H 10 -f 8 --cfg base.cfg --chr Chr1 --chr Chr3 --chr Chr4 --markers markers.bed
# R4
plotsr --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.5 -o R4 -W 7 -H 10 -f 8 --cfg base.cfg -v --markers markers.bed --chr Chr1 --chr Chr3 --chr Chr4
# R5
plotsr --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.5 -o R5 -W 7 -H 10 -f 8 --cfg base.cfg --reg col-0:Chr3:6600000-10000000 --markers markers.bed
# R6
plotsr --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.5 -o R6 -W 7 -H 10 -f 8 --cfg base.cfg -v --markers markers.bed --reg col-0:Chr3:6600000-10000000

# ITX
# R7
plotsr --itx --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.7 -o R7 -W 7 -H 10 -f 8 --cfg base.cfg --markers markers.bed
# R8
plotsr --itx --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.7 -o R8 -W 7 -H 10 -f 8 --cfg base.cfg -v --markers markers.bed
# R9
plotsr --itx --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.7 -o R9 -W 7 -H 10 -f 8 --cfg base.cfg --chr Chr1 --chr Chr3 --chr Chr4 --markers markers.bed
# R10
plotsr --itx --sr col_lersyri.filtered.out --sr ler_cvisyri.filtered.out --sr cvi_erisyri.filtered.out  --genomes genomes.txt --tracks tracks.txt -S 0.7 -o R10 -W 7 -H 10 -f 8 --cfg base.cfg -v --markers markers.bed --chr Chr1 --chr Chr3 --chr Chr4

# Analysing Col genomes
## Get input assemblies from 1001genome page

## Using assemblies in the example folder
gzip -d TAIR10.filtered.fa.gz
gzip -d cvi.filtered.fa.gz
gzip -d ler.filtered.fa.gz

## Run align genomes
minimap2 -ax asm5 -t 4 --eqx TAIR10.filtered.fa ler.filtered.fa \
 | samtools sort -O BAM - > col_ler.bam
minimap2 -ax asm5 -t 4 --eqx ler.filtered.fa cvi.filtered.fa \
 | samtools sort -O BAM - > ler_cvi.bam
samtools index col_ler.bam
samtools index ler_cvi.bam

## Run syri (V5)
syri -c col_ler.bam -r TAIR10.filtered.fa -q ler.filtered.fa -F B --prefix col_ler &
syri -c ler_cvi.bam -r ler.filtered.fa -q cvi.filtered.fa -F B --prefix ler_cvi &

## Get tracks
## SNPs and Indels
### Downloaded 1001 genome VCF from https://1001genomes.org/data/GMI-MPI/releases/v3.1/
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
### Split VCF into SNPs and indels
vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip > 1001genomes.snps.vcf.gz &
vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz --keep-only-indels  --recode --recode-INFO-all --stdout | gzip >  1001genomes.indels.vcf.gz &

zcat 1001genomes.snps.vcf.gz | vcf2bed --do-not-sort < /dev/stdin | cut -f1,2,3,6,7 | sed 's/^/Chr/g' > 1001genomes.snps.sorted.bed
zcat 1001genomes.indels.vcf.gz | vcf2bed < /dev/stdin | cut -f1,2,3,6,7 | sed 's/^/Chr/g' > 1001genomes.indels.sorted.bed &

### Genes and TE files are available in the example folder
#cat /srv/biodata/dep_mercier/grp_schneeberger/example/Athal/TAIR10/TAIR10_geneCoords.txt | awk '{print $1"\t"$2-1"\t"$3}'> TAIR10_GFF3_genes.bed
#cut -f1,4,5 /srv/biodata/dep_mercier/grp_schneeberger/example/Athal/TAIR10/TAIR10_GFF3_genes_transposons_onlyTE.gff | awk '{print $1"\t"$2-1"\t"$3}'> TAIR10_GFF3_transposons.bed

## Run plotsr
plotsr \
    --sr col_lersyri.out \
    --sr ler_cvisyri.out \
    --sr cvi_erisyri.out \
    --genomes genomes.txt \
    --tracks tracks.txt \
    --markers markers.bed \
    -S 0.65 \
    -o ampril_horizon.pdf \
    -W 10 \
    -H 14 \
    -f 8 \
    --cfg base.cfg

plotsr \
    --sr col_lersyri.out \
    --sr ler_cvisyri.out \
    --sr cvi_erisyri.out \
    --genomes genomes.txt \
    --tracks tracks.txt \
    --markers markers.bed \
    -S 0.65 \
    -o ampril_horizon.png \
    -W 10 \
    -H 14 \
    -f 8 \
    --cfg base.cfg

plotsr \
    --sr col_lersyri.out \
    --sr ler_cvisyri.out \
    --sr cvi_erisyri.out \
    --genomes genomes.txt \
    --tracks tracks.txt \
    --markers markers.bed \
    --reg col-0:Chr3:6600000-10000000 \
    -S 0.65 \
    -o ampril_col0_chr3_6600000_10000000.pdf \
    -W 10 \
    -H 6 \
    -f 8 \
    --cfg /biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/plotsr/config/base.cfg

plotsr \
    --sr col_lersyri.out \
    --sr ler_cvisyri.out \
    --sr cvi_erisyri.out \
    --genomes genomes.txt \
    --tracks tracks.txt \
    --markers markers.bed \
    --reg col-0:Chr3:6600000-10000000 \
    -S 0.65 \
    -o ampril_col0_chr3_6600000_10000000.png \
    -W 10 \
    -H 6 \
    -f 8 \
    --cfg base.cfg

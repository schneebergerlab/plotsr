# Analysing Col genomes
## Get input assemblies from 1001genome page
### Select chromosomes
hometools getchr --chrs chr1 chr2 chr3 chr4 chr5 -o an1.filtered.fa An-1.chr.all.v2.0.fasta.gz &
hometools getchr --chrs chr1 chr2 chr3 chr4 chr5 -o c24.filtered.fa C24.chr.all.v2.0.fasta.gz &
hometools getchr --chrs chr1 chr2 chr3 chr4 chr5 -o cvi.filtered.fa Cvi.chr.all.v2.0.fasta.gz &
hometools getchr --chrs chr1 chr2 chr3 chr4 chr5 -o eri.filtered.fa Eri.chr.all.v2.0.fasta.gz &
hometools getchr --chrs chr1 chr2 chr3 chr4 chr5 -o kyo.filtered.fa Kyo.chr.all.v2.0.fasta.gz &
hometools getchr --chrs chr1 chr2 chr3 chr4 chr5 -o ler.filtered.fa Ler.chr.all.v2.0.fasta.gz &
hometools getchr --chrs chr1 chr2 chr3 chr4 chr5 -o sha.filtered.fa Sha.chr.all.v2.0.fasta.gz &

sed -i 's/>chr/>Chr/g' an1.filtered.fa
sed -i 's/>chr/>Chr/g' c24.filtered.fa
sed -i 's/>chr/>Chr/g' cvi.filtered.fa
sed -i 's/>chr/>Chr/g' eri.filtered.fa
sed -i 's/>chr/>Chr/g' kyo.filtered.fa
sed -i 's/>chr/>Chr/g' ler.filtered.fa
sed -i 's/>chr/>Chr/g' sha.filtered.fa


## Run align genomes
minimap2 -ax asm5 -t 4 --eqx TAIR10_Filtered.fasta ler.filtered.fa \
 | samtools sort -O BAM - > col_ler.bam
minimap2 -ax asm5 -t 4 --eqx ler.filtered.fa cvi.filtered.fa \
 | samtools sort -O BAM - > ler_cvi.bam
minimap2 -ax asm5 -t 4 --eqx cvi.filtered.fa eri.filtered.fa  \
 | samtools sort -O BAM - > cvi_eri.bam &
minimap2 -ax asm5 -t 4 --eqx eri.filtered.fa sha.filtered.fa  \
 | samtools sort -O BAM - > eri_sha.bam &
minimap2 -ax asm5 -t 4 --eqx sha.filtered.fa kyo.filtered.fa  \
 | samtools sort -O BAM - > sha_kyo.bam &
minimap2 -ax asm5 -t 4 --eqx kyo.filtered.fa an1.filtered.fa  \
 | samtools sort -O BAM - > kyo_an1.bam &
minimap2 -ax asm5 -t 4 --eqx an1.filtered.fa c24.filtered.fa  \
 | samtools sort -O BAM - > an1_c24.bam &
samtools index col_ler.bam
samtools index ler_cvi.bam
samtools index cvi_eri.bam
samtools index eri_sha.bam
samtools index sha_kyo.bam
samtools index kyo_an1.bam
samtools index an1_c24.bam

## Run syri
### Run with syri3.8
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c col_ler.bam -r TAIR10_Filtered.fasta -q ler.filtered.fa -F B --prefix col_ler &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c ler_cvi.bam -r ler.filtered.fa -q cvi.filtered.fa -F B --prefix ler_cvi &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c cvi_eri.bam -r cvi.filtered.fa -q eri.filtered.fa -F B --prefix cvi_eri &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c eri_sha.bam -r eri.filtered.fa -q sha.filtered.fa -F B --prefix eri_sha &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c sha_kyo.bam -r sha.filtered.fa -q kyo.filtered.fa -F B --prefix sha_kyo &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c kyo_an1.bam -r kyo.filtered.fa -q an1.filtered.fa -F B --prefix kyo_an1 &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c an1_c24.bam -r an1.filtered.fa -q c24.filtered.fa -F B --prefix an1_c24 &

## Get tracks
## SNPs and Indels
### Downloaded 1001 genome VCF from https://1001genomes.org/data/GMI-MPI/releases/v3.1/
### Split VCF into SNPs and indels
vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip > 1001genomes.snps.vcf.gz &
vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz --keep-only-indels  --recode --recode-INFO-all --stdout | gzip >  1001genomes.indels.vcf.gz &

zcat 1001genomes.snps.vcf.gz | vcf2bed --do-not-sort < /dev/stdin | cut -f1,2,3,6,7 | sed 's/^/Chr/g' > 1001genomes.snps.sorted.bed
zcat 1001genomes.indels.vcf.gz | vcf2bed < /dev/stdin | cut -f1,2,3,6,7 | sed 's/^/Chr/g' > 1001genomes.indels.sorted.bed &

### Genes and TE
cat /srv/biodata/dep_mercier/grp_schneeberger/data/Athal/TAIR10/TAIR10_geneCoords.txt | awk '{print $1"\t"$2-1"\t"$3}'> TAIR10_GFF3_genes.bed
cut -f1,4,5 /srv/biodata/dep_mercier/grp_schneeberger/data/Athal/TAIR10/TAIR10_GFF3_genes_transposons_onlyTE.gff | awk '{print $1"\t"$2-1"\t"$3}'> TAIR10_GFF3_transposons.bed


## Run plotsr
/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/plotsr/plotsr/plotsr.py \
    --sr col_lersyri.out \
    --sr ler_cvisyri.out \
    --sr cvi_erisyri.out \
    --sr eri_shasyri.out \
    --sr sha_kyosyri.out \
    --sr kyo_an1syri.out \
    --sr an1_c24syri.out \
    --genomes genomes.txt \
    --tracks tracks.txt \
    --markers markers.bed \
    -S 0.6 \
    -o ampril_horizon.pdf \
    -W 10 \
    -H 12 \
    -f 8 \
    -R
/srv/biodata/dep_mercier/grp_schneeberger/data/Athal/TAIR10/TAIR10_GFF3_genes_transposons.gff

/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/plotsr/plotsr/plotsr.py \
    --sr col_lersyri.out \
    --sr ler_cvisyri.out \
    --sr cvi_erisyri.out \
    --sr eri_shasyri.out \
    --sr sha_kyosyri.out \
    --sr kyo_an1syri.out \
    --genomes genomes.txt \
    --tracks tracks.txt \
    --markers test_markers.bed \
    --reg col-0:Chr5:10000000-17000000 \
    -S 0.65 \
    -o ampril_col0_chr_10000000_17000000.pdf \
    -W 10 \
    -H 12 \
    -f 10 \
    -R

# Analysing potato genomes
CWD=/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/plotsr_align/otava/
gunzip -d ./*gz
for g in He1 He2 St1 St2; do
    hometools getchr --chrs Chr{1..12} -o ${g}_haplotype_genome.filtered.fasta ${g}_haplotype_genome.fasta
done
## Run align genomes
minimap2 -ax asm5 -t 4 --eqx He1_haplotype_genome.filtered.fasta He2_haplotype_genome.filtered.fasta \
 | samtools sort -O BAM - > he1_he2.bam &
minimap2 -ax asm5 -t 4 --eqx He2_haplotype_genome.filtered.fasta St1_haplotype_genome.filtered.fasta \
 | samtools sort -O BAM - > he2_st1.bam &
minimap2 -ax asm5 -t 4 --eqx St1_haplotype_genome.filtered.fasta St2_haplotype_genome.filtered.fasta  \
 | samtools sort -O BAM - > st1_st2.bam &
samtools index he1_he2.bam
samtools index he2_st1.bam
samtools index st1_st2.bam

## Run syri
### Run with syri3.8
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c he1_he2.bam -r He1_haplotype_genome.filtered.fasta -q He2_haplotype_genome.filtered.fasta -k -F B --prefix he1_he2 --nc 4 &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c he2_st1.bam -r He2_haplotype_genome.filtered.fasta -q St1_haplotype_genome.filtered.fasta -k -F B --prefix he2_st1 --nc 4 &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c st1_st2.bam -r St1_haplotype_genome.filtered.fasta -q St2_haplotype_genome.filtered.fasta -k -F B --prefix st1_st2 --nc 4 &

## Run plotsr
/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/plotsr/plotsr/plotsr.py \
    --sr he1_he2syri.out \
    --sr he2_st1syri.out \
    --sr st1_st2syri.out \
    --genomes genomes.txt \
    --tracks tracks.txt \
    --chr Chr1 \
    -S 0.85 \
    -o png \
    -v \
    -W 6 \
    -H 10 \
    -f 10 \
    -R

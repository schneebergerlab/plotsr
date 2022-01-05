# Analysing Col genomes
## Get input genomes
cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/plotsr_align/
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/TAIR10_Filtered.fasta .
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/ampril/v2/ler.filtered.fa .
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/ampril/v2/cvi.filtered.fa .
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/ampril/v2/eri.filtered.fa .
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/ampril/v2/sha.filtered.fa .
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/ampril/v2/kyo.filtered.fa .
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/ampril/v2/an1.filtered.fa .


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
samtools index col_ler.bam
samtools index ler_cvi.bam
samtools index cvi_eri.bam
samtools index eri_sha.bam
samtools index sha_kyo.bam
samtools index kyo_an1.bam

## Run syri
### Run with syri3.8
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c col_ler.bam -r TAIR10_Filtered.fasta -q ler.filtered.fa -k -F B --prefix col_ler &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c ler_cvi.bam -r ler.filtered.fa -q cvi.filtered.fa -k -F B --prefix ler_cvi &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c cvi_eri.bam -r cvi.filtered.fa -q eri.filtered.fa -k -F B --prefix cvi_eri &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c eri_sha.bam -r eri.filtered.fa -q sha.filtered.fa -k -F B --prefix eri_sha &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c sha_kyo.bam -r sha.filtered.fa -q kyo.filtered.fa -k -F B --prefix sha_kyo &
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c kyo_an1.bam -r kyo.filtered.fa -q an1.filtered.fa -k -F B --prefix kyo_an1 &

## Run plotsr
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
    -S 0.65 \
    -o ampril_horizon.pdf \
    -W 10 \
    -H 12 \
    -f 10 \
    -R

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
    --reg cvi:LR699761.1:12000000-13500000 \
    -S 0.65 \
    -o pdf \
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

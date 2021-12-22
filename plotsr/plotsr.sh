# Get input genomes
cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/results/plotsr_align/
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/TAIR10_Filtered.fasta .
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/ampril/v2/ler.filtered.fa .
cp /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/data/Arabidopsis-thaliana/ampril/v2/cvi.filtered.fa .

# Run align genomes
minimap2 -ax asm5 -t 4 --eqx TAIR10_Filtered.fasta ler.filtered.fa \
 | samtools sort -O BAM - > col_ler.bam
minimap2 -ax asm5 -t 4 --eqx ler.filtered.fa cvi.filtered.fa \
 | samtools sort -O BAM - > ler_cvi.bam
minimap2 -ax asm5 -t 4 --eqx cvi.filtered.fa TAIR10_Filtered.fasta  \
 | samtools sort -O BAM - > cvi_col.bam
samtools index col_ler.bam
samtools index ler_cvi.bam
samtools index cvi_col.bam

# Run syri
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c col_ler.bam -r TAIR10_Filtered.fasta -q ler.filtered.fa -k -F B --prefix col_ler
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c ler_cvi.bam -r ler.filtered.fa -q cvi.filtered.fa -k -F B --prefix ler_cvi
/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri -c cvi_col.bam -r cvi.filtered.fa -q TAIR10_Filtered.fasta -k -F B --prefix cvi_col


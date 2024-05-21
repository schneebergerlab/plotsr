# Converting syri.out to bedpe
# grep -P 'SYN\t|INV\t|TRA\t|INVTR\t|DUP\t|INVDP\t'  col_lersyri.out | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$11}' > col_ler.syri.bedpe

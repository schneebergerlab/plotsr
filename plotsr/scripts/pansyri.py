import pansyri.util as util
from pansyri.ordering import order_greedy
import os
os.chdir('/srv/netscratch/dep_mercier/grp_schneeberger/projects/pansr/data/ampril')

syns, alns = util.parse_input_tsv('pansr.tsv')
df = util.crosssyn_from_lists(syns, alns, SYNAL=False)
df2 = util.filter_multisyn_df_chr(df, "Chr3")
order = order_greedy(df2)
print(order)


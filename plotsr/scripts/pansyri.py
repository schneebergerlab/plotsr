import pansyri.util as util
from pansyri.ordering import order_greedy
import os
os.chdir('/srv/netscratch/dep_mercier/grp_schneeberger/projects/pansr/data/ampril')

syns, alns = util.parse_input_tsv('full.tsv')
df = util.crosssyn_from_lists(syns, alns, SYNAL=False)
df2 = util.filter_multisyn_df_chr(df, "Chr3")
order = order_greedy(df2)
print(order)

row = df.iloc[1][0]
row.ref.__len__()
Out[76]: 812
row.get_lens()
Out[77]: {'c24': 481, 'ler': 4737, 'sha': 1357}


synlen = defaultdict(int)
for row in df.itertuples(index=False):
    synlen[row[0].get_degree()] += row[0].ref.__len__()

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
collen = 119146348
x = np.array(list(range(1, 8)))
y = np.array([synlen[i] for i in range(1, 8)])
cm = np.cumsum(y)

# 300 represents number of points to make between x.min and x.max
xnew = np.linspace(x.min(), x.max(), 300)
spl = make_interp_spline(x, cm/collen, k=5)  # type: BSpline
cm_smooth = spl(xnew)
fig = plt.figure(figsize=[4, 3])
ax = fig.add_subplot()
ax.bar(x[:-1], y[:-1]/collen, zorder=1, color='#2da399', label='cross-synteny')
ax.bar(x[-1], y[-1]/collen, zorder=1, color='#006c66', label='core-synteny')
ax.plot(xnew, cm_smooth, zorder=2, color='#c6d325', lw=6, label='pansynteny')
ax.grid(which='major', axis='y', ls='--', zorder=0, lw=1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_axisbelow(True)
ax.set_xticks(range(1, 8))
ax.set_xlabel('Number of syntenic genomes')
ax.set_ylabel('Ratio of reference genome')
ax.legend()
plt.tight_layout()
plt.savefig("/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/pansyn_regions.png", dpi=300)
plt.close()
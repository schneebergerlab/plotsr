#!/usr/bin/env python3

## Constants
MARKERS = {".": "point",
           ",": "pixel",
           "o": "circle",
           "v": "triangle_down",
           "^": "triangle_up",
           "<": "triangle_left",
           ">": "triangle_right",
           "1": "tri_down",
           "2": "tri_up",
           "3": "tri_left",
           "4": "tri_right",
           "8": "octagon",
           "s": "square",
           "p": "pentagon",
           "P": "plus (filled)",
           "*": "star",
           "h": "hexagon1",
           "H": "hexagon2",
           "+": "plus",
           "x": "x",
           "X": "x (filled)",
           "D": "diamond",
           "d": "thin_diamond",
           "|": "vline",
           "_": "hline",
           "i0": "tickleft",
           "i1": "tickright",
           "i2": "tickup",
           "i3": "tickdown",
           "i4": "caretleft",
           "i5": "caretright",
           "i6": "caretup",
           "i7": "caretdown",
           "i8": "caretleft (centered at base)",
           "i9": "caretright (centered at base)",
           "i10": "caretup (centered at base)",
           "i11": "caretdown"}
VARS = ['SYN', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']
COLORS = ['#DEDEDE', '#FFA500', '#9ACD32', '#9ACD32', '#00BBFF', '#00BBFF', '#83AAFF', '#FF6A33']


"""
DEFINE READERS/PARSERS
"""
def readfasta(f):
    from gzip import open as gzopen
    from gzip import BadGzipFile
    from collections import deque
    import sys
    out = {}
    chrid = ''
    chrseq = deque()
    # Test if the file is Gzipped or not
    with gzopen(f, 'rb') as fin:
        try:
            fin.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzopen(f, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                            chrseq = deque()
                        else:
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                        if chrid in out.keys():
                            sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                    else:
                        chrseq.append(line.strip().decode())
        else:
            with open(f, 'r') as fin:
                for line in fin:
                    if '>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split('>')[1].split(' ')[0]
                            chrseq = deque()
                        else:
                            chrid = line.strip().split('>')[1].split(' ')[0]
                        if chrid in out.keys():
                            sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                    else:
                        chrseq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    if chrid != '':
        out[chrid] = ''.join(chrseq)
    # TODO: add check for the validation of input fasta files
    return out
# END


def readbed(path, v):
    from collections import deque
    mdata = deque()
    with open(path, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if len(line) not in [7, 10]:
                print(line)
                print("Incomplete data in BED file line: \n {}.\n Marker information (type, col, size) is necessary. Text annotation (text, col, size) is optional".format('\t'.join(line)))
                continue
            anno = bedanno(line[0], line[1], line[2], line[3], v)
            anno.setmarker(line[4], line[5], line[6])
            if len(line) > 7:
                anno.settext(line[7], line[8], line[9])
            mdata.append(anno)
    return mdata
# END


def readsyriout(f):
    from pandas import DataFrame
    import numpy as np
    from collections import deque, OrderedDict
    # Reads syri.out. Select: achr, astart, aend, bchr, bstart, bend, srtype
    syri_regs = deque()
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            # TODO: DECIDE WHETHER TO HAVE STATIC VARS OR FLEXIBLE ANNOTATION
            if l[10] in VARS:
                syri_regs.append(l)
    df = DataFrame(list(syri_regs))[[0, 1, 2, 5, 6, 7, 10]]
    df[[0, 5, 10]] = df[[0, 5, 10]].astype(str)
    df[[1, 2, 6, 7]] = df[[1, 2, 6, 7]].astype(int)
    # chr ID map
    chrid = []
    chrid_dict = OrderedDict()
    for i in np.unique(df[0]):
        chrid.append((i, np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]))
        chrid_dict[i] = np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]
    df.columns = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type']
    return df, chrid_dict
# END


def readbedout(f):
    from pandas import DataFrame
    import numpy as np
    from collections import deque, OrderedDict
    # BEDPE format: achr, astart, aend, bchr, bstart, bend, srtype
    bed_regs = deque()
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            # TODO: DECIDE WHETHER TO HAVE STATIC VARS OR FLEXIBLE ANNOTATION
            if l[6] in VARS:
                bed_regs.append(l)
    df = DataFrame(list(bed_regs))
    df[[0, 3, 6]] = df[[0, 3, 6]].astype(str)
    df[[1, 2, 4, 5]] = df[[1, 2, 4, 5]].astype(int)
    df[[1, 4]] = df[[1, 4]].astype(int) # Makes range closed as the terminal bases are also included
    # chr ID map
    chrid = []
    chrid_dict = OrderedDict()
    for i in np.unique(df[0]):
        chrid.append((i, np.unique(df.loc[(df[0] == i) & (df[6] == 'SYN'), 3])[0]))
        chrid_dict[i] = np.unique(df.loc[(df[0] == i) & (df[6] == 'SYN'), 3])[0]
    df.columns = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type']
    return df, chrid_dict
# END

"""
Validation and filtering
"""
def validalign2fasta(als, genf):
    """
    Check that the chromosome ID and length in the alignment file matches the
    genome fasta
    :param als: deque; each element is a dataframe containing the alignments to be plotted
    :param genf: path to file containing "genome_IDs:path_to_genome_fasta"
    :return: genome length dict.
    """
    from collections import OrderedDict
    import sys
    import numpy as np
    import os
    out = deque()
    with open(genf, 'r') as fin:
        i = 0
        for line in fin:
            line = line.strip().split(":")
            if len(line) != 2:
                raise SystemExit("Incorrect format for the genome path file.\nExpected format:\ngenome1_id:path_to_genome1\ngenome2_id:path_to_genome2")
                sys.exit()
            try:
                # length of individual chromosomes
                glen = {id: len(seq) for id, seq in readfasta(line[1]).items()}
            except Exception as e:
                raise SystemExit("Error in reading fasta: {}\n{}".format(line[1], e))
            # Check cases when the genome is the query-genome
            if i > 0:
                df = als[i-1][1]
                chrs = np.unique(df['bchr'])
                for c in chrs:
                    if c not in list(glen.keys()):
                        raise SystemExit('Chromosome ID: {} in structural annotation file: {} not present in genome fasta: {}. Exiting.'.format(c, als[i-1][0], os.path.basename(line[1])))
                    if np.max(np.max(df.loc[df['bchr'] == c, ['bstart', 'bend']])) > glen[c]:
                        raise SystemExit('For chromosome ID: {}, length in genome fasta: {} is less than the maximum coordinate in the structural annotation file: {}. Exiting.'.format(c, os.path.basename(fin), als[i-1][0]))
            # Check cases when the genome is the reference-genome
            if i < len(als):
                df = als[i][1]
                chrs = np.unique(df['achr'])
                for c in chrs:
                    if c not in list(glen.keys()):
                        raise SystemExit('Chromosome ID: {} in structural annotation file: {} not present in genome fasta: {}. Exiting.'.format(c, als[i][0], os.path.basename(line[1])))
                    if np.max(np.max(df.loc[df['achr'] == c, ['astart', 'aend']])) > glen[c]:
                        raise SystemExit('For chromosome ID: {}, length in genome fasta: {} is less than the maximum coordinate in the structural annotation file: {}. Exiting.'.format(c, os.path.basename(fin), als[i][0]))
            out.append((line[0], glen))
            i += 1
    return out
# END


def filterinput(args, df, chrid):
    # Get region length and filter out smaller SR
    df = df.loc[((df['aend'] - df['astart']) >= args.s) | ((df['bend'] - df['bstart']) >= args.s) | (df['type']=='SYN')]
    df = df.loc[df['bchr'] == [chrid[i] for i in df['achr']]]
    # Filter non-selected variations
    if args.nosyn:
        df = df.loc[df[10] != 'SYN']
    if args.noinv:
        df = df.loc[df[10] != 'INV']
    if args.notr:
        df = df.loc[~df[10].isin(['TRANS', 'INVTR'])]
    if args.nodup:
        df = df.loc[~df[10].isin(['DUP', 'INVDP'])]
    df.sort_values(['bchr', 'bstart', 'bend'], inplace=True)
    df.sort_values(['achr', 'astart', 'aend'], inplace=True)
    return df

def selectchrs(als, chrids, clens. chrs):
    """
    # TODO: Given a list of reference chromsome IDs, filter alignments to select
    alignments between chromosomes homologous to given reference chromosomes
    """
    # df = ''
    # from warnings import warn
    #
    # if args.chr is not None:
    #     for chr in args.chr:
    #         if chr not in list(reflenghts.keys()):
    #             warn('Selected chromosome ID: {} is not present in the reference genome. Only use reference chromsome ID for selecting chromsomes.'.format(chr))
    #     df = df.loc[df[0].isin(args.chr)]
    pass
# END


def createribbon(df):
    """
    Combine continuous syntenic regions to get larger ribbons for syntenic blocks
    """
    import numpy as np
    from pandas import DataFrame
    df.sort_values(['bchr', 'bstart', 'bend'], inplace=True)
    df['b'] = list(range(df.shape[0]))
    df.sort_values(['achr', 'astart', 'aend'], inplace=True)
    df['a'] = list(range(df.shape[0]))
    # Group syntenic regions
    groups = deque()
    cg = deque()
    ca = -1
    cb = -1
    issyn = list(df['type'] == 'SYN')
    a = list(df['a'])
    b = list(df['b'])
    chrchange = np.where((np.array(df['achr'][1:]) == np.array(df['achr'][0:-1])) == False)[0] + 1
    for i in range(df.shape[0]):
        if i in chrchange:
            if len(cg) > 0:
                groups.append(list(cg))
            cg = deque()
            if issyn[i]:
                cg.append(i)
                ca = i
                cb = i
        elif not issyn[i]:
            if len(cg) > 0:
                groups.append(list(cg))
            cg = deque()
        else:
            if len(cg) == 0:
                cg.append(i)
                ca = i
                cb = i
            elif i == ca+1 and i == cb+1:
                cg.append(i)
                ca = i
                cb = i
            else:
                if len(cg)>0:
                    groups.append(list(cg))
                cg = deque()
                cg.append(i)
                ca = i
                cb = i
    if len(cg)>0:
        groups.append(list(cg))
    newsyn = deque()
    for i in groups:
        tmpdf = df.iloc[i]
        newsyn.append([
            list(tmpdf['achr'])[0],
            list(tmpdf['astart'])[0],
            list(tmpdf['aend'])[-1],
            list(tmpdf['bchr'])[0],
            list(tmpdf['bstart'])[0],
            list(tmpdf['bend'])[-1]
        ])
    newsyn = DataFrame(list(newsyn), columns=['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend'])
    newsyn['type'] = 'SYN'

    df = df.drop(['a', 'b'], axis=1)
    df = df.loc[-(df['type'] == 'SYN')]
    df = df.append(newsyn)
    df.sort_values(['achr', 'astart', 'aend'], inplace=True)
    return df
# END


"""
Draw and plot
"""
def drawax(ax, chrgrps, chrlengths, V, S):
    import numpy as np
    nchr = len(chrgrps)
    # qchrs = [chrid_dict[k] for k in chrs]
    tick_pos = 1 - (S/2)
    ticklabels = list(chrgrps.keys())
    max_l = np.max([chrlengths[i][1][v[i]] for v in chrgrps.values() for i in range(len(v))])
    if not V:
        ax.set_ylim(0, nchr+0.2)
        ax.set_yticks([tick_pos+i for i in range(nchr)])
        ax.set_yticklabels(ticklabels[::-1])
        ax.tick_params(axis='y', right=False, left=False)
        ax.set_xlim(0, max_l)
        ax.ticklabel_format(axis='x', useOffset=False, style='plain')
        xticks = ax.get_xticks()
        if max_l >= 1000000000:
            xticks = xticks/1000000000
            ax.set_xlabel('chromosome position (in Gbp)')
        elif max_l >= 1000000:
            xticks = xticks/1000000
            ax.set_xlabel('chromosome position (in Mbp)')
        elif max_l >= 1000:
            xticks = xticks/1000
            ax.set_xlabel('chromosome position (in Kbp)')
        ax.set_xticklabels(xticks)
        ax.set_ylabel('reference chromosome id')
        ax.xaxis.grid(True, which='major', linestyle='--')
        ax.set_axisbelow(True)
    else:
        ax.set_xlim(0, nchr+0.2)
        ax.set_xticks([tick_pos+i for i in range(nchr)])
        ax.set_xticklabels(ticklabels)
        ax.tick_params(axis='x', top=False, bottom=False)
        ax.set_ylim(0, max_l)
        ax.ticklabel_format(axis='y', useOffset=False, style='plain')
        yticks = ax.get_yticks()
        if max_l >= 1000000000:
            yticks = yticks/1000000000
            ax.set_ylabel('chromosome position (in Gbp)')
        elif max_l >= 1000000:
            yticks = yticks/1000000
            ax.set_ylabel('chromosome position (in Mbp)')
        elif max_l >= 1000:
            yticks = yticks/1000
            ax.set_ylabel('chromosome position (in Kbp)')
        ax.set_yticklabels(yticks)
        ax.set_xlabel('reference chromosome id')
        ax.yaxis.grid(True, which='major', linestyle='--')
        ax.set_axisbelow(True)
    return ax, max_l
# END


def pltchrom(ax, chrs, chrgrps, chrlengths, V, S):
    import warnings
    chrlabs = [False]*len(chrlengths)
    max_l = np.max([chrlengths[i][1][v[i]] for v in chrgrps.values() for i in range(len(v))])
    CHRCOLS = plt.get_cmap('Dark2')                 # TODO: READ COLORS FROM CONFIG FILE
    if len(chrlengths) > 8:
        warnings.warn("More than 8 chromosomes are being analysed. This could result in different chromosomes having same color. Provide colors manually in config.")
    # Plot chromosomes
    pltchr = ax.axhline if not V else ax.axvline
    step = S/(len(chrlengths)-1)
    if not V:
        rend = len(chrs)-0.02
        indents = [rend - (i*step) for i in range(len(chrlengths))]
    elif V:
        rend = 1-S-0.02
        indents = [rend + (i*step) for i in range(len(chrlengths))]
    for s in range(len(chrlengths)):
        for i in range(len(chrs)):
            offset = i if not V else -i
            if not chrlabs[s]:
                pltchr(indents[s]-offset, 0, chrlengths[s][1][chrgrps[chrs[i]][s]]/max_l, color=CHRCOLS(s), linewidth=3, label=chrlengths[s][0])
                chrlabs[s] = True
            else:
                pltchr(indents[s]-offset, 0, chrlengths[s][1][chrgrps[chrs[i]][s]]/max_l, color=CHRCOLS(s), linewidth=3)
    return ax
#END

class bedanno():
    def __init__(self, chr, start, end, genome, V):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        if genome in ['R', 'Q']: self.genome = genome
        else: print("Incorrect genome for position {}:{}-{}. Setting marker to line".format(self.start, self.start, self.end))
        self.V = V

    def setmarker(self, mtype, mcol, msize):
        if mtype not in MARKERS.keys():
            print("Unrecongised marker used for position {}:{}-{}. Plotsr accepts markers defined in matplotlib (https://matplotlib.org/stable/api/markers_api.html) with some modifications.".format(self.start, self.start, self.end))
            for k, v in MARKERS.items():
                print("{} : {}".format(k, v))
            sys.exit()
        self.mtype = mtype if mtype[0] != i else int(mtype[1:])

        try:
            if mcol[0] == '#':
                matplotlib.colors.to_rgb(mcol)
            else:
                matplotlib.colors.to_hex(mcol)
        except ValueError as e:
            print("Error in using colour: {} for position {}:{}-{}. Use correct hexadecimal colours or named colours define in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html)".format(mcol, self.start, self.start, self.end))
            sys.exit()
        self.mcol = mcol
        self.msize = int(msize)
        if self.end - self.start > 1:
            print("Range selected for position {}:{}-{}. Setting marker to line".format(self.start, self.start, self.end))
            self.mtype = "|" if self.V else "_"

    def settext(self, text, col, size):
        self.text = text
        try:
            if col[0] == '#':
                matplotlib.colors.to_rgb(col)
            else:
                matplotlib.colors.to_hex(col)
        except ValueError as e:
            print("Error in using colour: {} for position {}:{}-{}. Use correct hexadecimal colours or named colours define in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html)".format(col, self.start, self.start, self.end))
            sys.exit()
        self.col = col
        self.size = int(size)
# END


def bezierpath(rs, re, qs, qe, ry, qy, V, col, alpha, label='', lw=0):
    import matplotlib.patches as patches
    from matplotlib.path import Path
    smid = (qs-rs)/2    # Start increment
    emid = (qe-re)/2    # End increment
    hmid = (qy-ry)/2    # Heinght increment
    if not V:
        verts = [(rs, ry),
                 (rs, ry+hmid),
                 (rs+2*smid, ry+hmid),
                 (rs+2*smid, ry+2*hmid),
                 (qe, qy),
                 (qe, qy-hmid),
                 (qe-2*emid, qy-hmid),
                 (qe-2*emid, qy-2*hmid),
                 (rs, ry),
                 ]
    else:
        verts = [(ry, rs),
                 (ry+hmid, rs),
                 (ry+hmid, rs+2*smid),
                 (ry+2*hmid, rs+2*smid),
                 (qy, qe),
                 (qy-hmid, qe),
                 (qy-hmid, qe-2*emid),
                 (qy-2*hmid, qe-2*emid),
                 (ry, rs),
                 ]
    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)
    # fig, ax = plt.subplots()
    patch = patches.PathPatch(path, facecolor=col, lw=lw, alpha=alpha, label=label, edgecolor=col)
    # patch = patches.PathPatch(path) #, facecolor=col, lw=lw, alpha=alpha, label=label, linecolor=c)
    # ax.set_xlim(-0.1, 2.5)
    # ax.set_ylim(-0.1, 2991784)
    # ax.add_patch(patch)
    return patch
# END


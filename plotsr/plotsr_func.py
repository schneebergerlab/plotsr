#!/usr/bin/env python3

import matplotlib.font_manager

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
FONT_NAMES = []
for fn in matplotlib.font_manager.findSystemFonts():
    try: FONT_NAMES.append(matplotlib.font_manager.get_font(fn).family_name)
    except RuntimeError: pass

for fn in FONT_NAMES:
    if re.findall('Comic', f, re.IGNORECASE) != []:
        print(fn)

"""
################################################################################
DEFINE READERS/PARSERS
################################################################################
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


class bedAnno():
    def __init__(self, c, start, end, genome, v):
        import matplotlib
        self.chr = c
        self.start = int(start)
        self.end = int(end)
        self.genome = genome
        # if genome in ['R', 'Q']: self.genome = genome
        # else: print("Incorrect genome for position {}:{}-{}. Setting marker to line".format(self.start, self.start, self.end))
        self.v = v
        self.mcol='black'
        self.msize=1
        self.mtype='.'
        self.text=''
        self.tcol='black'
        self.tsize=matplotlib.rcParamsDefault['font.size']
        self.tfont='Arial'
        self.tpos='0.05'
        self.logger = logging.getLogger("bedAnno")

# Define marker type
    def setmarker(self, marker):
        import matplotlib
        import sys
        marker = marker.split(":")
        if marker[0][:2] != 'm=' or marker[1][:2] != 'c=' or marker[2][:2] != 's=':
            raise ValueError("Incomplete information {} for marker at {}:{}-{}. Require: type, color, and color.".format(':'.join(marker), self.chr, self.start, self.end))
        # Set marker type
        m = marker[0].split("=")
        if m[1] not in MARKERS.keys():
            self.logger.error("Unrecongised marker used for position {}:{}-{}. Plotsr accepts markers defined in matplotlib (https://matplotlib.org/stable/api/markers_api.html) with some modifications.".format(self.start, self.start, self.end))
            for k, v in MARKERS.items():
                print("{} : {}".format(k, v))
            sys.exit()
        self.mtype = m[1] if m[1][0] != 'i' else int(m[1][1:])
        # Set marker color
        m = marker[1].split("=")
        try:
            if m[1][0] == '#': matplotlib.colors.to_rgb(m[1])
            else: matplotlib.colors.to_hex(m[1])
        except ValueError:
            self.logger.error("Error in using colour: {} for position {}:{}-{}. Use correct hexadecimal colours or named colours define in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html)".format(m[1], self.chr, self.start, self.end))
            sys.exit()
        self.mcol = m[1]
        # Set marker size
        m = marker[2].split("=")
        try: self.msize = int(m[1])
        except ValueError:
            self.logger.error("Non-integer size ({}) for marker at {}:{}-{}".format(m[1], self.start, self.start, self.end))
            sys.exit()
        # Set line marker when marker length > 1
        if self.end - self.start > 1:
            self.logger.warning("Range selected for position {}:{}-{}. Setting marker to line".format(self.start, self.start, self.end))
            self.mtype = "|" if self.v else "_"
        return

    # Define text type
    def settext(self, text):
        import warnings
        import matplotlib
        import sys
        text=text.split(":")
        if text[0][:2] != 't=' or text[1][:2] != 'c=' or text[2][:2] != 's=' or text[3][:2] != 'f=' or text[4][:2] != 'p=':
            raise ValueError("Incomplete information {} for text at {}:{}-{}. Require: text, color, size, and font_name.".format(':'.join(text), self.chr, self.start, self.end))
        # Set text
        t = text[0].split("=")
        self.text = t[1]
        # Set text color
        t = text[1].split("=")
        try:
            if t[1][0] == '#': matplotlib.colors.to_rgb(t[1])
            else: matplotlib.colors.to_hex(t[1])
        except ValueError:
            warnings.warn("Error in using colour: {} for position {}:{}-{}. Use correct hexadecimal colours or named colours define in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html)".format(t[1], self.chr, self.start, self.end))
            sys.exit()
        self.tcol = t[1]
        # Set text size
        t = text[2].split("=")
        try: self.tsize = int(t[1])
        except ValueError:
            raise ValueError("Non-integer size ({}) for marker at {}:{}-{}".format(t[1], self.chr, self.start, self.end))
        # Set text font
        t = text[3].split("=")
        if t[1] in FONT_NAMES: self.tfont = t[1]
        else:
            with open("plotsr_available_font_names.txt", 'w') as fout:
                fout.write("\n".join(FONT_NAMES))
            raise ValueError("Font ({}) for marker at {}:{}-{} is not available. Check plotsr_available_font_names.txt for list of available system markers".format(t[1], self.chr, self.start, self.end))
        # Set text position:
        t = text[4].split("=")
        self.tpos=float(t[1])
        return
# END


def readannobed(path, v, chrlengths):
    import warnings
    from collections import deque
    import matplotlib
    logger = logging.getLogger('readannobed')
    mdata = deque()
    with open(path, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                logger.warning("Skipping line\n{}".format(line.strip()))
                continue
            line = line.strip().split("\t")
            if len(line) not in [5, 6]:
                raise ValueError("Incomplete data in BED file line:\n{}\nMarker information (type:col:size) is necessary. Text annotation (text:col:size:font) is optional".format('\t'.join(line)))
                sys.exit()
            found = False
            for i in chrlengths:
                if i[0] == line[3]:
                    if line[0] in list(i[1].keys()):
                        found = True
            if not found:
                raise ValueError("Incorrect marker information. Chromosome: {} is not present in genome {}.".format(line[0], line[3]))
            anno = bedAnno(line[0], line[1], line[2], line[3], v)
            anno.setmarker(line[4])
            if len(line) == 6:
                anno.settext(line[5])
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

def gettrack(f, bw, chrlengths):
    """
    Reads BED file and creates a histogram for the number of bases per bin
    :param f:
    :param bw: binwidth
    :param chrlengths:
    :return:
    """
    from collections import deque
    import pandas as pd
    import sys
    bed = deque()
    logger = logging.getLogger("Reading track BED")
    # Read the BED file
    chrs = list(chrlengths[0][1].keys())
    with open(f, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if len(line) < 3:
                logger.error("Incomplete information in BED file at line: {}".format("\t".join(line)))
                sys.exit()
            if line[0] not in chrs:
                logger.error("Chromosome in BED is not present in FASTA at line: {}".format("\t".join(line)))
                sys.exit()
            try:
                bed.append([line[0], int(line[1]), int(line[2])])
            except ValueError:
                logger.error("Invalid values for line: {}".format("\t".join(line)))
                sys.exit()
    bed = pd.DataFrame(bed)

    # Create bins
    bins = {}
    for k,v in chrlengths[0][1].items():
        s = np.array(range(0, v, bw))
        e = np.concatenate((s[1:], [v]))
        bins[k] = np.array(list(zip(s, e)))

    bincnt = {}
    for k, v in bins.items():
        for r in v:
            garb = bed.loc[(bed[0]==k) & (bed[1]>=r[0]) & (bed[2]<r[1])]
            bincnt[(k, r[0], r[1])] = (garb[2] - garb[1]).sum()/bw


# END
"""
################################################################################
Validation and filtering
################################################################################
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

################################################################################
################################################################################
################################################################################
"""
################################################################################
Draw and plot
################################################################################
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
        ax.xaxis.grid(True, which='major', linestyle='--')
        ax.ticklabel_format(axis='x', useOffset=False, style='plain')
        xticks = ax.get_xticks()
        if max_l >= 1000000000:
            xticksl = xticks/1000000000
            ax.set_xlabel('chromosome position (in Gbp)')
        elif max_l >= 1000000:
            xticksl = xticks/1000000
            ax.set_xlabel('chromosome position (in Mbp)')
        elif max_l >= 1000:
            xticksl = xticks/1000
            ax.set_xlabel('chromosome position (in Kbp)')
        ax.set_xticks(xticks[:-1])
        ax.set_xticklabels(xticksl[:-1])
        ax.set_ylabel('reference chromosome id')
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
            yticksl = yticks/1000000000
            ax.set_ylabel('chromosome position (in Gbp)')
        elif max_l >= 1000000:
            yticksl = yticks/1000000
            ax.set_ylabel('chromosome position (in Mbp)')
        elif max_l >= 1000:
            yticksl = yticks/1000
            ax.set_ylabel('chromosome position (in Kbp)')
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticksl)
        ax.set_xlabel('reference chromosome id')
        ax.yaxis.grid(True, which='major', linestyle='--')
        ax.set_axisbelow(True)
    return ax, max_l
# END


def pltchrom(ax, chrs, chrgrps, chrlengths, V, S):
    import warnings
    import numpy as np
    chrlabs = [False]*len(chrlengths)
    max_l = np.max([chrlengths[i][1][v[i]] for v in chrgrps.values() for i in range(len(v))])
    CHRCOLS = plt.get_cmap('Dark2')                 # TODO: READ COLORS FROM CONFIG FILE
    if len(chrlengths) > 8:
        warnings.warn("More than 8 chromosomes are being analysed. This could result in different chromosomes having same color. Provide colors manually in config.")
    # Set chromosome direction
    pltchr = ax.axhline if not V else ax.axvline
    # Define indents # TODO: Check how the indents would change when plotting tracks
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
    return ax, indents
#END


def pltsv(ax, alignments, chrs, V, chrgrps, indents):
    adSynLab = False
    adInvLab = False
    adTraLab = False
    adDupLab = False

    alpha = 0.8 # TODO: Set alpha as a parameter
    for s in range(len(alignments)):
        df = alignments[s][1]
        for i in range(len(chrs)):
            offset = i if not V else -i
            # Plot syntenic regions
            for row in df.loc[(df['achr'] == chrgrps[chrs[i]][s]) & (df['type'] == 'SYN')].itertuples(index=False):
                if not adSynLab:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, V, col=COLORS[0], alpha=alpha, label='Syntenic')
                    adSynLab = True
                else:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, V, col=COLORS[0], alpha=alpha)
                ax.add_patch(p)

            # Plot Inversions
            for row in df.loc[(df['achr'] == chrgrps[chrs[i]][s]) & (df['type'] == 'INV')].itertuples(index=False):
                if not adInvLab:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, V, col=COLORS[1], alpha=alpha, label='Inversion', lw=0.1)
                    adInvLab=True
                else:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, V, col=COLORS[1], alpha=alpha, lw=0.1)
                ax.add_patch(p)
            # Plot Translocations
            for row in df.loc[(df['achr'] == chrgrps[chrs[i]][s]) & (df['type'].isin(['TRANS', 'INVTR']))].itertuples(index=False):
                if not adTraLab:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, V, col=COLORS[2], alpha=alpha, label='Translocation', lw=0.1)
                    adTraLab = True
                else:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, V, col=COLORS[2], alpha=alpha, lw=0.1)
                ax.add_patch(p)
            # Plot Duplications
            for row in df.loc[(df['achr'] == chrgrps[chrs[i]][s]) & (df['type'].isin(['DUP', 'INVDP']))].itertuples(index=False):
                if not adDupLab:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, V, col=COLORS[4], alpha=alpha, label='Duplication', lw=0.1)
                    adDupLab=True
                else:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, V, col=COLORS[4], alpha=alpha, lw=0.1)
                ax.add_patch(p)
    return ax
# END


def bezierpath(rs, re, qs, qe, ry, qy, v, col, alpha, label='', lw=0):
    import matplotlib.patches as patches
    from matplotlib.path import Path
    smid = (qs-rs)/2    # Start increment
    emid = (qe-re)/2    # End increment
    hmid = (qy-ry)/2    # Heinght increment
    if not v:
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


def drawmarkers(ax, b, v, chrlengths, indents, chrs, chrgrps):
    mdata = readannobed(b, v, chrlengths)
    for m in mdata:
        ind = [i for i in range(len(chrlengths)) if chrlengths[i][0] == m.genome][0]
        indent = indents[ind]
        offset = chrs.index([k for k, v in chrgrps.items() if v[ind] == m.chr][0])
        if not v:
            ax.plot(m.start, indent-offset, marker=m.mtype, color=m.mcol, markersize=m.msize)
            if m.text != '':
                ax.text(m.start, indent-offset+m.tpos, m.text, color=m.tcol, fontsize=m.tsize, fontfamily=m.tfont, ha='center', va='bottom')
        elif v:
            ax.plot(indent+offset, m.start, marker=m.mtype, color=m.mcol, markersize=m.msize)
            if m.text != '':
                ax.text(indent+offset-m.tpos, m.start, m.text, color=m.tcol, fontsize=m.tsize, fontfamily=m.tfont, ha='left', va='center', rotation='vertical')
    return ax
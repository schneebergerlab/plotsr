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

def filterinput(args, df, reflenghts, chrid_dict):
    # Get region length and filter out smaller SR
    df = df.loc[((df[2] - df[1]) >= args.s) | ((df[7] - df[6]) >= args.s) | (df[10]=='SYN')]
    df = df.loc[df[5] == [chrid_dict[i] for i in df[0]]]

    # Filter non-selected variations
    from warnings import warn
    if args.nosyn:
        df = df.loc[df[10] != 'SYN']
    if args.noinv:
        df = df.loc[df[10] != 'INV']
    if args.notr:
        df = df.loc[~df[10].isin(['TRANS', 'INVTR'])]
    if args.nodup:
        df = df.loc[~df[10].isin(['DUP', 'INVDP'])]
    if args.chr is not None:
        for chr in args.chr:
            if chr not in list(reflenghts.keys()):
                warn('Selected chromosome ID: {} is not present in the reference genome. Only use reference chromsome ID for selecting chromsomes.'.format(chr))
        df = df.loc[df[0].isin(args.chr)]
    return df
# END

def createribbon(df):
    """
    Combine continuous syntenic regions to get larger ribbons for syntenic blocks
    """
    df.sort_values([5, 6, 7], inplace=True)
    df['b'] = list(range(df.shape[0]))
    df.sort_values([0, 1, 2], inplace=True)
    df['a'] = list(range(df.shape[0]))
    # Group syntenic regions
    groups = deque()
    cg = deque()
    ca = -1
    cb = -1
    issyn = list(df[10] == 'SYN')
    a = list(df['a'])
    b = list(df['b'])
    chrchange = np.where((np.array(df[0][1:]) == np.array(df[0][0:-1])) == False)[0] + 1
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
            list(tmpdf[0])[0],
            list(tmpdf[1])[0],
            list(tmpdf[2])[-1],
            list(tmpdf[5])[0],
            list(tmpdf[6])[0],
            list(tmpdf[7])[-1]
        ])
    newsyn = DataFrame(list(newsyn), columns=[0, 1, 2, 5, 6, 7])
    newsyn[10] = 'SYN'

    df = df.drop(['a', 'b'], axis=1)
    df = df.loc[-(df[10] == 'SYN')]
    df = df.append(newsyn)
    df.sort_values([0, 1, 2], inplace=True)
    return df
# END


def drawax(ax, chrid_dict, reflenghts, qrylenghts, V, S, chrs):
    nchr = len(chrs)
    qchrs = [chrid_dict[k] for k in chrs]
    tick_pos = 1 - (S/2)
    ticklabels = chrs
    max_l = np.max([v for k, v in reflenghts.items() if k in chrs] + [v for k,v in qrylenghts.items() if k in qchrs])
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

def readbed(path, V):
    mdata = deque()
    with open(path, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if len(line) not in [7, 10]:
                print(line)
                print("Incomplete data in BED file line: \n {}.\n Marker information (type, col, size) is necessary. Text annotation (text, col, size) is optional".format('\t'.join(line)))
                continue
            anno = bedanno(line[0], line[1], line[2], line[3], V)
            anno.setmarker(line[4], line[5], line[6])
            if len(line) > 7:
                anno.settext(line[7], line[8], line[9])
            mdata.append(anno)
    return mdata
# END


def bezierpath(rs, re, qs, qe, ry, qy, V, col, alpha, label='', lw=0):
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
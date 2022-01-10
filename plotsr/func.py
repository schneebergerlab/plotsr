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
COLORS = ['#DEDEDE', '#FFA500', '#9ACD32', '#00BBFF']

FONT_NAMES = []
for fn in matplotlib.font_manager.findSystemFonts():
    try: FONT_NAMES.append(matplotlib.font_manager.get_font(fn).family_name)
    except RuntimeError: pass


"""
################################################################################
SUPPORT FUNCTIONS
################################################################################
"""

def mergeRanges(ranges):
    """
    Take a 2D numpy array, with each row as a range and return merged ranges
    i.e. ranges which are overlapping would be combined.
    :param ranges:
    :return:
    """
    from collections import deque
    import numpy as np
    if len(ranges) < 2:
        return ranges
    for i in ranges:
        if i[0] > i[1]:
            i[1], i[0] = i[0], i[1]
    ranges = ranges[ranges[:, 0].argsort()]
    min_value = ranges[0, 0]
    max_value = ranges[0, 1]
    out_range = deque()
    for i in ranges[1:]:
        if i[0] > max_value:
            out_range.append([min_value, max_value])
            min_value = i[0]
            max_value = i[1]
        elif i[1] > max_value:
            max_value = i[1]
    out_range.append([min_value, max_value])
    return np.array(out_range)
#END
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
    NUCBP = 'ACGTRYSWKMBDHVN'
    NUCBP += NUCBP.lower()
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
    ## Validate genomic fasta sequence: Count that only genomic BP are part of the sequence
    for k,v in out.items():
        if len(v) != sum([v.count(c) for c in NUCBP]):
            raise ValueError("Incorrect sequence for chromosome: {}. Genomic sequence can have only the following characters: {}".format(k, NUCBP))
    return out
# END


class bedAnno():
    def __init__(self, c, start, end, genome, v):
        import matplotlib
        import logging
        self.chr = c
        self.start = int(start)
        self.end = int(end)
        self.genome = genome
        self.v = v
        self.mt='.'
        self.mc='black'
        self.ms=8
        self.tt=''
        self.tc='black'
        self.ts=matplotlib.rcParams['font.size']
        self.tf='Arial'
        self.tp='0.05'
        self.logger = logging.getLogger("bedAnno")
        if self.start >= self.end:
            raise ValueError("Incorrect coordinates for marker position {}:{}:{}-{}. Start position must be less than end position. Skipping this marker.".format(self.genome, self.chr, self.start, self.end))

    # Add tags
    def addtags(self, tags):
        import matplotlib
        import sys
        for i in tags.split(";"):
            n, v = i.split(":")
            if hasattr(self, n):
                if n in ['mc', 'tc']:
                    try:
                        if v[0] == '#': matplotlib.colors.to_rgb(v)
                        else: matplotlib.colors.to_hex(v)
                    except ValueError:
                        self.logger.error("Error in using colour: {} for position {}:{}-{}. Use correct hexadecimal colours or named colours define in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html)".format(v, self.chr, self.start, self.end))
                        sys.exit()
                    setattr(self, n, v)
                elif n in ['ms', 'ts', 'tp']:
                    try: float(v)
                    except ValueError:
                        self.logger.error("Non-numerical value {} for {} at marker position: {}:{}-{}".format(v, n, self.chr, self.start, self.end))
                        sys.exit()
                    setattr(self, n, float(v))
                elif n in ['tf']:
                    if v not in FONT_NAMES:
                        with open("plotsr_available_font_names.txt", 'w') as fout:
                            fout.write("\n".join(FONT_NAMES))
                        raise ValueError("Selected font {} at marker position {}:{}-{} is not available. Check plotsr_available_font_names.txt for list of available system markers".format(v, self.chr, self.start, self.end))
                    setattr(self, n, v)
                elif n in ['mt']:
                    if v not in MARKERS.keys():
                        self.logger.error("Unrecongised marker used for position {}:{}-{}. Plotsr accepts markers defined in matplotlib (https://matplotlib.org/stable/api/markers_api.html) with some modifications.".format(self.chr, self.start, self.end))
                        for k, v in MARKERS.items():
                            print("{} : {}".format(k, v))
                        sys.exit()
                    if v[0] == 'i':
                        setattr(self, n, int(v[1:]))
                    else:
                        setattr(self, n, v)
                elif n in ['tt']:
                    setattr(self, n, v)
            else:
                raise ValueError("{} is not a valid tag".format(n))
        return
    # END
# END


def readannobed(path, v, chrlengths):
    from collections import deque
    import logging
    logger = logging.getLogger('readannobed')
    mdata = deque()
    with open(path, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                logger.debug("Skipping line\n{}".format(line.strip()))
                continue
            line = line.strip().split("\t")
            if len(line) < 4:
                logger.warning("Incomplete data in BED file line:\n{}\nGenome coordinate (chromosome, start,end) and genome name are required columns, while tags is an optional column".format('\t'.join(line)))
                continue
            found = False
            for i in chrlengths:
                if i[0] == line[3]:
                    if line[0] in list(i[1].keys()):
                        found = True
            if not found:
                logger.warning("Incorrect marker information. Chromosome: {} is not present in genome {}. Skipping it".format(line[0], line[3]))
                continue
            try:
                anno = bedAnno(line[0], line[1], line[2], line[3], v)
            except ValueError as e:
                print(e)
                continue
            if len(line) == 5:
                anno.addtags(line[4])
            # if len(line) == 6:
            #     anno.settext(line[5])
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
    try:
        df = DataFrame(list(syri_regs))[[0, 1, 2, 5, 6, 7, 10]]
    except KeyError:
        raise ImportError("Incomplete input file {}, syri.out file should have 11 columns.".format(f))
    df[[0, 5, 10]] = df[[0, 5, 10]].astype(str)
    try:
        df[[1, 2, 6, 7]] = df[[1, 2, 6, 7]].astype(int)
    except ValueError:
        raise ValueError("Non-numerical values used as genome coordinates in {}. Exiting".format(f))
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
    try:
        df[[0, 3, 6]] = df[[0, 3, 6]].astype(str)
    except KeyError:
        raise ImportError("Incomplete input file {}, BEDPE file should have 7 columns.".format(f))
    try:
        df[[1, 2, 4, 5]] = df[[1, 2, 4, 5]].astype(int)
    except ValueError:
        raise ValueError("Non-numerical values used as genome coordinates in {}. Exiting".format(f))
    df[[1, 2, 4, 5]] = df[[1, 2, 4, 5]].astype(int)
    df[[1, 4]] = df[[1, 4]] + 1 # Makes range closed as the terminal bases are also included
    # chr ID map
    chrid = []
    chrid_dict = OrderedDict()
    for i in np.unique(df[0]):
        chrid.append((i, np.unique(df.loc[(df[0] == i) & (df[6] == 'SYN'), 3])[0]))
        chrid_dict[i] = np.unique(df.loc[(df[0] == i) & (df[6] == 'SYN'), 3])[0]
    df.columns = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type']
    return df, chrid_dict
# END


class track():
    def __init__(self, f, n):
        import matplotlib
        import logging
        self.f = f
        self.n = n
        self.ft = 'bed'
        self.bw = 100000
        self.nc = 'black'
        self.ns = matplotlib.rcParams['font.size']
        self.nf = 'Arial'
        self.lc = 'black'
        self.lw = 1
        self.bc = 'lightgrey'
        self.ba = 0.7
        self.logger = logging.getLogger("track")
        self.bincnt = {}

    # Add tags
    def addtags(self, tags):
        import matplotlib
        import sys
        for i in tags.split(";"):
            n, v = i.split(":")
            if hasattr(self, n):
                ## Colour parameters
                if n in ['nc', 'lc', 'bc']:
                    try:
                        # print(n, v)
                        if v[0] == '#':
                            matplotlib.colors.to_rgb(v)
                        else:
                            matplotlib.colors.to_hex(v)
                    except ValueError:
                        self.logger.error("Error in using colour: {} for track {}. Use correct hexadecimal colours or named colours define in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html)".format(v, self.n))
                        # continue
                        sys.exit()
                    setattr(self, n, v)
                ## Numerical parameters
                elif n in ['ns', 'bw', 'lw', 'ba']:
                    try: float(v)
                    except ValueError:
                        self.logger.error("Non-numerical value {} for {} in track{}".format(v, n, self.n))
                        sys.exit()
                    setattr(self, n, float(v))
                ## Font parameters
                elif n in ['nf']:
                    if v not in FONT_NAMES:
                        with open("plotsr_available_font_names.txt", 'w') as fout:
                            fout.write("\n".join(FONT_NAMES))
                        raise ValueError("Font {} in track {} is not available. Check plotsr_available_font_names.txt for list of available system markers".format(v, self.n))
                    setattr(self, n, v)
                ## Input type parameters
                elif n == 'ft':
                    if v not in ['bed', 'bedgraph']:
                        raise ValueError("Unknown filetype specified for track {}. Please provide input in bed or bedgraph. Use bed for providing genome coordinate ranges or bedgraph for providing histogram.".format(v, self.n))
                    setattr(self, n, v)
            else:
                raise ValueError("{} is not a valid tag".format(n))
        return
    #END

    # Read input bed file and get histogram with binwidths==self.bw
    def _readbed(self, chrlengths):
        from collections import deque, defaultdict
        import numpy as np
        bw = int(self.bw)
        # Read the BED file
        chrpos = {k: np.zeros(v, dtype=np.int0) for k, v in chrlengths[0][1].items()}
        skipchrs = []
        with open(self.f, 'r') as fin:
            # with open('snps.sample.bed', 'r') as fin: # TODO: Delete line
            for line in fin:
                line = line.strip().split()
                if len(line) < 3:
                    self.logger.warning("Incomplete information in BED file at line: {}. Skipping it.".format("\t".join(line)))
                    continue
                if line[0] not in chrpos.keys():
                    if line[0] not in skipchrs:
                        self.logger.warning("Chromosome in BED is not present in FASTA or not selected for plotting. Skipping it. BED line: {}".format("\t".join(line)))
                        skipchrs.append(line[0])
                    continue
                try:
                    chrpos[line[0]][int(line[1]):int(line[2])] = 1
                except ValueError:
                    self.logger.warning("Invalid values for line: {}. Skipping it.".format("\t".join(line)))
        # Create bins
        bins = {}
        for k, v in chrlengths[0][1].items():
            s = np.array(range(0, v, bw))
            e = np.concatenate((s[1:], [v]))
            bins[k] = np.array(list(zip(s, e)))
        bincnt = defaultdict(deque)
        for k, v in bins.items():
            for r in v:
                bincnt[k].append(((r[0]+r[1])/2, np.sum(chrpos[k][r[0]:r[1]])/bw))
        self.bincnt = bincnt
        return
    # END

    ## Read input bedgraph file
    def _readbedgraph(self, chrlengths):
        from collections import deque, defaultdict
        import numpy as np
        bw = int(self.bw)
        chrpos = {k: np.zeros(v, dtype=np.int0) for k, v in chrlengths[0][1].items()}
        skipchrs = []
        with open(self.f, 'r') as fin:
            for line in fin:
                if line[0] == '#': continue
                line = line.strip().split()
                if line[0] == 'track': continue
                if len(line) < 4:
                    self.logger.warning("Incomplete information in bedgraph file at line: {}. Skipping it.".format("\t".join(line)))
                    continue
                if line[0] not in chrpos.keys():
                    if line[0] not in skipchrs:
                        self.logger.warning("Chromosome in BED is not present in FASTA or not selected for plotting. Skipping it. BED line: {}".format("\t".join(line)))
                        skipchrs.append(line[0])
                    continue
                try:
                    chrpos[line[0]][int(line[1]):int(line[2])] += int(line[3])
                except ValueError:
                    self.logger.warning("Invalid values for line: {}. Skipping it.".format("\t".join(line)))
        ## Create bins
        bins = {}
        for k, v in chrlengths[0][1].items():
            s = np.array(range(0, v, bw))
            e = np.concatenate((s[1:], [v]))
            bins[k] = np.array(list(zip(s, e)))
        bincnt = defaultdict(deque)
        maxv = 0
        for k, v in bins.items():
            for r in v:
                if np.sum(chrpos[k][r[0]:r[1]]) > maxv:
                    maxv = np.sum(chrpos[k][r[0]:r[1]])
        for k, v in bins.items():
            for r in v:
                bincnt[k].append(((r[0]+r[1])/2, np.sum(chrpos[k][r[0]:r[1]])/maxv))
        self.bincnt = bincnt
        return
    # END

    def readdata(self, chrlengths):
        if self.ft == 'bed':
            self._readbed(chrlengths)
        elif self.ft == 'bedgraph':
            self._readbedgraph(chrlengths)
        return
    # END
# END


def readtrack(f, chrlengths):
    """
    Reads BED file and creates a histogram for the number of bases per bin
    :param f:
    :param bw: binwidth
    :param chrlengths:
    :return:
    """
    from collections import deque
    tdata = deque()
    with open(f, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                # logger.warning("Skipping line\n{}".format(line.strip()))
                continue
            line = line.strip().split("\t")
            if len(line) == 0: continue
            if len(line) < 2:
                raise ValueError("Incomplete data in track file line:\n{}\nFile path and name are necessary columns, tags is optional columns".format('\t'.join(line)))
                sys.exit()
            t = track(line[0], line[1])
            # Read tags
            if len(line) == 3:
                t.addtags(line[2])
            # Reads BED
            t.readdata(chrlengths)
            # print(t.bincnt)
            tdata.append(t)
    return tdata
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
    from collections import deque
    import numpy as np
    import os
    from matplotlib.pyplot import get_cmap
    out = deque()
    errmess1 = 'Chromosome ID: {} in structural annotation file: {} not present in genome fasta: {}. Exiting.'
    errmess2 = 'For chromosome ID: {}, length in genome fasta: {} is less than the maximum coordinate in the structural annotation file: {}. Exiting.'
    tags = {'lc': {}, 'lw': {}}
    taglist = set(tags.keys())
    # Count number of genomes
    count = 0
    with open(genf, 'r') as fin:
        i = 0
        for line in fin:
            if line[0] =='#':
                continue
            count += 1
    if count <= 10:
        CHRCOLS = [matplotlib.colors.to_hex(get_cmap('tab10')(i)) for i in range(count)]
    else:
        CHRCOLS = [matplotlib.colors.to_hex(get_cmap('gist_rainbow')(int(255/count) * i)) for i in range(0, count)]
        if count % 2 != 0:
            m = CHRCOLS[int((count/2) + 1)]
        CHRCOLS = [j for i in range(int(count/2)) for j in [CHRCOLS[i]] + [CHRCOLS[int(i +count/2)]]] + [m]
    with open(genf, 'r') as fin:
        i = 0
        for line in fin:
            if line[0] =='#':
                continue
            line = line.strip().split("\t")
            if len(line) < 2:
                raise ImportError("Incomplete genomic information.\nExpected format for the genome file:\npath_to_genome1\tgenome1_id\ttags\npath_to_genome2\tgenome2_id\ttags")
            try:
                # length of individual chromosomes
                glen = {c: len(seq) for c, seq in readfasta(line[0]).items()}
            except Exception as e:
                raise ImportError("Error in reading fasta: {}\n{}".format(line[1], e))
            # Check cases when the genome is the query-genome
            if i > 0:
                try:
                    df = als[i-1][1]
                except IndexError:
                    raise ImportError("Cannot read data for genome: {}. Make sure that structural annotation data for all genomes is provided in correct order. Exiting.".format(line[0]))
                chrs = np.unique(df['bchr'])
                for c in chrs:
                    if c not in list(glen.keys()):
                        raise ImportError(errmess1.format(c, als[i-1][0], os.path.basename(line[1])))
                    if np.max(np.max(df.loc[df['bchr'] == c, ['bstart', 'bend']])) > glen[c]:
                        raise ImportError(errmess2.format(c, os.path.basename(fin), als[i-1][0]))
            # Check cases when the genome is the reference-genome
            if i < len(als):
                try:
                    df = als[i][1]
                except IndexError:
                    raise ImportError("Cannot read data for genome: {}. Make sure that structural annotation data for all genomes is provided in correct order. Exiting.".format(line[0]))
                chrs = np.unique(df['achr'])
                for c in chrs:
                    if c not in list(glen.keys()):
                        raise ImportError(errmess1.format(c, als[i][0], os.path.basename(line[1])))
                    if np.max(np.max(df.loc[df['achr'] == c, ['astart', 'aend']])) > glen[c]:
                        raise ImportError(errmess2.format(c, os.path.basename(fin), als[i][0]))
            out.append((line[1], glen))

            # TODO: create a class for genome
            # Set default tag values
            tags['lc'][line[1]] = CHRCOLS[i]
            tags['lw'][line[1]] = 1
            # Update tag values
            if len(line) > 2:
                tg = line[2].split(';')
                for t in tg:
                    t = t.split(':')
                    if t[0] in taglist:
                        if t[0] in ['lw']:
                            tags[t[0]][line[1]] = float(t[1])
                        else:
                            tags[t[0]][line[1]] = t[1]
            i += 1
    return out, tags
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
# END


def selectregion(reg, chrlengths, al, chrids):
    import pandas as pd
    import copy
    from collections import OrderedDict, deque
    alignments = copy.deepcopy(al)
    genids = [i[0] for i in chrlengths]
    if len(reg) != 3:
        raise ValueError("Incorrect values parsed to --reg. Provide values in the following format: GenomeID:ChromosomeID:Start-End. GenomeID and ChromosomeID cannot have ':' character.")
    # reg = ['cvi', 'LR699761.1', '12000000-13500000'] #TODO: Delete this line
    reg = reg[:2] + list(map(int, reg[2].split("-")))
    # Check that the region is in the given genomes
    if reg[0] not in genids:
        raise ValueError("Genome ID in --reg do not match with any genome ID provided with --genomes. Exiting.")
    if reg[1] not in chrlengths[genids.index(reg[0])][1].keys():
        raise ValueError("Chromosome ID in --reg do not match with any chromosome ID for genome {}. Exiting.".format(reg[0]))
    l = chrlengths[genids.index(reg[0])][1][reg[1]]
    if any([reg[2] < 1, reg[2] > l, reg[3] < 1, reg[3] > l, reg[2] > reg[3]]):
        raise ValueError("Incorrect chromosome coordinates provided for --reg. Exiting.")
    # Select alignments that are within the selected region in the focal genome. For other genomes, maximal syntenic region coordinate would be used as boundaries.
    ind = genids.index(reg[0])
    newal = [0]*len(alignments)
    # Alignments for genomes after the focal genome
    c, s, e = reg[1:]
    for i in range(0, ind).__reversed__():
        df = alignments[i][1].copy()
        syncrd = df.loc[(df['bchr'] == c) & (df['bstart'] <= e) & (df['bend'] >= s) & (df['type'] == 'SYN')].copy()

        if min(syncrd['bstart']) < s:
            lower = syncrd['bstart'] < s
            a = (syncrd['bend'][lower] - s)/(syncrd['bend'][lower] - syncrd['bstart'][lower])
            syncrd.loc[lower, 'astart'] = (syncrd['aend'][lower] - (syncrd['aend'][lower] - syncrd['astart'][lower])*a).astype(int)
            syncrd.loc[lower, 'bstart'] = s

        if max(syncrd['bend']) > e:
            higher = syncrd['bend'] > e
            a = (e-syncrd['bstart'][higher])/(syncrd['bend'][higher] - syncrd['bstart'][higher])
            syncrd.loc[higher, 'aend'] = (syncrd['astart'][higher] + (syncrd['aend'][higher] - syncrd['astart'][higher])*a).astype(int)
            syncrd.loc[higher, 'bend'] = e
        bc = syncrd['achr'].tolist()[0]
        bs = min(syncrd['astart'])
        be = max(syncrd['aend'])
        srcrd = df.loc[(df['achr'] == bc) &
                       (df['astart'] >= bs) &
                       (df['aend'] <= be) &
                       (df['bchr'] == c) &
                       (df['bstart'] >= s) &
                       (df['bend'] <= e) &
                       (df['type'] != 'SYN')].reset_index(drop=True).copy()
        garb = pd.concat([syncrd, srcrd])
        garb.sort_values(['achr', 'astart', 'aend'], inplace=True)
        garb.reset_index(inplace=True, drop=True)
        newal[i] = garb.copy()
        # print(newal)
        c, s, e = bc, bs, be
    # Alignments for genome after the focal genome
    c, s, e = reg[1:]
    for i in range(ind, len(alignments)):
        df = alignments[i][1].copy()
        syncrd = df.loc[(df['achr'] == c) & (df['astart'] <= e) & (df['aend'] >= s) & (df['type'] == 'SYN')].copy()
        # Alter alignments with coordinate less than the start position
        if min(syncrd['astart']) < s:
            lower = syncrd['astart'] < s
            a = (syncrd['aend'][lower] - s)/(syncrd['aend'][lower] - syncrd['astart'][lower])
            syncrd.loc[lower, 'bstart'] = (syncrd['bend'][lower] - (syncrd['bend'][lower] - syncrd['bstart'][lower])*a).astype(int)
            syncrd.loc[lower, 'astart'] = s
        # Alter alignments with coordinate more than the end position
        if max(syncrd['aend']) > e:
            higher = syncrd['aend'] > e
            a = (e - syncrd['astart'][higher])/(syncrd['aend'][higher] - syncrd['astart'][higher])
            syncrd.loc[higher, 'bend'] = (syncrd['bstart'][higher] + (syncrd['bend'][higher] - syncrd['bstart'][higher])*a).astype(int)
            syncrd.loc[higher, 'aend'] = e
        bc = syncrd['bchr'].tolist()[0]
        bs = min(syncrd['bstart'])
        be = max(syncrd['bend'])
        srcrd = df.loc[(df['achr'] == c) &
                       (df['astart'] >= s) &
                       (df['aend'] <= e) &
                       (df['bchr'] == bc) &
                       (df['bstart'] >= bs) &
                       (df['bend'] <= be) &
                       (df['type'] != 'SYN')].reset_index(drop=True).copy()
        garb = pd.concat([syncrd, srcrd])
        garb.sort_values(['achr', 'astart', 'aend'], inplace=True)
        garb.reset_index(inplace=False, drop=True)
        newal[i] = garb
        c, s, e = bc, bs, be
    # Update alignments, chrs and chrgrps
    for i in range(len(alignments)):
        alignments[i][1] = newal[i].copy()
    chrs = [k for k in chrids[0][1].keys() if k in alignments[0][1]['achr'].unique()]
    chrgrps = OrderedDict()
    for c in chrs:
        cg = deque([c])
        cur = c
        for i in range(len(chrids)):
            n = chrids[i][1][cur]
            cg.append(n)
            cur = n
        chrgrps[c] = cg
    return alignments, chrs, chrgrps



def createribbon(df):
    """
    Combine continuous syntenic regions to get larger ribbons for syntenic blocks
    """
    import numpy as np
    from pandas import DataFrame
    from collections import deque
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


def drawax(ax, chrgrps, chrlengths, v, s, minl=0, maxl=-1):
    import numpy as np
    nchr = len(chrgrps)
    # qchrs = [chrid_dict[k] for k in chrs]
    # tick_pos = 1 - (S/2)
    bottom_limit = -0.1
    upper_limit = 0.1
    ticklabels = list(chrgrps.keys())
    if maxl == -1:
        maxl = np.max([chrlengths[i][1][c[i]] for c in chrgrps.values() for i in range(len(c))])
    if not v:
        tick_pos = s/2 + 0.1
        ax.set_ylim(bottom_limit, nchr+upper_limit)
        ax.set_yticks([tick_pos+i for i in range(nchr)])
        ax.set_yticklabels(ticklabels[::-1])
        ax.tick_params(axis='y', right=False, left=False)
        ax.set_xlim(minl, maxl)
        ax.xaxis.grid(True, which='both', linestyle='--')
        ax.ticklabel_format(axis='x', useOffset=False, style='plain')
        xticks = ax.get_xticks()
        if maxl >= 1000000000:
            xticksl = xticks/1000000000
            ax.set_xlabel('Chromosome position (in Gbp)')
        elif maxl >= 1000000:
            xticksl = xticks/1000000
            ax.set_xlabel('Chromosome position (in Mbp)')
        elif maxl >= 1000:
            xticksl = xticks/1000
            ax.set_xlabel('Chromosome position (in Kbp)')
        ax.set_xticks(xticks[:-1])
        ax.set_xticklabels(xticksl[:-1])
        ax.set_ylabel('Reference Chromosome ID')
        ax.set_axisbelow(True)
    else:
        tick_pos = 1 - 0.1 - s/2
        ax.set_xlim(bottom_limit, nchr+upper_limit)
        ax.set_xticks([tick_pos+i for i in range(nchr)])
        ax.set_xticklabels(ticklabels)
        ax.tick_params(axis='x', top=False, bottom=False)
        ax.set_ylim(minl, maxl)
        ax.ticklabel_format(axis='y', useOffset=False, style='plain')
        yticks = ax.get_yticks()
        if maxl >= 1000000000:
            yticksl = yticks/1000000000
            ax.set_ylabel('Chromosome position (in Gbp)')
        elif maxl >= 1000000:
            yticksl = yticks/1000000
            ax.set_ylabel('Chromosome position (in Mbp)')
        elif maxl >= 1000:
            yticksl = yticks/1000
            ax.set_ylabel('Chromosome position (in Kbp)')
        ax.set_yticks(yticks[:-1])
        ax.set_yticklabels(yticksl[:-1])
        ax.set_xlabel('Reference Chromosome ID')
        ax.yaxis.grid(True, which='both', linestyle='--')
        ax.set_axisbelow(True)
    return ax, maxl
# END


def pltchrom(ax, chrs, chrgrps, chrlengths, v, S, chrtags, minl=0, maxl=-1):
    import warnings
    import numpy as np
    chrlabs = [False]*len(chrlengths)
    # maxl = np.max([chrlengths[i][1][v[i]] for v in chrgrps.values() for i in range(len(v))])
    # if len(chrlengths) > 8:
    #     warnings.warn("More than 8 chromosomes are being analysed. This could result in different chromosomes having same color. Provide colors manually in config.")
    # Set chromosome direction
    # pltchr = ax.axhline if not V else ax.axvline
    pltchr = ax.hlines if not v else ax.vlines
    # Define indents
    step = S/(len(chrlengths)-1)
    chrlabels = []
    if not v:
        rend = len(chrs)-1+S+0.1
        indents = [rend - (i*step) for i in range(len(chrlengths))]
    elif v:
        rend = 1-S-0.1
        indents = [rend + (i*step) for i in range(len(chrlengths))]
    for s in range(len(chrlengths)):
        for i in range(len(chrs)):
            offset = i if not v else -i
            if maxl == -1:
                maxcoord = chrlengths[s][1][chrgrps[chrs[i]][s]]
            else:
                maxcoord = maxl
            if not chrlabs[s]:
                chrlabels.append(pltchr(indents[s]-offset, minl, maxcoord,
                                        color=chrtags['lc'][chrlengths[s][0]],
                                        linewidth=chrtags['lw'][chrlengths[s][0]], label=chrlengths[s][0]))
                chrlabs[s] = True
            else:
                pltchr(indents[s]-offset, minl, maxcoord, color=chrtags['lc'][chrlengths[s][0]], linewidth=chrtags['lw'][chrlengths[s][0]])
    return ax, indents, chrlabels
# END


def pltsv(ax, alignments, chrs, v, chrgrps, indents):
    from collections import deque
    adSynLab = False
    adInvLab = False
    adTraLab = False
    adDupLab = False

    alpha = 0.8 # TODO: Set alpha as a parameter
    svlabels = deque()
    for s in range(len(alignments)):
        df = alignments[s][1]
        for i in range(len(chrs)):
            offset = i if not v else -i
            # Plot syntenic regions
            for row in df.loc[(df['achr'] == chrgrps[chrs[i]][s]) & (df['type'] == 'SYN')].itertuples(index=False):
                if not adSynLab:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, v, col=COLORS[0], alpha=alpha, label='Syntenic')
                    adSynLab = True
                    syn = ax.add_patch(p)
                    svlabels.append(syn)
                else:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, v, col=COLORS[0], alpha=alpha)
                    ax.add_patch(p)
            # Plot Inversions
            for row in df.loc[(df['achr'] == chrgrps[chrs[i]][s]) & (df['type'] == 'INV')].itertuples(index=False):
                if not adInvLab:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, v, col=COLORS[1], alpha=alpha, label='Inversion', lw=0.1)
                    adInvLab=True
                    inv = ax.add_patch(p)
                    svlabels.append(inv)
                else:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, v, col=COLORS[1], alpha=alpha, lw=0.1)
                    ax.add_patch(p)
            # Plot Translocations
            for row in df.loc[(df['achr'] == chrgrps[chrs[i]][s]) & (df['type'].isin(['TRANS', 'INVTR']))].itertuples(index=False):
                if not adTraLab:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, v, col=COLORS[2], alpha=alpha, label='Translocation', lw=0.1)
                    adTraLab = True
                    tra = ax.add_patch(p)
                    svlabels.append(tra)
                else:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, v, col=COLORS[2], alpha=alpha, lw=0.1)
                    ax.add_patch(p)
            # Plot Duplications
            for row in df.loc[(df['achr'] == chrgrps[chrs[i]][s]) & (df['type'].isin(['DUP', 'INVDP']))].itertuples(index=False):
                if not adDupLab:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, v, col=COLORS[3], alpha=alpha, label='Duplication', lw=0.1)
                    adDupLab=True
                    dup = ax.add_patch(p)
                    svlabels.append(dup)
                else:
                    p = bezierpath(row[1], row[2], row[4], row[5], indents[s]-offset, indents[s+1]-offset, v, col=COLORS[3], alpha=alpha, lw=0.1)
                    ax.add_patch(p)
    return ax, svlabels
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


def drawmarkers(ax, b, v, chrlengths, indents, chrs, chrgrps, minl=0, maxl=-1):
    import logging
    logger = logging.getLogger('drawmarkers')
    mdata = readannobed(b, v, chrlengths)
    for m in mdata:
        ind = [i for i in range(len(chrlengths)) if chrlengths[i][0] == m.genome][0]
        indent = indents[ind]
        chrid = [k for k, c in chrgrps.items() if c[ind] == m.chr]
        if chrid != []:
            offset = chrs.index(chrid[0])
        else:
            logger.warning("Cannot draw marker at {}:{}-{} on genome {} because the chromosome is not selected for plotting. Skipping it.".format(m.chr, m.start, m.end, m.genome))
            continue
        if maxl != -1:
            if m.start < minl or m.end > maxl:
                logger.warning("Cannot draw marker at {}:{}-{} on genome {} because the marker position is out of the selected range. Skipping it.".format(m.chr, m.start, m.end, m.genome))
                continue
        if not v:
            ax.plot(m.start, indent-offset, marker=m.mt, color=m.mc, markersize=m.ms)
            if m.tt != '':
                ax.text(m.start, indent-offset+m.tp, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='center', va='bottom')
        elif v:
            ax.plot(indent+offset, m.start, marker=m.mt, color=m.mc, markersize=m.ms)
            if m.tt != '':
                ax.text(indent+offset-m.tp, m.start, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='left', va='center', rotation='vertical')
    return ax
# END


def drawtracks(ax, tracks, s, chrgrps, chrlengths, v, minl, maxl):
    from matplotlib.patches import Rectangle
    import numpy as np
    # TODO: Read gap values from base.cfg
    th = (1 - s - 2*0.1)/len(tracks)
    if th < 0.01:
        raise RuntimeError("Decrease the value of -S to plot tracks correctly. Exiting.")
    cl = len(chrgrps.keys())
    chrs = list(chrgrps.keys())
    # Define space between track and track label
    # TODO: read margin spacing from base config
    if maxl == -1:
        margin = np.max([chrlengths[i][1][v[i]] for v in chrgrps.values() for i in range(len(v))])/500
    else:
        margin = maxl/500
    for i in range(len(tracks)):
        bedbin = tracks[i].bincnt
        # print(bedbin)
        for j in range(cl):
            if maxl != -1:
                chrpos = [k[0] for k in bedbin[chrs[j]] if minl <= k[0] <= maxl]
                tpos = [k[1] for k in bedbin[chrs[j]] if minl <= k[0] <= maxl]
            else:
                chrpos = [k[0] for k in bedbin[chrs[j]]]
                tpos = [k[1] for k in bedbin[chrs[j]]]
            tposmax = max(tpos)
            diff = 0.7*th
            if not v:
                y0 = cl - j - th*(i+1)
                ypos = [(t*diff/tposmax)+y0 for t in tpos]
                # TODO: parameterise colour
                ax.add_patch(Rectangle((0, y0), chrlengths[0][1][chrs[j]], diff,  linewidth=0, facecolor=tracks[i].bc, alpha=tracks[i].ba))
                # ax.plot(chrpos, ypos, color=tracks[i].lc, lw=tracks[i].lw)
                ax.fill_between(chrpos, ypos, y0, color=tracks[i].lc, lw=tracks[i].lw)
                if maxl == -1:
                    ax.text(chrlengths[0][1][chrs[j]] + margin, y0 + diff/2, tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='left', va='center', rotation='horizontal')
                else:
                    ax.text(maxl + margin, y0 + diff/2, tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='left', va='center', rotation='horizontal')
            else:
                x0 = j + (i+1)*th - diff
                xpos = [x0 + diff - (t*diff/tposmax) for t in tpos]
                ax.add_patch(Rectangle((x0, 0), diff, chrlengths[0][1][chrs[j]], linewidth=0, facecolor=tracks[i].bc, alpha=tracks[i].ba))
                # ax.plot(xpos, chrpos, color=tracks[i].lc, lw=tracks[i].lw)
                ax.fill_betweenx(chrpos, xpos, x0+diff, color=tracks[i].lc, lw=tracks[i].lw)
                if maxl == -1:
                    ax.text(x0 + diff/2, chrlengths[0][1][chrs[j]] + margin,tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='center', va='bottom', rotation='vertical')
                else:
                    ax.text(x0 + diff/2, maxl + margin, tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='center', va='bottom', rotation='vertical')
    return ax
# END

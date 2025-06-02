import sys
import logging
import matplotlib.font_manager
# Constants
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
    try:
        FONT_NAMES.append(matplotlib.font_manager.get_font(fn).family_name)
    except RuntimeError:
        pass
FONT_NAMES = sorted(set(FONT_NAMES))


'''
################################################################################
# SUPPORT FUNCTIONS
################################################################################
'''


def setlogconfig(lg):
    import logging.config
    logging.config.dictConfig({
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'log_file': {
                'format': "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s",
            },
            'stdout': {
                'format': "%(name)s - %(levelname)s - %(message)s",
            },
        },
        'handlers': {
            'stdout': {
                'class': 'logging.StreamHandler',
                'formatter': 'stdout',
                'level': 'WARNING',
            },
        },
        'loggers': {
            '': {
                'level': lg,
                'handlers': ['stdout'],
                # 'handlers': ['stdout', 'log_file'],
            },
        },
    })
# END


class CustomFormatter(logging.Formatter):
    """
    https://betterstack.com/community/questions/how-to-color-python-logging-output/
    """
    grey = "\x1b[0;49;90m"
    green = "\x1b[0;49;32m"
    yellow = "\x1b[0;49;93m"
    red = "\x1b[0;49;31m"
    bold_red = "\x1b[0;49;31;21m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: green + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
# END


def getfilehandler(fin, level):
    logformat = "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s"

    def filehandler():
        fhandle = logging.FileHandler(fin)
        fhandle.setLevel(level)
        fhandle.setFormatter(logging.Formatter(logformat))
        return fhandle
    return filehandler
# END


def streamhandler():
    shandler = logging.StreamHandler()
    shandler.setLevel(logging.INFO)
    shandler.setFormatter(CustomFormatter())
    return shandler
# END


def definelogger(fhandler):
    def getlogger(name):
        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)
        logger.addHandler(fhandler())
        logger.addHandler(streamhandler())
        logger.propagate = False
        return logger
    return getlogger
# END


def mergeranges(ranges):
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
# END


"""
################################################################################
DEFINE READERS/PARSERS
################################################################################
"""


def readbasecfg(f, v):
    import logging
    import matplotlib
    logger = logging.getLogger('readbasecfg')
    cfg = {}
    # Set alignment parameters
    cfg['syncol'] = '#DEDEDE'
    cfg['invcol'] = '#FFA500'
    cfg['tracol'] = '#9ACD32'
    cfg['dupcol'] = '#00BBFF'
    cfg['synlwd'] = 0
    cfg['invlwd'] = 0.1
    cfg['tralwd'] = 0.1
    cfg['duplwd'] = 0.1
    cfg['alpha'] = 0.8

    # Set chromosome margins
    cfg['chrmar'] = 0.1
    cfg['exmar'] = 0.1
    cfg['marginchr'] = 0.01     # Set ITX margin

    # Set legend properties
    cfg['legend'] = True
    cfg['genlegcol'] = -1
    cfg['bbox'] = [0, 1.01, 0.5, 0.3]
    cfg['bbox_v'] = [0, 1.1, 0.5, 0.3]
    cfg['bboxmar'] = 0.5

    # track properties
    cfg['norm'] = 'T'

    # axis properties
    cfg['maxl'] = -1
    cfg['genname'] = 'T'
    cfg['itxalign'] = 'D'

    if f == '':
        return cfg
    with open(f, 'r') as fin:
        cfgk = set(cfg.keys())
        for line in fin:
            if line[:2] == '##':
                continue
            line = line.strip().split()
            if len(line) == 0:
                continue
            line = line[0].split(':')
            if line[0] not in cfgk:
                logger.error("{} is not a valid config parameter. Using default value.".format(line[0]))
                continue
            if line[0] in ['syncol', 'invcol', 'tracol', 'dupcol']:
                try:
                    if line[1] == '#': matplotlib.colors.to_rgb(line[1])
                    else: matplotlib.colors.to_hex(line[1])
                except ValueError:
                    logger.error("Error in using colour: {} for {}. Use correct hexadecimal colours or named colours defined in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html). Using default value.".format(line[1], line[0]))
                    continue
                cfg[line[0]] = line[1]
            elif line[0] in ['synlwd', 'invlwd', 'tralwd', 'duplwd', 'alpha', 'chrmar', 'exmar', 'bboxmar', 'genlegcol', 'marginchr', 'maxl']:
                try:
                    float(line[1])
                except ValueError:
                    logger.error("Non-numerical value {} provided for {}. Using default value.".format(line[1], line[0]))
                    continue
                cfg[line[0]] = float(line[1])
            elif line[0] in ['bbox', 'bbox_v']:
                if line[0] == 'bbox' and v:
                    continue
                if line[0] == 'bbox_v' and not v:
                    continue
                line[1] = line[1].split(',')
                if len(line[1]) != 4:
                    logger.error("BBOX requires four values ({} provided: {}). Using default values.".format(len(line[1]), line[1]))
                    continue
                try:
                    cfg[line[0]] = [float(i) for i in line[1]]
                    cfg['bbox'] = cfg['bbox_v'] if v else cfg['bbox']
                except ValueError:
                    logger.error("Non-numerical values {} provided for {}. Using default value.".format(line[1], line[0]))
                    continue
                cfg['bboxmar'] = [float(i) for i in line[1]]
            elif line[0] in {'legend', 'norm', 'genname'}:
                if line[1] not in ['T', 'F']:
                    logger.warning("Invalid value {} for {} in base.cfg. Valid values: T/F".format(line[1], line[0]))
                    continue
                cfg[line[0]] = line[1] == 'T'
            elif line[0] == 'itxalign':
                if line[1] not in ['L', 'R', 'C', 'D']:
                    logger.error("Invalid value {} for {} in base.cfg. Valid values: L/R/C/D (left/right/center/default). Defaulting to equidistant.".format(line[1], line[0]))
                    cfg[line[0]] = 'D'
                    continue
                cfg[line[0]] = line[1]
    return cfg
# END



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
        self.mt = '.'
        self.mc = 'black'
        self.ms = 8
        self.mz = 1
        self.tt = ''
        self.tc = 'black'
        self.ts = matplotlib.rcParams['font.size']
        self.tf = 'Arial'
        self.tp = '0.05'
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
                elif n in ['ms', 'ts', 'tp', 'mz']:
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
                logger.warning("Incomplete example in BED file line:\n{}\nGenome coordinate (chromosome, start, end) and genome name are required columns; tags can be in optional column.\nPlease ensure separators are tabs, not spaces.".format('\t'.join(line)))
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
            mdata.append(anno)
    return mdata
# END


def readsyriout(f):
    from pandas import DataFrame
    import numpy as np
    from collections import deque, OrderedDict
    import logging
    # Reads syri.out. Select: achr, astart, aend, bchr, bstart, bend, srtype
    logger = logging.getLogger("readsyriout")
    syri_regs = deque()
    skipvartype = ['CPG', 'CPL', 'DEL', 'DUPAL', 'HDR', 'INS', 'INVAL', 'INVDPAL', 'INVTRAL', 'NOTAL', 'SNP', 'SYNAL', 'TDM', 'TRANSAL']
    logger.info('Reading input files generated by syri.')
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            # TODO: DECIDE WHETHER TO HAVE STATIC VARS OR FLEXIBLE ANNOTATION
            if l[10] in VARS:
                syri_regs.append(l)
            else:
                if l[10] not in skipvartype:
                    skipvartype.append(l[10])
                    logger.warning("{} is not a valid annotation for alignments in file {}. Alignments should belong to the following classes {}. Skipping alignment.".format(l[10], f, VARS))
    try:
        df = DataFrame(list(syri_regs))
    except KeyError:
        raise ImportError("Incomplete input file {}, syri.out file should have 11 columns.".format(f))
    if df.shape[1] == 13:
        logger.warning(f'{f} has an extra annotation column dedicated for alignment customisation. Using it.')
    else:
        logger.info(f'{f} does not have an extra annotation column dedicated for alignment customisation. Using default values.')
        df[12] = '-'
    df = df[[0, 1, 2, 5, 6, 7, 10, 12]]
    colnames = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type', 'anno']
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
    df.columns = colnames
    return df, chrid_dict
# END


def readbedout(f):
    from pandas import DataFrame
    import numpy as np
    from collections import deque, OrderedDict
    import logging
    # BEDPE format: achr, astart, aend, bchr, bstart, bend, srtype
    logger = logging.getLogger('readbedout')
    bed_regs = deque()
    skipvartype = []
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            # TODO: DECIDE WHETHER TO HAVE STATIC VARS OR FLEXIBLE ANNOTATION
            if l[6] in VARS:
                bed_regs.append(l)
            else:
                if l[6] not in skipvartype:
                    skipvartype.append(l[6])
                    logger.warning("{} is not a valid annotation for alignments in file {}. Alignments should belong to the following classes {}. Skipping alignment.".format(l[6], f, VARS))

    df = DataFrame(list(bed_regs))
    if df.shape[1] == 8:
        logger.warning(f'{f} have extra annotation column dedicated for alignment customisation. Using it.')
    else:
        logger.warning(f'{f} does have extra annotation column dedicated for alignment customisation. Using default values.')
        df[7] = '-'
    colnames = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type', 'anno']
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
    df.columns = colnames
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
        self.nm = 0
        self.lc = 'black'
        self.lw = 1
        self.bc = 'lightgrey'
        self.ba = 0.7
        self.logger = logging.getLogger("track")
        self.bincnt = {}
        self.gff = ''
        self.ti = -1
        self.tt = 'f'
        self.ta = 1

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
                        if v[0] == '#':
                            matplotlib.colors.to_rgb(v)
                        else:
                            matplotlib.colors.to_hex(v)
                    except ValueError:
                        self.logger.error("Error in using colour: {} for track {}. Use correct hexadecimal colours or named colours define in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html)".format(v, self.n))
                        # continue
                        sys.exit()
                    setattr(self, n, v)
                ## Integer parameters
                elif n in ['ti']:
                    try: int(v)
                    except ValueError:
                        self.logger.error("Non-integer value {} for {} in track{}".format(v, n, self.n))
                        sys.exit()
                    setattr(self, n, int(v))
                ## float parameters
                elif n in ['ns', 'bw', 'lw', 'ba', 'nm', 'ta']:
                    try: float(v)
                    except ValueError:
                        self.logger.error("Non-numerical value {} for {} in track{}".format(v, n, self.n))
                        sys.exit()
                    setattr(self, n, float(v))
                ## text parameters
                elif n in ['tt']:
                    if n == 'tt':
                        if v not in ['f', 'l']:
                            self.logger.error("Invalid value {} for {} in track{}. Allowed values: 'f' and 'l'.".format(v, n, self.n))
                            sys.exit()
                    setattr(self, n, v)
                ## Font parameters
                elif n in ['nf']:
                    if v not in FONT_NAMES:
                        with open("plotsr_available_font_names.txt", 'w') as fout:
                            fout.write("\n".join(FONT_NAMES))
                        raise ValueError("Font {} in track {} is not available. Check plotsr_available_font_names.txt for list of available system markers".format(v, self.n))
                    setattr(self, n, v)
                ## Input type parameters
                elif n == 'ft':
                    if v not in ['bed', 'bedgraph', 'gff']:
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
        _chrs = set([c for c in chrlengths[0][1].keys()])
        bincnt = defaultdict(deque)
        skipchrs = []
        curchr = ''
        pos = deque()
        added_chrs = list()
        with open(self.f, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if len(line) < 3:
                    self.logger.warning("Incomplete information in BED file at line: {}. Skipping it.".format("\t".join(line)))
                    continue
                if line[0] not in _chrs:
                    if line[0] not in skipchrs:
                        self.logger.info("Chromosome in BED is not present in FASTA or not selected for plotting. Skipping it. BED line: {}".format("\t".join(line)))
                        skipchrs.append(line[0])
                    continue
                if curchr == '':
                    curchr = line[0]
                    pos.append([int(line[1]), int(line[2])])
                elif curchr == line[0]:
                    pos.append([int(line[1]), int(line[2])])
                else:
                    if line[0] in added_chrs:
                        self.logger.error("BED file: {} is not sorted. For plotting tracks, sorted bed file is required. Exiting.".format(self.f))
                        sys.exit()
                    if len(pos) > 1:
                        rngs = mergeranges(np.array(pos))
                    else:
                        rngs = pos
                    chrpos = np.array(list(set([i for r in rngs for i in range(r[0], r[1])])))
                    # Get bin breakpoints for the chromosome
                    bins = np.concatenate((np.arange(0, chrlengths[0][1][curchr], bw), np.array([chrlengths[0][1][curchr]])), axis=0)
                    binval = np.histogram(chrpos, bins)[0]
                    bincnt[curchr] = deque([((bins[i] + bins[i+1])/2, binval[i]/bw) for i in range(len(binval))])
                    added_chrs.append(curchr)
                    # Set the new chromosome
                    curchr = line[0]
                    pos = deque([[int(line[1]), int(line[2])]])
            if len(pos) > 1:
                rngs = mergeranges(np.array(pos))
            else:
                rngs = pos
            chrpos = np.array(list(set([i for r in rngs for i in range(r[0], r[1])])))
            # Get bin breakpoints for the chromosome
            bins = np.concatenate((np.arange(0, chrlengths[0][1][curchr], bw), np.array([chrlengths[0][1][curchr]])), axis=0)
            binval = np.histogram(chrpos, bins)[0]
            bincnt[curchr] = deque([((bins[i] + bins[i+1])/2, binval[i]/bw) for i in range(len(binval))])
        self.bincnt = bincnt
        return
    # END

    # Read input bedgraph file
    def _readbedgraph(self, chrlengths):
        from collections import deque, defaultdict
        import numpy as np
        from math import ceil

        bw = int(self.bw)
        _chrs = set([c for c in chrlengths[0][1].keys()])
        bincnt = defaultdict(deque)
        skipchrs = []
        curchr = ''
        added_chrs = list()
        with open(self.f, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                try:
                    val = float(line[3])
                except ValueError:
                    if len(line) < 4:
                        self.logger.warning("Incomplete information in bedgraph file at line: {}. Skipping it.".format("\t".join(line)))
                        continue
                if curchr == '':
                    curchr = line[0]
                    binv = np.zeros(ceil(chrlengths[0][1][curchr] / bw), dtype=int)
                    s = deque([int(line[1])])
                    e = deque([int(line[2])])
                    v = deque([int(val)])
                elif curchr == line[0]:
                    s.append(int(line[1]))
                    e.append(int(line[2]))
                    v.append(val)
                elif line[0] not in _chrs:
                    if line[0] == '#': continue
                    if line[0] == 'track': continue
                    if line[0] not in skipchrs:
                        self.logger.warning("Chromosome in BEDGRAPH is not present in FASTA or not selected for plotting. Skipping it. BED line: {}".format("\t".join(line)))
                        skipchrs.append(line[0])
                    continue
                else:
                    if line[0] in added_chrs:
                        self.logger.error("BedGraph file: {} is not sorted. For plotting tracks, sorted BedGraph file is required. Exiting.".format(self.f))
                        sys.exit()

                    s = np.array(list(s))
                    e = np.array(list(e))
                    v = np.array(list(v))
                    spos = s//bw
                    epos = e//bw
                    ## Update bin value for rows which are in one BW
                    idx = np.where(np.equal(spos, epos))[0]
                    values = (e-s)*v
                    np.add.at(binv, spos[idx], values[idx])
                    ## Update bin value for rows which are in more than one BW
                    idx = np.where(np.not_equal(spos, epos))[0]
                    values_s = (bw-(s % bw))*v
                    values_e = (e % bw)*v
                    np.add.at(binv, spos[idx], values_s[idx])
                    np.add.at(binv, epos[idx], values_e[idx])
                    ### Update bin value for bins overlapped by rows
                    idx2 = np.where(np.not_equal(spos[idx]+1, epos[idx]))[0]
                    idx_to_add = [i for j in idx2 for i in range(spos[idx[j]]+1, epos[idx[j]])]
                    value_to_add = [v[idx[j]]*bw for j in idx2 for i in range(spos[idx[j]]+1, epos[idx[j]])]
                    np.add.at(binv, idx_to_add, value_to_add)

                    bins = np.concatenate((np.arange(0, chrlengths[0][1][curchr], bw), np.array([chrlengths[0][1][curchr]])), axis=0)
                    bins = [(bins[i] + bins[i+1])/2 for i in range(len(bins) - 1)]
                    bincnt[curchr] = deque([(bins[i], binv[i]) for i in range(len(bins))])
                    added_chrs.append(curchr)
                    # Set the new chromosome
                    curchr = line[0]
                    binv = np.zeros(ceil(chrlengths[0][1][curchr]/bw), dtype=int)
                    s = deque([int(line[1])])
                    e = deque([int(line[2])])
                    v = deque([int(val)])

        s = np.array(list(s))
        e = np.array(list(e))
        v = np.array(list(v))
        spos = s//bw
        epos = e//bw
        ## Update bin value for rows which are in one BW
        idx = np.where(np.equal(spos, epos))[0]
        values = (e-s)*v
        np.add.at(binv, spos[idx], values[idx])
        ## Update bin value for rows which are in more than one BW
        idx = np.where(np.not_equal(spos, epos))[0]
        values_s = (bw-(s % bw))*v
        values_e = (e % bw)*v
        np.add.at(binv, spos[idx], values_s[idx])
        np.add.at(binv, epos[idx], values_e[idx])
        ### Update bin value for bins overlapped by rows
        idx2 = np.where(np.not_equal(spos[idx]+1, epos[idx]))[0]
        idx_to_add = [i for j in idx2 for i in range(spos[idx[j]]+1, epos[idx[j]])]
        value_to_add = [v[idx[j]]*bw for j in idx2 for i in range(spos[idx[j]]+1, epos[idx[j]])]
        np.add.at(binv, idx_to_add, value_to_add)

        bins = np.concatenate((np.arange(0, chrlengths[0][1][curchr], bw), np.array([chrlengths[0][1][curchr]])), axis=0)
        bins = [(bins[i] + bins[i+1])/2 for i in range(len(bins) - 1)]
        bincnt[curchr] = deque([(bins[i], binv[i]) for i in range(len(bins))])
        ## Scale count values
        maxv = 0
        for k, v in bincnt.items():
            for r in v:
                if r[1] > maxv:
                    maxv = r[1]
        for k, v in bincnt.items():
            bincnt[k] = deque([(r[0], r[1]/maxv) for r in v])
        self.bincnt = bincnt
        return
    # END

    def _readgff(self, chrlengths):
        from collections import defaultdict, deque
        import sys
        annos = defaultdict(dict)
        acceptlist = {'mrna', 'cds'}
        featwarn = True
        skipchr = []
        chrs = set(chrlengths[0][1].keys())
        with open(self.f, 'r') as fin:
            self.logger.warning("Reading GFF file {}. Overlapping transcripts would be plotted as such without any filtering.".format(self.f))
            model = None
            for line in fin:
                if line.strip() == '': continue
                line = line.strip().split()
                if line[0] not in chrs:
                    if line[0] not in skipchr:
                        skipchr.append(line[0])
                        self.logger.info("Chromosome in GFF is not present in FASTA or not selected for plotting. Skipping it. GFF line: {}".format("\t".join(line)))
                    continue
                t = line[2].lower()
                if t not in acceptlist:
                    if featwarn:
                        self.logger.info("GFF feature: {} is not usable for plotting. Only {} are used for plotting. Will skipping all other features.".format(line[2], acceptlist))
                        featwarn = False
                else:
                    try:
                        int(line[3]), int(line[4])
                    except ValueError:
                        self.logger.error("Non-numerical value in GFF {} at line {}. Exiting.".format(self.f, "\t".join(line)))
                        sys.exit()

                    if model is not None:
                        if t == 'mrna':
                            k = list(model.keys())[0]
                            annos[k[0]][(k[1], k[2])] = model.values()
                            model = defaultdict()
                            cse = (line[0], int(line[3]), int(line[4])) #Chromosome, start, end
                            model[cse] = deque()
                        if t == 'cds':
                            model[cse].append((int(line[3]), int(line[4])))
                    elif model is None:
                        if t == 'cds':
                            self.logger.error("In GFF: {} at line {}, CDS feature is before the mRNA feature. Expected mRNA before CDS. Exiting.".format(self.f, line))
                            sys.exit()
                        else:
                            model = defaultdict()
                            cse = (line[0], int(line[3]), int(line[4]))
                            model[cse] = deque()
            k = list(model.keys())[0]
            annos[k[0]][(k[1], k[2])] = model.values()
        self.gff = annos
    # END

    def readdata(self, chrlengths):
        if self.ft == 'bed':
            self._readbed(chrlengths)
        elif self.ft == 'bedgraph':
            self._readbedgraph(chrlengths)
        elif self.ft == 'gff':
            self._readgff(chrlengths)
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
                raise ValueError("Incomplete line in track file.\n{}\nFile path and name are necessary columns, tags is optional column".format('\t'.join(line)))
            t = track(line[0], line[1])
            # Read tags
            if len(line) == 3:
                t.addtags(line[2])
            # Reads BED
            t.readdata(chrlengths)
            tdata.append(t)
    # add track index for tracks without index
    ## Get list of used index
    tiused = set([t.ti for t in tdata])
    tiused = list(tiused.difference(set([-1])))
    ## Set tracks names
    for ti in tiused:
        tnames = [t.n for t in tdata if t.ti == ti]
        tn = ','.join(tnames)
        first = True
        for t in tdata:
            if t.ti == ti:
                if first:
                    setattr(t, 'n', tn)
                    first = False
                else:
                    setattr(t, 'n', '')
    timax = max(tiused) if tiused != [] else 0
    if any([i not in tiused for i in range(1, timax+1)]):
        raise ValueError("Used track indices are not sequential. The input track indices should be a consecutive list of numbers.")
    for t in tdata:
        if t.ti == -1:
            setattr(t, 'ti', timax+1)
            timax += 1
    return tdata
# END


def readagp(f):
    """
    Reads an .agp file and returns the contig-breakpoints in the scaffolds
    :param f: Input file name
    :return: Dictionary with scaffold ids as keys and breakpoints as values
    """
    from collections import defaultdict
    outdict = defaultdict(list)
    with open(f, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            line = line.strip().split()
            if line[4] not in 'W':
                continue
            outdict[line[0]].append([int(line[1]), int(line[2])])
    return outdict
# END


def readbedasdict(f):
    """
    Reads a bed file and return a dict where the keys are chromosomes in bed and values are ranges
    :param f: Filename
    :return: Bed as dict
    """
    from collections import defaultdict, deque
    beddict = defaultdict(deque)
    with open(f, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            #TODO: ADD checks to ensure the correctness of BED file
            beddict[line[0]].append([line[1], line[2]])
    return beddict
# END

# </editor-fold>

"""
################################################################################
Validation and filtering
################################################################################
"""

class genome():
    def __init__(self, f, n, c):
        import logging
        self.f = f
        self.n = n
        self.lc = c
        self.lw = 1
        self.ft = 'fa'          # Input file type: "fa"=fasta, "cl"="chromosome length
        self.glen = None
        self.logger = logging.getLogger("genome")

    # Add tags
    def addtags(self, tags):
        import matplotlib
        import sys
        for i in tags.split(";"):
            n, v = i.split(":")
            if hasattr(self, n):
                ## Colour parameters
                if n in ['lc']:
                    try:
                        if v[0] == '#':
                            matplotlib.colors.to_rgb(v)
                        else:
                            matplotlib.colors.to_hex(v)
                    except ValueError:
                        self.logger.error("Error in using colour: {} for genome {}. Use correct hexadecimal colours or named colours define in matplotlib (https://matplotlib.org/stable/gallery/color/named_colors.html)".format(v, self.n))
                        sys.exit()
                    setattr(self, n, v)
                ## Numerical parameters
                elif n in ['lw']:
                    try:
                        float(v)
                    except ValueError:
                        self.logger.error("Non-numerical value {} for {} for genome {}".format(v, n, self.n))
                        sys.exit()
                    setattr(self, n, float(v))
                ## Input type parameters
                elif n == 'ft':
                    if v not in ['fa', 'cl']:
                        self.logger.error("Unknown filetype specified for genome {}. Please provide input type as 'fa' for fasta or 'cl' for chromosome-lengh in tab-separated file. 'cl' files are faster to read.".format(v, self.n))
                        sys.exit()
                    setattr(self, n, v)
            else:
                self.logger.error("{} is not a valid tag".format(n))
                sys.exit()
        return
    #END

    # Read chromosome lengths
    def readdata(self):
        import sys
        if self.ft == 'fa':
            try:
                self.logger.debug("Reading fasta file: {}".format(self.f))
                self.glen = {c: len(seq) for c, seq in readfasta(self.f).items()}
            except Exception as e:
                raise ImportError("Error in reading fasta: {}\n{}".format(self.n, e))
        elif self.ft == 'cl':
            glen = {}
            self.logger.debug("Reading chromosome length file: {}".format(self.f))
            with open(self.f, 'r') as fin:
                for line in fin:
                    line = line.strip().split()
                    if len(line) == 0: continue
                    try:
                        glen[line[0]] = int(line[1])
                    except IndexError:
                        self.logger.error("Incomplete input. Genome:{} chromosome:{} has no chromosome length. Exiting.".format(self.n, line[0]))
                        sys.exit()
                    except ValueError:
                        self.logger.error("Incorrect input. Genome:{} chromosome:{} has no-numerical chromosome length. Exiting.".format(self.n, line[0]))
                        sys.exit()
            self.glen = glen
# END


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
    errmess0 = "Cannot read annotations for genome: {}. Make sure that structural annotations for all genomes are provided in the same order as genomes. Exiting."
    errmess1 = 'Chromosome ID: {} in structural annotation file: {} not present in genome fasta: {}. Exiting.'
    errmess2 = 'For chromosome ID: {}, length in genome fasta: {} is less than the maximum coordinate in the structural annotation file: {}. Exiting.'
    tags = {'lc': {}, 'lw': {}}

    # Count number of genomes and set automatic colors
    count = 0
    with open(genf, 'r') as fin:
        i = 0
        for line in fin:
            if line[0] == '#':
                continue
            count += 1
    if count <= 10:
        CHRCOLS = [matplotlib.colors.to_hex(get_cmap('tab10')(i)) for i in range(count)]
    else:
        CHRCOLS = [matplotlib.colors.to_hex(get_cmap('gist_rainbow')(int(255/count) * i)) for i in range(0, count)]
        if count % 2 != 0:
            m = CHRCOLS[int((count/2) + 1)]
            CHRCOLS = [j for i in range(int(count/2)) for j in [CHRCOLS[i]] + [CHRCOLS[int(i +count/2)]]] + [m]
        else:
            CHRCOLS = [j for i in range(int(count/2)) for j in [CHRCOLS[i]] + [CHRCOLS[int(i +count/2)]]]
    # Read genomes and validate chrlengths
    with open(genf, 'r') as fin:
        i = 0
        genomes = deque()
        for line in fin:
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            line = line.strip().split("\t")
            if len(line) == 1 and not line[0]:
                continue
            if len(line) < 2:
                raise ImportError("Incomplete genomic information.\nExpected format for the genome file:\npath_to_genome1\tgenome1_id\ttags\npath_to_genome2\tgenome2_id\ttags\n\nMake sure that the columns are separated by tabs (and not spaces).")
            gen = genome(line[0], line[1], CHRCOLS[i])
            # Read tags
            if len(line) == 3:
                gen.addtags(line[2])
            # Reads fasta
            gen.readdata()
            # Validate chromosome lengths
            glen = gen.glen
            achr = set()
            bchr = set()
            if i > 0:
                try:
                    df = als[i-1][1]
                except IndexError:
                    raise ImportError(errmess0.format(line[0]))
                bchr = np.unique(df['bchr'])
                for c in bchr:
                    if c not in list(glen.keys()):
                        raise ImportError(errmess1.format(c, als[i-1][0], os.path.basename(line[1])))
                    if df.loc[df['bchr'] == c, ['bstart', 'bend']].max().max() > glen[c]:
                        raise ImportError(errmess2.format(c, os.path.basename(genf), als[i-1][0]))
            # Check cases when the genome is the reference-genome
            if i < len(als):
                try:
                    df = als[i][1]
                except IndexError:
                    raise ImportError(errmess0.format(line[0]))
                achr = np.unique(df['achr'])
                for c in achr:
                    if c not in list(glen.keys()):
                        raise ImportError(errmess1.format(c, als[i][0], os.path.basename(line[1])))
                    if df.loc[df['achr'] == c, ['astart', 'aend']].max().max() > glen[c]:
                        raise ImportError(errmess2.format(c, os.path.basename(genf), als[i][0]))
            out.append((line[1], {c: glen[c] for c in set(achr).union(set(bchr))}))
            i += 1
            genomes.append(gen)
    return out, genomes
# END


def filterinput(args, df, chrid, itx=False):
    # Get region length and filter out smaller SR
    df = df.loc[((df['aend'] - df['astart']) >= args.s) | ((df['bend'] - df['bstart']) >= args.s) | (df['type'] == 'SYN')].copy()
    if not itx:
        df = df.loc[df['bchr'] == [chrid[i] for i in df['achr']]]
    # Filter non-selected variations
    if args.nosyn:
        df = df.loc[df['type'] != 'SYN']
    if args.noinv:
        df = df.loc[df['type'] != 'INV']
    if args.notr:
        df = df.loc[~df['type'].isin(['TRANS', 'INVTR'])]
    if args.nodup:
        df = df.loc[~df['type'].isin(['DUP', 'INVDP'])]
    df.sort_values(['bchr', 'bstart', 'bend'], inplace=True)
    df.sort_values(['achr', 'astart', 'aend'], inplace=True)
    return df
# END


def selectchrom(CHRS, cs, chrgrps, alignments, chrlengths, chrids):
    import logging
    from collections import deque, OrderedDict
    homchrs = deque()
    chrs = deque()
    logger = logging.getLogger("Selecting chromosomes")
    for c in CHRS:
        if c not in cs:
            logger.warning("Selected chromosome: {} is not in reference genome. Skipping it.".format(c))
            continue
        homchrs.append(chrgrps[c])
        chrs.append(c)
    for i in range(len(alignments)):
        alignments[i][1] = alignments[i][1].loc[(alignments[i][1]['achr'].isin([h[i] for h in homchrs])) &
                                                (alignments[i][1]['bchr'].isin([h[i+1] for h in homchrs]))]
    ## Remove chromosomes that are not homologous to selected reference chromosomes
    if CHRS is not None:
        for i in range(len(chrlengths)):
            ks = list(chrlengths[i][1].keys())
            homs = [h[i] for h in homchrs]
            for k in ks:
                if k not in homs:
                    chrlengths[i][1].pop(k)
        # Update groups of homologous chromosomes
        # chrs = [k for k in chrids[0][1].keys() if k in alignments[0][1]['achr'].unique()]
        chrgrps = OrderedDict()
        for c in chrs:
            cg = deque([c])
            cur = c
            for i in range(len(chrids)):
                n = chrids[i][1][cur]
                cg.append(n)
                cur = n
            chrgrps[c] = cg
    return alignments, chrs, chrgrps, chrlengths
# END


def selectregion(reg, rtr, chrlengths, al, chrids):
    import pandas as pd
    import copy
    from collections import OrderedDict, deque
    alignments = copy.deepcopy(al)
    genids = [i[0] for i in chrlengths]
    if len(reg) != 3:
        raise ValueError("Incorrect values parsed to --reg. Provide values in the following format: GenomeID:ChromosomeID:Start-End. GenomeID and ChromosomeID cannot have ':' character.")
    reg = reg[:2] + list(map(int, reg[2].split("-")))
    # Check that the region is in the given genomes
    if reg[0] not in genids:
        raise ValueError("Genome ID in --reg do not match with any genome ID provided with --genomes. Exiting.")
    if reg[1] not in chrlengths[genids.index(reg[0])][1].keys():
        raise ValueError("Chromosome ID in --reg do not match with any chromosome ID for genome {}. Exiting.".format(reg[0]))
    # Get length of thetarget chromosome
    ind = genids.index(reg[0])
    l = chrlengths[ind][1][reg[1]]
    if any([reg[2] < 1, reg[2] > l, reg[3] < 1, reg[3] > l, reg[2] > reg[3]]):
        raise ValueError("Incorrect chromosome coordinates provided for --reg. Exiting.")
    # Select alignments that are within the selected region in the focal genome. For other genomes, maximal syntenic region coordinate would be used as boundaries.
    newal = [0]*len(alignments)
    # Alignments for genomes before the focal genome
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
        newal[i] = syncrd.copy()
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
        newal[i] = syncrd.copy()
        c, s, e = bc, bs, be
    # Update alignments, chrs and chrgrps
    if not rtr:
        for i in range(len(alignments)):
            achr = newal[i]['achr'].tolist()[0]
            astart = min(newal[i]['astart'])
            aend = max(newal[i]['aend'])
            bchr = newal[i]['bchr'].tolist()[0]
            bstart = min(newal[i]['bstart'])
            bend = max(newal[i]['bend'])
            srcrd = alignments[i][1].loc[(alignments[i][1]['achr'] == achr) &
                                         (alignments[i][1]['astart'] >= astart) &
                                         (alignments[i][1]['aend'] <= aend) &
                                         (alignments[i][1]['bchr'] == bchr) &
                                         (alignments[i][1]['bstart'] >= bstart) &
                                         (alignments[i][1]['bend'] <= bend) &
                                         (alignments[i][1]['type'] != 'SYN')].reset_index(drop=True).copy()
            tmp = pd.concat([newal[i], srcrd])
            tmp.sort_values(['achr', 'astart', 'aend'], inplace=True)
            tmp.reset_index(inplace=False, drop=True)
            alignments[i][1] = tmp
    else:
        allal = pd.concat(newal)[['astart', 'aend', 'bstart', 'bend']]
        minx = min([reg[2]] + allal.apply(min).tolist())
        maxx = max([reg[3]] + allal.apply(max).tolist())
        for i in range(len(alignments)):
            ac = list(set(newal[i]['achr']))[0]
            bc = list(set(newal[i]['bchr']))[0]
            # Get SR alignments within the view
            srcrd = alignments[i][1].loc[(alignments[i][1]['achr'] == ac) &
                                         (alignments[i][1]['astart'] >= minx) &
                                         (alignments[i][1]['aend'] <= maxx) &
                                         (alignments[i][1]['bchr'] == bc) &
                                         (alignments[i][1]['bstart'] >= minx) &
                                         (alignments[i][1]['bend'] <= maxx) &
                                         (alignments[i][1]['type'] != 'SYN')].reset_index(drop=True).copy()
            # Get syntenic regions within the view
            syncrd = alignments[i][1].loc[(alignments[i][1]['achr'] == ac) &
                                          (alignments[i][1]['astart'] <= maxx) &
                                          (alignments[i][1]['aend'] >= minx) &
                                          (alignments[i][1]['bchr'] == bc) &
                                          (alignments[i][1]['bstart'] <= maxx) &
                                          (alignments[i][1]['bend'] >= minx) &
                                          (alignments[i][1]['type'] == 'SYN')].reset_index(drop=True).copy()
            # Alter syn alignments to fit within the selected view
            ## Trim astart
            if min(syncrd['astart']) < minx:
                lower = syncrd['astart'] < minx
                a = (syncrd['aend'][lower] - minx)/(syncrd['aend'][lower] - syncrd['astart'][lower])
                syncrd.loc[lower, 'bstart'] = (syncrd['bend'][lower] - (syncrd['bend'][lower] - syncrd['bstart'][lower])*a).astype(int)
                syncrd.loc[lower, 'astart'] = minx
            ## Trim aend
            if max(syncrd['aend']) > maxx:
                higher = syncrd['aend'] > maxx
                a = (maxx - syncrd['astart'][higher])/(syncrd['aend'][higher] - syncrd['astart'][higher])
                syncrd.loc[higher, 'bend'] = (syncrd['bstart'][higher] + (syncrd['bend'][higher] - syncrd['bstart'][higher])*a).astype(int)
                syncrd.loc[higher, 'aend'] = maxx
            ## Trim bstart
            if min(syncrd['bstart']) < minx:
                lower = syncrd['bstart'] < minx
                a = (syncrd['bend'][lower] - minx)/(syncrd['bend'][lower] - syncrd['bstart'][lower])
                syncrd.loc[lower, 'astart'] = (syncrd['aend'][lower] - (syncrd['aend'][lower] - syncrd['astart'][lower])*a).astype(int)
                syncrd.loc[lower, 'bstart'] = minx
            ## Trim bend
            if max(syncrd['bend']) > maxx:
                higher = syncrd['bend'] > maxx
                a = (maxx-syncrd['bstart'][higher])/(syncrd['bend'][higher] - syncrd['bstart'][higher])
                syncrd.loc[higher, 'aend'] = (syncrd['astart'][higher] + (syncrd['aend'][higher] - syncrd['astart'][higher])*a).astype(int)
                syncrd.loc[higher, 'bend'] = maxx


            tmp = pd.concat([syncrd, srcrd])
            tmp.sort_values(['achr', 'astart', 'aend'], inplace=True)
            tmp.reset_index(inplace=False, drop=True)
            alignments[i][1] = tmp
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
# END


def createribbon(df):
    """
    Combine continuous syntenic regions to get larger ribbons for syntenic blocks
    """
    import numpy as np
    from pandas import DataFrame, concat
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
    chrchange = np.where((np.array(df['achr'].iloc[1:]) == np.array(df['achr'].iloc[0:-1])) == False)[0] + 1
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
    if len(cg) > 0:
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
    newsyn['anno'] = '-'
    newsyn['type'] = 'SYN'

    df = df.drop(['a', 'b'], axis=1)
    df = df.loc[-(df['type'] == 'SYN')]
    df = concat([df, newsyn])
    df.sort_values(['achr', 'astart', 'aend'], inplace=True)
    return df
# END

################################################################################
################################################################################
"""
################################################################################
Draw and plot
################################################################################
"""


# for itx-mode, generate dict of starting coordinates of chromosome boxes per genome
def genbuff(chrlengths, chrgrps, chrs, maxl, v, cfg):
    # read from config file
    mchr = cfg['marginchr']
    # declare output dict
    chr_start_coord = {}
    for s, chrl in enumerate(chrlengths):
        # add genome_name as key
        chr_start_coord[chrl[0]] = {}
        # initial start coordinate = 0
        start = 0
        if not v:
            chr_iterate_order = list(chrs)
        else:#elif v
            chr_iterate_order = list(chrs)[::-1]
        for cid in chr_iterate_order:
            longer_chrom_length = max([c[1][chrgrps[cid][j]] for j, c in enumerate(chrlengths)])
            self_length = chrl[1][chrgrps[cid][s]]
            # four itx align modes
            if cfg['itxalign'] == 'L':
                end = start + longer_chrom_length
                chr_start_coord[chrl[0]][cid] = start
            elif cfg['itxalign'] == 'R':
                padding = longer_chrom_length - self_length
                end = start + longer_chrom_length
                chr_start_coord[chrl[0]][cid] = start + padding
            elif cfg['itxalign'] == 'C':
                padding = (1/2) * (longer_chrom_length - self_length)
                end = start + longer_chrom_length
                chr_start_coord[chrl[0]][cid] = start + padding
            else:# defaulting to cfg['itxalign'] == 'D'
                end = start + self_length
                chr_start_coord[chrl[0]][cid] = start
            # update start
            start = end + mchr * maxl
    return chr_start_coord
# END


def drawax(ax, chrgrps, chrlengths, v, s, cfg, itx, chr_start_coord, minl=0, maxl=-1, chrname=None):
    import numpy as np
    from collections import deque
    bottom_limit = -cfg['exmar']
    upper_limit = cfg['exmar']
    if chrname is not None:
        chrnamedict = {}
        with open(chrname, 'r') as fin:
            for line in fin:
                line = line.strip().split('\t')
                if len(line) != 2:
                    raise ImportError("Error in reading --chrname. Each row needs to have two columns. Error in line: " + " ".join(line) + ". Exiting.")
                chrnamedict[line[0]] = line[1]
        try:
            _ = [chrnamedict[k] for k in chrgrps.keys()]
        except KeyError:
            raise ImportError("Error in reading --chrname. Chromosome IDs do not match. Exiting.")
    else:
        chrnamedict = {k: k for k in chrgrps.keys()}
    if not itx:
        ticklabels = [chrnamedict[k] for k in chrgrps.keys()]
        nchr = len(chrgrps)
        if maxl == -1:
            maxl = np.max([chrlengths[i][1][c[i]] for c in chrgrps.values() for i in range(len(c))])
        if not v:
            tick_pos = s/2 + cfg['chrmar']
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
            else:
                xticksl = xticks
                ax.set_xlabel('Chromosome position')
            ax.set_xticks(xticks[:-1])
            ax.set_xticklabels(xticksl[:-1])
            ax.set_ylabel('Reference Chromosome ID')
            ax.set_axisbelow(True)
        else:
            tick_pos = 1 - cfg['chrmar'] - s/2
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
            else:
                yticksl = yticks
                ax.set_ylabel('Chromosome position')
            ax.set_yticks(yticks[:-1])
            ax.set_yticklabels(yticksl[:-1])
            ax.set_xlabel('Reference Chromosome ID')
            ax.yaxis.grid(True, which='both', linestyle='--')
            ax.set_axisbelow(True)
    elif itx:
        step = s/(len(chrlengths)-1)
        if not v:
            ax.set_xlim(0, maxl)
            ax.set_ylim(bottom_limit, 1+upper_limit)
            tick_pos = deque()
            tick_lab = deque()
            # get chrbox start coord of the last genome (plotted at bottom, closest to ticks)
            offset = chr_start_coord[chrlengths[-1][0]]
            for chr in chrgrps:
                # tick position is offset plus half of last genome chrbox length
                tick_pos.append(offset[chr] + (chrlengths[-1][1][chr])/2)
                tick_lab.append(chrnamedict[chr])
            ax.set_xticks(tick_pos)
            ax.set_xticklabels(tick_lab)
            step = s/(len(chrlengths)-1)
            ax.set_yticks([round(step*i, 2) for i in range(len(chrlengths))])
            ax.set_yticklabels([i[0] for i in chrlengths][::-1])
            ax.tick_params(axis='y', right=False, left=False)
            ax.set_axisbelow(True)
            ax.set_xlabel('Reference Chromosome ID')
            ax.set_ylabel('Genome')
        else:
            ax.set_ylim(0, maxl)
            ax.set_xlim(bottom_limit, 1+upper_limit)
            tick_pos = deque()
            tick_lab = deque()
            # get chrbox start coord of the first genome (plotted leftmost, closest to ticks)
            offset = chr_start_coord[chrlengths[0][0]]
            for chr in chrgrps.__reversed__():
                # tick position is offset plus half of first/leftmost genome chrbox length
                tick_pos.append(offset[chr] + (chrlengths[0][1][chr])/2)
                tick_lab.append(chrnamedict[chr])
            ax.set_yticks(tick_pos)
            ax.set_yticklabels(tick_lab)
            ax.set_xticks([round((1-s) + step*i, 2) for i in range(len(chrlengths))])
            ax.set_xticklabels([i[0] for i in chrlengths])
            ax.tick_params(axis='x', right=False, left=False)
            ax.set_axisbelow(True)
            ax.set_ylabel('Reference Chromosome ID')
            ax.set_xlabel('Genome')
    return ax
# END


def pltchrbox(ax, v, loc, col, cw, maxl, minl=0, lab=None):
    from matplotlib.patches import Rectangle
    from matplotlib import colors as mc
    chrlabel = None
    col = mc.rgb_to_hsv(mc.to_rgb(col))
    col1 = mc.to_hex(mc.hsv_to_rgb((col[0], 1, col[2])))
    col2 = mc.to_hex(mc.hsv_to_rgb((col[0], 0.75, col[2])))
    if lab is not None:
        chrlabel = ax.hlines(0, 0, 0, color=col1, linewidth=1, label=lab, zorder=2)
    if not v:
        ax.add_patch(Rectangle((minl+1, loc-cw/2), width=maxl, height=cw,  linewidth=0, facecolor=col1, alpha=1, zorder=2))
    else:
        ax.add_patch(Rectangle((loc-cw/2, minl+1), width=cw, height=maxl, linewidth=0, facecolor=col1, alpha=1, zorder=2))
    '''
    # Add track for AGP. Currently overdrawn on the chromsome
    ax.add_patch(Rectangle((minl+1, loc-(cw/2)), maxv, cw,  linewidth=0, facecolor=col1, alpha=1, zorder=2))
    for i, v in enumerate(regions[1::2]):
        ax.add_patch(Rectangle((int(v[0])+minl, loc-cw), int(v[1]) - int(v[0]), cw*2,  linewidth=0, facecolor=col2, alpha=1, zorder=2))
    # Add centromere track
    for v in centroms:
        ax.add_patch(Rectangle((int(v[0])+minl, loc-(cw/4)), int(v[1]) - int(v[0]), cw/2,  linewidth=0, facecolor='black', alpha=1, zorder=2))
    '''
    return ax, chrlabel


def pltchrom(ax, chrs, chrgrps, chrlengths, v, S, genomes, cfg, itx, chr_start_coord, minl=0):
    chrlabs = [False]*len(chrlengths)
    # Set chromosome direction
    pltchr = ax.hlines if not v else ax.vlines
    chrlabels = []
    indents = []
    chrwidth = 0.025
    '''
    # read agp data
    agpdata = dict()
    with open('agps.txt', 'r') as fin:
        for line in fin:
            line = line.strip().split()
            agpdata[line[1]] = readagp(line[0])
    # Read centromere data
    centrodata = dict()
    with open('centromeres.txt', 'r') as fin:
        for line in fin:
            line = line.strip().split()
            centrodata[line[1]] = readbedasdict(line[0])
    '''
    if not itx:
        # define indents
        step = S/(len(chrlengths)-1)
        if not v:
            rend = len(chrs)-1+S+cfg['chrmar']
            indents = [rend - (i*step) for i in range(len(chrlengths))]
        elif v:
            rend = 1-S-cfg['chrmar']
            indents = [rend + (i*step) for i in range(len(chrlengths))]
        for s in range(len(chrlengths)):
            for i in range(len(chrs)):
                offset = i if not v else -i
                genome = [gen for gen in genomes if gen.n == chrlengths[s][0]][0]
                chrlength = chrlengths[s][1][chrgrps[chrs[i]][s]]
                if not chrlabs[s]:
                    #ax, c = pltchrbox(ax, agpdata[genome.n][chrgrps[chrs[i]][s]], centrodata[genome.n][chrgrps[chrs[i]][s]], indents[s]-offset, genome.lc, chrwidth, minl=minl, lab=chrlengths[s][0])
                    ax, c = pltchrbox(ax, v, indents[s] - offset, genome.lc, chrwidth, maxl=chrlength, minl=minl, lab=chrlengths[s][0])
                    chrlabels.append(c)
                    chrlabs[s] = True
                else:
                    #ax, c = pltchrbox(ax, agpdata[genome.n][chrgrps[chrs[i]][s]], centrodata[genome.n][chrgrps[chrs[i]][s]], indents[s]-offset, genome.lc, chrwidth, minl=minl)
                    ax, c = pltchrbox(ax, v, indents[s] - offset, genome.lc, chrwidth, maxl=chrlength, minl=minl)
    elif itx:
        step = S/(len(chrlengths)-1)
        for s, chrl in enumerate(chrlengths):
            # format of chrl: ('genome_name', {'chr1': 000000, 'chr2': 000000, ...})
            genome = [gen for gen in genomes if gen.n == chrl[0]][0]
            # starting loc of the chrbox being plotted: used for x-coord if not v, for y-coord if v; is fixed when itx mode because chrboxes are plotted in one line
            fixed = S - (step*s) if not v else 1 - S + (step*s)
            for cid in chr_start_coord[chrl[0]]:
                # chrl[0] gets the genome name; chr_start_coord[chrl[0]] gets the dict of chromosome start coords, where the key/cid is chromosome name
                chrlength = chrl[1][chrgrps[cid][s]]
                # chrlength is the length of the chromosome being plotted by the function below
                # minl gets the chromosome starting coord: {dict of chromosome start coords}[cid]
                ax, c = pltchrbox(ax, v, fixed, genome.lc, chrwidth, maxl=chrlength, minl=chr_start_coord[chrl[0]][cid])
    return ax, indents, chrlabels
# END


def pltsv(ax, alignments, chrs, v, chrgrps, chrlengths, indents, S, cfg, itx, chr_start_coord):
    from collections import deque
    from copy import deepcopy
    def annotodict(anno):
        return {a.split(':')[0]: a.split(':')[1] for a in anno.split(';')}
    #END
    alpha = cfg['alpha']
    adsynlab = False
    adinvlab = False
    adtralab = False
    adduplab = False
    svlabels = dict()
    legenddict = {'SYN': adsynlab, 'INV': adinvlab, 'TRANS': adtralab, 'DUP': adduplab}
    for s in range(len(alignments)):
        df = deepcopy(alignments[s][1])
        df.loc[df['type'] == 'INVTR', 'type'] = 'TRANS'
        df.loc[df['type'] == 'INVDP', 'type'] = 'DUP'
        coldict = {'SYN': cfg['syncol'],
                   'INV': cfg['invcol'],
                   'TRANS': cfg['tracol'],
                   'DUP': cfg['dupcol']}
        lwddict = {'SYN': cfg['synlwd'],
                   'INV': cfg['invlwd'],
                   'TRANS': cfg['tralwd'],
                   'DUP': cfg['duplwd']}
        # df['col'] = [coldict[c] for c in df['type']]
        labdict = {'SYN': 'Syntenic', 'INV': 'Inversion', 'TRANS': 'Translocation', 'DUP': 'Duplication'}
        df['lab'] = [labdict[c] for c in df['type']]
        df.loc[df.duplicated(['lab']), 'lab'] = ''
        newcol = deque()
        newlw = deque()
        newz = deque()
        for row in df.itertuples():
            if row.anno == '-':
                newcol.append(coldict[row.type])
                newlw.append(lwddict[row.type])
                newz.append(0 if row.type == 'SYN' else 1)
            else:
                anno = annotodict(row.anno)
                try:
                    newcol.append(anno['cl'])
                except KeyError:
                    newcol.append(coldict[row.type])
                try:
                    newlw.append(float(anno['lw']))
                except KeyError:
                    newlw.append(0 if row.type == 'SYN' else 0.1)
                try:
                    newz.append(float(anno['z']))
                except KeyError:
                    newz.append(0 if row.type == 'SYN' else 1)
        df['col'] = list(newcol)
        df['lw'] = list(newlw)
        df['zorder'] = list(newz)
        # df['zorder'] = 0
        # df.loc[df['type'] != 'SYN', 'zorder'] = 1
        if not itx:
            df['ry'] = indents[s]
            df['qy'] = indents[s+1]
            for i in range(len(chrs)):
                offset = i if not v else -i
                df.loc[df['achr'] == chrgrps[chrs[i]][s], 'ry'] -= offset
                df.loc[df['achr'] == chrgrps[chrs[i]][s], 'qy'] -= offset
            for row in df.itertuples(index=False):
                p = bezierpath(row.astart, row.aend, row.bstart, row.bend, row.ry, row.qy, v, row.col, alpha=alpha, label=row.lab, lw=row.lw, zorder=row.zorder)
                l = ax.add_patch(p)
                if row.lab != '':
                    if not legenddict[row.type]:
                        svlabels[row.type] = l
                        legenddict[row.type] = True
        elif itx:
            step = S/(len(chrlengths)-1)
            S - (step*s) if not v else 1 - S + (step*s)
            rbuff = chr_start_coord[chrlengths[s][0]]
            rbuff = [rbuff[c] for c in df['achr']]
            df['astart'] += rbuff
            df['aend'] += rbuff
            qbuff = chr_start_coord[chrlengths[s+1][0]]
            qbuff = [qbuff[c] for c in df['bchr']]
            df['bstart'] += qbuff
            df['bend'] += qbuff
            df['ry'] = S - (step*s) if not v else 1 - S + (step*s) # len(chrlengths) - s - 0.5 if not v else s + 0.5
            df['qy'] = S - (step*(s+1)) if not v else 1 - S + (step*(s+1)) # len(chrlengths) - s - 1 - 0.5 if not v else s + 1 + 0.5
            for row in df.itertuples(index=False):
                p = bezierpath(row.astart, row.aend, row.bstart, row.bend, row.ry, row.qy, v, row.col, alpha=alpha, label=row.lab, lw=row.lw, zorder=row.zorder)
                l = ax.add_patch(p)
                if row.lab != '':
                    if not legenddict[row.type]:
                        svlabels[row.type] = l
                        legenddict[row.type] = True
    return ax, [svlabels[i] for i in ['SYN', 'INV', 'TRANS', 'DUP'] if i in svlabels]
# END


def bezierpath(rs, re, qs, qe, ry, qy, v, col, alpha, label='', lw=0, zorder=0):
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
    patch = patches.PathPatch(path, facecolor=col, lw=lw, alpha=alpha, label=label, edgecolor=col, zorder=zorder)
    return patch
# END


def drawmarkers(ax, b, v, chrlengths, indents, chrs, chrgrps, S, cfg, itx, chr_start_coord, minl=0, maxl=-1):
    import logging
    logger = logging.getLogger('drawmarkers')
    mdata = readannobed(b, v, chrlengths)
    for m in mdata:
        # get genome index for genomes in "markers.bed" from the order in chrlengths (parallels the order in "genomes.txt")
        ind = [i for i in range(len(chrlengths)) if chrlengths[i][0] == m.genome][0]
        # get chromosome id (e.g. "Chr1") and in the order of chrord if given, or in the order in "syri.out"
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

        if not itx:
            indent = indents[ind]
            if not v:
                if m.mt != '_':
                    ax.plot(m.start, indent-offset, marker=m.mt, color=m.mc, markersize=m.ms)
                    ax.text(m.start, indent-offset+m.tp, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='center', va='bottom')
                else:
                    ax.axhline(y=indent-offset, xmin=m.start/ax.get_xlim()[1], xmax=m.end/ax.get_xlim()[1], color=m.mc, linewidth=m.ms)
                    ax.text((m.start + m.end)/2, indent-offset+m.tp, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='center', va='bottom')
                # if m.tt != '':
                #     ax.text(m.start, indent-offset+m.tp, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='center', va='bottom')
            elif v:
                if m.mt != '|':
                    ax.plot(indent+offset, m.start, marker=m.mt, color=m.mc, markersize=m.ms)
                    ax.text(indent+offset-m.tp, m.start, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='left', va='center', rotation='vertical')
                else:
                    ax.axvline(x=indent+offset, ymin=m.start/ax.get_ylim()[1], ymax=m.end/ax.get_ylim()[1], color=m.mc, linewidth=m.ms)
                    ax.text(indent+offset-m.tp, (m.start+m.end)/2, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='left', va='center', rotation='vertical')
                # if m.tt != '':
                #     ax.text(indent+offset-m.tp, m.start, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='left', va='center', rotation='vertical')
        elif itx:
            # chr_start_coord is dict with genomes as keys and dict of chromosome plotting start coordinates per chromosome as values
            # buff selects the dict by genome name (m.genome)
            buff = chr_start_coord[m.genome]
            chrid = chrid[0]
            step = S/(len(chrlengths)-1)
            if not v:
                if m.mt != '_':
                    # buff[chrid] selects the chromosome plotting start coordinate by chromosome name (chrid)
                    ax.plot(m.start + buff[chrid], S - (step * ind), marker=m.mt, color=m.mc, markersize=m.ms)
                    ax.text(m.start + buff[chrid], S - (step * ind) + m.tp, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='center', va='bottom')
                else:
                    ax.axhline(y=S - (step * ind), xmin=(m.start + buff[chrid]) / ax.get_xlim()[1], xmax=(m.end + buff[chrid]) / ax.get_xlim()[1], color=m.mc, linewidth=m.ms)
                    ax.text(((m.start + m.end) / 2) + buff[chrid], S - (step * ind) + m.tp, m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='center', va='bottom')
            else:
                if m.mt != '|':
                    ax.plot(1 - S + (step * ind), m.start + buff[chrid], marker=m.mt, color=m.mc, markersize=m.ms)
                    ax.text(1 - S + (step * ind) - m.tp, m.start + buff[chrid], m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='left', va='center', rotation='vertical')
                else:
                    ax.axvline(x=1 - S + (step * ind), ymin=(m.start + buff[chrid]) / ax.get_ylim()[1], ymax=(m.end + buff[chrid]) / ax.get_ylim()[1], color=m.mc, linewidth=m.ms)
                    ax.text(1 - S + (step * ind) - m.tp, ((m.start + m.end) / 2) + buff[chrid], m.tt, color=m.tc, fontsize=m.ts, fontfamily=m.tf, ha='left', va='center', rotation='vertical')
    return ax
# END


def drawtracks(ax, tracks, s, chrgrps, chrlengths, v, itx, cfg, chr_start_coord, minl=0, maxl=-1):
    from matplotlib.patches import Rectangle
    import numpy as np
    from collections import deque
    from functools import partial
    import pandas as pd
    th = (1 - s - 2*cfg['chrmar'])/max([t.ti for t in tracks])
    if th < 0.01:
        raise RuntimeError("Decrease the value of -S to plot tracks correctly. Exiting.")
    cl = len(chrgrps.keys())
    chrs = list(chrgrps.keys())
    diff = 0.7*th
    # Define space between track and track label
    # TODO: read margin spacing from base config
    if maxl == -1:
        margin = np.max([chrlengths[i][1][v[i]] for v in chrgrps.values() for i in range(len(v))])/300
    else:
        margin = (maxl-minl)/300
    if itx:
        if maxl < 1:
            raise ValueError("Incorrect value for maxl. This is a bug. Contact developers.")
        rbuff = chr_start_coord[chrlengths[-1][0]] if not v else chr_start_coord[chrlengths[0][0]]# chr_start_coord is dict of chromosome plotting starting coordinates
        # chrlengths[-1][0] gets last genome name; chr_start_coord[chrlengths[0][0]] gets the first for v mode
    for i in range(len(tracks)):
        # Plot background rectangles for the tracks
        ti = tracks[i].ti   # track index
        for j in range(cl):
            if not v:
                x0 = 0 if not itx else 0 + rbuff[chrs[j]]# no modification because the way it is now works
                y0 = cl - j - th*(ti) if not itx else 1 - th*(ti)
                ax.add_patch(Rectangle((x0, y0), chrlengths[0][1][chrs[j]], diff,  linewidth=0, facecolor=tracks[i].bc, alpha=tracks[i].ba, zorder=1))
            else:
                x0 = j + (ti)*th - diff if not itx else (ti)*th - diff
                y0 = 0 if not itx else 0 + rbuff[chrs[j]]
                ax.add_patch(Rectangle((x0, y0), diff, chrlengths[0][1][chrs[j]], linewidth=0, facecolor=tracks[i].bc, alpha=tracks[i].ba, zorder=1))

        if tracks[i].ft in ['bed', 'bedgraph']:
            bedbin = tracks[i].bincnt
            globaltposmax = None
            if cfg['norm'] == 'F':
                for j in range(cl):
                    # Select positions that are within the limits
                    if maxl != -1:
                        tpos = [k[1] for k in bedbin[chrs[j]] if minl <= k[0] <= maxl]
                    else:
                        tpos = [k[1] for k in bedbin[chrs[j]]]
                    globaltposmax = (max(tpos) if globaltposmax is None else max(tpos + [globaltposmax])) if len(tpos) > 0 else globaltposmax

            for j in range(cl):
                # Select positions that are within the limits
                if maxl != -1:
                    chrpos = [k[0] if not itx else k[0] + rbuff[chrs[j]] for k in bedbin[chrs[j]] if minl <= k[0] <= maxl]
                    tpos = [k[1] for k in bedbin[chrs[j]] if minl <= k[0] <= maxl]
                else:
                    chrpos = [k[0] if not itx else k[0] + rbuff[chrs[j]] for k in bedbin[chrs[j]]]
                    tpos = [k[1] for k in bedbin[chrs[j]]]
                tposmax = (max(tpos) if len(tpos) > 0 else 0) if cfg['norm'] == 'T' else globaltposmax
                tvars = {'color': tracks[i].lc, 'lw': tracks[i].lw, 'zorder': 2, 'alpha': tracks[i].ta}
                if not v:
                    y0 = cl - j - th*(ti) if not itx else 1 - th*(ti)
                    ypos = [(t*diff/tposmax)+y0 for t in tpos] if len(tpos) > 0 else y0
                    if tracks[i].tt == 'f':
                        ax.fill_between(chrpos, ypos, y0, **tvars)
                    elif tracks[i].tt == 'l':
                        ax.plot(chrpos, ypos, **tvars)
                    else:
                        raise ValueError(f"Invalid value for track type for track {tracks[i].n}")
                    if not itx:
                        xpos = chrlengths[0][1][chrs[j]] + margin if maxl == -1 else maxl + margin
                        ax.text(xpos, y0 + diff/2, tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='left', va='center', rotation='horizontal')
                else:
                    x0 = j + (ti)*th - diff if not itx else (ti)*th - diff
                    xpos = [x0 + diff - (t*diff/tposmax) for t in tpos] if len(tpos) > 0 else x0
                    if tracks[i].tt == 'f':
                        ax.fill_betweenx(chrpos, xpos, x0+diff, **tvars)
                    elif tracks[i].tt == 'l':
                        ax.plot(xpos, chrpos, **tvars)
                    else:
                        raise ValueError(f"Invalid value for track type for track {tracks[i].n}")

                    if not itx:
                        ypos = chrlengths[0][1][chrs[j]] + margin if maxl == -1 else maxl + margin
                        ax.text(x0 + diff/2, ypos, tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='center', va='bottom', rotation='vertical')
            if itx:
                if not v:
                    xpos = chrlengths[0][1][chrs[j]] + margin if maxl == -1 else maxl + margin
                    ax.text(xpos, y0 + diff/2, tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='left', va='center', rotation='horizontal')
                else:
                    ypos = chrlengths[0][1][chrs[j]] + margin if maxl == -1 else maxl + margin
                    ax.text(x0 + diff/2, ypos, tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='center', va='bottom', rotation='vertical')

        elif tracks[i].ft in ['gff']:
            from matplotlib import collections as mc
            annos = tracks[i].gff
            l = minl
            r = maxl if maxl != -1 else np.inf
            anno = deque()
            for j in range(cl):
                for mloc, cloc in annos[chrs[j]].items():
                    if l <= mloc[0] and mloc[1] <= r:
                        anno.append([chrs[j], mloc[0], mloc[1], 'mrna'])
                        for c in list(cloc)[0]:
                            anno.append([chrs[j], c[0], c[1], 'cds'])
            anno = pd.DataFrame(anno)
            anno.columns = ['chr', 'xstart', 'xend', 'type']
            anno['zorder'] = 2
            anno.loc[anno['type'] == 'cds', 'zorder'] = 3
            anno['lw'] = tracks[i].lw
            anno.loc[anno['type'] == 'cds', 'lw'] = 2*tracks[i].lw
            anno['colour'] = tracks[i].lc
            anno['fixed'] = round(1 - th*(ti) + diff/2, 4) if not v else round((ti)*th - diff/2, 4)
            if not itx:
                # Update fixed coordinate when not using ITX mode
                anno['fixed'] = [round(cl-chrs.index(c)-th*(ti) + diff/2, 4) if not v else round(chrs.index(c) + (ti)*th - diff/2, 4) for c in anno['chr']]
            if itx:
                # Update genome coordinate when using the ITX mode
                for c in chrs:
                    anno.loc[anno['chr'] == c, 'xstart'] += rbuff[c]
                    anno.loc[anno['chr'] == c, 'xend'] += rbuff[c]

            for grp in anno.groupby(['chr', 'type']):
                mcl = partial(mc.LineCollection, colors=tracks[i].lc, linewidths=float(pd.unique(grp[1]['lw'])), zorder=float(pd.unique(grp[1]['zorder'])))
                lc = mcl([[(row.xstart, row.fixed), (row.xend, row.fixed)] if not v else [(row.fixed, row.xstart), (row.fixed, row.xend)] for row in grp[1].itertuples(index=False)])
                ax.add_collection(lc)
            if not v:
                axt = partial(ax.text, s=tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='left', va='center', rotation='horizontal')
            else:
                axt = partial(ax.text, s=tracks[i].n, color=tracks[i].nc, fontsize=tracks[i].ns, fontfamily=tracks[i].nf, ha='center', va='bottom', rotation='vertical')
            if not itx:
                for j in range(cl):
                    pos = chrlengths[0][1][chrs[j]] + margin + (tracks[i].nm*chrlengths[0][1][chrs[j]]) if maxl == -1 else maxl + margin + (tracks[i].nm*maxl)
                    if not v:
                        axt(pos, float(pd.unique(anno.loc[anno['chr'] == chrs[j], 'fixed'])))
                    else:
                        axt(float(pd.unique(anno.loc[anno['chr'] == chrs[j], 'fixed'])), pos)
            else:
                pos = maxl + margin + (tracks[i].nm*maxl)
                if not v:
                    axt(pos, float(pd.unique(anno['fixed'])))
                else:
                    axt(float(pd.unique(anno['fixed'])), pos)
    return ax
# END

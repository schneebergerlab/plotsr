import sys
import logging.config
import logging
from pandas import concat as pdconcat
from pandas import unique
from plotsr.func import *
from collections import deque, OrderedDict
import os
from math import ceil
import matplotlib

def plotsr(args):
    ## Define logger
    setlogconfig(args.log)
    logger = logging.getLogger("Plotsr")

    ###################################################################
    # Check python and pandas version. Add other checks (if required!!)
    ###################################################################
    try:
        assert sys.version_info.major == 3
        assert sys.version_info.minor >= 8
    except AssertionError:
        logger.warning('\nPlotsr is tested for Python >=3.8. Currently using Python {}.{}. This may result in errors.'.format(sys.version_info.major, sys.version_info.minor))
    except KeyboardInterrupt:
        raise()
    except Exception as E:
        sys.exit(E)

    ## Validate input
    if args.sr is None and args.bp is None:
        logger.error("No structural annotations provided. Use --sr or -bp to provide path to input files")
        sys.exit()

    if args.sr is not None and args.bp is not None:
        logger.error("Both --sr and --bp cannot be used. Use single file type for all input structural annotations files. User converter to reformat BEDPE/syri.out files")
        sys.exit()

    # Check if both --chr and --reg are defined
    if args.chr is not None and args.reg is not None:
        logger.error("Both --chr and --reg are provided. Only one parameter can be provided at a time. Exiting.")
        sys.exit()

    # Check if both --chr and --chrord are defined
    if args.chr is not None and args.chrord is not None:
        logger.error("Both --chr and --chrord are provided. Only one parameter can be provided at a time. Exiting.")
        sys.exit()

    # Check if --rtr is used without --reg
    if args.rtr and args.reg is None:
        logger.error("Cannot use --rtr without --reg. Exiting.")
        sys.exit()


    ###################################################################
    # Declare variable using argument values
    ###################################################################

    # Set Figure height and width. Change later based on chromosome number and size
    FS = args.f             # Font size
    H = args.H              # Height
    W = args.W              # Width
    O = args.o              # Output file name
    D = args.d              # Output file DPI
    R = args.R              # Create ribbons
    V = args.v              # Vertical chromosomes
    S = args.S              # Space between homologous chromosomes
    B = None if args.markers is None else args.markers.name              # Annotation bed file
    TRACKS = None if args.tracks is None else args.tracks.name
    REG = None if args.reg is None else args.reg.strip().split(":")
    RTR = args.rtr
    CHRS = args.chr
    ITX = args.itx


    ## Get config
    cfg = readbasecfg('', V) if args.cfg is None else readbasecfg(args.cfg.name, V)
    if S < 0.1 or S > 0.75:
        sys.exit('Out of range value for S. Please provide a value in the range 0.1-0.75')

    ## Check output file extension
    if len(O.split('.')) == 1:
        logger.warning("Output filename has no extension. Plot would be saved as a pdf")
        O = O + ".pdf"
    elif O.split('.')[-1] not in ['pdf', 'png', 'svg']:
        logger.warning("Output file extension is not in {'pdf','png', 'svg'}. Plot would be saved as a pdf")
        O = O.rsplit(".", 1)[0] + ".pdf"

    ## Set matplotlib backend

    try :
        # matplotlib.use(args.b)
        matplotlib.use('Qt5Agg')    # TODO: Delete this line
    except :
        sys.exit('Matplotlib backend cannot be selected')

    fins = ['col_lersyri.out', 'ler_cvisyri.out', 'cvi_erisyri.out', 'eri_shasyri.out', 'sha_kyosyri.out', 'kyo_an1syri.out', 'an1_c24syri.out'] #TODO: Delete this line
    # fins = ['ler_cvisyri.out', 'cvi_erisyri.out', 'eri_shasyri.out', 'sha_kyosyri.out'] #TODO: Delete this line
    # Read alignment coords
    alignments = deque()
    chrids = deque()
    if args.sr is not None:
        for f in args.sr:
            fin = f.name
        for fin in fins: #TODO: Delete this line
            al, cid = readsyriout(fin) #TODO: Delete this line
            alignments.append([os.path.basename(fin), al])
            chrids.append((os.path.basename(fin), cid))
    elif args.bp is not None:
        for f in args.bp:
            fin = f.name
            al, cid = readbedout(fin)
            alignments.append([os.path.basename(fin), al])
            chrids.append((os.path.basename(fin), cid))

    # Get groups of homologous chromosomes. Use the order from the user if provided.
    cs = set(unique(alignments[0][1]['achr']))
    if args.chrord is None:
        chrs = [k for k in chrids[0][1].keys() if k in alignments[0][1]['achr'].unique()]
    else:
        chrs = deque()
        with open(args.chrord.name, 'r') as fin:
            for line in fin:
                c = line.strip()
                if c not in cs:
                    logger.error("Chromosome {} in {} is not a chromosome in alignment file {}. Exiting.".format(c, args.chrord.name, alignments[0][0]))
                    sys.exit()
                chrs.append(c)
        chrs = list(chrs)
        # Check that the chrorder file contains all chromosomes
        if len(chrs) != len(cs):
            logger.error("Number of chromsomes in {} is less than the number of chromsomes in the alignment file {}. Either list the order of all chromosomes or use --chr if chromosome selection is requires. Exiting.".format(args.chrord.name, alignments[0][0]))
            sys.exit()

    chrgrps = OrderedDict()
    for c in chrs:
        cg = deque([c])
        cur = c
        for i in range(len(chrids)):
            n = chrids[i][1][cur]
            cg.append(n)
            cur = n
        chrgrps[c] = cg

    # Filter alignments to select long alignments between homologous chromosomes
    for i in range(len(alignments)):
        alignments[i][1] = filterinput(args, alignments[i][1], chrids[i][1], ITX)

        # Check chromsome IDs and sizes
    chrlengths, genomes = validalign2fasta(alignments, args.genomes.name)
    # chrlengths, genomes = validalign2fasta(alignments, 'genomes.txt') # TODO: Delete this line


    # Select only chromosomes selected by --chr
    if CHRS is not None:
        alignments, chrs, chrgrps, chrlengths = selectchrom(CHRS, cs, chrgrps, alignments, chrlengths, chrids)


    if REG is not None:
        alignments, chrs, chrgrps = selectregion(REG, RTR, chrlengths, alignments, chrids)

    # Combine Ribbon is selected than combine rows
    if R:
        for i in range(len(alignments)):
            alignments[i][1] = createribbon(alignments[i][1])

    # invert coord for inverted query genome
    for i in range(len(alignments)):
        df = alignments[i][1].copy()
        invindex = ['INV' in i for i in df['type']]
        g = set(df.loc[invindex, 'bstart'] < df.loc[invindex, 'bend'])
        if len(g) == 2:
            logger.error("Inconsistent coordinates in input file {}. For INV, INVTR, INVDUP annotations, either bstart < bend for all annotations or bstart > bend for all annotations. Mixing is not permitted.".format(alignments[i][0]))
            sys.exit()
        elif False in g:
            continue
        df.loc[invindex, 'bstart'] = df.loc[invindex, 'bstart'] + df.loc[invindex, 'bend']
        df.loc[invindex, 'bend'] = df.loc[invindex, 'bstart'] - df.loc[invindex, 'bend']
        df.loc[invindex, 'bstart'] = df.loc[invindex, 'bstart'] - df.loc[invindex, 'bend']
        alignments[i][1] = df.copy()


    # from matplotlib import pyplot as plt
    plt = matplotlib.pyplot
    plt.rcParams['font.size'] = FS
    try:
        if H is None and W is None:
            H = len(chrs)
            W = 3
            fig = plt.figure(figsize=[W, H])
        elif H is not None and W is None:
            fig = plt.figure(figsize=[H, H])
        elif H is None and W is not None:
            fig = plt.figure(figsize=[W, W])
        else:
            fig = plt.figure(figsize=[W, H])
    except Exception as e:
        logger.error("Error in initiliazing figure. Try using a different backend.\n{}".format(e.with_traceback()))
        sys.exit()
    ax = fig.add_subplot(111, frameon=False)

    allal = pdconcat([alignments[i][1] for i in range(len(alignments))])
    if ITX:
        minl = 0
        MCHR = 0.01     # TODO : read spacing between neighbouring chromosome from config file
        maxchr = max([sum(chrlengths[i][1].values()) for i in range(len(chrlengths))])
        maxl = int(maxchr/(MCHR + 1 - (MCHR*len(chrgrps))))
    elif REG is None:
        minl, maxl = 0, -1
    else:
        minl = min(allal[['astart', 'bstart']].apply(min))
        maxl = max(allal[['aend', 'bend']].apply(max))
    labelcnt = 0
    if 'SYN' in allal['type'].array:
        labelcnt += 1
    if 'INV' in allal['type'].array:
        labelcnt += 1
    if 'TRA' in allal['type'].array or 'INVTR' in allal['type'].array:
        labelcnt += 1
    if 'DUP' in allal['type'].array or 'INVDP' in allal['type'].array:
        labelcnt += 1

    ## Draw Axes
    ax, max_l = drawax(ax, chrgrps, chrlengths, V, S, cfg, ITX, minl=minl, maxl=maxl)

    ## Draw Chromosomes
    ax, indents, chrlabels = pltchrom(ax, chrs, chrgrps, chrlengths, V, S, genomes, cfg, ITX, minl=minl, maxl=maxl)

    if cfg['genlegcol'] < 1:
        ncol = ceil(len(chrlengths)/labelcnt)
    else:
        ncol = int(cfg['genlegcol'])

    # Get Genome legend
    if cfg['legend']:
        bbox_to_anchor = cfg['bbox']
        l1 = plt.legend(handles=chrlabels, loc='lower left', bbox_to_anchor=bbox_to_anchor, ncol=ncol, mode=None, borderaxespad=0., frameon=False, title='Genomes')
        l1._legend_box.align = "left"
        plt.gca().add_artist(l1)

    # Plot structural annotations
    ax, svlabels = pltsv(ax, alignments, chrs, V, chrgrps, chrlengths, indents, S, cfg, ITX)

    # TODO: Uncomment this
    # if cfg['legend']:
    #     bbox_to_anchor[0] += cfg['bboxmar']
    #     plt.legend(handles=svlabels, loc='lower left', bbox_to_anchor=bbox_to_anchor, ncol=1, mode='expand', borderaxespad=0., frameon=False, title='Annotations')._legend_box.align = "left"

    # Plot markers
    if B is not None:
        ax = drawmarkers(ax, B, V, chrlengths, indents, chrs, chrgrps, S, ITX, minl=minl, maxl=maxl)

    # Draw tracks
    if TRACKS is not None:
        tracks = readtrack(TRACKS, chrlengths)
        # tracks = readtrack(f, chrlengths) #TODO: delete this
        ax = drawtracks(ax, tracks, S, chrgrps, chrlengths, V, cfg, minl=minl, maxl=maxl)

    # Save the plot
    try:
        fig.savefig(O, dpi=D, bbox_inches='tight', pad_inches=0.01)
    except Exception as e:
        sys.exit('Error in saving the figure. Try using a different backend.' + '\n' + e.with_traceback())

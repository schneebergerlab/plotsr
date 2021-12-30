#!/usr/bin/env python3

"""
Author: Manish Goel
Date: 30.12.2021
Description: Plotting multi genome structural annotations 
"""

from __init__ import __version__
import argparse


if __name__ == '__main__':
    from matplotlib.rcsetup import non_interactive_bk as bklist
    parser = argparse.ArgumentParser("Plotting structural rearrangements between genomes", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sr', help='Structural annotation mappings (syri.out) identified by SyRI', action='append', type=argparse.FileType('r'))
    parser.add_argument('--bp', help='Structural annotation mappings in BEDPE format', action='append', type=argparse.FileType('r'))
    parser.add_argument('--genomes', help='Path to genome files', type=argparse.FileType('r'), required=True)
    parser.add_argument('--markers', help='Path to markers (bed format)', type=argparse.FileType('r'))
    parser.add_argument('--tracks', help='File listing paths and details for all tracks to be plotted', type=argparse.FileType('r'))
    # parser.add_argument('-B', help='Annotation bed file for marking specific positions on genome', type=argparse.FileType('r'))
    parser.add_argument('--chr', help='Select specific chromosomes on reference for plotting.', type=str, action='append')
    parser.add_argument('--nosyn', help='Do not plot syntenic regions', default=False, action='store_true')
    parser.add_argument('--noinv', help='Do not plot inversions', default=False, action='store_true')
    parser.add_argument('--notr', help='Do not plot translocations regions', default=False, action='store_true')
    parser.add_argument('--nodup', help='Do not plot duplications regions', default=False, action='store_true')
    parser.add_argument('-s', help='minimum size of a SR to be plotted', type=int, default=10000)
    parser.add_argument('-R', help='Create ribbons', default=False, action="store_true")
    parser.add_argument('-f', help='font size', type=int, default=6)
    parser.add_argument('-H', help='height of the plot', type=int)
    parser.add_argument('-W', help='width of the plot', type=int)
    parser.add_argument('-S', help='Space between homologous chromosome (0.1-0.85). Adjust this to make more space for annotation marker/text.', default=0.7, type=float)
    parser.add_argument('-o', help='output file format (pdf, png, svg)', default="pdf", choices=['pdf', 'png', 'svg'])
    parser.add_argument('-d', help='DPI for the final image', default="300", type=int)
    parser.add_argument('-b', help='Matplotlib backend to use', default="agg", type=str, choices=bklist)
    parser.add_argument('-v', help='Plot vertical chromosome', default=False, action='store_true')
    parser.add_argument('--log', help='Log-level', choices=['DEBUG', 'INFO', 'WARN'], default='WARN', type=str)
    parser.add_argument('--version', action='version', version='{version}'.format(version=__version__))

    # args = parser.parse_args([]) # TODO: Delete this line
    args = parser.parse_args()
    ## Define logger
    import logging
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
            # 'log_file': {
            #     'class': 'logging.FileHandler',
            #     'filename': args.log_fin.name,
            #     'mode': 'a',
            #     'formatter': 'log_file',
            #     'level': args.log,
            # },
        },
        'loggers': {
            '': {
                'level': args.log,
                'handlers': ['stdout'],
                # 'handlers': ['stdout', 'log_file'],
            },
        },
    })
    logger = logging.getLogger("Plotsr")

    ## Validate input
    import sys
    if len(args.sr) == 0 and len(args.bp) == 0:
        logger.error("No structural annotations provided. Use --sr or -bp to provide path to input files")
        sys.exit()
    try:
        if len(args.sr) > 0 and len(args.bp) > 0:
            logger.error("Both --sr and --bp cannot be used. Use single file type for all input structural annotations files. User converter to reformat BEDPE/syri.out files")
            sys.exit()
    except TypeError:
        pass

    # Set Figure height and width. Change later based on chromosome number and size
    FS = args.f             # Font size
    H = args.H              # Height
    W = args.W              # Width
    O = args.o              # Output file format
    D = args.d              # Output file DPI
    R = args.R              # Create ribbons
    V = args.v              # Vertical chromosomes
    S = args.S              # Space between homologous chromosomes
    MARKERS = None if args.markers is None else args.markers.name              # Annotation bed file
    TRACKS = None if args.tracks is None else args.tracks.name
    if S < 0.1 or S > 0.85:
        sys.exit('Out of range value for S. Please provide a value in the range 0.1-0.85')

    from plotsr.func import VARS, readfasta
    from collections import deque, OrderedDict

    ## Set matplotlib backend
    import matplotlib
    try:
        # matplotlib.use(args.b)
        matplotlib.use('Qt5Agg')
    except:
        sys.exit('Matplotlib backend cannot be selected')

    # fins = ['col_lersyri.out', 'ler_cvisyri.out', 'cvi_colsyri.out'] #TODO: Delete this line
    # Read alignment coords
    alignments = deque()
    chrids = deque()
    if len(args.sr) > 0:
        for fin in args.sr:
            al, cid = readsyriout(fin.name)
        # for fin in fins: #TODO: Delete this line
        #     al, cid = readsyriout(fin) #TODO: Delete this line
            alignments.append([os.path.basename(fin), al])
            chrids.append((os.path.basename(fin), cid))
    elif len(args.bp) > 0:
        for fin in args.bp:
            al, cid = readbedout(fin.name)
            alignments.append([os.path.basename(fin), al])
            chrids.append((os.path.basename(fin), cid))

    # Get groups of homologous chromosomes
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

    # Filter alignments to select long alignments between homologous chromosomes
    for i in range(len(alignments)):
        alignments[i][1] = filterinput(args, alignments[i][1], chrids[i][1])

    # Select only chromosomes selected by --chr
    if args.chr is not None:
        homchrs = deque()
        for c in args.chr:
            if c not in chrs:
                logger.warning("Selected chromosome: {} is not in reference genome. Skipping it.".format(c))
                continue
            homchrs.append(chrgrps[c])
        for i in range(len(alignments)):
            alignments[i][1] = alignments[i][1].loc[alignments[i][1]['achr'].isin([h[i] for h in homchrs])]

    # Check chromsome IDs and sizes
    chrlengths = validalign2fasta(alignments, args.genomes)
    # chrlengths, chrtags = validalign2fasta(alignments, 'genomes.txt') # TODO: Delete this line
    ## Remove chromosomes that are not homologous to selected reference chromosomes
    for i in range(len(chrlengths)):
        ks = list(chrlengths[i][1].keys())
        homs = [h[i] for h in homchrs]
        for k in ks:
            if k not in homs:
                chrlengths[i][1].pop(k)

    # Update groups of homologous chromosomes
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


    # Combine Ribbon is selected than combine rows
    if R:
        for i in range(len(alignments)):
            alignments[i][1] = createribbon(alignments[i][1])

    # invert coord for inverted query genome

    for i in range(len(alignments)):
        df = alignments[i][1].copy()
        invindex = ['INV' in i for i in df['type']]
        df.loc[invindex, 'bstart'] = df.loc[invindex, 'bstart'] + df.loc[invindex, 'bend']
        df.loc[invindex, 'bend'] = df.loc[invindex, 'bstart'] - df.loc[invindex, 'bend']
        df.loc[invindex, 'bstart'] = df.loc[invindex, 'bstart'] - df.loc[invindex, 'bend']
        alignments[i][1] = df.copy()




    from matplotlib import pyplot as plt
    plt.rcParams['font.size'] = FS
    try:
        if H is None and W is None:
            H = len(chrs)
            W = 3
            fig = plt.figure(figsize=[W, H])
        if H is not None and W is None:
            fig = plt.figure(figsize=[H, H])
        if H is None and W is not None:
            fig = plt.figure(figsize=[W, W])
        if H is not None and W is not None:
            fig = plt.figure(figsize=[W, H])
    except Exception as e:
        sys.exit("Error in initiliazing figure. Try using a different backend." + '\n' + e.with_traceback())
    ax = fig.add_subplot(111, frameon=False)

    ## Draw Axes
    ax, max_l = drawax(ax, chrgrps, chrlengths, V, S)

    ## Draw Chromosomes
    ax, indents, chrlabels = pltchrom(ax, chrs, chrgrps, chrlengths, V, S, chrtags)

    # Get Chromosome legend
    bbox_to_anchor = [0, 1.01, 0.5, 0.3] if not V else [0, 1.1, 0.5, 0.3] # TODO: READ from base.cfg
    l1 = plt.legend(handles=chrlabels, loc='lower left', bbox_to_anchor=bbox_to_anchor, ncol=1, mode='expand', borderaxespad=0., frameon=False, title='Genome')
    l1._legend_box.align = "left"
    plt.gca().add_artist(l1)

    # Plot structural annotations
    # TODO: Parameterise: colors, alpha,
    ax, svlabels = pltsv(ax, alignments, chrs, V, chrgrps, indents)
    bbox_to_anchor[0] += 0.5
    plt.legend(handles=svlabels, loc='lower left', bbox_to_anchor=bbox_to_anchor, ncol=1, mode='expand', borderaxespad=0., frameon=False, title='Annotation')._legend_box.align = "left"

    # Plot markers
    if B is not None:
        ax = drawmarkers(ax, B, V, chrlengths, indents, chrs, chrgrps)

    # Draw tracks
    if args.tracks is not None:
        tracks = readtrack(args.tracks.name, chrlengths)
        # tracks = readtrack(f, chrlengths) #TODO: delete this
        ax = drawtracks(ax, tracks, S, chrgrps, chrlengths, V)

    # Save the plot
    try:
        fig.savefig('syri.'+O, dpi=D, bbox_inches='tight', pad_inches=0.01)
    except Exception as e:
        sys.exit('Error in saving the figure. Try using a different backend.' + '\n' + e.with_traceback())

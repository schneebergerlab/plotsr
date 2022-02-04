#!/usr/bin/env python3

"""
Author: Manish Goel
Date: 30.12.2021
Description: Plotting multi genome structural annotations
"""

import argparse
from plotsr import __version__
from matplotlib.rcsetup import non_interactive_bk as bklist
from plotsr.plotsr import plotsr

def main(cmd):
    parser = argparse.ArgumentParser("Plotting structural rearrangements between genomes", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    other = parser._action_groups.pop()
    inputfiles = parser.add_argument_group("Input/Output files")
    inputfiles.add_argument('--sr', help='Structural annotation mappings (syri.out) identified by SyRI', action='append', type=argparse.FileType('r'))
    inputfiles.add_argument('--bp', help='Structural annotation mappings in BEDPE format', action='append', type=argparse.FileType('r'))
    inputfiles.add_argument('--genomes', help='File containing path to genomes', type=argparse.FileType('r'), required=True)
    inputfiles.add_argument('--markers', help='File containing path to markers (bed format)', type=argparse.FileType('r'))
    inputfiles.add_argument('--tracks', help='File listing paths and details for all tracks to be plotted', type=argparse.FileType('r'))
    inputfiles.add_argument('--chrord', help='File containing reference (first genome) chromosome IDs in the order in which they are to be plotted. File requires one chromosome ID per line. Not compatible with --chr', type=argparse.FileType('r'))
    inputfiles.add_argument('-o', help='Output file name. Acceptable format: pdf, png, svg', default="plotsr.pdf")

    filtering = parser.add_argument_group("Data filtering")
    # filtering.add_argument('--itx', help='Plot inter-chromosomal SRs as well (experimental)', default='FALSE', action='store_true')
    filtering.add_argument('--chr', help='Select specific chromosome on reference (first genome) and plots them in the given order. Not compatible with --chrord. Can be used multiple time to select more than one chromosomes.', type=str, action='append')
    filtering.add_argument('--reg', help='Plots a specific region. Use as: GenomeID:ChromosomeID:Start-End. Not compatible with --chr and -R.', type=str)
    filtering.add_argument('--rtr', help='When using --reg, plot all SRs that are within the boundaries of the homologous regions. For highly zoomed regions, this could result in visually disconnected alignments.', default=False, action='store_true')
    filtering.add_argument('--nosyn', help='Do not plot syntenic regions', default=False, action='store_true')
    filtering.add_argument('--noinv', help='Do not plot inversions', default=False, action='store_true')
    filtering.add_argument('--notr', help='Do not plot translocations regions', default=False, action='store_true')
    filtering.add_argument('--nodup', help='Do not plot duplications regions', default=False, action='store_true')
    filtering.add_argument('-s', help='minimum size of a SR to be plotted', type=int, default=10000)

    plotting = parser.add_argument_group("Plot adjustment")
    plotting.add_argument('--cfg', help='Path to config file containing parameters to adjust plot.', type=argparse.FileType('r'))
    plotting.add_argument('-R', help='Join adjacent syntenic blocks if they are not interrupted by SRs. Using this would decrease gaps in the visualisation.', default=False, action="store_true")
    plotting.add_argument('-f', help='font size', type=int, default=6)
    plotting.add_argument('-H', help='height of the plot', type=int)
    plotting.add_argument('-W', help='width of the plot', type=int)
    plotting.add_argument('-S', help='Space between homologous chromosome (0.1-0.75). Adjust this to make more space for annotation marker/text and tracks.', default=0.7, type=float)
    plotting.add_argument('-d', help='DPI for the final image', default="300", type=int)
    plotting.add_argument('-b', help='Matplotlib backend to use', default="agg", type=str, choices=bklist)
    plotting.add_argument('-v', help='Plot vertical chromosome', default=False, action='store_true')

    other.add_argument('--log', help='Log-level', choices=['DEBUG', 'INFO', 'WARN'], default='WARN', type=str)
    other.add_argument('--version', action='version', version='{version}'.format(version=__version__))
    parser._action_groups.append(other)

    # args = parser.parse_args([]) # TODO: Delete this line
    args = parser.parse_args(cmd)
    plotsr(args)
# END

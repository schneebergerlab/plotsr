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
    parser.add_argument('--sr', help='Structural annotation mappings (syri.out) identified by SyRI', action='append', type=argparse.FileType('r'))
    parser.add_argument('--bp', help='Structural annotation mappings in BEDPE format', action='append', type=argparse.FileType('r'))
    parser.add_argument('--genomes', help='File containing path to genomes', type=argparse.FileType('r'), required=True)
    parser.add_argument('--markers', help='File containing path to markers (bed format)', type=argparse.FileType('r'))
    parser.add_argument('--tracks', help='File listing paths and details for all tracks to be plotted', type=argparse.FileType('r'))
    # parser.add_argument('-B', help='Annotation bed file for marking specific positions on genome', type=argparse.FileType('r'))
    parser.add_argument('--chr', help='Select specific chromosomes on reference for plotting.', type=str, action='append')
    parser.add_argument('--reg', help='Plots a specific region. Use as: GenomeID:ChromosomeID:Start-End. Not compatible with --chr and -R.', type=str)
    parser.add_argument('--cfg', help='Path to config file containing parameters to adjust plot.', type=argparse.FileType('r'))
    parser.add_argument('--nosyn', help='Do not plot syntenic regions', default=False, action='store_true')
    parser.add_argument('--noinv', help='Do not plot inversions', default=False, action='store_true')
    parser.add_argument('--notr', help='Do not plot translocations regions', default=False, action='store_true')
    parser.add_argument('--nodup', help='Do not plot duplications regions', default=False, action='store_true')
    parser.add_argument('-s', help='minimum size of a SR to be plotted', type=int, default=10000)
    parser.add_argument('-R', help='Create ribbons', default=False, action="store_true")
    parser.add_argument('-f', help='font size', type=int, default=6)
    parser.add_argument('-H', help='height of the plot', type=int)
    parser.add_argument('-W', help='width of the plot', type=int)
    parser.add_argument('-S', help='Space between homologous chromosome (0.1-0.75). Adjust this to make more space for annotation marker/text.', default=0.7, type=float)
    parser.add_argument('-o', help='Output file name. Acceptable format: pdf, png, svg', default="plotsr.pdf")
    parser.add_argument('-d', help='DPI for the final image', default="300", type=int)
    parser.add_argument('-b', help='Matplotlib backend to use', default="agg", type=str, choices=bklist)
    parser.add_argument('-v', help='Plot vertical chromosome', default=False, action='store_true')
    parser.add_argument('--log', help='Log-level', choices=['DEBUG', 'INFO', 'WARN'], default='WARN', type=str)
    parser.add_argument('--version', action='version', version='{version}'.format(version=__version__))
    # args = parser.parse_args([]) # TODO: Delete this line
    args = parser.parse_args(cmd)
    plotsr(args)
# END

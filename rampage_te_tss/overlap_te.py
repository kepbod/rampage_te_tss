#!/usr/bin/env python

'''
Usage: overlap_te.py [options] -t te <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -t te                          TE annotations (BED format).
    --extend length                TE extended length. [default: 50]
    --span span                    RAMPAGE span cutoff. [default: 1000]
    --entropy entropy              RAMPAGE entropy cutoff. [default: 2.5]
'''

import os.path
import tempfile
from pybedtools import BedTool
from collections import defaultdict

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '1.0.1'


def overlap(options):
    '''
    Overlap peaks with TE
    '''
    # check options
    extend = int(options['--extend'])
    te_file = options['-t']
    span = int(options['--span'])
    entropy = float(options['--entropy'])
    # prepare tempfile
    temp_dir = tempfile.mkdtemp()
    te = os.path.join(temp_dir,  'te.bed')
    # parse TE annotations
    with open(te_file, 'r') as f, open(te, 'w') as out:
        for line in f:
            chrom, start, end, name, _, strand = line.rstrip().split()
            start = int(start)
            end = int(end)
            if strand == '+':
                start -= extend
                if start < 0:
                    continue
            else:
                end += extend
            out.write('%s\t%d\t%d\t%s\t0\t%s\n' % (chrom, start, end, name,
                                                   strand))
    # fetch TE peaks
    rampage_bed = BedTool(os.path.join(options['<rampagedir>'],
                                       'rampage_entropy.txt'))
    te_bed = BedTool(te)
    te_peak = rampage_bed.intersect(te_bed, wa=True, wb=True)
    # filter TE peaks
    peak_lst = defaultdict(dict)
    for p in te_peak:
        if float(p[13]) < entropy:
            continue
        if not check_span(p, span):
            continue
        peak_site = int(p[6])
        te_chr = p[16]
        te_start = int(p[17])
        te_end = int(p[18])
        te_name = p[19]
        te_strand = p[21]
        if peak_site < te_start or peak_site > te_end:
            continue
        if te_strand == '+':
            te_start += extend
            te_tss = te_start
        else:
            te_end -= extend
            te_tss = te_end
        te_info = '%s\t%d\t%d\t%s\t0\t%s' % (te_chr, te_start, te_end, te_name,
                                             te_strand)
        peak_strand = p[5]
        if peak_strand == te_strand:
            if peak_strand == '+':
                peak_dis = peak_site - te_tss
            else:
                peak_dis = te_tss - peak_site
            te_info += '\tsense|%d' % peak_dis
            peak_dis = abs(peak_dis)
        else:
            if te_start <= peak_site <= te_end:
                peak_dis = 10000
                te_info += '\tantisense'
            else:
                continue
        peak_info = '\t'.join(p[:16])
        peak_lst[peak_info][peak_dis] = te_info
    outfile = os.path.join(options['<rampagedir>'], 'rampage_TE.txt')
    with open(outfile, 'w') as out:
        for peak in peak_lst:
            dis = min(peak_lst[peak])
            out.write('%s\t%s\n' % (peak, peak_lst[peak][dis]))


def check_span(peak, span_cutoff):
    sites = [int(x) for x in peak[14].split('|')]
    if peak.strand == '+':
        span = max(sites) - peak.start
    else:
        span = peak.end - min(sites)
    if span > span_cutoff:
        return True
    else:
        return False


if __name__ == '__main__':
    from docopt import docopt
    overlap(docopt(__doc__, version=__version__))

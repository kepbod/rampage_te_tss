#!/usr/bin/env python

'''
Usage: assign_peak.py [options] -a ANNO -g GTF <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -a ANNO                        Known gene annotation BED12 file.
    -g ANNO                        Assembled gene annotation BED12 file.
    -p THREAD                      Threads. [default: 10]
'''

import os.path
from pybedtools import BedTool
import tempfile
import re
from joblib import Parallel, delayed
from collections import defaultdict

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def assign_peak(options):
    '''
    Assign peaks to genes
    '''
    #  parse options
    known_anno_f = options['-a']
    anno_f = options['-g']
    te_peak = os.path.join(options['<rampagedir>'], 'rampage_TE.txt')
    threads = int(options['-p'])
    # construct known gene info
    known_anno_set = set()
    with open(known_anno_f, 'r') as f:
        for line in f:
            items = line.rstrip().split()
            known_anno_set.add('\t'.join(items[:3] + [items[5]] + items[9:12]))
    # filter out novel gene annotations
    temp_dir = tempfile.mkdtemp()
    novel_anno_f = os.path.join(temp_dir, 'novel_iso.bed')
    novel_anno = BedTool(anno_f)
    novel_anno.filter(filter_bed, anno=known_anno_set).saveas(novel_anno_f)
    # assign rampage peaks to genes
    with open(te_peak, 'r') as f:
        result = Parallel(n_jobs=threads)(delayed(assign)(line.rstrip(),
                                                          known_anno_f,
                                                          novel_anno_f)
                                          for line in f)
    outfile = os.path.join(options['<rampagedir>'], 'rampage_TE_gene.txt')
    with open(outfile, 'w') as out:
        for r in result:
            if r is not None:
                out.write(r)


def assign(peak_info, known_anno_f, novel_anno_f):
    known_anno = BedTool(known_anno_f).sort()
    novel_anno = BedTool(novel_anno_f).sort()
    items = peak_info.split()
    info_str = '%s\t%d\t%d\tpeak\t0\t%s\n'
    chrom = items[0]
    strand = items[5]
    peak_site = int(items[6])
    total = int(items[9])
    sites = [int(x) for x in items[14].split('|')]
    counts = [int(x) for x in items[15].split('|')]
    temp_dir = tempfile.mkdtemp()
    peak_f = os.path.join(temp_dir, 'sites.bed')
    with open(peak_f, 'w') as f:
        for s, n in zip(sites, counts):
            if strand == '+':
                a = s - 1
                b = s
            else:
                a = s
                b = s + 1
            for i in range(n):
                f.write(info_str % (chrom, a, b, strand))
    peak = BedTool(peak_f)
    known_ratio, known_dis, known_gene = overlap_bed(peak, known_anno, 'known',
                                                     peak_site, peak_info,
                                                     total)
    novel_ratio, novel_dis, novel_gene = overlap_bed(peak, novel_anno, 'novel',
                                                     peak_site, peak_info,
                                                     total)
    if known_ratio > novel_ratio:
        return known_gene
    elif known_ratio < novel_ratio:
        return novel_gene
    else:
        return known_gene if known_dis < novel_dis else novel_gene


def overlap_bed(peak, anno, flag, peak_site, peak_info, total):
    isoform = defaultdict(int)
    tss_dis = {}
    abs_tss_dis = {}
    for o in peak.intersect(anno, wa=True, wb=True, s=True, split=True):
        iso_info = '\t'.join(o[6:])
        isoform[iso_info] += 1
        if iso_info in tss_dis:
            continue
        strand = o[11]
        if strand == '+':
            dis = peak_site - int(o[7])
        else:
            dis = int(o[8]) - peak_site
        tss_dis[iso_info] = dis
        abs_dis = abs(dis)
        abs_tss_dis[iso_info] = abs_dis
    peak_gene = ''
    max_ratio = 0
    min_dis = None
    for iso in isoform:
        ratio = isoform[iso] * 1.0 / total
        abs_dis = abs_tss_dis[iso]
        dis = tss_dis[iso]
        gene_info = '%s\t%s\t%s\t%d\t%f\n' % (peak_info, iso, flag, dis,
                                              ratio)
        if ratio == max_ratio:  # same ratio
            if abs_dis < min_dis:  # shorter distance
                max_ratio = ratio
                min_dis = abs_dis
                peak_gene = gene_info
        elif ratio > max_ratio:  # more ratio
            max_ratio = ratio
            min_dis = abs_dis
            peak_gene = gene_info
    if peak_gene:
        return max_ratio, min_dis, peak_gene
    else:
        return None, None, None


def filter_bed(bed, anno):
    if not re.match(r'chr(\d{1,2}|X|Y)$', bed.chrom):  # incorrect chrom
        return False
    bed_info = '\t'.join(bed[:3] + [bed[5]] + bed[9:12])
    if bed_info in anno:  # known isoform
        return False
    return True


if __name__ == '__main__':
    from docopt import docopt
    assign_peak(docopt(__doc__, version=__version__))

#!/usr/bin/env python

'''
Usage: assemble.py [options] -g GTF <bam>...

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -g GTF --gtf=GTF               Gene annotation GTF file.
    -p THREAD --thread=THREAD      Threads. [default: 10]
    --dir=DIR                      Output directory. [default: ./]
    --prefix=PREFIX                Prefix for merged GTF.
                                   [default: merged_stringtie]
'''

import sys
import os
import os.path
from seqlib.path import which
from seqlib.helper import run_command

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '1.0.0'


def assemble(options):
    '''
    Assemble RNA-seq with StringTie
    '''
    # parse options
    if not which('stringtie'):
        sys.exit('Error: No StringTie installed!')
    if not which('gtfToGenePred'):
        sys.exit('Error: No gtfToGenePred installed!')
    if not which('genePredToBed'):
        sys.exit('Error: No genePredToBed installed!')
    gtf = options['--gtf']
    thread = options['--thread']
    bamf = options['<bam>']
    dir_path = options['--dir']
    prefix = options['--prefix']
    if len(bamf) == 1:  # no replicates
        # run StringTie
        out_gtf = run_stringtie(bamf[0], gtf, dir_path, thread)
        # convert GTF to GenePred
        convert_gtf(out_gtf, dir_path, prefix)
    else:  # have replicates
        # run StringTie
        gtf_list = []
        for f in bamf:
            gtf_list.append(run_stringtie(f, gtf, dir_path, thread))
        out_gtf = merge_stringtie(gtf_list, gtf, dir_path, prefix,
                                  thread)
        # convert GTF to GenePred
        convert_gtf(out_gtf, dir_path, prefix)


def run_stringtie(bam, gtf, dir_path, thread):
    if not os.path.isfile(bam):
        sys.exit('No BAM file: %s!' % bam)
    if not os.path.isfile(gtf):
        sys.exit('No GTF file: %s' % gtf)
    fname = os.path.basename(bam)
    prefix = os.path.splitext(fname)[0]
    outf = os.path.join(dir_path, prefix + '_stringtie.gtf')
    command = 'stringtie -G %s -p %s -o %s --rf %s' % (gtf, thread, outf, bam)
    run_command(command, 'Error in StringTie!')
    return outf


def merge_stringtie(gtf_list, gtf, dir_path, prefix, thread):
    outf = os.path.join(dir_path, prefix + '.gtf')
    command = 'stringtie --merge -G %s -p %s -o %s %s' % (gtf, thread, outf,
                                                          '\t'.join(gtf_list))
    run_command(command, 'Error in StringTie --merge!')
    return outf


def convert_gtf(out_gtf, dir_path, prefix):
    # convert GTF to BED 12
    out_genepred = os.path.join(dir_path, prefix + '.txt')
    out_bed = os.path.join(dir_path, prefix + '.bed')
    command = 'gtfToGenePred %s %s' % (out_gtf, out_genepred)
    run_command(command, 'Error in gtfToGenePred!')
    command = 'genePredToBed %s %s' % (out_genepred, out_bed)
    run_command(command, 'Error in genePredToBed!')


if __name__ == '__main__':
    from docopt import docopt
    assemble(docopt(__doc__, version=__version__))

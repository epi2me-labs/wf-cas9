import argparse
import sys

import argparse
from functools import reduce
from pathlib import Path
import math
import sys
from pybedtools import BedTool
import pandas as pd
import numpy as np
import pysam
import seaborn as sns
import matplotlib.pyplot as pl
from natsort import natsort_keygen


def prepare_data(targets, aln, chrom_sizes, stats):
    targets = BedTool(targets)
    aln = BedTool(aln)

    df_chrom_sizes = pd.read_csv(chrom_sizes, sep='\t', index_col=0)
    # chr_used = np.unique(targets.to_dataframe())
    # df_chrom_sizes = df_chrom_sizes[df_chrom_sizes.index.isin(chr_used)]

    # Coverage on tiles
    # tiles = BedTool().makewindows(g=chrom_sizes, w=100, i='winnum')
    tiles = BedTool('/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/window_name_strand.bed')
    # dft = tiles.to_dataframe()
    # dft = pd.read_csv("/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/windows_small", sep='\t')
    # dft[['0']] = '0'
    # dft[['+']] = '+'
    # BedTool().from_dataframe(dft).saveas('/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/window_name_strand.bed')


    # Each window

    # header= [chrom, tile_start, tile_end, name, score, strand, cov, #coverged_bases, tile_len, fracCovBases]
    cov_f = tiles.coverage(aln, s=True)
    cov_r = tiles.coverage(aln, S=True)
    # f = cov_f.to_dataframe()
    # r = cov_r.to_dataframe()

    ct = target_coverage_table(cov_f, cov_r, targets)
    ct.to_csv('target_coverage')

    # ts = target_summary(aln, cov_f, cov_r, targets, stats, tiles)



def target_summary(aln, cov_f, cov_r, targets, stats_file, tiles) -> pd.DataFrame:
    """Create table with summary stats for each target"""

    stats = pd.read_csv(stats_file, index_col=0, sep='\t')
    aln_cov = aln.coverage(targets)
    aln_cov_df = aln_cov.to_dataframe().sort_values(['chrom', 'start'])

    b = targets.intersect(tiles, wb=True, wa=True)
    bdf = b.to_dataframe()

    # Get errors about end being before start so do it in pandas
    # c = b.groupby(g=[1, 2], c=10, o=['count'])

    read_per_tile = bdf.groupby([bdf.columns[0], bdf.columns[1]]).count().iloc[:, 0]

    e = c.to_dataframe()

    on_target_depth = targets.coverage(aln).to_dataframe().sort_index()

    f_hits = targets.coverage(aln, s=True).to_dataframe()
    r_hits = targets.coverage(aln, S=True).to_dataframe()

    on_target_depth.drop(columns=['score'], inplace=True)
    on_target_depth.header = ['chrom', 'start', 'end', 'basesCovered', 'tsize', 'fracTAln']

    # target summaries table

    on = aln_cov_df[aln_cov_df['frac_overlap'] > 0]


    f = on[on['strand'] == '+'].groupby(
        ['target']).count()['chrom']
    r = on[on['strand'] == '-'].groupby(
        ['target']).count()['chrom']




    bias = (f - r) / (f + r)
    bias.columns = ['strand_bias']
    mean_read_len = on.groupby(['target']).mean()[['read_len']]
    result_df['all_overlaps'] = result_df.overlaps_f + result_df.overlaps_r
    median_coverage = result_df.groupby(['target']).median()[['all_overlaps']]
    median_coverage.columns = ['median_coverage']
    # on_target_depth['kbases'] = on.bes + result_df.overlaps_r
    kbases = on.groupby(['target']).sum()[['bases_aligning']] / 1000
    kbases.columns = ['kbases']
    mean_quality = on.groupby(['target']).mean()[['mean_quality']]

    frames = [on_target_depth, kbases, median_coverage, mean_quality,
              mean_read_len, bias]
    on_target_depth = reduce(lambda left, right: pd.merge(
        left, right, on='target', how='left'), frames)
    on_target_depth = on_target_depth.rename(
        columns={'read_len': 'mean_read_len', 'chrom_x': 'chrom'}) \
        .sort_values(by=['chrom', 'start']) \
        .drop(columns=['chrom_y'])
    on_target_depth.reset_index(inplace=True, drop=False)
    on_target_depth.sort_values(['chrom', 'start'], key=natsort_keygen(),
                                inplace=True, ascending=True)
    on_target_depth.to_csv('target_summary.csv')


def target_coverage_table(cov_f: BedTool, cov_r:BedTool, targets:BedTool) \
        -> pd.DataFrame:
    """Create table with stranded, tiled coverage."""


def coverage_summary():
    ...


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", help="bed file of single target region")
    parser.add_argument("alignment_bed", help="alignment")
    parser.add_argument("ref_genome", help='Reference genome fasta')
    parser.add_argument("sample_id", help="alignment")
    parser.add_argument("seq_stats", help="Sequence summary stats")
    parser.add_argument("--chrom_sizes", help="chromosome sizes csv")
    parser.add_argument("--proximal_pad", type=int, default=1000)
    args = parser.parse_args()
    prepare_data(args.targets,
                 args.alignment_bed,
                 args.chrom_sizes,
                 args.seq_stats)


if __name__ == '__main__':
    target_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/targets.bed"
    aln_bed = "/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/fastq_pass.bed"
    genome_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/grch38.fasta.gz"
    stats_file = "/Users/Neil.Horner/work/workflow_outputs/cas9/seqstats.csv"
    chrome_sizes = "/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/chrom_small.txt"
    sys.argv.extend([target_file, aln_bed, genome_file, 'test', stats_file,
                     '--chrom_sizes', '/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/chrom.sizes'])
    main()


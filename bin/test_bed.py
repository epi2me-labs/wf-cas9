from pathlib import Path
import sys
from pybedtools import BedTool
import pandas as pd
import pysam
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Tiles

# targets = pd.read_csv('/Users/Neil.Horner/work/testing/cas9/bg/targets.csv', index_col=0)
# t_read_map = pd.read_csv('/Users/Neil.Horner/work/testing/cas9/bg/target_read_map.csv', index_col=0)
# tile_reads_slop = pd.read_csv('/Users/Neil.Horner/work/testing/cas9/bg/tile_reads_slop.csv', index_col=0)
# # tile_reads['cat'] = 'background'
#
tile_reads = pd.read_csv("/Users/Neil.Horner/work/testing/cas9/bg/tile_cov.csv")
slop_targets = pd.read_csv('/Users/Neil.Horner/work/testing/cas9/bg/slop_targets')


tiles = BedTool(tile_reads)
targets = BedTool(slop_targets)

cov = tiles.intersect(targets)



print('p')





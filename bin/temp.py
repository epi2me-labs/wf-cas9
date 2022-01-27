import pandas as pd

df = pd.read_csv('/Users/Neil.Horner/work/testing/cas9/output/target_coverage.csv', index_col=0)

df['total'] = df.overlaps_f + df.overlaps_r
p = df.groupby(['target']).median()[['total']]

k = df.groupby(['target']).sum()[['total']]
print(k)



s_f = df.groupby(['target']).sum()['overlaps_f']
s_r = df.groupby(['target']).sum()['overlaps_r']
bias = 0 - (s_f / s_r)
print('sb')

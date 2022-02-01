import pandas as pd

df_on = pd.read_csv("/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/on_intervals.bed", sep='\t',
                    names=['chr', 'start', 'end'])

df_on['int_len'] = df_on.end - df_on.start
total_on = df_on['int_len'].sum()
print(total_on)


df_off = pd.read_csv("/Users/Neil.Horner/work/workflow_outputs/cas9/grch38/off_intervals.bed", sep='\t',
                    names=['chr', 'start', 'end'])

df_off['int_len'] = df_off.end - df_off.start
total_off = df_off['int_len'].sum()
print(total_off)

import sys
import pandas as pd
infile = sys.argv[1]
df = pd.read_csv(infile, sep='\t')
temp = df.columns.tolist()
temp[0] = ''
df.columns = temp
df.to_csv(infile, sep='\t', index=False)

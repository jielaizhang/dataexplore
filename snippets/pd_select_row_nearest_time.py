#!/usr/bin/env python


import pandas as pd

# write text file to use in code snippet
txt = '''blah,10-9-2022@0828-1732,5,90
blah,10-9-2022@0853-0332,46,98
blah,11-9-2022@0853-0332,7,65
blah,12-9-2022@0853-0332,48,3
'''
with open('example.csv','w') as f:
    f.write(txt)
    
# Read in file, specify DatetimeStr column as index. 
a   = pd.read_csv('example.csv', 
                  names=['NA','DatetimeStr','val1','val2'], 
                  index_col='DatetimeStr', 
                  parse_dates=['DatetimeStr'], 
                  date_parser=lambda x: pd.to_datetime(x, format='%d-%m-%Y@%H%M-%S%f'))

# Date time you want to find row with closest, closest before, closest after datetime
xxx = datetime(2022, 9, 10, 7, 50, 39)

# Get the index (unique only) with closest date time occuring after xxx 
# (bfill = closest after, ffill = closest before, nearest = closest)
index = a.index.get_indexer([xxx],method="bfill")

# Get the row data only
row_data_only = a.iloc[index[0]]
print(row_data_only)

# get the row as a dataframe
row_df = a.iloc[index]
print('----')
print(row_df)
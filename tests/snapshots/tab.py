#!/usr/bin/env python

# Print the csv file as readable table:
# ./tab.py *.csv

import pandas as pd
import sys

myfile = sys.argv[1]

print(myfile)

df = pd.read_csv(myfile, index_col = 0, na_values = None)
print(df.where(df.notnull(), None).to_string(index = False, justify = 'left'))


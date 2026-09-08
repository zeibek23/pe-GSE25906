import pandas as pd
import numpy as np
import scipy 

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 250)

with open ("data/GSE25906_series_matrix.txt") as f:
    meta_lines = [line.rstrip("\n") for line in f][:28]

gsm_line = [l for l in meta_lines if l.startswith("!Sample_geo_accession")][0]

expression = pd.read_csv("results/01_expression.csv", index_col=0)

print(len(meta_lines))
print(expression.shape)
print(expression.index.dtype)
print(expression.head(10))
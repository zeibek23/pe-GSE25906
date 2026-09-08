import pandas as pd
import numpy as np
import scipy 

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 250)

with open ("data/GSE25906_series_matrix.txt") as f:
    meta_lines = [line.rstrip("\n") for line in f if line.startswith("!")]

gsm_line = [l for l in meta_lines if l.startswith("!Sample_geo_accession")][0]

expression = pd.read_csv("results/01_expression.csv", index_col=0)

print(len(meta_lines))
print(expression.shape)
print(expression.index.dtype)
print(expression.head(10))

def get_characteristic(meta_lines, key):
    matches = [l for l in meta_lines
               if l.startswith("!Sample_characteristics_ch1") and f'"{key}' in l]
    assert len(matches) == 1, f"{key}: found {len(matches)} lines"
    values = [p.strip('"').split(": ", 1)[1] for p in matches[0].split("\t")[1:]]
    return values

gsm_ids = [p.strip ('"') for p in gsm_line.split("\t")[1:]]
meta = pd.DataFrame({
    "gsm": gsm_ids,
    "diagnosis": get_characteristic(meta_lines, "classification:"),
    "gest_age": get_characteristic(meta_lines, "gestational age:"),
    "batch": get_characteristic(meta_lines, "batch:"),
    "gender": get_characteristic(meta_lines, "gender:"),
    "labour": get_characteristic(meta_lines, "induction of labor:"),
}).set_index("gsm")

assert all (meta.index == expression.columns)
print(meta.dtypes)
print(pd.crosstab(meta["batch"], meta["diagnosis"]))
print(pd.crosstab(meta["labour"], meta["diagnosis"]))
meta.to_csv("results/meta.csv")
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

expression = pd.read_csv("results/01_expression.csv", index_col=0)      
meta = pd.read_csv("results/meta.csv", index_col=0)
meta["gest_age"] = meta["gest_age"].astype(int)
annot = pd.read_csv("data/GPL6102.annot", sep="\t", skiprows=28,
                    skipfooter=1, engine="python")[["ID", "Gene symbol"]]

assert all(meta.index == expression.columns)

panel = ["CTH", "CBS", "MPST", "FLT1", "ENG", "SLC5A2"]
panel_probes = annot[annot["Gene symbol"].isin(panel)].set_index("ID")
print(panel_probes)              # expect ~7 rows: CTH has two probes

is_pe = (meta["diagnosis"] == "preeclamptic").values
FLOOR = 6.0                      # the threshold you chose in script 03

rows = []
for probe, symbol in panel_probes["Gene symbol"].items():
    x = expression.loc[probe]
    pe, ctrl = x[is_pe], x[~is_pe]
    t, p = stats.ttest_ind(pe, ctrl, equal_var=False)
    rows.append({
        "probe": probe, "symbol": symbol,
        "mean_PE": pe.mean(), "mean_ctrl": ctrl.mean(),
        "log2FC": pe.mean() - ctrl.mean(),                       # you know this one
        "pval": p,
        "detected": x.mean() > FLOOR,
    })

panel_df = pd.DataFrame(rows).set_index("probe").sort_values("symbol")
print(panel_df.round(3))
panel_df.to_csv("results/06_gene_panel.csv")

# --- 4. Strip plot grid ---
fig, axes = plt.subplots(2, 4, figsize=(14, 7))
colors = {"induced": "tab:orange", "spontaneous": "tab:blue"}

for ax, (probe, row) in zip(axes.flat, panel_df.iterrows()):
    x = expression.loc[probe]
    for i, (label, mask) in enumerate([("Control", ~is_pe), ("PE", is_pe)]):
        jitter = np.random.normal(i, 0.06, mask.sum())
        ax.scatter(jitter, x[mask],
                   c=meta.loc[mask, "labour"].map(colors), s=18, alpha=0.8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Control", "PE"])
    ax.set_title(f'{row["symbol"]}  ({probe})\nlog2FC={row["log2FC"]:.2f}  p={row["pval"]:.2g}',
                 fontsize=9)
    ax.axhline(FLOOR, ls=":", lw=0.8, c="grey")       # detection floor
    if not row["detected"]:
        ax.set_facecolor("#f4f4f4")                    # shade undetected panels

for ax in axes.flat[len(panel_df):]:
    ax.set_visible(False)                              # hide unused panels

fig.suptitle("GSE25906 gene panel — orange = induced, blue = spontaneous labour")
fig.tight_layout()
fig.savefig("results/06_gene_panel.png", dpi=300)
print("Saved panel figure")
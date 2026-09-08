import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

final = pd.read_csv("results/03_stats_annotated.csv")

padj_safe = np.maximum(final["padj"], 1e-300)
final["neglog10padj"] = -np.log10(padj_safe)

sig = (final["padj"] < 0.05) & (final["log2FC"].abs()>0.5)

plt.figure(figsize=(10,8))
plt.scatter(final["log2FC"], final["neglog10padj"], s=4, c="grey", alpha=0.4)
plt.scatter(final.loc[sig, "log2FC"], final.loc[sig, "neglog10padj"], s=6, c="red")

label_genes = {"LEP", "FLT1", "ENG", "INHA", "PAPPA2", "FSTL3", "HTRA4", "CTH"}
to_label = final[sig & final["Gene symbol"].isin(label_genes)]
for _, row in to_label.iterrows():
    plt.annotate(row["Gene symbol"], (row["log2FC"], row["neglog10padj"]), fontsize=7)
    
plt.axhline(-np.log10(0.05), ls="--", lw=0.8)        
plt.axvline(0.5, ls="--", lw=0.8)
plt.axvline(-0.5, ls="--", lw=0.8)
plt.xlabel("log2 fold change (Preeclamptic vs Control)")
plt.ylabel("-log10 adjusted p-value")
plt.title("GSE25906 Preeclamptic vs Control Analysis")
plt.tight_layout()
plt.savefig("results/04_volcano.png", dpi=1200)
print("Figure saved")
print(sig.sum())
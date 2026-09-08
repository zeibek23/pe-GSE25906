import pandas as pd

expression = pd.read_csv("results/01_expression.csv", index_col=0)
diagnosis = pd.read_csv("results/01_diagnosis.csv", index_col=0).squeeze()

print( expression.shape)
print(expression.head(5))
print(type(diagnosis))
print(len(diagnosis))
print(all(diagnosis.index == expression.columns))

overall_mean = expression.mean(axis=1)
expression_filter = expression[overall_mean > 6]
assert "ILMN_1705002" in expression_filter.index, "FLT1 was filtered out"
print(expression_filter)
expression = expression[~expression.index.str.startswith("ILMN_13432")]


from scipy import stats

is_pe = diagnosis == "preeclamptic"
group_pe = expression_filter.loc[:, is_pe.values]
group_non = expression_filter.loc[:, (~is_pe).values]

print(group_pe.shape, group_non.shape)

tstat, pval = stats.ttest_ind(group_pe, group_non, axis=1, equal_var=False)
print(len(pval))
print(pval[:5])
print(pd.isna(pval).sum())

pval = pd.Series(pval, index=expression_filter.index)
from statsmodels.stats.multitest import multipletests

reject, padj, _, _ = multipletests(pval.values, alpha=0.05, method="fdr_bh")
padj = pd.Series(padj, index=expression_filter.index)

print ((pval < 0.05). sum())
print ((padj < 0.05).sum())

mean_pe = group_pe.mean(axis=1)
mean_non = group_non.mean(axis=1)

stats_df = pd.DataFrame({
    "mean_PE": mean_pe,
    "mean_nonPE": mean_non,
    "log2FC": mean_pe-mean_non,
    "pval": pval,
    "padj": padj,
})

stats_df.to_csv("results/03_stats.csv")
print(stats_df.shape)
print(stats_df.head(10))
print(stats_df.sort_values("padj").head(20))


annotation = pd.read_csv("data/GPL6102.annot", sep="\t", skiprows=28, skipfooter=1, engine="python")
columns_to_keep = ["ID", "Gene symbol"]
annotation = annotation[columns_to_keep]
final = stats_df.merge(annotation, left_index=True, right_on="ID", how="left")
print(final.shape)
print(final.head(5))

top20 = final.sort_values("padj").head(20)
print(top20[["ID", "Gene symbol", "log2FC", "pval", "padj"]])
final.to_csv("results/03_stats_annotated.csv", index=False)
print("stats_annotated have been saved inside the results folder!!")


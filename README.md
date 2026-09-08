# pe-GSE25906 — Preeclampsia placental transcriptomics (Illumina HumanWG-6 v2)

Independent re-analysis of **GSE25906** (Tsai *et al.*, 2011): genome-wide expression profiling of
60 human term placentas, 23 preeclamptic (PE) and 37 control, on the Illumina Human-6 v2 BeadChip
(GPL6102). The analysis is written from scratch in Python/pandas as a learning project, and serves
as the second dataset in a two-cohort validation (the first was GSE75010, Affymetrix HuGene 1.0 ST).

**Biological question.** Which genes change in the PE placenta, does this dataset reproduce the
canonical PE signature (sFlt-1, sEng, LEP, INHA, PAPPA2 …), and what does bulk placental expression
say about the hydrogen sulfide–producing enzymes CTH (CSE), CBS and MPST?

---

## Data

| Item | Source |
|---|---|
| Expression matrix | `GSE25906_series_matrix.txt` (log2, quantile-normalised, 48,701 probes × 60 samples) |
| Probe annotation | `GPL6102.annot` (GEO curated annotation, 48,702 rows) |
| Sample metadata | parsed from the `!Sample_characteristics_ch1` lines of the series matrix |

`data/` and `results/` are git-ignored (raw files are ~33 MB each). To regenerate:

```bash
mkdir -p data results
cd data
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25906/matrix/GSE25906_series_matrix.txt.gz
gunzip GSE25906_series_matrix.txt.gz
curl -O https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6102/annot/GPL6102.annot.gz
gunzip GPL6102.annot.gz
cd ..
```

Sample metadata available per sample: classification (PE / Control), gestational age (weeks),
fetal sex, induction of labour (induced / spontaneous), hybridisation batch (A / B).

---

## Pipeline

All scripts are run **from the project root** (`python scripts/NN_name.py`); every path is relative
to the root. Each script reads the outputs of the previous one from `results/`.

| Script | Reads | Writes | What it does |
|---|---|---|---|
| `01_load_data.py` | series matrix | `01_expression.csv`, `01_diagnosis.csv`, `01_log2fc.csv` | Parses the 67-line GEO header, extracts GSM IDs and PE/Control labels, verifies label↔column alignment, computes group means and log2 fold change (log2 scale → subtraction) |
| `02_annotate.py` | `01_log2fc.csv`, GPL6102 annotation | `02_annotated_log2fc.csv` | Left-merges gene symbols onto probes (checks key dtypes and duplicates first); inspects H₂S-pathway and marker genes |
| `03_statistics.py` | `01_expression.csv`, `01_diagnosis.csv` | `03_stats.csv`, `03_stats_annotated.csv` | Removes Illumina housekeeping control probes (`ILMN_13432xx`), filters probes with mean log2 intensity ≤ 6, Welch's t-test per probe, Benjamini–Hochberg FDR, annotation |
| `04_valcano.py` | `03_stats_annotated.csv` | `04_volcano.png` | Volcano plot; significance = FDR < 0.05 and \|log2FC\| > 0.5; curated gene labels |
| `05_subgroups.py` | series matrix, `01_expression.csv` | `meta.csv` | `get_characteristic()` parser → sample metadata table (diagnosis, gestational age, batch, sex, labour); batch × diagnosis and labour × diagnosis crosstabs |
| `06_gene_panel.py` | `01_expression.csv`, `meta.csv`, annotation | `06_gene_panel.csv`, `06_gene_panel.png` | Targeted test of CTH, CBS, MPST, FLT1, ENG, SLC5A2 on the unfiltered matrix, with a detection-floor flag and per-sample strip plots coloured by labour onset |

---

## Key results

**Cohort.** 23 PE vs 37 Control, matching the published design. Metadata and expression columns
verified to be in the same order (`all(meta.index == expression.columns)`).

**Differential expression** (13,423 probes tested after control-probe removal and the mean > 6
filter): 3,304 probes at nominal *p* < 0.05, **1,354 at FDR < 0.05**, and **98** at FDR < 0.05 with
\|log2FC\| > 0.5.

**The canonical PE signature is reproduced.** Largest effects (all FDR < 0.05):

| Gene | log2FC | Fold | Note |
|---|---|---|---|
| LEP | +2.26 | 4.8× | Leptin — strongest signal in the dataset |
| FSTL3 | +1.28 | 2.4× | |
| PAPPA2 | +1.04 / +0.87 | ~2× | Two independent probes agree |
| ENG | +0.91 | 1.9× | Source of soluble endoglin (sEng) |
| INHA | +0.85 | 1.8× | Inhibin A — maternal serum marker |
| HTRA4 | +0.84 | 1.8× | |
| FLT1 | +0.55 | 1.5× | Source of sFlt-1; FDR ≈ 0.03 |
| SIAE | +0.44 | 1.4× | Most significant probe (FDR 8 × 10⁻⁶) — the gene in the original paper's title |

Genes lower in PE include HBD, DEFA1B and LYZ (haemoglobin, neutrophil defensin, lysozyme) —
a maternal-blood / inflammatory signature more plausibly linked to mode of delivery than to disease
(see confounders below).

**Hydrogen sulfide enzymes are below the array's detection floor in bulk placenta.**

| Gene | Probe(s) | Mean log2 (PE / Ctrl) | log2FC |
|---|---|---|---|
| CTH (CSE) | ILMN_1777060, ILMN_1784112 | 5.50 / 5.53, 5.72 / 5.75 | −0.03, −0.02 |
| CBS | ILMN_1804735 | 6.52 / 6.52 | 0.00 |
| MPST | ILMN_1800096 | 6.35 / 6.24 | +0.11 |
| SLC5A2 (SGLT2) | ILMN_1666972 | 5.68 / 5.70 | −0.02 |

The array-wide median probe mean is 5.69 and 72% of probes fall below 5.7; CTH sits inside that
noise cluster. The correct statement is *undetected*, not *unchanged* — see Interpretation.

**Confounders.**

| | Control | PE |
|---|---|---|
| Batch A | 8 | 7 |
| Batch B | 29 | 16 |
| Induced labour | 8 | 16 |
| Spontaneous labour | 29 | 7 |

Batch is balanced (the authors randomised samples across chips), so batch is noise rather than a
confound. **Mode of delivery is strongly confounded with diagnosis** (70% of PE inductions vs 22% of
controls), as is gestational age (PE deliveries range 27–38 weeks; controls cluster at 37–40).
Neither is adjusted for in the current t-test; this is the main limitation.

---

## Interpretation and caveats

* Fold changes are computed on log2, quantile-normalised data (difference of means).
* The Illumina housekeeping control probes (RPS9, UBC, EEF1A1, TUBB2A, TXN, ACTB, GAPDH) are floored
  to a constant value in every batch A sample and carry real values in batch B. They were removed
  before testing. Some additional probes (e.g. ILMN_20387xx) show the same batch-A floor pattern and
  should be removed by pattern rather than by ID prefix (open item).
* The detection floor of 6.0 is approximately the 75th percentile of probe means, just above the
  modal noise cluster; it is a pragmatic choice, not derived from Illumina detection p-values (which
  are only in the non-normalised supplementary file). A sensitivity check at 5.5 / 6 / 6.5 is planned.
* "Not detected" for CTH means the bulk signal is indistinguishable from background on this
  platform. CSE protein is readily demonstrated in placental vasculature by IHC; vascular smooth
  muscle and endothelium are a small fraction of villous tissue mass, so a cell-type-restricted
  transcript is diluted below array sensitivity by trophoblast RNA. The same pattern was seen in
  GSE75010 (Affymetrix), i.e. it is not a platform artefact. Bulk arrays are the wrong instrument for
  this question; targeted qPCR, IHC or single-cell datasets are the right ones.
* Results are probe-level, not gene-level; several genes have multiple probes.
* Illumina probes are single 50-mers near the canonical 3′ end. The FLT1 probe may not capture the
  intronic-polyadenylation sFlt-1 isoforms (sFLT1-i13, sFLT1-e15a), which is one reason the FLT1
  effect is modest here compared with Affymetrix exon-tiled arrays.

---

## Environment

Conda env `pe-learning`: Python 3.12, pandas, numpy, scipy, statsmodels, matplotlib.

```bash
conda activate pe-learning
python scripts/01_load_data.py
python scripts/02_annotate.py
python scripts/03_statistics.py
python scripts/04_valcano.py
python scripts/05_subgroups.py
python scripts/06_gene_panel.py
```

---

## Lessons logged (Protocol 1 — silent-failure catalogue)

* Run from the project root; relative paths resolve against the shell, not the script.
* A string filter that matches nothing (missing `!`, misspelt key) returns empty and does not complain — assert on the match count.
* A filter applied to the wrong variable, or after the derived table was built, runs clean and changes nothing — verify by checking the thing that should have disappeared.
* An expression filter can delete your positive control (FLT1 at `mean > 7`); assert it survived.
* Housekeeping genes are not batch-proof: ACTB and GAPDH showed a −0.7 log2FC purely from batch-A flooring.
* Everything parsed from GEO metadata is a string — `astype(int)` before any numeric comparison.
* CSV round-trips: `index_col=` on read, check merge-key dtypes and duplicates before merging.
* Correct output does not prove correct code.

---

## Next steps

1. Pattern-based removal of all batch-A-floored probes (`expression.loc[:, meta.batch=="A"].nunique(axis=1) == 1`).
2. Fetal-sex sanity check (XIST vs RPS4Y1/DDX3Y against the `gender` column).
3. Per-gene linear model `expr ~ diagnosis + gest_age + labour + batch` (statsmodels) to replace the unadjusted t-test.
4. Early- (< 34 wk) vs late-onset PE subgroup for FLT1 and ENG.
5. Detection-floor sensitivity analysis; optionally pull Illumina detection p-values from `GSE25906_non-normalized.txt.gz`.
6. GO enrichment on the FDR < 0.05 set after collapsing probes to genes.

---

## References

* Tsai S, Hardison NE, James AH, *et al.* Transcriptional profiling of human placentas from pregnancies complicated by preeclampsia reveals disregulation of sialic acid acetylesterase and immune signalling pathways. *Placenta* 2011. PMID 21183218. GEO: GSE25906.
* Leavey K, *et al.* Unsupervised placental gene expression profiling identifies clinically relevant subclasses of human preeclampsia. *Hypertension* 2016. GEO: GSE75010 (re-analysed the GSE25906 samples).

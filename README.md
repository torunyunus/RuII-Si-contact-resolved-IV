# Contact-resolved Ru(II)/p-Si Schottky diode I–V measurements

This repository provides the dataset and analysis code for the contact-resolved comparison of reverse-bias photocurrent and photovoltaic power in Al/Ru(II)-complex/p-Si/Al diodes. It contains **50 dark/illuminated measurement pairs (100 instrument-export curves)**: 30 contacts from the 1000 rpm specimen and 20 from the 2000 rpm specimen.

## Contents

- `raw_iv/1000/` and `raw_iv/2000/`: the signed dark and illuminated I–V records used in the article.
- `PAIRING_MAP.csv`: specimen identifier, contact identifier, paired source files and original archive-member names.
- `MANIFEST_SHA256.csv`: file paths, sizes and SHA-256 checksums.
- `contact_metrics.csv` and `rank_summary.csv`: contact-level and specimen-level rank comparisons.
- `derived/`: supplementary calculations and figures, including the per-contact leave-one-out Spearman analysis (`LOCO_spearman.csv`).
- `reproduce_rank_analysis.py`, `reproduce_step2.py`, `reproduce_step3.py`, `reproduce_superposition.py`: analysis scripts.
- `requirements.txt`: Python dependencies.

The original instrument exports are preserved byte for byte. Each `.dat` file contains recorded voltage (V), current magnitude (A), and **signed current (A)** in its first three numeric columns. Analyses use the signed-current column. The pairing map links the repository file names to the retained instrument-export names. Contact identifiers label the measured diode contacts; they are not a verified spatial map.

## Reproduce the main and supplementary calculations

With Python 3.10+ installed, run:

```bash
python -m pip install -r requirements.txt
python reproduce_rank_analysis.py
python reproduce_step2.py . --out derived/step2
python reproduce_step3.py . --out derived/step3
python reproduce_superposition.py
```

`reproduce_rank_analysis.py` also generates `derived/LOCO_spearman.csv` with one leave-one-contact-out Spearman coefficient per contact. The remaining scripts reproduce the voltage-dependent controls, contact-specific quality assessments, and the normalized photocurrent/superposition analysis. See the accompanying paper and electronic supplementary material for the definitions and interpretation of each derived quantity.

## Data scope and citation

The dataset includes the measurements actually analyzed in the article; unrelated preparation speeds and other characterization records are not deposited. For citation and reproducibility, refer to the **specific tagged GitHub Release** and its commit SHA, rather than the mutable default branch. The paper DOI and Release URL can be added to `CITATION.cff` when assigned. GitHub is not a guaranteed permanent data archive.

The data-sharing permissions and licence should be agreed by the experimental data rights holders and co-authors before public release.

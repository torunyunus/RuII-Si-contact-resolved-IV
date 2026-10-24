# Contact-resolved Ru(II)/p-Si Schottky diode I–V dataset

This repository accompanies the study of reverse-bias photocurrent and photovoltaic power from contact-matched measurements of Al/Ru(II)-complex/p-Si/Al Schottky diodes.

## Dataset scope

The prepared deposition package contains **50 matched dark/illuminated current–voltage pairs**: 30 contacts from the 1000 rpm specimen and 20 from the 2000 rpm specimen, for a total of 100 signed I–V export files. It also includes a contact pairing map, SHA-256 checksums, derived metrics, supplementary figures and Python reproduction scripts. No 3000 rpm records form part of this dataset.

## Repository status

**Dataset upload pending.** This README does not mean that the measurement files, code or a tagged release have been uploaded. Please do not cite a release until the files and version are publicly available and validated.

## Expected contents

- `raw_iv/1000/`, `raw_iv/2000/` — measured signed dark/light I–V exports
- `PAIRING_MAP.csv` and `MANIFEST_SHA256.csv` — source pairing and integrity manifest
- `contact_metrics.csv`, `rank_summary.csv`, `derived/` — contact-level and sensitivity results
- `reproduce_rank_analysis.py`, `reproduce_step2.py`, `reproduce_step3.py`, `reproduce_superposition.py` — reproducibility scripts
- `requirements.txt` and `CITATION.cff`

Once uploaded, verify that all 100 I–V files and accompanying outputs are present, run the four scripts, and create a versioned release (e.g. `v1.0.0`). Data-sharing permission should be confirmed by the relevant authors before public deposition.

## Related manuscript

*Contact-resolved comparison of reverse-bias photocurrent and photovoltaic power in Ru(II)-interlayer/p-Si Schottky diodes.*

This repository is the intended code/data access point; the article DOI will be linked after publication.

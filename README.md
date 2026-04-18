## Feature planning helper

This repo now includes `feature_pipeline.py`, which prepares a clean descriptor list for your peptide ML workflow.

### What it does
- Loads sequences from Excel and uppercases them.
- Starts with your requested BBPpred-style descriptors:
  - AAC, BIT188, BIT12, BIT21, CTD, GDPC, GTPC, IT, OLP, DPC, ASDC, GGAP, PAAC, CTF, BIT20, QSO
- Excludes descriptors you already have from Modlamp/iFeature/pFeature.
- Optionally forces PAAC back in with dimension formula: **20 + 2λ**.

### Run
```bash
python feature_pipeline.py
```

If `combined.xlsx` is in the working directory and contains a `sequence` column, it prints a preview and the final feature plan.

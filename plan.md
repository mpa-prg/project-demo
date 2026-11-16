# CP1 plan: demo (replication origin from GC skew)

An example of a plan. A real CP1 plan is longer, because the core of a project topic has more parts.

## 1. Pseudocode of the core

```
FIND-ORIGIN(s)                       # s = circular chromosome, n = |s|
    c[0] ← 0
    for i ← 1 to n:
        c[i] ← c[i-1] + [s_i = G] − [s_i = C]
    T ← c[n]                          # total skew, rarely exactly 0
    for i ← 0 to n − 1:
        d[i] ← c[i] − i · T / n       # remove the trend so the curve is circular
    ori ← argmin_i d[i] + 1
    ter ← argmax_i d[i] + 1
    return ori, ter
```

The key difficulty: the chromosome is circular, so the answer must not depend on where the sequence file starts. Without the trend correction, the minimum can jump to the end of the file.

Time O(n), memory O(n).

## 2. Validation plan

- Simulation: a circular genome of 500 kb; from the origin to the terminus, G is more frequent than C by the factor (1 + skew) / (1 − skew), on the other half the opposite. The true origin is random and known.
- Varied parameter: skew 0.002, 0.005, 0.01, 0.02, 0.05; 30 genomes (seeds) for each.
- Metric: circular distance between the found and the true origin (bp); median and spread per skew.
- Later, if time allows: GC content, unequal replication arms.

## 3. Real data

- U00096.3 (*E. coli* K-12 MG1655) from ENA, downloaded by `download_data.py`, MD5 checked.
- One record, 4 641 652 bp, GC content 50.8 %.
- Comparison: the annotated oriC (3 925 744–3 925 975) from the annotation of the same ENA record.

## 4. Split of work

| Member | Functions and modules | Experiments |
|---|---|---|
| A | `skew.py`, `find_origin.py`, tests of the core | real genome |
| B | `simulate.py` and its tests | skew experiment |
| C | `download_data.py`, `run_all.py`, `experiments.py` (figures) | GC-content experiment |

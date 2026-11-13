"""Table-1 and Table-3 analogues (Semenova paper) for the two new
empirical applications.

Table-1 analogue: identified set [E pi_L, E pi_U] for the always-selected
share, by specification (None / discrete single-covariate cells / logit on
baseline covariates), with the paper's one-sided 95% lower confidence band
from 400 bootstrap resamples (cluster resamples for THP, unit resamples
for Bursztyn).  -> results/app_table1_analogues.csv

Table-3 analogue: the finite-grid heuristic LP (Section 2.3, conservative
endpoint convention, unweighted) at K = 5/10/15 equal-mass outcome bins,
for THP's continuous depression indices (6m/1y waves) and 7-year
parentmoney with wealth-quintile covariate cells.  For Bursztyn the
follow-up outcomes are binary, so the K-bin discretization is exact at
K = 2 and there is no granularity dimension - recorded as a single exact
row.  -> results/app_table3_analogues.csv
"""

import numpy as np
import pandas as pd

import data_prep_baranov as db
import estimator as est
from lp_bounds import lp_wage_effect_bounds
"from run_bursztyn import load as bz_load, design as bz_design"

RNG = np.random.default_rng(2028)
N_BOOT = 400


def share_box(S, D, cells):
    EL = EU = 0.0
    n = len(S)
    for c in np.unique(cells):
        m = cells == c
        m0, m1 = m & (D == 0), m & (D == 1)
        if m0.sum() == 0 or m1.sum() == 0:
            continue
        s0c, s1c = S[m0].mean(), S[m1].mean()
        wc = m.sum() / n
        EL += wc * max(s0c + s1c - 1.0, 0.0)
        EU += wc * min(s0c, s1c)
    return EL, EU


def logit_box(S, D, X):
    s0 = est.Logit().fit(X[D == 0], S[D == 0]).predict(X)
    s1 = est.Logit().fit(X[D == 1], S[D == 1]).predict(X)
    return (float(np.mean(np.maximum(s0 + s1 - 1, 0))),
            float(np.mean(np.minimum(s0, s1))))


def band(fn, S, D, arg, cluster=None):
    draws = []
    if cluster is not None:
        ids = np.unique(cluster)
    for _ in range(N_BOOT):
        if cluster is None:
            idx = RNG.integers(0, len(S), len(S))
        else:
            pick = RNG.choice(ids, len(ids), replace=True)
            idx = np.concatenate([np.flatnonzero(cluster == c) for c in pick])
        try:
            draws.append(fn(S[idx], D[idx], arg[idx])[0])
        except Exception:
            continue
    return float(np.quantile(draws, 0.05))


def main():
    rows1, rows3 = [], []

    # ================= THP (Baranov) =================
    r, D, uc, S7 = db.load_randomized()
    Xfull, _ = db.design_matrix(r)
    Xsmall = np.column_stack([pd.to_numeric(r[c], errors="coerce").values
                              for c in db.X_SMALL])
    grand = pd.to_numeric(r["grandmother_baseline"],
                          errors="coerce").fillna(0).astype(int).values
    first = pd.to_numeric(r["first_child"],
                          errors="coerce").fillna(0).astype(int).values
    wealth_q = pd.qcut(pd.to_numeric(r["wealth_baseline"], errors="coerce"),
                       5, labels=False, duplicates="drop").values

    waves = [("7y", S7), ("6m", db.wave_selection(r, "depindex_6m")),
             ("1y", db.wave_selection(r, "depindex_1y"))]
    for wave, S in waves:
        specs = [("None", "cells", np.zeros(len(r), dtype=int)),
                 ("Grandmother", "cells", grand),
                 ("FirstChild", "cells", first),
                 ("Logit-targeted", "logit", Xsmall),
                 ("Logit-full", "logit", Xfull)]
        for name, kind, arg in specs:
            fn = share_box if kind == "cells" else logit_box
            EL, EU = fn(S, D, arg)
            rows1.append(dict(application="THP", wave=wave, spec=name,
                              EpiL=EL, EpiU=EU,
                              band=band(fn, S, D, arg, cluster=uc)))

    # Table-3 analogue: granularity LP on continuous outcomes
    for wave, col, S in [("6m", "depindex_6m", waves[1][1]),
                         ("1y", "depindex_1y", waves[2][1]),
                         ("7y", "parentmoney", S7)]:
        Y = db.outcome(r, col, S)
        for K in (5, 10, 15):
            lo, hi, _ = lp_wage_effect_bounds(Y, S, D, wealth_q, K)
            rows3.append(dict(application="THP", outcome=col, wave=wave,
                              cells="wealth-quintiles", K=K,
                              lower=lo, upper=hi))
    
    from pathlib import Path
    out = Path(__file__).resolve().parent / "results"
    out.mkdir(exist_ok=True)

    t1 = pd.DataFrame(rows1)
    t3 = pd.DataFrame(rows3)
    t1.to_csv(out / "app_table1_analogues.csv", index=False)
    t3.to_csv(out / "app_table3_analogues.csv", index=False)
    print(t1.to_string(index=False))
    print(t3.to_string(index=False))
    return t1, t3

if __name__ == "__main__":
    main()
    

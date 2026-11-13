"""Main estimator sweep: Miguel-Kremer PSDP deworming test-score effect
among always-tested pupils, with and without monotone-attendance.

Columns per outcome (only the ICS composite here):
  (1) ITT on tested       -- naive selected-sample difference (the null)
  (2) Lee, monotone       -- monotone test-day attendance (S(1) >= S(0))
  (3) no mon., marginal   -- exact fractional-atom bounds, no covariates
  (4) no mon., covariates -- Section 5 debiased estimator, pupil+school X
Cluster = school (75); cluster multiplier bootstrap / IM regions.
Also a Table-1 analogue: E[pi_L], E[pi_U] by specification with bands.
"""

import numpy as np
import pandas as pd

import data_prep_deworming as dd
import estimator as est
from run_baranov import lee_bidirectional, itt_respondents, fmt

SEED = 42
N_BOOT = 400


def cluster_folds(schid, D, n_folds=3, seed=SEED):
    rng = np.random.default_rng(seed)
    fold = {}
    for arm in (0, 1):
        ids = np.unique(schid[D == arm])
        perm = rng.permutation(len(ids))
        for j, c in enumerate(ids):
            fold[c] = perm[j] % n_folds
    return np.array([fold[c] for c in schid])


def main():
    df, D, schid, S, Y = dd.load_frame()
    w = np.ones(len(df))
    folds = cluster_folds(schid, D)
    Xp, _ = dd.design(df, dd.X_PUPIL + dd.X_SCHOOL)
    s1, s0 = S[D == 1].mean(), S[D == 0].mean()
    print(f"Deworming test-frame N={len(df)}; schools={len(np.unique(schid))}"
          f" ({(np.unique(schid[D==1]).size)}T/{np.unique(schid[D==0]).size}C)")
    print(f"selection: s0={s0:.4f} s1={s1:.4f} sum={s0+s1:.4f} -> "
          f"box [{max(s0+s1-1,0):.3f}, {min(s0,s1):.3f}]\n")

    rows = []
    beta, se = itt_respondents(Y, S, D, w, schid)
    print(f"  (1) ITT on tested   : {beta:+.3f} (cluster se {se:.3f})  "
          f"[M&K null]")
    r2 = lee_bidirectional(Y, S, D, w, schid)
    print(f"  (2) Lee (monotone)  : {fmt(r2)}")
    r3 = est.exact_marginal_bounds(Y, S, D, w, cluster=schid, n_boot=N_BOOT)
    print(f"  (3) no mon, marginal: {fmt(r3)}  [exact]")
    cache = est.first_stage_cache(Xp, Y, S, D, w, folds, qr_alpha=1e-2)
    r4 = est.bound_interval(Xp, Y, S, D, w, cache, cluster=schid,
                            n_boot=N_BOOT)
    ex = ""
    if r4.get("status") == "ok":
        rs = r4["regime_shares"]
        ex = f"  regimes {rs['int']:.2f}/{rs['flr']:.2f}/{rs['cei']:.2f}"
    print(f"  (4) no mon, X       : {fmt(r4)}{ex}")

    for lab, r in [("itt", dict(lower=beta, upper=beta, se_lower=se,
                                se_upper=se, IM_lo=beta - 1.96 * se,
                                IM_hi=beta + 1.96 * se, status="ok")),
                   ("lee_monotone", r2), ("nomon_marginal", r3),
                   ("nomon_X", r4)]:
        base = dict(outcome="ics98", spec=lab, s0=s0, s1=s1)
        if isinstance(r, dict) and r.get("status") == "ok":
            base.update(lower=r["lower"], upper=r["upper"],
                        IM_lo=r["IM_lo"], IM_hi=r["IM_hi"],
                        se_lower=r.get("se_boot_lower", r.get("se_lower")),
                        se_upper=r.get("se_boot_upper", r.get("se_upper")),
                        status="ok")
        else:
            base.update(status=r.get("status"))
        rows.append(base)

    # Table-1 analogue: does the always-tested stratum exist / where?
    print("\nTable-1 analogue: E[pi_L], E[pi_U] by specification")
    t1 = []

    def cell_epi(cells):
        EL = EU = 0.0
        n = len(cells)
        for cv in np.unique(cells):
            m = cells == cv
            m0, m1 = m & (D == 0), m & (D == 1)
            if m0.sum() < 5 or m1.sum() < 5:
                continue
            s0c, s1c = S[m0].mean(), S[m1].mean()
            wc = m.sum() / n
            EL += wc * max(s0c + s1c - 1, 0)
            EU += wc * min(s0c, s1c)
        return EL, EU

    def logit_epi(X):
        p0 = est.Logit().fit(X[D == 0], S[D == 0]).predict(X)
        p1 = est.Logit().fit(X[D == 1], S[D == 1]).predict(X)
        return (float(np.mean(np.maximum(p0 + p1 - 1, 0))),
                float(np.mean(np.minimum(p0, p1))))

    specs = [("None", "cells", np.zeros(len(df), dtype=int)),
             ("Sex", "cells", df.sex.fillna(-1).astype(int).values),
             ("Grade", "cells", df.grade.fillna(-1).astype(int).values),
             ("Logit pupil+school", "logit", Xp)]
    for name, kind, arg in specs:
        EL, EU = (cell_epi(arg) if kind == "cells" else logit_epi(arg))
        t1.append(dict(spec=name, EpiL=EL, EpiU=EU))
        print(f"  {name:20s}: E[piL]={EL:.4f} E[piU]={EU:.4f}")

    from pathlib import Path
    path = Path(__file__).resolve().parent / "results"
    path.mkdir(exist_ok=True)
    
    pd.DataFrame(rows).to_csv(path / "deworming_bounds.csv", index=False)
    pd.DataFrame(t1).to_csv(path / "deworming_table1.csv", index=False)
    print("\nSaved deworming_bounds.csv, deworming_table1.csv")


if __name__ == "__main__":
    main()

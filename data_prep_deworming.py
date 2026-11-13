"""Data mapping for the Miguel-Kremer (2004) PSDP deworming application.

Setting -> paper notation:
  - Sample: the 1998 pupil-questionnaire / ICS-test frame (N = 9,905), the
    pupils sampled for academic testing.  Sampling into this frame is
    balanced across arms (comparison 0.344 / treated 0.327 of eligibles),
    so it is a valid pre-selection frame.  The 20k eligibles never in the
    test frame are a sampling-design exclusion, not attrition.
  - Treatment D = treat98 = 1 if the school's deworming group is 1 or 2
    (assigned to the program), 0 if group 3 (1998 comparison).  School-
    level (cluster) randomization; 50 vs 25 schools.
  - Selection S = tested98 = present/scored on the 1998 ICS test day.
    s0 = 0.914, s1 = 0.922, s0+s1 = 1.835 -> THP-type (floor 0.835 > 0):
    the always-tested stratum exists marginally; deworming raised test-day
    presence only +0.8pp, so Lee monotonicity here is S(1) >= S(0).
  - Outcome Y = ics98, the authors' combined English/Math/Science ICS
    score (table10), standardized to comparison-group mean 0 / sd 1;
    missing when not tested.
  - Covariates X: baseline pupil (sex, age, grade, SES index) + school
    (baseline infection z_inf98, school population) -- all complete within
    the test frame.
  - Estimand beta = E[Y(1)-Y(0) | tested under both arms]: the deworming
    effect on academic scores among always-tested pupils, WITHOUT the
    monotone-attendance assumption.  The paper reports a null (naive ITT
    -0.092) and notes its Manski bounds were "uninformative" (fn 34); no
    Lee/covariate-assisted bounds were reported.
"""

import numpy as np
import pandas as pd

CSV = "/Users/yaroslavkorobka/Documents/PhD/Thesis/Code/Dinkelbach/code/deworming_clean.csv"
X_PUPIL = ["sex", "age", "grade", "ses_index"]
X_SCHOOL = ["z_inf98", "pup_pop"]


def load_frame():
    df = pd.read_csv(CSV)
    df = df[df["sex"].notna()].reset_index(drop=True)   # test frame
    D = df["treat98"].astype(int).values
    schid = df["schid"].astype(int).values
    S = df["tested98"].astype(int).values
    Y = np.where(S == 1, pd.to_numeric(df["testscore98"], errors="coerce"),
                 np.nan).astype(float)
    return df, D, schid, S, Y


def design(df, cols):
    out, names = [], []
    for c in cols:
        v = pd.to_numeric(df[c], errors="coerce").values.astype(float)
        v = np.where(np.isnan(v), np.nanmedian(v), v)
        if len(np.unique(v)) > 2:
            lo, hi = np.nanpercentile(v, [1, 99])
            v = np.clip(v, lo, hi)
        out.append(v)
        names.append(c)
    X = np.column_stack(out)
    keep = X.std(axis=0) > 0
    return X[:, keep], [n for n, k in zip(names, keep) if k]

"""Data mapping for the Baranov et al. (AER 2020) THP application.

Setting -> paper notation:
  - Sample: the 903 cluster-randomized prenatally depressed mothers
    (THP_sample == 1); the 300 non-randomized non-depressed 7-year
    comparison mothers (sample == 3) are excluded.
  - Treatment D = `treat` (Thinking Healthy Program CBT offer), randomized
    at the Union Council level: 40 clusters (`uc`), 20 treated / 20
    control, treatment constant within cluster.
  - Selection S = 1 - attrit2 at the 7-year wave (re-interviewed: 585 of
    903; s1 = 289/463 = 0.624, s0 = 296/440 = 0.673, so the marginal
    Frechet floor s0+s1-1 = 0.297 > 0: the always-responder stratum exists
    WITHOUT covariates, unlike JobCorps).  Earlier waves via outcome
    non-missingness: 6m (s0+s1 = 1.81), 12m (1.75).
  - Outcome Y: 7-year Anderson indices standardized to control mean 0/sd 1
    (values blanked for non-respondents, matching the replication
    package's own `replace var = . if sample != 1` before leebounds).
  - Covariates X: the paper's own attrition-IPW control list
    $X_control_vars MINUS `doi0` (baseline interview date - missing for
    exactly the 119 post-1y attriters, i.e. selection-revealed) and never
    `girl` (child gender - missing for exactly the 142 attriters without a
    live birth observed).  All retained covariates are complete for the
    903.
  - Weights: none (w = 1); estimand is the ATE among always-responders,
    beta = E[Y(1) - Y(0) | S(0) = S(1) = 1].
"""

import numpy as np
import pandas as pd

CSV = ("/Users/yaroslavkorobka/Documents/PhD/Thesis/Code/Dinkelbach/empirical app/maternalDepression-1.2/dataClean/THP_clean.csv")

# paper's $X_control_vars minus doi0 (selection-revealed missingness);
# age_baseline_sq is re-derived below
X_CONTROLS = ["age_baseline", "employed_mo_baseline", "mo_emp",
              "grandmother_baseline", "MIL", "wealth_baseline",
              "edu_lvl_mo_1", "edu_lvl_mo_2", "edu_lvl_mo_3",
              "edu_mo_baseline", "edu_fa_baseline", "kids_no",
              "first_child", "hamd_baseline", "mspss_baseline"]

# targeted covariates for the conditional column: the baseline variables
# the paper's own attrition analysis (Table A1) finds predictive of
# attrition (empowerment, employment, social support) plus wealth,
# baseline severity, and education
X_SMALL = ["mo_emp", "employed_mo_baseline", "mspss_baseline",
           "wealth_baseline", "hamd_baseline", "edu_mo_baseline"]

# 7-year outcomes (label, column, higher-is)
OUTCOMES_7Y = [
    ("Mother financial empowerment", "motherfinancial", "better"),
    ("Mother financial empowerment (7y only)", "motherfinancial_7y", "better"),
    ("Parental investment (monetary)", "parentmoney", "better"),
    ("Parental investment (time)", "parenttime", "better"),
    ("Parenting style", "parentstyle", "better"),
    ("Child physical development", "healthindex", "better"),
    ("Child cognitive development", "cognindex", "better"),
    ("Child socio-emotional development", "emoindex", "better"),
    ("Depression index (7y)", "depindex_7y", "worse"),
]
OUTCOMES_EARLY = [
    ("Depression index (6m)", "depindex_6m", "worse"),
    ("Depression index (1y)", "depindex_1y", "worse"),
]


def load_randomized():
    df = pd.read_csv(CSV, low_memory=False)
    r = df[df.THP_sample == 1].reset_index(drop=True)
    assert len(r) == 903
    D = r["treat"].astype(int).values
    uc = r["uc"].astype(int).values
    S7 = (r["attrit2"] == 0).astype(int).values
    return r, D, uc, S7


def design_matrix(r):
    """Winsorized, complete-case (verified complete) baseline design."""
    cols, names = [], []
    for c in X_CONTROLS:
        v = pd.to_numeric(r[c], errors="coerce").values.astype(float)
        assert not np.isnan(v).any(), f"{c} unexpectedly missing"
        if len(np.unique(v)) > 2:      # winsorize continuous at 1%/99%
            lo, hi = np.nanpercentile(v, [1, 99])
            v = np.clip(v, lo, hi)
        cols.append(v)
        names.append(c)
    age = pd.to_numeric(r["age_baseline"], errors="coerce").values
    cols.append((age ** 2).astype(float))
    names.append("age_baseline_sq")
    X = np.column_stack(cols)
    keep = X.std(axis=0) > 0
    return X[:, keep], [n for n, k in zip(names, keep) if k]


def outcome(r, col, S):
    """Outcome observed iff S = 1 (blank any attriter values, per the
    replication package's own convention before leebounds)."""
    Y = pd.to_numeric(r[col], errors="coerce").values.astype(float)
    Y = np.where(S == 1, Y, np.nan)
    return Y


def wave_selection(r, col):
    """Selection for the 6m/1y waves = outcome non-missingness."""
    return pd.to_numeric(r[col], errors="coerce").notna().astype(int).values

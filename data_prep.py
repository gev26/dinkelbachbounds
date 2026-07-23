"""Data preparation for the JobCorps application.

Variable mapping (validated against the paper's reported moments):
  - N = 9,145 applicants (Lee 2009 sample).
  - Treatment  D = TREATMNT.y  (random assignment to the program group).
  - Selection  S_w = 1{EARNH_w > 0 and HWH_w > 0}: employed with a recorded
    wage and hours in week w.  This reproduces the pooled selection rates
    reported in Section 3.1: s0+s1 = 0.77 / 0.92 / 0.96 at weeks 45/90/100.
  - Outcome    Y_w = log(EARNH_w / HWH_w) on S_w = 1 (log hourly wage =
    weekly earnings over weekly hours).
  - Design weights = DSGN_WGT.x (study design weights; used wherever the paper
    says "computations use study design weights"; the heuristic LP table is
    computed unweighted per its table notes).

Covariate sets:
  - "Earnings"  : annual baseline earnings EARN_YR, quantile-binned into five
                  cells (the point mass at zero collapses the five quantile
                  edges to four distinct cells; missing EARN_YR forms its own
                  cell).  Used by the discrete-heuristic LP table.
  - "Income"    : household income HH_INC on its native brackets (five recorded
                  brackets + a missing bracket = six cells).  Used by the
                  discrete-heuristic LP table.
  - "Reasons"   : the eight "Reasons for joining" baseline items
                  (R_HOME, R_COMM, R_TRAIN, R_CRGOAL, R_GETGED, R_NOWORK,
                  R_OTHER and the most-important-reason code MOSTIMPR).
  - "Lee28"     : 28 baseline covariates in the spirit of Lee (2009) /
                  Semenova (2024): demographics, marital status, children,
                  education, work history, earnings, household and personal
                  income brackets.

Missing covariate values are median-imputed with missingness indicators
appended (an implementation choice the paper does not specify; flagged in the
run reports).
"""

import numpy as np
import pandas as pd

DATA_PATH = "/Users/gevorgkhandamiryan/Desktop/DinkelbachBounds/dataLee2009.csv"

REASONS_ITEMS = ["R_HOME", "R_COMM", "R_TRAIN", "R_CRGOAL", "R_GETGED",
                 "R_NOWORK", "R_OTHER", "MOSTIMPR"]

LEE28 = ["FEMALE", "AGE", "BLACK", "HISP", "OTHERRAC",
         "NEVERMARRIED", "MARRIED", "TOGETHER", "SEPARATED",
         "HASCHLD", "NCHLD", "HGC", "HGC_MOTH", "HGC_FATH",
         "EVERWORK", "YR_WORK", "EARN_YR", "WKEARNR", "WELF_KID",
         "HH_INC1", "HH_INC2", "HH_INC3", "HH_INC4", "HH_INC5",
         "PERS_INC1", "PERS_INC2", "PERS_INC3", "PERS_INC4"]


def load_data():
    df = pd.read_csv(DATA_PATH, low_memory=False)
    assert len(df) == 9145
    return df


def outcome_selection(df, week):
    """Return (S, Y) at horizon `week`.

    S = 1{EARNH_w > 0 and HWH_w > 0} (employed with a recorded wage);
    Y = log hourly wage = log(weekly earnings EARNH_w / weekly hours HWH_w).
    Validated: weighted arm-means of Y at weeks 45/90 reproduce the basic-Lee
    column of the paper's Table 5 (+0.027 point at week 45, [+0.046,+0.049]
    at week 90), and the pooled selection rates match Section 3.1.
    """
    earnh = pd.to_numeric(df[f"EARNH{week}"], errors="coerce").values
    hwh = pd.to_numeric(df[f"HWH{week}"], errors="coerce").values
    S = ((earnh > 0) & (hwh > 0)).astype(int)
    with np.errstate(divide="ignore", invalid="ignore"):
        Y = np.where(S == 1, np.log(earnh / hwh), np.nan)
    return S, Y


def treatment(df):
    return df["TREATMNT.y"].astype(int).values


def weights(df):
    return pd.to_numeric(df["DSGN_WGT.x"], errors="coerce").values


def earnings_cells(df):
    """Baseline annual earnings, quantile-binned into five cells.

    The mass point at EARN_YR = 0 (~38% of the sample) collapses the five
    quantile edges to four distinct cells; missing EARN_YR is its own cell.
    """
    earn = pd.to_numeric(df["EARN_YR"], errors="coerce")
    cells = pd.qcut(earn, 5, labels=False, duplicates="drop")
    out = cells.copy()
    out[earn.isna()] = -1  # missing cell
    return out.astype(int).values


def income_cells(df):
    """Household income on its native brackets (5 recorded + missing = 6)."""
    inc = pd.to_numeric(df["HH_INC"], errors="coerce")
    out = inc.fillna(-1).astype(int)
    return out.values


def design_matrix(df, columns, add_missing_dummies=True):
    """Numeric design matrix: median-impute + missingness indicators.

    MOSTIMPR (a categorical code 1-8) is expanded into dummies.
    """
    blocks, names = [], []
    for c in columns:
        v = pd.to_numeric(df[c], errors="coerce")
        if c == "MOSTIMPR":
            for lev in sorted(v.dropna().unique())[1:]:  # drop first level
                blocks.append((v == lev).astype(float).values)
                names.append(f"{c}_{int(lev)}")
        else:
            med = v.median()
            vv = v.fillna(med).values.astype(float)
            if v.nunique() > 2:  # winsorize continuous covariates at 1%/99%
                lo, hi = np.nanpercentile(vv, [1, 99])
                vv = np.clip(vv, lo, hi)
            blocks.append(vv)
            names.append(c)
        if add_missing_dummies and v.isna().any():
            blocks.append(v.isna().astype(float).values)
            names.append(f"{c}_miss")
    X = np.column_stack(blocks)
    # drop constant columns
    keep = X.std(axis=0) > 0
    return X[:, keep], [n for n, k in zip(names, keep) if k]

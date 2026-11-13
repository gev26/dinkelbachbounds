"""Additional empirical application: Baranov et al. (AER 2020) THP trial.

Bounds on treatment effects among ALWAYS-RESPONDERS (mothers re-interviewed
at follow-up under either assignment), without attrition monotonicity:

  (1) ITT on respondents (naive reference; cluster SE)
  (2) basic Lee bounds, monotone attrition (trims the larger-response arm,
      the direction the data selects: here control responds MORE, so
      monotonicity is S(1) <= S(0)) - validation anchor vs the paper's own
      Table 8 leebounds column
  (3) no-monotonicity Dinkelbach bounds, MARGINAL (no covariates): here the
      marginal floor s0+s1-1 = 0.297 > 0, so the bound exists without
      covariates - the regime JobCorps could not deliver
  (4) no-monotonicity Dinkelbach bounds with baseline covariates

Inference: cluster multiplier bootstrap (one Exp(1) weight per Union
Council, Remark 6.3's statistic) and cluster-summed IF variance; cross-fit
folds are drawn at the CLUSTER level so fold independence respects the
design.  95% Imbens-Manski regions with q = Phi^-1(sqrt(0.95)).
"""

import numpy as np
import pandas as pd

import data_prep_baranov as db
import estimator as est

SEED = 42
N_BOOT = 400


def cluster_folds(uc, D, n_folds=3, seed=SEED):
    """Assign whole clusters to folds, STRATIFIED by treatment arm.

    Treatment is cluster-constant (20 T / 20 C), so an unstratified
    permutation can produce badly imbalanced fold-level treatment shares
    (seed 42 unstratified: 7/4/9 treated clusters per fold, fold D-means
    0.31-0.72, cross-fitted D/p-hat averaging 2.3 instead of 1 - caught by
    the adversarial referee).  Stratifying preserves fold independence
    (assignment still ignores outcomes and selection) and pins every
    fold's treatment share near the design's 1/2.
    """
    rng = np.random.default_rng(seed)
    fold_of_cluster = {}
    for arm in (0, 1):
        ids = np.unique(uc[D == arm])
        perm = rng.permutation(len(ids))
        for j, c in enumerate(ids):
            fold_of_cluster[c] = perm[j] % n_folds
    return np.array([fold_of_cluster[c] for c in uc])


def _frac_tail_mean(y, w, share, top):
    """Weighted mean of the top/bottom `share` fraction, fractional at the
    boundary (Neyman-Pearson boundary weight, eq. (4.2) convention)."""
    o = np.argsort(-y if top else y, kind="mergesort")
    target = share * w.sum()
    tot = got = 0.0
    for k in o:
        take = min(target - got, w[k])
        tot += take * y[k]
        got += take
        if got >= target - 1e-12:
            break
    return tot / target


def lee_bidirectional(Y, S, D, w, uc, n_boot=N_BOOT, seed=11):
    """Lee (2009) bounds trimming the larger-response arm (leebounds
    convention), fractional boundary mass, cluster bootstrap for the SEs."""
    def point(Yb, Sb, Db, wb):
        i1 = (Db == 1) & (Sb == 1)
        i0 = (Db == 0) & (Sb == 1)
        s1 = np.average(Sb[Db == 1], weights=wb[Db == 1])
        s0 = np.average(Sb[Db == 0], weights=wb[Db == 0])
        m1 = np.average(Yb[i1], weights=wb[i1])
        m0 = np.average(Yb[i0], weights=wb[i0])
        if s1 <= s0:                          # trim control down to s1/s0
            p = s1 / s0
            lo = m1 - _frac_tail_mean(Yb[i0], wb[i0], p, top=True)
            hi = m1 - _frac_tail_mean(Yb[i0], wb[i0], p, top=False)
        else:                                 # trim treated down to s0/s1
            p = s0 / s1
            lo = _frac_tail_mean(Yb[i1], wb[i1], p, top=False) - m0
            hi = _frac_tail_mean(Yb[i1], wb[i1], p, top=True) - m0
        return lo, hi
    lo, hi = point(Y, S, D, w)
    rng = np.random.default_rng(seed)
    ids = np.unique(uc)
    los, his = [], []
    for _ in range(n_boot):
        pick = rng.choice(ids, len(ids), replace=True)
        idx = np.concatenate([np.flatnonzero(uc == c) for c in pick])
        try:
            l, h = point(Y[idx], S[idx], D[idx], w[idx])
            los.append(l)
            his.append(h)
        except Exception:
            continue
    q = est.norm.ppf(np.sqrt(0.95))
    seL, seU = float(np.std(los)), float(np.std(his))
    return dict(lower=lo, upper=hi, se_lower=seL, se_upper=seU,
                IM_lo=lo - q * seL, IM_hi=hi + q * seU, status="ok")


def itt_respondents(Y, S, D, w, uc):
    """Naive ITT on respondents with cluster-robust SE (reference row)."""
    i = S == 1
    m1 = np.average(Y[i & (D == 1)], weights=w[i & (D == 1)])
    m0 = np.average(Y[i & (D == 0)], weights=w[i & (D == 0)])
    beta = m1 - m0
    resid = np.where(D == 1, Y - m1, Y - m0)
    contrib = np.where(i, resid * (D / np.mean(D[i]) -
                                   (1 - D) / np.mean(1 - D[i])), 0.0)
    cl = np.array([contrib[uc == c].sum() for c in np.unique(uc[i])])
    G, n, k = len(cl), int(i.sum()), 2
    corr = (G / (G - 1)) * ((n - 1) / (n - k))   # CR1 small-sample factor
    se = float(np.sqrt(corr * np.sum(cl ** 2)) / i.sum())
    return beta, se


def fmt(res):
    if res.get("status") != "ok":
        return res.get("status", "failed")
    return (f"[{res['lower']:+.3f}, {res['upper']:+.3f}]  "
            f"CR=({res['IM_lo']:+.3f}, {res['IM_hi']:+.3f})")


def main():
    r, D, uc, S7 = db.load_randomized()
    w = np.ones(len(r))
    folds = cluster_folds(uc, D)
    # conditional column: the targeted attrition-relevant baseline set with
    # stronger quantile regularization -- with fold-level arm-selected
    # samples of ~95 obs, the full 16-covariate basis overfits the first
    # stages and inflates the estimated bounds beyond the marginal ones
    # (which Prop. 2.2's refinement ordering forbids in population)
    X = np.column_stack([pd.to_numeric(r[c], errors="coerce").values
                         for c in db.X_SMALL])
    Xnone = np.zeros((len(r), 1))             # marginal: single-cell design

    s1 = S7[D == 1].mean()
    s0 = S7[D == 0].mean()
    print(f"7-year response: s1={s1:.4f} (T), s0={s0:.4f} (C), "
          f"s0+s1={s0+s1:.4f} -> marginal Frechet box "
          f"[{max(s0+s1-1,0):.3f}, {min(s0,s1):.3f}]")
    print(f"clusters: {len(np.unique(uc))}, folds by cluster: "
          f"{np.bincount(folds).tolist()} obs per fold\n")

    rows = []
    specs = [(lab, col, hib, S7) for lab, col, hib in db.OUTCOMES_7Y]
    specs += [(lab, col, hib, db.wave_selection(r, col))
              for lab, col, hib in db.OUTCOMES_EARLY]

    for lab, col, hib, S in specs:
        Y = db.outcome(r, col, S)
        print("=" * 72)
        s1w, s0w = S[D == 1].mean(), S[D == 0].mean()
        print(f"{lab}  [{col}; higher = {hib}; s0+s1 = {s0w+s1w:.3f}]")

        beta, se = itt_respondents(Y, S, D, w, uc)
        print(f"  (1) ITT respondents  : {beta:+.3f} (cluster se {se:.3f})")

        r2 = lee_bidirectional(Y, S, D, w, uc)
        print(f"  (2) Lee (monotone)   : {fmt(r2)}")

        r3 = est.exact_marginal_bounds(Y, S, D, w, cluster=uc,
                                       n_boot=N_BOOT)
        print(f"  (3) no mon., marginal: {fmt(r3)}  [exact]")

        cacheX = est.first_stage_cache(X, Y, S, D, w, folds, qr_alpha=1e-2)
        r4 = est.bound_interval(X, Y, S, D, w, cacheX, cluster=uc,
                                n_boot=N_BOOT)
        extra = ""
        if r4.get("status") == "ok":
            rs = r4["regime_shares"]
            extra = (f"  regimes int/flr/cei = "
                     f"{rs['int']:.2f}/{rs['flr']:.2f}/{rs['cei']:.2f}")
        print(f"  (4) no mon., X       : {fmt(r4)}{extra}")

        def rec(spec, res):
            base = dict(outcome=col, label=lab, spec=spec,
                        s0=s0w, s1=s1w)
            if isinstance(res, dict) and res.get("status") == "ok":
                base.update(lower=res["lower"], upper=res["upper"],
                            IM_lo=res["IM_lo"], IM_hi=res["IM_hi"],
                            se_lower=res.get("se_boot_lower",
                                             res.get("se_lower")),
                            se_upper=res.get("se_boot_upper",
                                             res.get("se_upper")),
                            status="ok")
            elif isinstance(res, dict):
                base.update(status=res.get("status"))
            return base
        rows.append(dict(outcome=col, label=lab, spec="itt", s0=s0w, s1=s1w,
                         lower=beta, upper=beta, se_lower=se, se_upper=se,
                         IM_lo=beta - 1.96 * se, IM_hi=beta + 1.96 * se,
                         status="ok"))
        rows += [rec("lee_monotone", r2), rec("nomon_marginal", r3),
                 rec("nomon_X", r4)]

    from pathlib import Path
    path = Path(__file__).resolve().parent / "results"
    path.mkdir(exist_ok=True)
    
    out = pd.DataFrame(rows)
    out.to_csv(path / "baranov_bounds.csv", index=False)
    print("\nSaved: results/baranov_bounds.csv")


if __name__ == "__main__":
    main()

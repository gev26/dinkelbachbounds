"""Section 8 empirical application: JobCorps.

Produces:
  * Table 1 validation: E[pi_L(X)] under the logit specifications
    (baseline earnings + household income; Lee's 28) at weeks 45/90/100.
  * Table 5: bounds on the always-takers' wage effect by identification
    regime at weeks 45, 90, 110, 150, 208 --
      (1) basic Lee, no covariates (monotonicity);
      (2) conditional monotonicity, Lee's 28 covariates (generalized Lee /
          Appendix A: box collapsed to the ceiling);
      (3) no monotonicity, Reasons-for-joining covariates (Section 5
          estimator);
      (4) no monotonicity + Masten-Poirier pointwise share floor,
          pibar at the midpoint of [E pi_L, E pi_U] (see ASSUMPTIONS.md),
          with the full breakdown frontier saved separately.
  * Table 7: Table 5 divided by the one-sided compliance rate.

All computations use study design weights.  95% Imbens-Manski regions with
q = Phi^{-1}(sqrt(0.95)); endpoint SEs from the multiplier bootstrap
(Remark 6.2).
"""

import numpy as np
import pandas as pd

import data_prep as dp
import estimator as est

WEEKS_T5 = [45, 90, 110, 150, 208]
SEED = 42


def fmt(res):
    if res.get("status") != "ok":
        return res.get("status", "failed")
    return (f"[{res['lower']:+.3f}, {res['upper']:+.3f}]  "
            f"CR=({res['IM_lo']:+.3f}, {res['IM_hi']:+.3f})")


def main():
    df = dp.load_data()
    D = dp.treatment(df)
    w = dp.weights(df)
    N = len(df)
    rng = np.random.default_rng(SEED)
    folds = rng.permutation(np.arange(N) % 3)

    X_reasons, _ = dp.design_matrix(df, dp.REASONS_ITEMS)
    X_lee28, _ = dp.design_matrix(df, dp.LEE28)
    X_earninc, _ = dp.design_matrix(
        df, ["EARN_YR", "HH_INC1", "HH_INC2", "HH_INC3", "HH_INC4", "HH_INC5"])

    # ---- Table 1: share bounds, all columns, with the paper's one-sided
    # 95% lower confidence band for E[pi_L] from 400 nonparametric
    # bootstrap resamples (Section 3.1 / Table 1 notes) ----
    T1_BOOT = 400

    def cell_epi(S_, D_, w_, cells_):
        """Cell-averaged Frechet endpoints (2.19), design-weighted."""
        EL = EU = 0.0
        tot = w_.sum()
        for c in np.unique(cells_):
            m = cells_ == c
            m0, m1 = m & (D_ == 0), m & (D_ == 1)
            if m0.sum() == 0 or m1.sum() == 0:
                continue
            s0c = np.average(S_[m0], weights=w_[m0])
            s1c = np.average(S_[m1], weights=w_[m1])
            wc = w_[m].sum() / tot
            EL += wc * max(s0c + s1c - 1.0, 0.0)
            EU += wc * min(s0c, s1c)
        return EL, EU

    def logit_epi(S_, D_, w_, X_):
        s0 = est.Logit().fit(X_[D_ == 0], S_[D_ == 0],
                             w_[D_ == 0]).predict(X_)
        s1 = est.Logit().fit(X_[D_ == 1], S_[D_ == 1],
                             w_[D_ == 1]).predict(X_)
        return (float(np.average(np.maximum(s0 + s1 - 1, 0), weights=w_)),
                float(np.average(np.minimum(s0, s1), weights=w_)))

    def cat_cells(col):
        v = pd.to_numeric(df[col], errors="coerce")
        return v.fillna(-1).astype(int).values

    discrete_specs = {"None": np.zeros(N, dtype=int),
                      "Hispanic": cat_cells("HISP"),
                      "Kidcount": cat_cells("NCHLD"),
                      "Married": cat_cells("MARRIED")}
    logit_specs = {"EarnInc": X_earninc, "Lee28": X_lee28,
                   "Reasons": X_reasons}

    print("=" * 70)
    print("Table 1: share bounds with one-sided 95% lower bands "
          f"({T1_BOOT} bootstrap resamples)")
    t1 = []
    rng_t1 = np.random.default_rng(2027)
    for week in (45, 90, 100):
        S, _ = dp.outcome_selection(df, week)
        row = {"week": week}
        boot_idx = [rng_t1.integers(0, N, N) for _ in range(T1_BOOT)]
        for name, cells in discrete_specs.items():
            row[f"EpiL_{name}"], row[f"EpiU_{name}"] = cell_epi(S, D, w, cells)
            draws = [cell_epi(S[i], D[i], w[i], cells[i])[0]
                     for i in boot_idx]
            row[f"band_{name}"] = float(np.quantile(draws, 0.05))
        for name, X in logit_specs.items():
            row[f"EpiL_{name}"], row[f"EpiU_{name}"] = logit_epi(S, D, w, X)
            draws = [logit_epi(S[i], D[i], w[i], X[i])[0]
                     for i in boot_idx]
            row[f"band_{name}"] = float(np.quantile(draws, 0.05))
        t1.append(row)
        print({k: (round(v, 4) if isinstance(v, float) else v)
               for k, v in row.items()})
    pd.DataFrame(t1).to_csv("../results/table1_validation.csv", index=False)

    # ---- compliance rate for Table 7 ----
    enroll = pd.to_numeric(df["ENROLL"], errors="coerce").values
    ok = ~np.isnan(enroll)
    compl_w = np.average(enroll[ok & (D == 1)], weights=w[ok & (D == 1)])
    print(f"\ncompliance Pr(ENROLL=1|Z=1): weighted={compl_w:.4f} "
          f"(paper: 0.725)")

    # ---- Table 5 ----
    rows5, frontier_rows = [], []
    for week in WEEKS_T5:
        S, Y = dp.outcome_selection(df, week)
        print("\n" + "=" * 70)
        print(f"WEEK {week}")

        r1 = est.basic_lee(Y, S, D, w)
        print("  (1) basic Lee        :", fmt({**r1, "status": "ok"}))

        cache28 = est.first_stage_cache(X_lee28, Y, S, D, w, folds)
        r2 = est.bound_interval(X_lee28, Y, S, D, w, cache28,
                                collapse_ceiling=True)
        print("  (2) cond. mon. Lee28 :", fmt(r2))

        cacheR = est.first_stage_cache(X_reasons, Y, S, D, w, folds)
        r3 = est.bound_interval(X_reasons, Y, S, D, w, cacheR)
        print("  (3) no mon. Reasons  :", fmt(r3),
              "regimes:", r3.get("regime_shares"))

        # Masten-Poirier pointwise floor: frontier and midpoint choice
        s0f = est.Logit().fit(X_reasons[D == 0], S[D == 0],
                              w[D == 0]).predict(X_reasons)
        s1f = est.Logit().fit(X_reasons[D == 1], S[D == 1],
                              w[D == 1]).predict(X_reasons)
        EpiL = np.average(np.maximum(s0f + s1f - 1, 0), weights=w)
        EpiU = np.average(np.minimum(s0f, s1f), weights=w)
        grid = np.linspace(EpiL, EpiU, 7)
        pibar_mid = 0.5 * (EpiL + EpiU)
        r4 = None
        for pibar in list(grid) + [pibar_mid]:
            rr = est.bound_interval(X_reasons, Y, S, D, w, cacheR,
                                    floor_bar=pibar)
            frontier_rows.append(dict(week=week, pibar=pibar,
                                      lower=rr.get("lower"),
                                      upper=rr.get("upper"),
                                      status=rr.get("status")))
            if pibar == pibar_mid:
                r4 = rr
        print(f"  (4) no mon. + MP (pibar={pibar_mid:.3f}):", fmt(r4))

        def rec(col, r):
            if r.get("status", "ok") != "ok":
                return dict(week=week, col=col, status=r.get("status"))
            return dict(week=week, col=col, status="ok",
                        lower=r["lower"], upper=r["upper"],
                        IM_lo=r["IM_lo"], IM_hi=r["IM_hi"],
                        se_lower=r.get("se_boot_lower", r.get("se_lower")),
                        se_upper=r.get("se_boot_upper", r.get("se_upper")))
        rows5 += [rec("basic_lee", {**r1, "status": "ok"}),
                  rec("cond_mon_lee28", r2),
                  rec("no_mon_reasons", r3),
                  rec("no_mon_MP", r4)]

    t5 = pd.DataFrame(rows5)
    t5.to_csv("../results/table5_bounds.csv", index=False)
    pd.DataFrame(frontier_rows).to_csv("../results/mp_frontier.csv",
                                       index=False)

    # ---- Table 7: complier rescaling ----
    # The paper's Table 7 notes fix Pr(D=1|Z=1) = 0.725; the data-derived
    # weighted rate (printed above) is reported alongside for transparency.
    compl_paper = 0.725
    t7 = t5.copy()
    for c in ("lower", "upper", "IM_lo", "IM_hi", "se_lower", "se_upper"):
        if c in t7:
            t7[c] = t7[c] / compl_paper
    t7["compliance_used"] = compl_paper
    t7["compliance_data"] = compl_w
    t7.to_csv("../results/table7_compliers.csv", index=False)
    print("\nSaved: table1_validation.csv, table5_bounds.csv, "
          "mp_frontier.csv, table7_compliers.csv")


if __name__ == "__main__":
    main()

# Implementation choices and assumptions

Everything below is a decision the paper leaves open (or a data-driven
necessity), stated explicitly so nothing silently deviates from the
methodology. Items marked **[validated]** are corroborated by reproducing a
number reported in the paper.

## Data mapping

1. **Outcome** `Y_w = log(EARNH_w / HWH_w)` (weekly earnings over weekly
   hours = log hourly wage) and **selection** `S_w = 1{EARNH_w > 0 & HWH_w > 0}`.
   **[validated]** Pooled selection rates s0+s1 = 0.77/0.92/0.96 at weeks
   45/90/100 (Section 3.1) and the weighted arm-difference in mean log wage
   (+0.027 at week 45, +0.048 at week 90) reproduce Table 5, column (1).
2. **Treatment** = `TREATMNT.y` (random assignment), **weights** = `DSGN_WGT.x`.
3. **Compliance rate** for the Table 7 rescaling is computed from the data
   (`ENROLL` among Z=1, design-weighted). The paper states 0.725; the raw
   unweighted rate is 0.751. The runner reports the computed value next to
   the paper's.

## Covariate sets

4. **"Reasons for joining" (8 items)** = the seven binary items `R_HOME,
   R_COMM, R_TRAIN, R_CRGOAL, R_GETGED, R_NOWORK, R_OTHER` plus the
   most-important-reason code `MOSTIMPR` (8 categories, entered as dummies).
   The paper does not list the item names. **[validated]** E[pi_L] ≈ 0.05 at
   week 90 vs. the paper's "about 0.06".
5. **"Lee's 28"**: the 28 baseline covariates are not listed in the paper;
   we use the standard Lee (2009)/Semenova (2024) set (demographics, marital
   status, children, education incl. parents', work history, baseline
   earnings, household & personal income brackets) — see `LEE28` in
   `data_prep.py`. **[validated]** E[pi_L] compared against Table 1, col (6).
6. **Missing covariates**: median imputation + missingness indicator dummies
   (the paper is silent on missing-data handling; 36% of `HH_INC` and 27% of
   `R_GETGED` are missing, so dropping rows is not an option).
7. **LP cells**: "annual baseline earnings quantile-binned into five cells" —
   the ~38% point mass at EARN_YR = 0 collapses five quantile edges into
   four distinct cells; missing EARN_YR forms its own cell. "Household income
   on its six native brackets" — HH_INC records five brackets; the sixth
   cell is the missing-income category.

## First-stage estimation

8. **Penalty levels** (the paper specifies l1-penalized fits but not the
   tuning constants): logistic C = 1.0 (liblinear), quantile regression
   alpha = 1e-3, hinge Lasso alpha = 1e-3, all on standardized covariates.
9. **Quantile model**: linear-in-X quantile regression on the u-grid
   {0.005, 0.02, ..., 0.98, 0.995}, monotonically rearranged (row-wise sort,
   as prescribed in §5.2), linearly interpolated between grid points, flat
   extrapolation beyond the grid ends, and clipped to the observed outcome
   range of the fitting arm.
10. **Trim integrals** psi(t,x) are computed as s_d * int_0^{t/s_d} F_d^{-1}(u) du
    by trapezoid on the fitted (rearranged) quantile curve.
11. **Overlap clipping**: fitted selection probabilities s0, s1 clipped to
    [0.01, 0.99] and the fitted propensity p to [0.05, 0.95]
    (strong-overlap Assumption 6.1 made operational). Continuous covariates
    are winsorized at their 1st/99th percentiles inside the basis b(X):
    without this, heavy-tailed baseline earnings give a few observations
    extreme leverage in the propensity logit (fitted p up to 0.999) and the
    control-arm AIPW factor 1/(1-p) explodes.

## Estimator structure

12. **Type signal for the difference parameter**: psi_t is the treated-arm
    signal (5.6) with threshold tau_1, exactly as printed for the level case.
    The paper does not print a difference-specific psi_t; the treated-arm
    signal remains conditionally unbiased for pi* in all regimes.
13. **Combined moment for the difference parameter** (the paper prints
    (5.16) only for the level target):
    Psi = eta1 - eta0 - lambda*psi_pi + tau1*(psi_pi - psi_t1)
          - tau0*(psi_pi - psi_t0),
    with one type signal per revealing arm (psi_t0 uses the control arm and
    the upper-tail indicator 1{Y >= tau0}).  This is the unique assembly
    that (a) reduces verbatim to (5.16) when arm 0 is absent, (b) has
    conditional mean equal to the per-cell Dinkelbach value at the truth,
    and (c) is Neyman-orthogonal in both threshold channels: d/dtau0 = 0
    identically by Rockafellar-Uryasev stationarity and d/dtau1 =
    (tau1 - tau0 - lambda) f1(tau1), zero on the interior by the QTE FOC
    (5.22) and cancelled on binding faces by the correction terms.  (An
    earlier draft of this code used the prefactor (tau1 - tau0) on a single
    correction term, which fails threshold orthogonality; caught by the
    adversarial referee pass and verified analytically.)
13b. **Masten-Poirier active floor face**: with the pointwise floor
    pi_L^new(x) = min(max(pi_L(x), pibar), pi_U(x)), a floor-binding cell's
    active facet is the Frechet floor (signal psi_LB) where pi_L(x) >= pibar,
    the *known constant* pibar (influence = the constant itself) where
    pi_L(x) < pibar <= pi_U(x), and the ceiling (psi_UB) where
    pibar > pi_U(x).  The paper marks MP inference conjectural (Sec. 7
    Remark 7.1(c)); point bounds are unaffected.
14. **Upper endpoint** computed by running the identical lower-endpoint
    machinery on -Y and negating (the remark after (5.9)); the fitted
    quantile curves are mirrored, F_{-Y}^{-1}(u) = -F_Y^{-1}(1-u), rather
    than refit.
15. **lambda_init** uses sign-bisection (the variant Proposition 5.1 covers)
    rather than the fixed-point iteration; both converge to the same root.
16. **Basic Lee (col. 1)**: trimming proportion 1 - s0/s1 clamped at zero
    when s1 < s0 (no arm flipping). **[validated]** reproduces the week-45
    point interval [+0.027, +0.027].
17. **Masten-Poirier floor (col. 4)**: enforced pointwise,
    pi_L^new(x) = min(max(pi_L(x), pibar), pi_U(x)), per Section 8 point 3.
    The paper reports "the regret-chosen / median pibar" without a formula;
    we report the bound at the **midpoint of [E pi_L, E pi_U]** and save the
    full frontier over a 7-point grid in `results/mp_frontier.csv` so any
    other pibar convention can be read off.
18. **Folds and aggregation**: three equal folds by seeded permutation
    (seed 42), three cyclic rotations per Definition 5.1.  The point
    estimate is the **pooled cross-fit ratio** of (5.17)'s 1/N sums (each
    observation enters exactly one J3, so the pooled averages run over the
    full sample); per-rotation ratios are retained as diagnostics.  This is
    the reading of the (5.17)-vs-Definition-5.1 tension under which the
    point estimate coincides with the statistic the influence-function
    variance and the Remark 6.3 bootstrap describe.  A rotation whose
    debiased J3 denominator is non-positive (possible where the stratum
    share is tiny, e.g. week 45) makes the cell **degenerate**: the run is
    surfaced with status "degenerate-denominator" instead of averaging a
    sign-flipped ratio (Theorem 6.1 divides by piAT > 0; Algorithm 5.1
    gates the degenerate cases).  A no-sign-change bisection reports the
    degenerate floor value 0 directly per Algorithm 5.1 step (0), without
    running the one-step at an arbitrary lambda.

## Inference

19. Endpoint SEs: both the plug-in IF variance (Thm 6.1(iii)) and the
    **Remark 6.3 multiplier bootstrap** are computed.  The bootstrap holds
    the cross-fitted Psi_i and piAT fixed and draws the *linear* statistic
    (weighted mean of e_i Psi_i)/piAT with e_i ~ Exp(1) -- the denominator
    is never re-randomized, so the draws have all moments (an earlier
    ratio-form bootstrap with a re-randomized denominator, patched by a
    MAD-based scale, was removed after the conformance audit: neither
    appears in the paper).  Draws are shared across endpoints so the joint
    covariance Sigma_LU of Remark 6.2 is estimable.  The IM region (6.3)
    uses the bootstrap SDs with q = Phi^{-1}(sqrt(0.95)); the Remark 6.3
    quantile critical values are stored alongside as diagnostics.
19b. **Type-signal centering in the grid-cap windows**: the IPW residuals
    of psi_t1/psi_t0 are centered at the kept mass the fitted (grid-capped,
    saturated) quantile surrogate actually delivers, s_d * u_eff with
    u_eff = clip(a, 0.005, 0.995) (or 1/0 under saturation), instead of at
    pi_star.  On interior cells this equals pi_star; it differs only in the
    narrow cap windows on binding cells, where centering at pi_star would
    leak a first-order propensity term.
19c. **Table 1 inference**: the one-sided 95% lower confidence band for
    E[pi_L(X)] from 400 nonparametric bootstrap resamples (Section 3.1 /
    Table 1 notes), for the discrete-cell columns (None/Hispanic/Kidcount/
    Married, cell-averaged endpoints per (2.19)) and the logit columns.
20. Design weights enter all moment averages and the variance as a weighted
    sandwich; the weighted ratio SE treats weights as fixed.

---

## Items 21-34: added after the four-referee code audit (see
## AUDIT_CODE_ASSUMPTIONS.md in the project root for the full findings)

### Corrections to items above
- Item 8 is incomplete: qr_alpha = 1e-3 holds only for JobCorps; the
  Baranov, Bursztyn and deworming runs use qr_alpha = 1e-2 (small
  fold-arm samples).  The OLS sensitivity check (olsversions/) shows no
  signed conclusion depends on the penalty.
- Item 11's blanket winsorization/clipping claim is inaccurate in two
  places: Baranov's X_SMALL covariates enter raw (unwinsorized; verified
  complete, no NaN), and the MP-frontier / Table-1 logit share fits are
  full-sample, non-cross-fitted and UNclipped (a different object from
  the clipped cross-fitted shares in the estimator's box).
- Item 13b describes the MP face split as pi_L(x) < pibar <= pi_U(x);
  the code assigns the constant-pibar face only when pibar < pi_U(x)
  STRICTLY (at equality the cell gets the ceiling face).
- Item 16 documents basic Lee with no arm flipping; the applications'
  lee_bidirectional (Baranov/Bursztyn/deworming) instead picks the
  trimming direction from the data (trims the higher-response arm),
  re-choosing it inside every bootstrap resample.  Two Lee conventions
  coexist deliberately.

### Core estimator numerical conventions (not in the paper)
21. dinkelbach_root accepts a bracket ENDPOINT as the root whenever
    |G(endpoint)| <= tol * max(mean pi*, 1e-12) before declaring the
    degenerate 0 that Algorithm 5.1 step (0) prescribes; the relative
    tolerance denominator is floored at 1e-12; 200-iteration bisection cap.
22. The QTE FOC (5.22) is inverted on a fixed 101-point t-grid per cell
    with cummax monotonization and linear interpolation (O((piU-piL)/100)
    discretization); on flat QTE segments the RIGHT edge (largest root)
    is taken.
23. Saturation cutoffs a1 >= 0.9995 / a0 <= 0.0005 (relative to grid ends
    0.995/0.005) decide whole-arm-kept saturation in both the signal
    indicators and the Lasso hinge targets.
24. Tie-breaks: ceiling binding arm d_dag = 1{s1 < s0} picks arm 0 on
    exact ties; classify() breaks equidistant floor/ceiling toward floor.
25. The difference-parameter bisection bracket is
    [min_i(Q1(0.005)-Q0(0.995)), max_i(Q1(0.995)-Q0(0.005))] (capped
    grid ends standing in for F^{-1}(0)/F^{-1}(1)).
26. Logit falls back to the weighted mean on single-class outcomes or
    constant designs (avoids liblinear's penalized-intercept shrinkage).

### Cluster inference (paper is iid-only; all of this is our extension)
27. Cluster designs (Baranov 40 UCs, deworming 75 schools, Bursztyn
    _sessclust): one Exp(1) multiplier per cluster, cluster-summed IF
    variance with G/(G-1), bootstrap draws rescaled by sqrt(G/(G-1)),
    and cross-fit folds assigned at cluster level stratified by arm.
    None of this appears in the paper (Remark 6.3 is per-observation).
    buffer delta_N and tolerances use N observations, not G clusters.
28. exact_marginal_bounds (marginal columns, all applications) bypasses
    the cross-fitted one-step entirely: exact atom-boundary enumeration
    with fractional trimming, NONPARAMETRIC (cluster) bootstrap -- not
    the Remark 6.3 multiplier scheme (runner docstrings saying otherwise
    are wrong for that row).  Bootstrap resamples whose Frechet floor is
    <= 0 are silently dropped: the SE is conditional on non-degeneracy.

### Application-level conventions
29. Baranov: (a) X_SMALL (6 covariates) and qr_alpha=1e-2 were selected
    PARTLY AFTER observing that the full 16-covariate basis produced
    conditional bounds wider than the marginal ones (violating Prop 2.2's
    refinement ordering) -- a data-informed specification choice that
    should be disclosed wherever column (4) is reported; (b) the pooled
    motherfinancial index blanks 235 genuinely observed values held by
    7-year attriters (selection is attrit2-based, not observability-
    based); (c) 6m/1y wave selection = depindex item non-missingness,
    conflating unit attrition with item nonresponse; (d) depression
    outcomes are reported on the raw index ('higher = worse' is a label
    only; welfare direction is reversed vs other rows).
30. Bursztyn: (a) two respondent definitions coexist -- Exercise A uses
    applied_out_fl non-missingness (N=290, incl. 6 matched==0 subjects
    Kwon-Roth exclude), Exercise B uses matched==1 (N=284); outcome-
    specific Sa quietly reverts outside_others_fl to the 284 sample;
    (b) binary outcomes run through machinery whose CLT assumes a
    continuous outcome density at the trim threshold (Assumption
    6.3(v)/Remark 6.4) -- exact fractional-atom handling mitigates but
    the inference approximation is unproven for 0/1 outcomes;
    (c) covariate imputation/winsorization is full-sample (mild leak vs
    Def 5.1's fold discipline) and median imputation has NO missingness
    dummies (contradicting item 6's convention); (d) the 23 pseudo-
    session ids > 1000 (all non-respondents) are merged into ONE cluster;
    (e) Exercise B additionally conditions on matched==1 (response
    ignorability stated in the docstring); estimand alignment with
    Kwon-Roth's never-signer bound verified sign-correct.
31. JobCorps: (a) THE headline divergence -- the paper's printed Table 5
    columns (2)-(4) are NOT reproduced (week 90 col (3): ours
    [-0.92,+1.12] vs printed [-0.40,+0.26]; col (2) [+0.015,+0.060] vs
    [-0.067,+0.002]; col (4) [-0.41,+0.51] vs [-0.245,+0.218]); only
    column (1) and Table 1 match.  COMPARISON.md's conclusion that the
    printed column (3) is not reproducible belongs here on the record.
    Week 45 col (3) is gated degenerate here vs finite in print.
    (b) The Table-3 LP is a two-arm common-share difference program --
    our extension; the paper prints only the one-arm (2.16)-(2.17), and
    its own Table 3 appears to be an arm-wise outer/Minkowski box.
    (c) Wage bins: pooled-sample equal-mass edges, clamped/de-duplicated.
    (d) Table 1 discrete cells give missing HISP/MARRIED their own cell.
    (e) cell_epi skips cells empty in one arm (contributes 0 to both
    endpoints while keeping the weight in the denominator).
    (f) Table 1 band = 5th percentile of 400 draws (percentile method),
    shared resample indices (seed 2027) across columns.
    (g) MP frontier endpoints/midpoint from a full-sample unclipped
    logit (see item 11 correction).  (h) 400 draws and seeds 11/42/2027
    are code choices; the paper fixes 400 only for Table 1.
32. LP granularity table computed unweighted per Table 3 notes, silently
    resolving the paper's internal contradiction with Section 2.3's
    'study-weighted' text.
33. Column (4) (MP) confidence regions are exported unqualified although
    the paper marks the MP influence function conjectural.
34. ITT reference rows use CR1 ((G/(G-1))((n-1)/(n-k)), k=2) with z=1.96
    stored in the IM_lo/IM_hi columns; bound rows use q=1.9545 -- two
    critical-value conventions share one column name in the CSVs.

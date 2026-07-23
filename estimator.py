"""The debiased Dinkelbach estimator of Section 5 (difference parameter).

Implements, for the always-takers' wage effect beta = E[Y(1)-Y(0) | AT]:

  * first-stage nuisances (Section 5.1): l1-penalised logistic for the
    propensity p(x) and arm-specific selection probabilities s_d(x);
    l1-penalised quantile regression for the conditional outcome quantiles
    F_d^{-1}(u, x), monotonically rearranged (Section 5.2);
  * the preliminary Dinkelbach root lambda_init (Algorithm 5.1), adapted to
    the difference parameter through the QTE first-order condition (5.22):
    per cell, QTE(t,x) = F_1^{-1}(t/s1,x) - F_0^{-1}(1-t/s0,x) = lambda,
    solved by a monotone root-find and clipped to the Frechet box;
  * the localized one-step (Algorithm 5.2 + Section 5.4.3): freeze
    (pi*, tau_1, tau_0) at lambda_init, classify cells by the buffer rule
    (5.20), fit the hinge regressions, and form the debiased ratio (5.17)
    with the binding-face correction, here with prefactor
    (tau_1 - tau_0 - lambda) -- the difference-parameter analogue of
    (tau - lambda) in (5.16), vanishing on interior cells by the QTE FOC;
  * three-way cross-fitting (Definition 5.1) with cyclic rotation;
  * the upper endpoint by running the identical machinery on -Y and negating
    (the upper trim is the lower trim of -Y, cf. the remark after (5.9));
  * endpoint sorting (Algorithm 5.2, step 4);
  * inference: plug-in influence-function variance (Theorem 6.1(iii)),
    multiplier-bootstrap endpoint covariance (Remark 6.2), and the
    Imbens-Manski confidence region (6.3) with q = Phi^{-1}(sqrt(1-alpha));
  * the Masten-Poirier pointwise share floor of Section 8, point 3:
    pi_L^new(x) = min{max(pi_L(x), pibar), pi_U(x)};
  * the conditional-monotonicity (generalized Lee) special case of
    Appendix A: box collapsed to pi*(x) = pi_U(x), ceiling signal everywhere.

Implementation choices the paper leaves open are flagged in ASSUMPTIONS.md.
"""

import numpy as np
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression, Lasso, QuantileRegressor


# ----------------------------------------------------------------------
# first stages
# ----------------------------------------------------------------------

U_GRID = np.concatenate(([0.005], np.linspace(0.02, 0.98, 49), [0.995]))


def _standardize(X):
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd, mu, sd


class Logit:
    """l1-penalised logistic regression (5.3), weighted."""

    def __init__(self, C=1.0):
        self.C = C

    def fit(self, X, y, w=None):
        self.Xs_, self.mu_, self.sd_ = None, None, None
        Xs, self.mu_, self.sd_ = _standardize(X)
        y = np.asarray(y)
        if len(np.unique(y)) < 2:
            self.const_ = float(np.average(y, weights=w))
            self.m_ = None
            return self
        self.const_ = None
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.m_ = LogisticRegression(penalty="l1", C=self.C,
                                         solver="liblinear", max_iter=2000)
            self.m_.fit(Xs, y, sample_weight=w)
        return self

    def predict(self, X):
        if self.m_ is None:
            return np.full(len(X), self.const_)
        Xs = (X - self.mu_) / self.sd_
        return self.m_.predict_proba(Xs)[:, 1]


class QuantileCurves:
    """l1-penalised quantile regression (5.4) on a u-grid, rearranged.

    predict(Xnew) returns an (n, U) array Q[i, j] = F^{-1}(U_GRID[j], x_i),
    monotone in j by rearrangement (row-wise sort), clipped to the observed
    outcome range of the fitting arm.
    """

    def __init__(self, alpha=1e-3):
        self.alpha = alpha

    def fit(self, X, y, w=None):
        Xs, self.mu_, self.sd_ = _standardize(X)
        self.ymin_, self.ymax_ = y.min(), y.max()
        self.models_ = []
        for u in U_GRID:
            m = QuantileRegressor(quantile=u, alpha=self.alpha,
                                  solver="highs")
            m.fit(Xs, y, sample_weight=w)
            self.models_.append(m)
        return self

    def predict(self, X):
        Xs = (X - self.mu_) / self.sd_
        Q = np.column_stack([m.predict(Xs) for m in self.models_])
        Q = np.sort(Q, axis=1)                      # monotone rearrangement
        return np.clip(Q, self.ymin_, self.ymax_)


def quantile_at(Q, u):
    """Row-wise interpolation of the fitted quantile curve at levels u.

    Q: (n, U) rearranged curves on U_GRID; u: scalar or (n,) array in [0,1].
    Flat extrapolation beyond the grid ends.
    """
    u = np.broadcast_to(np.asarray(u, dtype=float), (Q.shape[0],))
    uu = np.clip(u, U_GRID[0], U_GRID[-1])
    j = np.clip(np.searchsorted(U_GRID, uu) - 1, 0, len(U_GRID) - 2)
    u0, u1 = U_GRID[j], U_GRID[j + 1]
    t = (uu - u0) / (u1 - u0)
    rows = np.arange(Q.shape[0])
    return (1 - t) * Q[rows, j] + t * Q[rows, j + 1]


def lower_partial_int(Q, a):
    """int_0^a F^{-1}(u) du, row-wise, trapezoid on the u-grid.

    Constant extrapolation of F^{-1} below U_GRID[0] and above U_GRID[-1].
    a: (n,) array in [0, 1].
    """
    n, U = Q.shape
    a = np.clip(np.asarray(a, dtype=float), 0.0, 1.0)
    # cumulative trapezoid over the grid, with flat left tail [0, U_GRID[0]]
    left = Q[:, 0] * U_GRID[0]
    du = np.diff(U_GRID)
    mid = np.cumsum(0.5 * (Q[:, 1:] + Q[:, :-1]) * du, axis=1)
    cum = np.column_stack([left, left[:, None] + mid])  # int_0^{U_GRID[j]}
    j = np.clip(np.searchsorted(U_GRID, a) - 1, 0, U - 2)
    rows = np.arange(n)
    below = a <= U_GRID[0]
    u0 = U_GRID[j]
    q0, q1 = Q[rows, j], Q[rows, j + 1]
    frac = np.where(a > u0, (a - u0) / (U_GRID[j + 1] - u0), 0.0)
    qa = (1 - frac) * q0 + frac * q1
    extra = 0.5 * (q0 + qa) * (a - u0)
    out = cum[rows, j] + np.where(a > u0, extra, 0.0)
    out = np.where(below, Q[:, 0] * a, out)
    over = a >= U_GRID[-1]
    out = np.where(over, cum[:, -1] + Q[:, -1] * (a - U_GRID[-1]), out)
    return out


# ----------------------------------------------------------------------
# per-cell program: QTE inversion and trims (Section 5.4.3)
# ----------------------------------------------------------------------

class CellProgram:
    """Vectorised per-cell objects at held-out points.

    Given fitted quantile curves Q1, Q0 (n, U) and selection rates s1, s0,
    provides, for the *lower* endpoint of the difference target:
      QTE(t)   = F1^{-1}(t/s1) - F0^{-1}(1 - t/s0)   (5.22), increasing in t
      n(t)     = psi^(1)(t) - psibar^(0)(t)          (5.21)
      pi*(lam) = clip_[piL, piU]( t : QTE(t) = lam )
    """

    def __init__(self, Q1, Q0, s1, s0, piL, piU, n_t=101):
        self.Q1, self.Q0 = Q1, Q0
        self.s1, self.s0 = s1, s0
        self.piL, self.piU = piL, piU
        n = len(s1)
        grid01 = np.linspace(0.0, 1.0, n_t)
        # per-point t-grid spanning the Frechet box
        self.T = piL[:, None] + (piU - piL)[:, None] * grid01[None, :]
        Tc = np.clip(self.T, 1e-12, None)
        self.QTE = np.column_stack([
            quantile_at(Q1, Tc[:, k] / s1) - quantile_at(Q0, 1.0 - Tc[:, k] / s0)
            for k in range(n_t)])
        self.QTE = np.maximum.accumulate(self.QTE, axis=1)  # enforce monotone

    def pi_star(self, lam):
        """clip to box of the root of QTE(t) = lam (monotone interpolation)."""
        n, K = self.QTE.shape
        out = np.empty(n)
        for i in range(n):
            out[i] = np.interp(lam, self.QTE[i], self.T[i])
        return out

    def trim_value(self, t):
        """n(t) = psi^(1)(t) - psibar^(0)(t) per point."""
        t = np.asarray(t)
        a1 = np.clip(t / self.s1, 0, 1)
        psi1 = self.s1 * lower_partial_int(self.Q1, a1)
        # upper trim of arm 0 at mass t: s0 * int_{1 - t/s0}^1 F0^{-1}
        a0 = np.clip(1.0 - t / self.s0, 0, 1)
        full0 = self.s0 * lower_partial_int(self.Q0, np.ones(len(t)))
        psibar0 = full0 - self.s0 * lower_partial_int(self.Q0, a0)
        return psi1 - psibar0

    def thresholds(self, t):
        """Trim thresholds and kept fractions.

        a1 = t/s1 is the kept (bottom) fraction of arm 1; a0 = 1 - t/s0 is
        the discarded (bottom) fraction of arm 0.  The fractions are
        returned so that downstream indicators can saturate exactly when a
        trim keeps the whole arm (a1 = 1, or a0 = 0), where the fitted
        end-of-grid quantile would otherwise wrongly exclude the tail.
        """
        a1 = np.clip(t / self.s1, 0, 1)
        a0 = np.clip(1.0 - t / self.s0, 0, 1)
        tau1 = quantile_at(self.Q1, a1)
        tau0 = quantile_at(self.Q0, a0)
        return tau1, tau0, a1, a0


def dinkelbach_root(cp, wts, N_total=None, tol=None):
    """Algorithm 5.1 (sign-bisection variant) for the plug-in moment (5.5),

    G(lam) = wmean[ n(pi*(lam)) - lam * pi*(lam) ],
    on the natural bracket spanned by the full quantile support.  Returns
    (lambda_init, status); on no sign change the paper prescribes the
    degenerate floor value 0 (Algorithm 5.1, step 0).
    """
    N = N_total if N_total is not None else len(wts)
    tol = tol if tol is not None else N ** (-0.25) / np.log(N)
    if np.average(cp.piL, weights=wts) <= 0:
        return np.nan, "undefined"          # 0/0 gate

    # natural bracket (Algorithm 5.1, step 0): the full quantile support,
    # [min_i (F1^{-1}(0) - F0^{-1}(1)), max_i (F1^{-1}(1) - F0^{-1}(0))],
    # so every cell-wise QTE transition lies inside it even when the floor
    # is raised (Masten-Poirier).
    lam_lo = (cp.Q1[:, 0] - cp.Q0[:, -1]).min()
    lam_hi = (cp.Q1[:, -1] - cp.Q0[:, 0]).max()

    def G(lam):
        ps = cp.pi_star(lam)
        return np.average(cp.trim_value(ps) - lam * ps, weights=wts), ps

    g_lo, _ = G(lam_lo)
    g_hi, _ = G(lam_hi)
    if not (g_lo > 0 > g_hi):
        return 0.0, "degenerate"            # degenerate floor (Alg. 5.1)
    for _ in range(200):
        lam = 0.5 * (lam_lo + lam_hi)
        g, ps = G(lam)
        denom = np.average(ps, weights=wts)
        if abs(g) <= tol * max(denom, 1e-12) or lam_hi - lam_lo < 1e-10:
            return lam, "ok"
        if g > 0:
            lam_lo = lam
        else:
            lam_hi = lam
    return lam, "ok"


# ----------------------------------------------------------------------
# orthogonal signals (Section 5.3)
# ----------------------------------------------------------------------

def buffer_delta(N):
    return N ** (-0.25) / np.log(N)         # (5.19)


def classify(pi_star, piL, piU, dN):
    """Buffer rule (5.20): 0 = interior, 1 = floor, 2 = ceiling."""
    reg = np.zeros(len(pi_star), dtype=int)
    flr = pi_star <= piL + dN
    cei = pi_star >= piU - dN
    both = flr & cei
    d_flr = np.abs(pi_star - piL)
    d_cei = np.abs(piU - pi_star)
    reg[flr] = 1
    reg[cei] = 2
    reg[both] = np.where(d_flr[both] <= d_cei[both], 1, 2)
    return reg


def signals(Y, S, D, w, p, s0, s1, piL, piU, pi_star, tau1, tau0,
            g1, g0, regime, floor_face=None, a1=None, a0=None):
    """The orthogonal signals of Section 5.3, difference parameter.

    The combined moment generalizing (5.16) to the difference target is

        Psi = eta1 - eta0 - lambda*psi_pi
              + tau1*(psi_pi - psi_t1) - tau0*(psi_pi - psi_t0),

    with a type signal per revealing arm (the per-configuration signals of
    Section 4.4, signed by the contrast):
        psi_t1 = pi* + D/p       * (S*1{Y <= tau1} - pi*)
        psi_t0 = pi* + (1-D)/(1-p) * (S*1{Y >= tau0} - pi*).
    Channel-by-channel: d/dtau1 E[Psi|X] = (tau1-tau0-lambda) f1(tau1) = 0 on
    the interior by the QTE FOC (5.22) and cancels exactly on binding faces;
    d/dtau0 E[Psi|X] = 0 identically (Rockafellar-Uryasev stationarity); with
    arm 0 absent the moment reduces verbatim to (5.16).  Solving the sample
    moment in lambda gives the one-step ratio (5.17) with numerator summand
        A_i = eta1 - eta0 + tau1*(psi_pi - psi_t1) - tau0*(psi_pi - psi_t0).

    floor_face (optional, Masten-Poirier pointwise floor): per-observation
    code for which face the *raised* floor pi_L^new = min(max(pi_L, pibar),
    pi_U) presents on floor-binding cells: 0 = Frechet floor (psi_LB),
    1 = the known constant pibar (influence function = the constant itself),
    2 = the ceiling (psi_UB).  None means the plain Frechet floor.
    """
    S = S.astype(float)
    Ysafe = np.where(np.isnan(Y), 0.0, Y)
    # saturate when the trim keeps the whole arm (kept fraction 1 in arm 1 /
    # discarded fraction 0 in arm 0): the whole-arm cell is kept in full
    keep_all_1 = a1 >= 0.9995 if a1 is not None else np.zeros(len(S), bool)
    keep_all_0 = a0 <= 0.0005 if a0 is not None else np.zeros(len(S), bool)
    ind1 = (S == 1) & ((Ysafe <= tau1) | keep_all_1)
    ind0 = (S == 1) & ((Ysafe >= tau0) | keep_all_0)

    # (5.8) numerator AIPW, treated lower trim
    eta1 = s1 * g1 + (D / p) * (S * Ysafe * ind1 - s1 * g1)
    # control upper trim
    eta0 = s0 * g0 + ((1 - D) / (1 - p)) * (S * Ysafe * ind0 - s0 * g0)

    # (5.6) type signals, one per revealing arm.  The IPW residuals are
    # centered at the kept mass the fitted (grid-capped) quantile surrogate
    # actually delivers -- s_d * u_eff -- rather than at pi_star, so that the
    # residual is exactly conditionally mean-zero under the fitted model in
    # every regime (on interior cells u_eff = a and s*u_eff = pi_star, so
    # nothing changes; the difference is confined to the narrow grid-cap
    # windows on binding cells).
    if a1 is not None:
        u1_eff = np.where(keep_all_1, 1.0,
                          np.clip(a1, U_GRID[0], U_GRID[-1]))
        center1 = s1 * u1_eff
    else:
        center1 = pi_star
    if a0 is not None:
        u0_eff = np.where(keep_all_0, 0.0,
                          np.clip(a0, U_GRID[0], U_GRID[-1]))
        center0 = s0 * (1.0 - u0_eff)
    else:
        center0 = pi_star
    psi_t1 = pi_star + (D / p) * (S * ind1 - center1)
    psi_t0 = pi_star + ((1 - D) / (1 - p)) * (S * ind0 - center0)

    # (5.13) floor feasibility signal
    resid = (D / p) * (S - s1) + ((1 - D) / (1 - p)) * (S - s0)
    psi_LB = (s0 + s1 - 1.0) * (s0 + s1 > 1.0) + resid

    # (5.14) ceiling feasibility signal
    d_dag = (s1 < s0).astype(int)             # binding arm argmin_d s_d
    s_dag = np.where(d_dag == 1, s1, s0)
    p_dag = np.where(d_dag == 1, p, 1 - p)
    psi_UB = s_dag + (D == d_dag) / p_dag * (S - s_dag)

    # active-facet selector (5.15); on the floor face the active share may
    # be the raised Masten-Poirier constant or the ceiling (clipped floor)
    if floor_face is None:
        psi_floor = psi_LB
    else:
        psi_floor = np.where(floor_face == 0, psi_LB,
                             np.where(floor_face == 1, pi_star, psi_UB))
    psi_pi = np.where(regime == 0, psi_t1,
                      np.where(regime == 1, psi_floor, psi_UB))

    A = (eta1 - eta0 + tau1 * (psi_pi - psi_t1)
         - tau0 * (psi_pi - psi_t0))
    return dict(eta1=eta1, eta0=eta0, psi_t=psi_t1, psi_t0=psi_t0,
                psi_pi=psi_pi, A=A)


# ----------------------------------------------------------------------
# the cross-fitted localized one-step (Definition 5.1 + Algorithm 5.2)
# ----------------------------------------------------------------------

def first_stage_cache(X, Y, S, D, w, folds, qr_alpha=1e-3, logit_C=1.0,
                      clip_s=0.01, clip_p=0.05):
    """Fit and pre-evaluate all parameter-free first stages, per rotation.

    Per Definition 5.1: for rotation r, (p, s0, s1) are fit on J2 and the
    conditional quantiles F_d^{-1} on J1.  Predictions are stored for the
    full sample; downstream code slices the relevant folds.
    """
    fold_ids = np.unique(folds)
    cache = []
    for r in range(3):
        J1 = folds == fold_ids[r]
        J2 = folds == fold_ids[(r + 1) % 3]
        ent = {"J1": J1, "J2": J2, "J3": folds == fold_ids[(r + 2) % 3]}
        ent["p"] = np.clip(
            Logit(logit_C).fit(X[J2], D[J2], w[J2]).predict(X),
            clip_p, 1 - clip_p)
        for d in (0, 1):
            m = J2 & (D == d)
            ent[f"s{d}"] = np.clip(
                Logit(logit_C).fit(X[m], S[m], w[m]).predict(X),
                clip_s, 1 - clip_s)
            mq = J1 & (D == d) & (S == 1)
            ent[f"Q{d}"] = QuantileCurves(qr_alpha).fit(
                X[mq], Y[mq], w[mq]).predict(X)
        cache.append(ent)
    return cache


def _box(piL_raw, piU_raw, floor_bar=None, collapse_ceiling=False):
    piL, piU = piL_raw.copy(), piU_raw.copy()
    if floor_bar is not None:                 # Masten-Poirier pointwise floor
        piL = np.minimum(np.maximum(piL, floor_bar), piU)
    if collapse_ceiling:                      # conditional monotonicity (App A)
        piL = piU.copy()
    return piL, piU


def one_sided_bound(X, Y, S, D, w, cache, negate=False,
                    lasso_alpha=1e-3, floor_bar=None, collapse_ceiling=False):
    """Cross-fitted localized one-step for ONE side (coded as the lower
    endpoint; `negate=True` runs the machinery on -Y for the other side --
    the caller negates the resulting estimate).

    Returns dict with the rotation-averaged ratio, pooled per-observation
    summands (A_i, psi_pi_i, psi_t_i, tau-diff), and diagnostics.
    """
    N = len(Y)
    dN = buffer_delta(N)
    rot_estimates, lam_inits = [], []
    A_pool = np.full(N, np.nan)
    pi_pool = np.full(N, np.nan)
    psit_pool = np.full(N, np.nan)
    taud_pool = np.full(N, np.nan)
    regime_pool = np.full(N, -1)

    if negate:
        Y = -Y

    for ent in cache:
        J1, J2, J3 = ent["J1"], ent["J2"], ent["J3"]
        s0_all, s1_all, p_all = ent["s0"], ent["s1"], ent["p"]
        if negate:
            # F_{-Y}^{-1}(u) = -F_Y^{-1}(1-u); U_GRID is symmetric.
            Q1_all, Q0_all = -ent["Q1"][:, ::-1], -ent["Q0"][:, ::-1]
        else:
            Q1_all, Q0_all = ent["Q1"], ent["Q0"]
        piL_all, piU_all = _box(np.maximum(s0_all + s1_all - 1.0, 0.0),
                                np.minimum(s0_all, s1_all),
                                floor_bar, collapse_ceiling)

        # --- lambda_init on J1 (Algorithm 5.1) ---
        cp1 = CellProgram(Q1_all[J1], Q0_all[J1], s1_all[J1], s0_all[J1],
                          piL_all[J1], piU_all[J1])
        if collapse_ceiling:
            lam_init, status = 0.0, "collapsed"  # box is a point; no root
        else:
            lam_init, status = dinkelbach_root(cp1, w[J1], N_total=N)
            if status == "undefined":        # the 0/0 case: no floor mass
                return dict(estimate=np.nan, status=status)
            if status == "degenerate":
                # Algorithm 5.1 step (0): no sign change on the bracket ->
                # report the degenerate floor value beta = 0 directly; the
                # one-step at an arbitrary lambda is undefined by the paper.
                return dict(estimate=0.0, status="degenerate")
        lam_inits.append(lam_init)

        # --- freeze pi*, tau at lambda_init; hinge regressions on J2 ---
        def frozen(idx):
            piLx, piUx = piL_all[idx], piU_all[idx]
            cp = CellProgram(Q1_all[idx], Q0_all[idx],
                             s1_all[idx], s0_all[idx], piLx, piUx)
            ps = np.clip(cp.pi_star(lam_init), piLx, piUx)
            t1, t0, a1, a0 = cp.thresholds(ps)
            return s0_all[idx], s1_all[idx], piLx, piUx, ps, t1, t0, a1, a0

        (s0_2, s1_2, piL_2, piU_2, ps_2,
         tau1_2, tau0_2, a1_2, a0_2) = frozen(J2)

        h1 = J2 & (D == 1) & (S == 1)
        h0 = J2 & (D == 0) & (S == 1)
        sel1 = (D[J2] == 1) & (S[J2] == 1)
        sel0 = (D[J2] == 0) & (S[J2] == 1)
        Xs1, mu1, sd1 = _standardize(X[h1])
        y1 = Y[h1] * ((Y[h1] <= tau1_2[sel1]) | (a1_2[sel1] >= 0.9995))
        las1 = Lasso(alpha=lasso_alpha, max_iter=50000).fit(
            Xs1, y1, sample_weight=w[h1])
        Xs0, mu0, sd0 = _standardize(X[h0])
        y0 = Y[h0] * ((Y[h0] >= tau0_2[sel0]) | (a0_2[sel0] <= 0.0005))
        las0 = Lasso(alpha=lasso_alpha, max_iter=50000).fit(
            Xs0, y0, sample_weight=w[h0])

        # --- debiased averages on J3 (Algorithm 5.2, step 3) ---
        (s0_3, s1_3, piL_3, piU_3, ps_3,
         tau1_3, tau0_3, a1_3, a0_3) = frozen(J3)
        p_3 = p_all[J3]
        reg_3 = classify(ps_3, piL_3, piU_3, dN)
        if collapse_ceiling:
            reg_3 = np.full(J3.sum(), 2)      # ceiling everywhere (App. A)
        g1_3 = las1.predict((X[J3] - mu1) / sd1)
        g0_3 = las0.predict((X[J3] - mu0) / sd0)

        # which face the (possibly raised) floor presents (Masten-Poirier)
        floor_face = None
        if floor_bar is not None:
            piL_raw = np.maximum(s0_3 + s1_3 - 1.0, 0.0)
            piU_raw = np.minimum(s0_3, s1_3)
            floor_face = np.where(piL_raw >= floor_bar, 0,
                                  np.where(floor_bar < piU_raw, 1, 2))

        sig = signals(Y[J3], S[J3], D[J3], w[J3], p_3, s0_3, s1_3,
                      piL_3, piU_3, ps_3, tau1_3, tau0_3, g1_3, g0_3, reg_3,
                      floor_face=floor_face, a1=a1_3, a0=a0_3)
        NL = np.average(sig["A"], weights=w[J3])
        piAT = np.average(sig["psi_pi"], weights=w[J3])
        if piAT <= 0:
            # No reading of the paper licenses a ratio with a non-positive
            # debiased fold denominator (Theorem 6.1 divides by piAT > 0;
            # Algorithm 5.1 gates the degenerate cases).  Surface the cell
            # as degenerate rather than averaging a sign-flipped ratio.
            return dict(estimate=np.nan, status="degenerate-denominator",
                        fold_denominator=float(piAT))
        rot_estimates.append(NL / piAT)

        A_pool[J3] = sig["A"]
        pi_pool[J3] = sig["psi_pi"]
        psit_pool[J3] = sig["psi_t"]
        taud_pool[J3] = tau1_3 - tau0_3
        regime_pool[J3] = reg_3

    # Point estimate: the pooled cross-fit ratio of (5.17)'s 1/N sums --
    # each observation enters exactly one J3, so the pooled averages run over
    # the full sample.  This is the same statistic the influence-function
    # variance and the Remark 6.3 multiplier bootstrap describe; per-rotation
    # ratios are kept as diagnostics.
    est = float(np.average(A_pool, weights=w) / np.average(pi_pool, weights=w))
    return dict(estimate=est, status="ok", lam_inits=lam_inits,
                A=A_pool, psi_pi=pi_pool, psi_t=psit_pool, tau_diff=taud_pool,
                regime=regime_pool, rot_estimates=rot_estimates)


def bound_interval(X, Y, S, D, w, cache, alpha=0.05, n_boot=400, seed=11,
                   **kwargs):
    """Both endpoints + inference.  Runs the machinery on Y (one side) and
    on -Y (other side, negated), sorts (Algorithm 5.2, step 4), and returns
    point bounds, IF standard errors, and the Imbens-Manski region (6.3).
    """
    lo_run = one_sided_bound(X, Y, S, D, w, cache, negate=False, **kwargs)
    hi_run = one_sided_bound(X, Y, S, D, w, cache, negate=True, **kwargs)
    if lo_run["status"] != "ok" or hi_run["status"] != "ok":
        return dict(status="degenerate/undefined",
                    detail=(lo_run["status"], hi_run["status"]))
    cand = [(lo_run["estimate"], lo_run, +1), (-hi_run["estimate"], hi_run, -1)]
    cand.sort(key=lambda z: z[0])             # Algorithm 5.2, step 4 (sort)
    (bL, runL, sgnL), (bU, runU, sgnU) = cand

    out = dict(status="ok", lower=bL, upper=bU)
    rng = np.random.default_rng(seed)
    # Remark 6.3 multiplier bootstrap: holding the cross-fitted Psi_i and
    # piAT fixed, draw e_i ~ Exp(1) and recompute the LINEAR statistic
    # (weighted mean of e_i * Psi_i) / piAT.  The denominator is NOT
    # re-randomized, so the draws have all moments; shared draws across the
    # two endpoints make the joint covariance Sigma_LU of Remark 6.2
    # estimable from the same summands.
    e = rng.exponential(1.0, size=(n_boot, len(w)))
    wsum = w.sum()
    all_draws = {}
    for name, b, run, sign in (("lower", bL, runL, sgnL),
                               ("upper", bU, runU, sgnU)):
        # Psi_i at the root: sign * (A_i - lambda_hat * psi_pi_i), with
        # lambda_hat = sign * b the root in the (possibly negated) run units.
        # With the pooled-ratio estimate the weighted mean of Psi is zero.
        Psi = sign * (run["A"] - (sign * b) * run["psi_pi"])
        piAT = np.average(run["psi_pi"], weights=w)
        c = w / w.sum()                        # weighted IF variance
        var = np.sum(c ** 2 * (Psi - np.sum(c * Psi)) ** 2) / piAT ** 2
        out[f"se_{name}"] = float(np.sqrt(var))
        draws = (e * w) @ Psi / wsum / piAT    # fixed-denominator, linear
        all_draws[name] = draws
        out[f"se_boot_{name}"] = float(np.std(draws))
        # Remark 6.3 quantile critical values (diagnostic alternative)
        out[f"boot_q05_{name}"] = float(np.quantile(draws, 0.05))
        out[f"boot_q95_{name}"] = float(np.quantile(draws, 0.95))
    out["cov_boot_LU"] = float(np.cov(all_draws["lower"],
                                      all_draws["upper"])[0, 1])
    q = norm.ppf(np.sqrt(1 - alpha))
    # IM region (6.3) with the multiplier-bootstrap endpoint SEs of
    # Remark 6.2/6.3 and the common critical value q = Phi^-1(sqrt(1-alpha)).
    out["IM_lo"] = bL - q * out["se_boot_lower"]
    out["IM_hi"] = bU + q * out["se_boot_upper"]
    out["q"] = q
    out["regime_shares"] = {k: float(np.mean(lo_run["regime"] == v))
                            for k, v in (("int", 0), ("flr", 1), ("cei", 2))}
    out["lam_inits"] = (lo_run["lam_inits"], hi_run["lam_inits"])
    out["piAT"] = float(np.average(lo_run["psi_pi"], weights=w))
    return out


# ----------------------------------------------------------------------
# classical Lee bounds, no covariates (Table 5, column 1)
# ----------------------------------------------------------------------

def _wquantile(y, w, q):
    o = np.argsort(y)
    cw = np.cumsum(w[o]) / w[o].sum()
    return np.interp(q, cw, y[o])


def basic_lee(Y, S, D, w, n_boot=400, alpha=0.05, seed=11):
    """Unconditional Lee (2009) trimming bounds on the wage effect.

    Maintains the classical monotonicity direction S(1) >= S(0): the treated
    arm is trimmed by the proportion 1 - s0/s1, clamped at zero when s1 < s0
    (then no trimming occurs and the interval collapses to the difference in
    mean log wages -- the convention that reproduces the week-45 point
    interval of the paper's Table 5, column 1).  Bootstrap SEs; IM region.
    """
    def endpoints(Y, S, D, w):
        i1 = (D == 1) & (S == 1)
        i0 = (D == 0) & (S == 1)
        s1 = np.average(S[D == 1], weights=w[D == 1])
        s0 = np.average(S[D == 0], weights=w[D == 0])
        y1, w1 = Y[i1], w[i1]
        y0, w0 = Y[i0], w[i0]
        frac = min(s0 / s1, 1.0)        # kept share of the treated arm
        qlo = _wquantile(y1, w1, frac)
        qhi = _wquantile(y1, w1, 1 - frac)
        m0 = np.average(y0, weights=w0)
        lo = np.average(y1[y1 <= qlo], weights=w1[y1 <= qlo]) - m0
        hi = np.average(y1[y1 >= qhi], weights=w1[y1 >= qhi]) - m0
        return lo, hi

    lo, hi = endpoints(Y, S, D, w)
    rng = np.random.default_rng(seed)
    N = len(Y)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, N, N)
        bs.append(endpoints(Y[i], S[i], D[i], w[i]))
    bs = np.array(bs)
    seL, seU = bs[:, 0].std(), bs[:, 1].std()
    q = norm.ppf(np.sqrt(1 - alpha))
    return dict(lower=lo, upper=hi, se_lower=seL, se_upper=seU,
                IM_lo=lo - q * seL, IM_hi=hi + q * seU)

"""Discrete-heuristic LP bounds on the always-takers' wage effect.

Implements Section 2.3 of the paper: partition X into cells C_j, discretize
the log wage onto a K-bin equal-mass grid, and solve the Charnes-Cooper LP
(2.17) for the fractional program (2.16), extended to the *difference* target
of Section 5.4.3: for the lower endpoint the treated-arm cell measure is
trimmed from below and the control-arm cell measure from above at the same
mass pi_j, which is a free variable in the cell-averaged Frechet box
[pi_L,j, pi_U,j] of (2.19).

Conservative endpoint convention of Proposition 2.2: mass relocated to the
bin's *lower* endpoint for a lower trim and to the bin's *upper* endpoint for
an upper trim, so the discrete program is a guaranteed outer bound and the
interval tightens as the grid refines.

The program (lower endpoint; the upper endpoint swaps the endpoint
conventions and the trim directions and maximizes):

    min  sum_j w_j [ sum_k ylo_k r1_{j,k}  -  sum_k yhi_k r0_{j,k} ]
         ------------------------------------------------------------
                       sum_j w_j pi_j

    s.t. 0 <= r1_{j,k} <= m1_{j,k},   0 <= r0_{j,k} <= m0_{j,k},
         sum_k r1_{j,k} = pi_j,       sum_k r0_{j,k} = pi_j,
         pi_L,j <= pi_j <= pi_U,j,

with m_{d,j,k} = Pr(Y in B_k, S=1 | D=d, C_j), w_j = Pr(X in C_j),
pi_L,j = (s0_j + s1_j - 1)_+, pi_U,j = min(s0_j, s1_j).

Charnes-Cooper: u = t*r, v_j = t*pi_j, t = 1/(sum_j w_j pi_j) turns the ratio
into an LP solved at standard cost (HiGHS).  Computed unweighted, per the
table notes.
"""

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix


def cell_moments(Y, S, D, cells, K):
    """Equal-mass wage grid (pooled selected sample) and cell moments.

    Returns dict with bin edges, per-cell weights w_j, selection rates s_{d,j},
    gridded masses m_{d,j,k}, and the Frechet box endpoints.
    """
    sel = S == 1
    y_sel = Y[sel]
    # K equal-mass bins of the pooled selected log wage
    edges = np.quantile(y_sel, np.linspace(0, 1, K + 1))
    edges[0], edges[-1] = y_sel.min(), y_sel.max()
    edges = np.unique(edges)
    Keff = len(edges) - 1

    cell_ids = np.unique(cells)
    J = len(cell_ids)
    N = len(Y)

    w = np.zeros(J)
    s = np.zeros((2, J))
    m = np.zeros((2, J, Keff))
    for j, cj in enumerate(cell_ids):
        in_cell = cells == cj
        w[j] = in_cell.mean()
        for d in (0, 1):
            arm = in_cell & (D == d)
            n_arm = arm.sum()
            s[d, j] = S[arm].mean()
            ya = Y[arm & sel]
            if len(ya) > 0:
                idx = np.clip(np.searchsorted(edges, ya, side="right") - 1,
                              0, Keff - 1)
                cnt = np.bincount(idx, minlength=Keff).astype(float)
                m[d, j] = cnt / n_arm  # Pr(Y in B_k, S=1 | D=d, C_j)

    piL = np.maximum(s[0] + s[1] - 1.0, 0.0)
    piU = np.minimum(s[0], s[1])
    return dict(edges=edges, K=Keff, J=J, w=w, s=s, m=m, piL=piL, piU=piU,
                cell_ids=cell_ids)


def _solve_cc(mom, lower=True):
    """Charnes-Cooper LP for one endpoint of the wage-effect bound.

    lower=True : minimize  (lower trim of D=1) - (upper trim of D=0);
                 treated mass valued at bin lower endpoints,
                 control mass valued at bin upper endpoints.
    lower=False: maximize  (upper trim of D=1) - (lower trim of D=0);
                 conventions swapped.
    """
    K, J = mom["K"], mom["J"]
    w, m, piL, piU = mom["w"], mom["m"], mom["piL"], mom["piU"]
    edges = mom["edges"]
    ylo, yhi = edges[:-1], edges[1:]

    if piL @ w <= 0:
        return np.nan  # 0/0 gate: no cell-averaged Frechet floor mass

    # variable order: u1 (J*K), u0 (J*K), v (J), t (1)
    nvar = 2 * J * K + J + 1
    iu1 = lambda j, k: j * K + k
    iu0 = lambda j, k: J * K + j * K + k
    iv = lambda j: 2 * J * K + j
    it = 2 * J * K + J

    c = np.zeros(nvar)
    for j in range(J):
        if lower:
            c[[iu1(j, k) for k in range(K)]] = w[j] * ylo
            c[[iu0(j, k) for k in range(K)]] = -w[j] * yhi
        else:  # maximize -> minimize the negative
            c[[iu1(j, k) for k in range(K)]] = -w[j] * yhi
            c[[iu0(j, k) for k in range(K)]] = w[j] * ylo

    # equalities: sum_j w_j v_j = 1 ; sum_k u_d(j,.) = v_j
    Aeq = lil_matrix((1 + 2 * J, nvar))
    beq = np.zeros(1 + 2 * J)
    for j in range(J):
        Aeq[0, iv(j)] = w[j]
    beq[0] = 1.0
    for j in range(J):
        for k in range(K):
            Aeq[1 + j, iu1(j, k)] = 1.0
            Aeq[1 + J + j, iu0(j, k)] = 1.0
        Aeq[1 + j, iv(j)] = -1.0
        Aeq[1 + J + j, iv(j)] = -1.0

    # inequalities: u_d <= m_d t ; piL t <= v <= piU t
    nin = 2 * J * K + 2 * J
    Ain = lil_matrix((nin, nvar))
    bin_ = np.zeros(nin)
    r = 0
    for j in range(J):
        for k in range(K):
            Ain[r, iu1(j, k)] = 1.0
            Ain[r, it] = -mom["m"][1, j, k]
            r += 1
            Ain[r, iu0(j, k)] = 1.0
            Ain[r, it] = -mom["m"][0, j, k]
            r += 1
    for j in range(J):
        Ain[r, iv(j)] = 1.0
        Ain[r, it] = -piU[j]
        r += 1
        Ain[r, iv(j)] = -1.0
        Ain[r, it] = piL[j]
        r += 1

    res = linprog(c, A_ub=Ain.tocsr(), b_ub=bin_,
                  A_eq=Aeq.tocsr(), b_eq=beq,
                  bounds=[(0, None)] * nvar, method="highs")
    if not res.success:
        return np.nan
    val = res.fun if lower else -res.fun
    return val


def lp_wage_effect_bounds(Y, S, D, cells, K):
    """Return (lower, upper) discrete-heuristic bounds at K wage bins."""
    mom = cell_moments(Y, S, D, cells, K)
    lo = _solve_cc(mom, lower=True)
    hi = _solve_cc(mom, lower=False)
    return lo, hi, mom

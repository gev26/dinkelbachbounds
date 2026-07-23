"""Produce the wage-granularity table (heuristic LP bounds, Section 2.3/3.2).

Rows: weeks 45, 90, 100.  Panels: 5 / 10 / 15 equal-mass wage bins.
Columns: (1) baseline annual earnings quantile-binned into five cells,
         (2) household income on its native brackets.
Computed unweighted with the conservative endpoint convention, per the
paper's table notes.
"""

import numpy as np
import pandas as pd

import data_prep as dp
from lp_bounds import lp_wage_effect_bounds

def main():
    df = dp.load_data()
    D = dp.treatment(df)
    specs = {"Earnings": dp.earnings_cells(df), "Income": dp.income_cells(df)}

    rows = []
    for K in (5, 10, 15):
        for week in (45, 90, 100):
            for name, cells in specs.items():
                S, Y = dp.outcome_selection(df, week)
                lo, hi, mom = lp_wage_effect_bounds(Y, S, D, cells, K)
                rows.append(dict(bins=K, week=week, covariate=name,
                                 lower=lo, upper=hi, width=hi - lo,
                                 EpiL=float(mom["piL"] @ mom["w"]),
                                 EpiU=float(mom["piU"] @ mom["w"])))
                print(f"K={K:2d} week={week:3d} {name:9s} "
                      f"[{lo:+.3f}, {hi:+.3f}]  width={hi-lo:.2f}  "
                      f"E[piL]={mom['piL'] @ mom['w']:.4f}")
    out = pd.DataFrame(rows)
    out.to_csv("../results/table2_lp_granularity.csv", index=False)
    print("\nSaved to results/table2_lp_granularity.csv")

if __name__ == "__main__":
    main()

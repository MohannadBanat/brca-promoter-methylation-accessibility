"""H3 proportions among silenced genes that have both a promoter methylation value and a promoter peak (Table 3; deviation 8).
Run: python3 -I scripts/05b_h3_both_measures.py results"""
import json, sys
import pandas as pd
from statsmodels.stats.proportion import proportion_confint
R = sys.argv[1]; out = {}
for lab in ['all', 'HM450']:
    c = pd.read_csv(f'{R}/H3_{lab}_classes.csv'); c = c[c.has_meth & c.has_peak]; out[lab] = {}
    for d, g in c.groupby('direction'):
        n = len(g); r = {'n': n}
        for k in ['methylation', 'accessibility', 'both', 'neither']:
            x = int((g['class'] == k).sum()); lo, hi = proportion_confint(x, n, method='wilson')
            r[k] = [x, round(x / n, 3), round(lo, 3), round(hi, 3)]
        out[lab][d] = r
json.dump(out, open(f'{R}/H3_both_measures_table.json', 'w'), indent=1)
print(json.dumps(out, indent=1))

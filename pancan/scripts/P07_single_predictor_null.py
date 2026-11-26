"""POST HOC, EXPLORATORY (deviations log entry 4). Single-predictor mismatched-promoter null for M and A, nine eligible types.
Imports the feature builder of P05 so definitions are identical. Run: python3 -I P07_single_predictor_null.py <processed> <ref> <results>"""
import json, sys, importlib.util
import numpy as np, pandas as pd
P, REF, R = sys.argv[1:4]
spec = importlib.util.spec_from_file_location('p05', __file__.replace('P07_single_predictor_null.py', 'P05_posthoc_robustness.py'))
src = open(spec.origin).read().split("out = {'note'")[0]          # helpers and features() only
sys.argv = [spec.origin, P, REF, R, '']; ns = {}; exec(compile(src, spec.origin, 'exec'), ns)
rng = np.random.default_rng(20261127); out = {}
for T in ns['ELIGIBLE']:
    gc, M, A, E, npr, cohort = ns['features'](T); v = ns['variable'](M, A, E); M, A, E = M[v], A[v], E[v]; G, n = E.shape
    m = ns['cv_r2'](ns['design']([ns['zs'](M)], G, n), E); a = ns['cv_r2'](ns['design']([ns['zs'](A)], G, n), E)
    dec = pd.qcut(E.mean(1), 10, labels=False, duplicates='drop'); nm, na = [], []
    for _ in range(100):
        perm = np.arange(G)
        for q in np.unique(dec):
            ii = np.where(dec == q)[0]; perm[ii] = rng.permutation(ii)
        nm.append(np.median(ns['cv_r2'](ns['design']([ns['zs'](M[perm])], G, n), E)))
        na.append(np.median(ns['cv_r2'](ns['design']([ns['zs'](A[perm])], G, n), E)))
    out[T] = {'n_tumours': n, 'n_genes': G, 'observed_median_M': float(np.median(m)), 'null_median_M_mean': float(np.mean(nm)),
              'null_median_M_95': [float(np.percentile(nm, 2.5)), float(np.percentile(nm, 97.5))],
              'observed_median_A': float(np.median(a)), 'null_median_A_mean': float(np.mean(na)),
              'null_median_A_95': [float(np.percentile(na, 2.5)), float(np.percentile(na, 97.5))]}
    print(T, n, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in out[T].items() if 'median' in k and not k.endswith('95')}, flush=True)
json.dump({'note': 'POST HOC, EXPLORATORY (deviation 4).', 'per_type': out}, open(f'{R}/posthoc_single_predictor_null.json', 'w'), indent=1)

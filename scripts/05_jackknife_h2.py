"""Delete-one-tumour jackknife CIs for H2 (deviation 7). Gene set fixed as selected on the full cohort.
Run: python3 -I scripts/05_jackknife_h2.py results"""
import json, sys
import numpy as np

R = sys.argv[1]
z = np.load(f'{R}/H2_matrices.npz', allow_pickle=True)
M, A, E, is450, pur, hcp = z['M'], z['A'], z['E'], z['is450'], z['purity'], z['hcp']

def cv_r2(Xs, y):
    XtX = Xs.transpose(0, 2, 1) @ Xs + 1e-9 * np.eye(Xs.shape[2]); Ai = np.linalg.inv(XtX)
    b = Ai @ (Xs.transpose(0, 2, 1) @ y[..., None]); res = y - (Xs @ b)[..., 0]
    h = np.einsum('gnk,gkl,gnl->gn', Xs, Ai, Xs); e = res / np.clip(1 - h, 1e-6, None)
    return 1 - (e ** 2).sum(1) / ((y - y.mean(1, keepdims=True)) ** 2).sum(1)
def zs(X): return (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-12)
def stats_(M, A, E, cov, gsel=None):
    G, n = E.shape; ex = [np.broadcast_to(c, (G, n)) for c in cov]
    d = lambda cols: np.stack([np.ones((G, n))] + cols, axis=2)
    m = cv_r2(d([zs(M)] + ex), E); a = cv_r2(d([zs(A)] + ex), E); ma = cv_r2(d([zs(M), zs(A)] + ex), E)
    if gsel is not None: m, a, ma = m[gsel], a[gsel], ma[gsel]
    return {'median_cvR2_M': np.median(m), 'median_cvR2_A': np.median(a), 'median_cvR2_MA': np.median(ma),
            'median_MA_minus_A': np.median(ma - a), 'median_MA_minus_M': np.median(ma - m), 'frac_A_gt_M': np.mean(a > m)}

def jack(M, A, E, covs, gsel=None):
    n = E.shape[1]; full = stats_(M, A, E, covs, gsel); reps = []
    for i in range(n):
        k = np.arange(n) != i
        reps.append(stats_(M[:, k], A[:, k], E[:, k], [c[k] for c in covs], gsel))
    out = {}
    for key in full:
        v = np.array([r[key] for r in reps]); se = np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))
        out[key] = {'estimate': float(full[key]), 'jackknife_SE': float(se),
                    'CI95': [float(full[key] - 1.96 * se), float(full[key] + 1.96 * se)]}
    return out

res = {'n_genes': int(E.shape[0])}
res['pooled'] = jack(M, A, E, [is450])
res['pooled_HCP'] = jack(M, A, E, [is450], hcp)
res['pooled_LCP'] = jack(M, A, E, [is450], ~hcp)
for lab, sel in [('HM450', is450 == 1), ('HM27', is450 == 0)]:
    res[lab] = jack(M[:, sel], A[:, sel], E[:, sel], [])
ok = ~np.isnan(pur)
res['purity_adjusted'] = jack(M[:, ok], A[:, ok], E[:, ok], [is450[ok], pur[ok]])
res['purity_adjusted']['n_tumours'] = int(ok.sum())
json.dump(res, open(f'{R}/H2_jackknife.json', 'w'), indent=1)
for k in ['pooled', 'HM450', 'HM27', 'purity_adjusted']:
    print(k, {kk: (round(v['estimate'], 3), [round(x, 3) for x in v['CI95']]) for kk, v in res[k].items() if isinstance(v, dict)})

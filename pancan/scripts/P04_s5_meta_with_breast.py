"""S5 (pre-specified, second part): random-effects meta-analysis of D and G including breast, reported beside R1-R3.
Breast enters with its re-run values through the pan-cancer code (masked probes; identical to unmasked, see breast_rerun.json).
Run: python3 -I P04_s5_meta_with_breast.py <results dir>"""
import json, sys
import numpy as np
R = sys.argv[1]
def dersimonian_laird(est, se):
    est, se = np.asarray(est, float), np.asarray(se, float); w = 1 / se ** 2; k = len(est)
    fe = np.sum(w * est) / np.sum(w); Q = float(np.sum(w * (est - fe) ** 2))
    tau2 = max(0.0, (Q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
    wr = 1 / (se ** 2 + tau2); re = np.sum(wr * est) / np.sum(wr); sre = np.sqrt(1 / np.sum(wr))
    I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    return {'estimate': float(re), 'SE': float(sre), 'CI95': [float(re - 1.96 * sre), float(re + 1.96 * sre)],
            'tau2': float(tau2), 'Q': Q, 'I2': float(I2), 'k': int(k)}
p = json.load(open(f'{R}/pancan_results.json')); b = json.load(open(f'{R}/breast_rerun.json'))['masked_M_general']
types = p['meta_types_included']; out = {'types': types + ['BRCA']}
for key in ('D_median_A_minus_M', 'G_median_MA_minus_A'):
    est = [p['per_type'][t]['H2'][key]['estimate'] for t in types] + [b[key]['estimate']]
    se = [p['per_type'][t]['H2'][key]['SE'] for t in types] + [b[key]['SE']]
    out[key] = dersimonian_laird(est, se)
    chk = dersimonian_laird(est[:-1], se[:-1])
    ref = p['R1_pooled_D' if key.startswith('D') else 'R3_pooled_G']
    assert abs(chk['estimate'] - ref['estimate']) < 1e-12, 'meta function differs from P02'
out['D_positive_including_breast'] = int(sum(p['per_type'][t]['H2']['D_median_A_minus_M']['estimate'] > 0 for t in types)) + int(b['D_median_A_minus_M']['estimate'] > 0)
json.dump(out, open(f'{R}/S5_meta_with_breast.json', 'w'), indent=1)
print(json.dumps(out, indent=1))

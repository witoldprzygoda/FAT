#!/usr/bin/env python3
# ============================================================================
# flex_tables.py -- package flex of report_all: tables of the flexible SMASH description
# ============================================================================
# Reads /hdd3/przygoda/HADES/FAT_out/all_report/flex/flex_results.txt (flex_fit.C) and flex_bias_<tag>.txt (flex_toys.C),
# writes docs/latex/tables_all/flex_scales.{txt,tex} (every scale: value, errors, constrained / degenerate / bias test,
# Part III equivalent, alternative model, variants), flex_grid.{txt,tex} (model ladder and the D11 choice),
# flex_variants.txt (every scale in every variant), flex_ratios.txt (ratios and sums).
# Run (from /home/przygoda/HADES/FAT/FAT): python3 plots/all_report/flex_tables.py
# ============================================================================
import math
import re

OUT = '/hdd3/przygoda/HADES/FAT_out/all_report/flex'
TAB = '/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all'


def read_results():
    fits, extra = {}, {'d11': [], 'grid': [], 'ratio': [], 'sumscale': [], 'p3equiv_err': {}, 'p3equiv': {}}
    lines = open(f'{OUT}/flex_results.txt').read().split('\n')
    i = 0
    while i < len(lines):
        w = lines[i].split()
        i += 1
        if not w:
            continue
        k = w[0]
        if k in ('names', 'values', 'err_stat', 'err_mc', 'err_total', 'fixed', 'mc_shift'):
            f = fits.setdefault(w[1], {})
            f[k] = w[2:] if k == 'names' else [float(x) for x in w[2:]]
        elif k == 'text':
            fits.setdefault(w[1], {})['text'] = ' '.join(w[2:])
        elif k == 'gof':
            f = fits.setdefault(w[1], {})
            for q in range(2, len(w) - 1):
                if w[q] in ('nbins', 'nfree', 'ndf', 'chi2P', 'chi2BC', 'chi2GLShyb', 'chi2GLSdata'):
                    f[w[q]] = float(w[q + 1])
        elif k == 'cor_total':
            n = int(w[2])
            m = []
            for a in range(n):
                m.append([float(x) for x in lines[i].split()])
                i += 1
            fits.setdefault(w[1], {})['cor'] = m
        elif k in ('cov_total', 'cov_stat'):
            i += int(w[2])
        elif k == 'nominal':
            extra['nominal'] = w[1]
        elif k == 'alternative':
            extra['alternative'] = w[1]
        elif k == 'd11':
            extra['d11'].append(' '.join(w[1:]))
        elif k == 'grid':
            extra['grid'].append(w[1:])
        elif k == 'ratio':
            extra['ratio'].append(lines[i - 1])
        elif k == 'sumscale':
            extra['sumscale'].append(lines[i - 1])
        elif k == 'p3equiv':
            for q in range(1, len(w) - 1, 2):
                extra['p3equiv'][w[q]] = float(w[q + 1])
        elif k == 'p3equiv_err':
            extra['p3equiv_err'][w[1]] = (float(w[2]), float(w[4]), float(w[6]), float(w[8]))
    return fits, extra


def read_bias(tag):
    out = {0: {}, 1: {}}
    scen = -1
    gof = {}
    for l in open(f'{OUT}/flex_bias_{tag}.txt'):
        if l.startswith('# scenario 0'):
            scen = 0
        elif l.startswith('# scenario 1'):
            scen = 1
        w = l.split()
        if not w:
            continue
        if w[0] == 'gof':
            gof[int(w[2])] = l.strip()
        if scen >= 0 and (w[0].startswith('s_') or w[0].startswith('k_') or w[0].startswith('b')) and len(w) >= 14 and w[1] != 'truth':
            out[scen][w[0]] = dict(truth=float(w[1]), sstat=float(w[2]), smc=float(w[3]), mcshift=float(w[4]), rf=float(w[5]), rfsig=float(w[6]),
                                   pullmean=float(w[7]), pullwid=float(w[8]), cover=float(w[9]), meth=float(w[10]), verdict=w[13])
    return out, gof


def val(f, n):
    if n not in f.get('names', []):
        return None
    k = f['names'].index(n)
    return f['values'][k], f['err_stat'][k], f['err_mc'][k], f['err_total'][k]


fits, ex = read_results()
NOM, ALT = ex['nominal'], ex['alternative']
bias, bgof = read_bias(NOM)
biasA, bgofA = read_bias(ALT)
F = fits[NOM]
names = F['names']
cor = F['cor']
p3map = {'s_C2': 'C2', 's_C3': 'C3C4', 's_C4': 'C3C4', 's_C5': 'C5', 's_Geta': 'Geta', 's_GD_om': 'GD', 's_GD_nom': 'GD', 's_GC': 'GC'}
labels = {'s_C2': 'C2 bump, eta->pi+pi-pi0, pi0 Dalitz', 's_C3': 'C3 bump, pi0 -> gamma gamma, conversion', 's_C4': 'C4 peak, eta->pi+pi-gamma, conversion',
          's_C5': 'C5 peak, direct eta->pi+pi-e+e-', 's_Geta': 'G_eta other eta (pair not C2-C5)', 's_GD_om': 'G_D non-eta Dalitz pair, omega/rho0 pions',
          's_GD_nom': 'G_D non-eta Dalitz pair, other pions', 's_GC': 'G_C non-eta conversion pair', 'k_E2': 'k_E2 CB scale (constraint 1 +- 0.05)'}


def cls(v, e):
    r = e / abs(v) if v != 0 else 99
    return ('constrained' if r < 0.3 else 'weak' if r < 0.5 else 'not constrained'), r


rows = []
for k, n in enumerate(names):
    v, es, em, et = F['values'][k], F['err_stat'][k], F['err_mc'][k], F['err_total'][k]
    c, r = cls(v, et)
    deg = [f"{names[q]} ({cor[k][q]:+.2f})" for q in range(len(names)) if q != k and abs(cor[k][q]) > 0.8]
    b0 = bias[0].get(n)
    b1 = bias[1].get(n)
    p3 = ex['p3equiv_err'].get(p3map.get(n, ''), None)
    if n == 'k_E2':
        p3 = (ex['p3equiv'].get('k_E2', 0.0), 0, 0, 0.068)
    alt = val(fits[ALT], n) or val(fits[ALT], 's_GD' if n.startswith('s_GD') else n)
    e2s = val(fits.get('V_E2split', {}), n)
    e1k = val(fits.get('V_E1k', {}), n)
    corr = None
    if b0 and n != 'k_E2' and abs(b0['truth']) > et:
        corr = v / (1 + b0['rf'] / 100.)
    rows.append(dict(n=n, v=v, es=es, em=em, et=et, cls=c, rel=r, deg=deg, b0=b0, b1=b1, p3=p3, alt=alt, e2s=e2s, e1k=e1k, corr=corr))


def fmt(t):
    return '-' if t is None else f"{t[0]:.3f} +- {t[3]:.3f}"


# ---------------------------------------------------------------- flex_scales.txt
with open(f'{TAB}/flex_scales.txt', 'w') as o:
    o.write('# flex_scales.txt -- plots/all_report/flex_tables.py (from flex_fit.C, flex_toys.C)\n')
    o.write(f'# flexible SMASH description of the data, nominal model {NOM}: N = k_E2 CB + sum_i s_i T_i, F10 E2 (M(ee) >= 10 MeV), M4 0.40-0.70,\n')
    o.write('# 5 M(pi+pi-) slices, 350 bins; S0, isBest >= 0, L = 6079.2 nb^-1; s = 1 is SMASH sigma x PDG BR x V3 (M4\' = 0.99615 M4).\n')
    o.write('# errors: stat = sandwich (event groups + CB replicas + k_E2 constraint), MC = spread of 100 simulation-replica refits.\n')
    o.write('# constrained: sigma_tot/|s| < 30 %; weak 30-50 %; degenerate: |corr(stat+MC)| > 0.8 with another scale.\n')
    o.write('# bias test (flex_toys.C, 500 pseudo-experiments each): truth = full simulation x the fitted scales (scen 0) or x 1 (scen 1);\n')
    o.write('# realistic fit with bootstrap-replica templates; PASS = |bias| < max(0.5 sigma_stat, 10 % truth) and 68 % coverage in 0.60-0.76.\n')
    o.write('# attenuation-corrected = s / (1 + RF bias of scen 0): INFORMATION only (WP4: the corrected estimator under-covers).\n')
    gf = F
    o.write(f"# goodness of fit: Pearson {gf['chi2P']:.1f}, hybrid GLS {gf['chi2GLShyb']:.1f}, ndf {int(gf['ndf'])}; pseudo-data (scen 0): {bgof.get(0, '')}\n\n")
    hdr = f"{'scale':10s} {'value':>7s} {'stat':>6s} {'MC':>6s} {'tot':>6s} {'rel':>5s}  {'verdict':16s} {'bias0[%]':>8s} {'cov0':>5s} {'test0':15s} {'test1':15s} {'att.corr':>8s}  {'Part III equiv.':16s} {'alt '+ALT:18s} {'E2 split':16s} {'+E1 k_E1':16s}  degenerate with"
    o.write(hdr + '\n')
    for r in rows:
        b0, b1 = r['b0'], r['b1']
        p3s = f"{r['p3'][0]:.3f} +- {r['p3'][3]:.3f}" if r['p3'] else '-'
        bias0 = b0['rf'] if b0 else float('nan')
        cov0 = b0['cover'] if b0 else float('nan')
        t0 = b0['verdict'] if b0 else 'nuisance'
        t1 = b1['verdict'] if b1 else 'nuisance'
        ac = r['corr'] if r['corr'] is not None else float('nan')
        dg = ', '.join(r['deg']) if r['deg'] else '-'
        o.write(f"{r['n']:10s} {r['v']:7.3f} {r['es']:6.3f} {r['em']:6.3f} {r['et']:6.3f} {100*r['rel']:4.0f}%  {r['cls']:16s} "
                f"{bias0:8.1f} {cov0:5.2f} {t0:15s} {t1:15s} {ac:8.3f}  {p3s:16s} {fmt(r['alt']):18s} {fmt(r['e2s']):16s} {fmt(r['e1k']):16s}  {dg}\n")
    o.write('\nmeaning of the scales:\n')
    for r in rows:
        o.write(f"  {r['n']:10s} {labels.get(r['n'], '')}\n")

# ---------------------------------------------------------------- flex_scales.tex
def tx(s):
    return s.replace('_', r'\_')


with open(f'{TAB}/flex_scales.tex', 'w') as o:
    o.write('% flex_scales.tex -- plots/all_report/flex_tables.py\n')
    o.write('\\begin{tabular}{lrrrrlrll rr}\n\\toprule\n')
    o.write('scale & value & stat & MC & $\\sigma/s$ & verdict & bias [\\%] & bias test & degenerate with & Part III & alt.\\ ' + tx(ALT) + ' \\\\\n\\midrule\n')
    for r in rows:
        b0 = r['b0']
        deg = ', '.join(x.split(' ')[0] for x in r['deg']) if r['deg'] else '--'
        bs = f"{b0['rf']:+.0f}" if b0 else '--'
        bv = b0['verdict'] if b0 else 'nuisance'
        p3s = f"{r['p3'][0]:.2f}" if r['p3'] else '--'
        als = f"{r['alt'][0]:.2f}" if r['alt'] else '--'
        o.write(f"\\texttt{{{tx(r['n'])}}} & {r['v']:.3f} & {r['es']:.3f} & {r['em']:.3f} & {100*r['rel']:.0f}\\% & {r['cls']} & "
                f"{bs} & {bv} & \\texttt{{{tx(deg)}}} & {p3s} & {als} \\\\\n")
    o.write('\\bottomrule\n\\end{tabular}\n')

# ---------------------------------------------------------------- flex_grid
with open(f'{TAB}/flex_grid.txt', 'w') as o:
    o.write('# flex_grid.txt -- plots/all_report/flex_tables.py: the model ladder (F10 E2, M4 0.40-0.70, 350 bins) and the D11 choice\n')
    o.write('# hybrid GLS chi2 = data covariance (event groups + CB) + BB-lite variance of the raw templates + replica covariance of the smoothed ones\n')
    o.write(f"{'model':8s} {'free':>4s} {'Pearson':>8s} {'hybGLS':>8s}  description\n")
    for t in ['ASIS', 'ASIS9', 'P3', 'F7', 'F7b', 'F7bj', 'F8', 'F8b', 'F8bj', 'F9', 'F9b', 'F9bj', 'F7cj']:
        f = fits.get(t, {})
        if 'chi2P' not in f:
            continue
        o.write(f"{t:8s} {int(f.get('nfree', 0)):4d} {f['chi2P']:8.1f} {f['chi2GLShyb']:8.1f}  {f.get('text', '')}\n")
    o.write('\nD11 steps (physics first; then the tilt-first order for the alternative):\n')
    for l in ex['d11']:
        o.write('  ' + l + '\n')
    o.write(f'\nnominal {NOM}, alternative {ALT}\n')
with open(f'{TAB}/flex_grid.tex', 'w') as o:
    o.write('% flex_grid.tex -- plots/all_report/flex_tables.py\n\\begin{tabular}{lrrrl}\n\\toprule\nmodel & free & Pearson $\\chi^2$ & hybrid GLS $\\chi^2$ & content \\\\\n\\midrule\n')
    for t in ['ASIS', 'P3', 'F7', 'F7b', 'F7bj', 'F8', 'F8b', 'F8bj', 'F9', 'F9b', 'F9bj', 'F7cj']:
        f = fits.get(t, {})
        if 'chi2P' not in f:
            continue
        o.write(f"{tx(t)} & {int(f.get('nfree', 0))} & {f['chi2P']:.0f} & {f['chi2GLShyb']:.0f} & {tx(f.get('text', ''))[:90]} \\\\\n")
    o.write('\\bottomrule\n\\end{tabular}\n')

# ---------------------------------------------------------------- variants
order = ['P3', NOM, ALT, 'F7', 'F7b', 'F8b', 'F8bj', 'F9', 'F9b', 'F9bj', 'F7cj', 'V_E2split', 'V_E1k', 'V_E1T3', 'V_E1kT3', 'V_E1k3', 'V_r42', 'V_r45', 'V_rawBg', 'V_T3all',
         'V_r100', 'V_r100raw', 'V_K7r100', 'V_K9r100', 'V_bjr100', 'A_E2split', 'A_E1k', 'A_E1T3', 'A_E1kT3', 'A_E1k3', 'A_r42', 'A_r45', 'A_rawBg', 'A_T3all', 'A_r100', 'A_r100raw']
allnames = ['s_C2', 's_C3', 's_C4', 's_C5', 's_Geta', 's_GD', 's_GD_om', 's_GD_nom', 's_GC', 's_GC_om', 's_GC_nom', 's_T3', 'k_E1', 'k_E2']
with open(f'{TAB}/flex_variants.txt', 'w') as o:
    o.write('# flex_variants.txt -- plots/all_report/flex_tables.py: every scale (value +- stat (+) MC) in every fit; tilts and c_j in flex_results.txt\n')
    o.write('# V_ = variants of the nominal model, A_ = of the alternative; layouts and texts below the table\n')
    o.write(f"{'fit':11s} {'ndf':>4s} {'Pearson':>8s} {'hybGLS':>8s} " + ' '.join(f'{n:>15s}' for n in allnames) + '\n')
    for t in order:
        f = fits.get(t)
        if not f or 'values' not in f:
            continue
        cells = []
        for n in allnames:
            x = val(f, n)
            cells.append(f'{x[0]:6.3f}+-{x[3]:6.3f}' if x else f"{'-':>15s}")
        o.write(f"{t:11s} {int(f.get('ndf', 0)):4d} {f.get('chi2P', 0):8.1f} {f.get('chi2GLShyb', 0):8.1f} " + ' '.join(f'{c:>15s}' for c in cells) + '\n')
    o.write('\n')
    for t in order:
        if t in fits and 'text' in fits[t]:
            o.write(f"{t:11s} {fits[t]['text']}\n")

with open(f'{TAB}/flex_ratios.txt', 'w') as o:
    o.write('# flex_ratios.txt -- plots/all_report/flex_tables.py: ratios of scales (Part III-like parameters) and content-weighted sums over the fit bins\n')
    for l in ex['ratio'] + ex['sumscale']:
        o.write(l + '\n')
    o.write('\n# Part III equivalents (flex_fit.C, from the reproduced Part III F-main): value stat mc total\n')
    for k, v in ex['p3equiv_err'].items():
        o.write(f'p3equiv {k:6s} {v[0]:.4f} stat {v[1]:.4f} mc {v[2]:.4f} tot {v[3]:.4f}\n')
print('flex_tables: written', NOM, ALT)
for r in rows:
    print(r['n'], round(r['v'], 3), round(r['et'], 3), r['cls'], r['deg'], r['b0']['verdict'] if r['b0'] else '', r['corr'])

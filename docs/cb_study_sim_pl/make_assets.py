#!/usr/bin/env python3
# make_assets.py -- figures (figs/), equations (eq/) and table numbers (data/deck_tables.json) of the deck.
# Input: data/deck_numbers.json from export_deck_data.C. Style: SLIDES_STYLE.md section 7, Polish labels and decimal commas.
import json, os, subprocess, tempfile, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'data/deck_numbers.json')))
x_ = json.load(open(os.path.join(HERE, 'data/deck_extra.json')))
X = lambda k: np.array(x_[k], dtype=float)
INK, INK2, GRID = '#111827', '#4B5563', '#E5E7EB'
C1, C2, C3, C4, C5 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#6250d6'
RED = '#D00000'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 15, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
    'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8, 'legend.frameon': False,
    'savefig.dpi': 200, 'savefig.bbox': 'tight'})
comma = FuncFormatter(lambda v, p: ('%g' % v).replace('.', ','))
ed = np.array(d['edges']); lo, hi = ed[:-1], ed[1:]; wid = hi - lo; cen = 0.5*(lo + hi)
A = lambda k: np.array(d[k], dtype=float)
MLAB = 'M$_{ee}$ [GeV/c$^2$]'
os.makedirs(os.path.join(HERE, 'figs'), exist_ok=True); os.makedirs(os.path.join(HERE, 'eq'), exist_ok=True)
def save(fig, name): fig.savefig(os.path.join(HERE, 'figs', name + '.png')); plt.close(fig)
def steps(ax, y, **kw): ax.stairs(y, ed, **kw, baseline=None)
def ax_mass(ax, xmax=0.8, log=False):
    ax.set_xlim(0, xmax); ax.xaxis.set_major_formatter(comma)
    if log: ax.set_yscale('log')
    else: ax.yaxis.set_major_formatter(comma)
def cat(v, oa, p, cats):
    s = np.zeros(len(lo))
    for c in cats: s += A(f'h_all_{v}_oa{oa}_{p}_cat{c}')
    return s
Sg = lambda v: A(f'Sgen_all_{v}_oa0'); N = lambda v: A(f'Npm_all_{v}_oa0'); B = lambda v: A(f'Btrue_all_{v}_oa0'); CB = lambda v: A(f'CBhat_all_{v}_oa0'); CBe = lambda v: A(f'CBhat_all_{v}_oa0_err'); Npp = lambda v: A(f'Npp_all_{v}_oa0'); Nmm = lambda v: A(f'Nmm_all_{v}_oa0')
coarse = [0, .005, .01, .02, .03, .05, .08, .11, .14, .16, .18, .20, .22, .24, .26, .28, .30, .34, .38, .42, .46, .50, .56, .60, .70, .80]
def coarse_sum(y):
    out = []
    for a, b in zip(coarse[:-1], coarse[1:]):
        m = (lo >= a - 1e-9) & (hi <= b + 1e-9); out.append(y[m].sum())
    return np.array(out)
cx = 0.5*(np.array(coarse[:-1]) + np.array(coarse[1:])); cw = 0.5*(np.array(coarse[1:]) - np.array(coarse[:-1]))
def rat(a, b, ea2, eb2, minb=30):
    with np.errstate(divide='ignore', invalid='ignore'):
        r = a/b; e = r*np.sqrt(ea2/a**2 + eb2/b**2)
    m = (b >= minb) & (a > 0)
    return r, e, m

# ---- channels: mass spectra of genuine pairs
fig, ax = plt.subplots(figsize=(7.6, 5.5))
chs = [([0], 'Dalitz π⁰', C1, 'o'), ([1], 'Dalitz η', C5, 's'), ([2], 'konwersja γ z π⁰', C4, '^'), ([3, 4], 'konwersja γ z η i inne', C3, 'v')]
for cats, lab, c, m in chs:
    y = cat('all', 0, 'pm', cats)
    ok = y > 0
    ax.errorbar(cen[ok], y[ok]/wid[ok], yerr=np.sqrt(y[ok])/wid[ok], fmt=m + '-', ms=5, color=c, lw=0.8, elinewidth=1.0, label=lab)
ax_mass(ax, 0.6, True); ax.set_ylim(1e2, 5e9); ax.set_xlabel(MLAB); ax.set_ylabel('pary prawdziwe / (GeV/c$^2$)'); ax.legend(loc='upper right', fontsize=14)
save(fig, 'channels')

# ---- PID: beta_new vs q*p with polygons, e- mirrored
nx, ny = map(int, d['pid_nx_ny']); z = A('pid_h2').reshape(ny, nx); xe = np.linspace(*d['pid_x'], nx + 1); ye = np.linspace(*d['pid_y'], ny + 1)
fig, ax = plt.subplots(figsize=(8.0, 5.3)); ax.grid(False)
m = ax.pcolormesh(xe, ye, np.ma.masked_where(z <= 0, z), norm=LogNorm(vmin=1, vmax=z.max()), cmap='jet', shading='flat', rasterized=True)
cb_ = fig.colorbar(m, ax=ax, pad=0.015); cb_.set_label('ślady / przedział'); cb_.outline.set_visible(False)
styles = {'1sig': (':', 1.8, '1σ'), '3sig': ('--', 1.8, '3σ'), '5sig': ('-', 2.4, '5σ (zastosowane)')}
for s, (ls, lw, lab) in styles.items():
    for q, sign in (('ep', 1), ('em', -1)):
        x = X(f'v6_{q}_{s}_x')*sign; y = X(f'v6_{q}_{s}_y')
        ax.plot(x, y, ls=ls, lw=lw, color='k', label=lab if q == 'ep' else None)
ax.axvline(0, color=INK2, lw=0.8)
ax.set_xlim(-1200, 1200); ax.set_ylim(0.92, 1.08); ax.set_xlabel('q · p [MeV/c]'); ax.set_ylabel('β$_{new}$'); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: ('%g' % v).replace('-', '−'))); ax.yaxis.set_major_formatter(comma); ax.grid(False)
ax.text(0.03, 0.05, 'e⁻', transform=ax.transAxes, fontsize=18, color=INK); ax.text(0.95, 0.05, 'e⁺', transform=ax.transAxes, fontsize=18, color=INK, ha='right')
ax.legend(loc='upper center', ncol=3, fontsize=13, handlelength=2.6, bbox_to_anchor=(0.5, 1.10))
save(fig, 'pid')

# ---- main spectrum with S/CB strip
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.4, 6.0), sharex=True, gridspec_kw={'height_ratios': [3.2, 1], 'hspace': 0.07})
n, s, b, cb = N('all'), Sg('all'), B('all'), CB('all')
s_true_cb = n - b; s_form = n - cb; eform = np.sqrt(n + CBe('all')**2)
a1.errorbar(cen, n/wid, yerr=np.sqrt(n)/wid, fmt='o', ms=4.5, color=INK, lw=1, label='wszystkie pary e⁺e⁻, N$_{+-}$', zorder=3)
a1.stairs(s/wid, ed, color=C1, lw=2.2, label='sygnał prawdziwy (z historii śladów)', zorder=2, baseline=None)
ok = s_true_cb > 0
a1.plot(cen[ok], s_true_cb[ok]/wid[ok], 'o', mfc='none', mec=C1, ms=8, mew=1.2, label='wszystkie − tło prawdziwe', zorder=4)
ok = s_form > 0
a1.errorbar(cen[ok], s_form[ok]/wid[ok], yerr=eform[ok]/wid[ok], fmt='^', ms=5, color=C5, lw=1, label='wszystkie − tło ze wzoru 2√(N$_{++}$N$_{--}$)', zorder=5)
ax_mass(a1, 0.8, True); a1.set_ylim(1e2, 5e9); a1.set_ylabel('pary / (GeV/c$^2$)'); h_, l_ = a1.get_legend_handles_labels(); o_ = [2, 0, 1, 3]; a1.legend([h_[i] for i in o_], [l_[i] for i in o_], loc='upper right', fontsize=13)
with np.errstate(divide='ignore', invalid='ignore'):
    r1 = s/b; r2 = s_form/cb
a2.stairs(np.where(b > 0, r1, 0), ed, color=C1, fill=True, alpha=0.30, label='sygnał / tło, prawdziwe'); a2.stairs(np.where(b > 0, r1, 0), ed, color=C1, lw=1.5, baseline=None)
ok = (cb > 0) & (s_form > 0); a2.plot(cen[ok], r2[ok], '^', ms=4.5, color=C5, label='sygnał / tło, ze wzoru')
a2.axhline(1, color=INK2, lw=1, ls='--'); ax_mass(a2, 0.8, True); a2.set_ylim(0.02, 1000); a2.set_xlabel(MLAB); a2.set_ylabel('sygnał / tło', fontsize=13); a2.legend(loc='upper right', fontsize=12, ncol=2)
a1.text(0.03, 0.06, 'isBest ≥ 0 (wszystkie kombinacje)', transform=a1.transAxes, fontsize=15, fontweight='bold', color=INK)
save(fig, 'spectrum')

# ---- the same spectrum and ratio on a linear scale, above 0.1 GeV/c2
fig, (a1, a2) = plt.subplots(2, 1, figsize=(10.4, 6.0), sharex=True, gridspec_kw={'height_ratios': [3.2, 1.25], 'hspace': 0.07})
sc = 1e6
a1.errorbar(cen, n/wid/sc, yerr=np.sqrt(n)/wid/sc, fmt='o', ms=4.5, color=INK, lw=1, label='wszystkie pary e⁺e⁻, N$_{+-}$', zorder=3)
a1.stairs(s/wid/sc, ed, color=C1, lw=2.2, label='sygnał prawdziwy (z historii śladów)', zorder=2, baseline=None)
a1.plot(cen, s_true_cb/wid/sc, 'o', mfc='none', mec=C1, ms=8, mew=1.2, label='wszystkie − tło prawdziwe', zorder=4)
a1.errorbar(cen, s_form/wid/sc, yerr=eform/wid/sc, fmt='^', ms=5, color=C5, lw=1, label='wszystkie − tło ze wzoru 2√(N$_{++}$N$_{--}$)', zorder=5)
a1.stairs(b/wid/sc, ed, color=RED, lw=1.8, label='tło prawdziwe', zorder=2, baseline=None)
a1.plot(cen, cb/wid/sc, 's', ms=4.5, color=RED, label='tło ze wzoru', zorder=5)
a1.set_xlim(0.10, 0.80); a1.set_ylim(0, 1.08*(n[cen > 0.10]/wid[cen > 0.10]).max()/sc); a1.xaxis.set_major_formatter(comma); a1.yaxis.set_major_formatter(comma); a1.set_ylabel('pary / (GeV/c$^2$)  ×10$^6$'); a1.legend(loc='upper right', fontsize=13)
a1.text(0.14, 0.90, 'isBest ≥ 0', transform=a1.transAxes, fontsize=15, fontweight='bold', color=INK)
with np.errstate(divide='ignore', invalid='ignore'):
    r1 = s/b; r2 = s_form/cb
a2.stairs(np.where(b > 0, r1, 0), ed, color=C1, fill=True, alpha=0.30, label='sygnał / tło, prawdziwe'); a2.stairs(np.where(b > 0, r1, 0), ed, color=C1, lw=1.5, baseline=None)
ok = cb > 0; a2.plot(cen[ok], r2[ok], '^', ms=4.5, color=C5, label='sygnał / tło, ze wzoru')
a2.axhline(1, color=INK2, lw=1, ls='--'); a2.set_xlim(0.10, 0.80); a2.set_ylim(0, 7); a2.xaxis.set_major_formatter(comma); a2.yaxis.set_major_formatter(comma); a2.set_xlabel(MLAB); a2.set_ylabel('sygnał / tło', fontsize=13); a2.legend(loc='upper right', fontsize=12, ncol=2)
save(fig, 'spectrum_lin')

# ---- ratio CB formula / true, isBest >= 0
def ratio_series(v):
    a = coarse_sum(CB(v)); e2 = coarse_sum(CBe(v)**2); t = coarse_sum(B(v)); return rat(a, t, e2, t)
fig, ax = plt.subplots(figsize=(7.2, 5.0))
r, e, m = ratio_series('all'); ax.errorbar(cx[m], r[m], xerr=cw[m], yerr=e[m], fmt='o', ms=6, color=RED, lw=1.2, label='isBest ≥ 0')
ax.axhline(1, color=INK2, lw=1, ls='--'); ax_mass(ax, 0.8); ax.set_ylim(0.75, 1.45); ax.set_xlabel(MLAB); ax.set_ylabel('tło ze wzoru / tło prawdziwe'); ax.legend(loc='upper right', fontsize=14)
save(fig, 'ratio_all')

# ---- like-sign sources
fig, axs = plt.subplots(1, 2, figsize=(13.2, 5.15), sharey=True, gridspec_kw={'wspace': 0.04})
src = [([6], 'Dalitz + Dalitz', C1), ([7], 'Dalitz + konwersja', C2), ([8], 'konwersja + konwersja', C5), ([9, 10], 'leptony wtórne i inne', C3), ([11, 12, 13], 'fałszywy ślad', INK2)]
for ax, p, t in zip(axs, ('pp', 'mm'), ('e⁺e⁺', 'e⁻e⁻')):
    for cats, lab, c in src:
        y = cat('all', 0, p, cats); ok = y > 0; ax.errorbar(cen[ok], y[ok]/wid[ok], yerr=np.sqrt(y[ok])/wid[ok], fmt='o', ms=4.5, color=c, lw=0.9, label=lab)
    ax_mass(ax, 0.7, True); ax.set_ylim(5e1, 5e6); ax.set_xlabel(MLAB); ax.set_title(t, loc='left', fontsize=17)
axs[1].set_xticks([.1, .2, .3, .4, .5, .6, .7]); axs[1].tick_params(labelleft=False)
axs[0].set_ylabel('pary / (GeV/c$^2$)'); h_, l_ = axs[1].get_legend_handles_labels(); fig.legend(h_, l_, loc='lower center', ncol=5, fontsize=14, bbox_to_anchor=(0.5, 0.965), handletextpad=0.1, columnspacing=1.4)
save(fig, 'like_sign')

# ---- fraction kept by isBest==1
def kept(a, b):
    A0 = coarse_sum(a); B0 = coarse_sum(b)
    with np.errstate(divide='ignore', invalid='ignore'):
        f = B0/A0; e = np.sqrt(np.maximum(f*(1 - f), 1/np.maximum(A0, 1))/np.maximum(A0, 1))
    return f, e, A0 >= 30
def kept_plot(series, fname, figsize, xmax=0.7, ylo=None, leg=13):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=figsize, sharex=True, gridspec_kw={'height_ratios': [2.5, 1.2], 'hspace': 0.07})
    for (y_all, y_best, lab, c, mk) in series:
        A0 = coarse_sum(y_all); B0 = coarse_sum(y_best); w2 = 2*cw; ok = A0 > 0
        a1.stairs(A0/w2, coarse, color=c, lw=1.8, label=lab + ', isBest ≥ 0', baseline=None)
        a1.errorbar(cx[ok], B0[ok]/w2[ok], xerr=cw[ok], yerr=np.sqrt(B0[ok])/w2[ok], fmt=mk, ms=5, color=c, lw=1, label=lab + ', isBest == 1')
        f, e, m = kept(y_all, y_best); a2.errorbar(cx[m], f[m], xerr=cw[m], yerr=e[m], fmt=mk + '-', ms=5, color=c, lw=0.9, elinewidth=1.1, label=lab)
    ax_mass(a1, xmax, True); a1.set_ylabel('pary / (GeV/c$^2$)'); a1.legend(loc='upper right', fontsize=leg, ncol=1)
    ax_mass(a2, xmax); a2.set_ylim(*ylo); a2.set_xlabel(MLAB); a2.set_ylabel('isBest == 1 /\nisBest ≥ 0', fontsize=12.5); a2.axhline(1, color=INK2, lw=1, ls='--')
    save(fig, fname)
kept_plot([(Sg('all'), Sg('best'), 'pary prawdziwe', C1, 'o'), (B('all'), B('best'), 'tło prawdziwe', RED, 's')], 'kept_pm', (7.8, 5.6), ylo=(0.6, 1.04))
kept_plot([(Npp('all'), Npp('best'), 'e⁺e⁺', C2, '^'), (Nmm('all'), Nmm('best'), 'e⁻e⁻', C5, 'v')], 'kept_ls', (6.0, 3.9), ylo=(0.85, 1.02), leg=11)

# ---- ratio formula / truth, both variants
fig, ax = plt.subplots(figsize=(7.2, 5.0))
for v, lab, c, mk in (('all', 'isBest ≥ 0', INK, 'o'), ('best', 'isBest == 1', RED, 's')):
    r, e, m = ratio_series(v); ax.errorbar(cx[m], r[m], xerr=cw[m], yerr=e[m], fmt=mk, ms=6, color=c, lw=1.2, label=lab)
ax.axhline(1, color=INK2, lw=1, ls='--'); ax_mass(ax, 0.8); ax.set_ylim(0.8, 1.7); ax.set_xlabel(MLAB); ax.set_ylabel('tło ze wzoru / tło prawdziwe'); ax.legend(loc='upper right', fontsize=14)
save(fig, 'ratio_best')

# ---- signal: isBest >= 0 vs isBest == 1
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.0, 5.9), sharex=True, gridspec_kw={'height_ratios': [2.6, 1], 'hspace': 0.07})
sa, sb = Sg('all'), Sg('best'); fa = N('all') - CB('all'); fb = N('best') - CB('best'); ea = np.sqrt(N('all') + CBe('all')**2); eb = np.sqrt(N('best') + CBe('best')**2)
a1.stairs(sa/wid, ed, color=C1, lw=2.2, label='sygnał prawdziwy, isBest ≥ 0', baseline=None)
ok = fa > 0; a1.errorbar(cen[ok], fa[ok]/wid[ok], yerr=ea[ok]/wid[ok], fmt='o', ms=4.5, color=INK, lw=1, label='wszystkie − wzór, isBest ≥ 0')
ok = fb > 0; a1.errorbar(cen[ok], fb[ok]/wid[ok], yerr=eb[ok]/wid[ok], fmt='s', ms=4.5, color=RED, lw=1, label='wszystkie − wzór, isBest == 1')
ax_mass(a1, 0.56, True); a1.set_ylim(3e3, 3e9); a1.set_ylabel('sygnał / (GeV/c$^2$)'); a1.legend(loc='upper right', fontsize=13)
for f, st, e_, c, mk, lab in ((fa, sa, ea, INK, 'o', 'isBest ≥ 0'), (fb, sb, eb, RED, 's', 'isBest == 1')):
    m = (st > 30) & (cen < 0.56); r_ = f[m]/st[m]; a2.errorbar(cen[m], r_, yerr=e_[m]/st[m], fmt=mk, ms=4.5, color=c, lw=1, label=lab)
a2.axhline(1, color=INK2, lw=1, ls='--'); ax_mass(a2, 0.56); a2.set_ylim(0.6, 1.2); a2.set_xlabel(MLAB); a2.set_ylabel('ze wzoru / prawdziwy', fontsize=12)
save(fig, 'signal_best')

# ---- all combinations against one pair per event: background and like-sign pairs, same Y range in both panels
fig, axs = plt.subplots(2, 2, figsize=(11.0, 6.3), sharex=True, gridspec_kw={'height_ratios': [3.0, 1.15], 'hspace': 0.07, 'wspace': 0.05})
for iax, v, t in zip((0, 1), ('all', 'best'), ('isBest ≥ 0', 'isBest == 1')):
    ax, axr = axs[0][iax], axs[1][iax]
    ax.stairs(B(v)/wid, ed, color=INK, lw=2.0, label='tło prawdziwe e⁺e⁻', zorder=2, baseline=None)
    y = CB(v); ok = y > 0
    ax.errorbar(cen[ok], y[ok]/wid[ok], yerr=CBe(v)[ok]/wid[ok], fmt='s', ms=4.5, color=RED, lw=1, label='tło ze wzoru 2√(N$_{++}$N$_{--}$)', zorder=4)
    for yy, mk, c, lab in ((Npp(v), '^', C2, 'pary e⁺e⁺'), (Nmm(v), 'v', C5, 'pary e⁻e⁻')):
        ok = yy > 0; ax.errorbar(cen[ok], yy[ok]/wid[ok], yerr=np.sqrt(yy[ok])/wid[ok], fmt=mk, ms=3.8, color=c, lw=0.8, label=lab, zorder=3)
    ax_mass(ax, 0.8, True); ax.set_ylim(3e2, 2e7); ax.set_title(t, loc='left', fontsize=17)
    t_ = coarse_sum(B(v)); a_ = coarse_sum(CB(v)); e2_ = coarse_sum(CBe(v)**2); r_, e_, m_ = rat(t_, a_, t_, e2_)
    axr.errorbar(cx[m_], r_[m_], xerr=cw[m_], yerr=e_[m_], fmt='o', ms=4.5, color=INK, lw=1); axr.axhline(1, color=INK2, lw=1, ls='--'); ax_mass(axr, 0.8); axr.set_ylim(0.55, 1.35); axr.set_xlabel(MLAB)
    if iax == 1: ax.tick_params(labelleft=False); axr.tick_params(labelleft=False); axr.set_xticks([.2, .4, .6, .8])
axs[0][0].set_ylabel('pary / (GeV/c$^2$)'); axs[1][0].set_ylabel('prawdziwe /\nze wzoru', fontsize=13)
h_, l_ = axs[0][1].get_legend_handles_labels(); fig.legend(h_, l_, loc='lower center', ncol=4, fontsize=12.5, bbox_to_anchor=(0.5, 0.955), handletextpad=0.6, columnspacing=1.2)
save(fig, 'cb_both')

# ---- k factor: simulation truth against event mixing with two pools; data; final curve
kc = np.array([0, .02, .05, .08, .11, .14, .16, .18, .20, .24, .28, .34, .42, .56, .80])
kx = 0.5*(kc[:-1] + kc[1:]); kw = 0.5*(kc[1:] - kc[:-1])
def inw(a, b): return (lo >= a - 1e-9) & (hi <= b + 1e-9)
def k_of_mix(pm, pp, mm, edges=kc):
    v, e = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = inw(a, b); g = (pp[m] > 0) & (mm[m] > 0); G = (2*np.sqrt(pp[m][g]*mm[m][g])).sum(); P = pm[m].sum()
        k_ = P/G if G > 0 else np.nan; v.append(k_); e.append(k_*np.sqrt(1/max(P, 1) + 0.25/max(pp[m].sum(), 1) + 0.25/max(mm[m].sum(), 1)))
    return np.array(v), np.array(e)
def k_true(st, edges=kc):
    bb, cc, e2 = X(f'sim_Btrue_{st}'), X(f'sim_CBhat_{st}'), X(f'sim_CBhat_{st}_err')**2; v, e = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = inw(a, b); B_, C_ = bb[m].sum(), cc[m].sum(); k_ = B_/C_; v.append(k_); e.append(k_*np.sqrt(1/B_ + e2[m].sum()/C_**2))
    return np.array(v), np.array(e)
pools = lambda kind, pool, st='PT3': (X(f'kpool_sim_{pool}_{st}_pm'), X(f'kpool_sim_{pool}_{st}_pp'), X(f'kpool_sim_{pool}_{st}_mm')) if kind == 'sim' else (X(f'kpool_data_{pool}_pm'), X(f'kpool_data_{pool}_pp'), X(f'kpool_data_{pool}_mm'))
def k_final(edges=kc):
    kt, ekt = k_true('PT3', edges); ksa, _ = k_of_mix(*pools('sim', 'all'), edges); ksn, _ = k_of_mix(*pools('sim', 'nge3'), edges); kda, _ = k_of_mix(*pools('data', 'all'), edges); kdn, _ = k_of_mix(*pools('data', 'nge3'), edges)
    k1 = kda*kt/ksa; k2 = kdn*kt/ksn; kf = 0.5*(k1 + k2); es = 0.5*np.abs(k1 - k2); et = kf*ekt/kt
    return dict(kt=kt, ekt=ekt, ksa=ksa, ksn=ksn, kda=kda, kdn=kdn, k1=k1, k2=k2, k=kf, syst=es, stat=et, tot=np.hypot(es, et))
KF = k_final()
# central value: classes n = 2 (pool c2) and n >= 3 (pool nge3), weighted with the like-sign population of the data; sim per class from cb_study_sim, data per class from the emulation
from scipy.optimize import curve_fit
two_ = lambda pre: (X(pre % 2), X(pre % 3) + X(pre % 4) + X(pre % 5))
gm0 = lambda p_, m_: np.where((p_ > 0) & (m_ > 0), 2*np.sqrt(np.maximum(p_, 0)*np.maximum(m_, 0)), 0.)
Sp2, Sp3 = two_('sc_pp_c%d'); Sm2, Sm3 = two_('sc_mm_c%d'); B2_, B3_ = two_('scB_c%d'); Dp2, Dp3 = two_('emu_pp_c%d'); Dm2, Dm3 = two_('emu_mm_c%d')
Fs2, Fs3, Fd2, Fd3, Fdt_ = gm0(Sp2, Sm2), gm0(Sp3, Sm3), gm0(Dp2, Dm2), gm0(Dp3, Dm3), gm0(X('data_all_pp'), X('data_all_mm'))
kcf = kc[kc >= 0.05 - 1e-9]; xcf = 0.5*(kcf[:-1] + kcf[1:]); wcf = 0.5*(kcf[1:] - kcf[:-1])
def k_two(edges=kcf):
    ks2_, _ = k_of_mix(X('kclass_sim_c2_pm'), X('kclass_sim_c2_pp'), X('kclass_sim_c2_mm'), edges); kd2_, _ = k_of_mix(X('kclass_data_c2_pm'), X('kclass_data_c2_pp'), X('kclass_data_c2_mm'), edges)
    ks3_, _ = k_of_mix(*pools('sim', 'nge3'), edges); kd3_, _ = k_of_mix(*pools('data', 'nge3'), edges); v, e_ = [], []
    for i, (a_, b_) in enumerate(zip(edges[:-1], edges[1:])):
        m = inw(a_, b_); k2_ = B2_[m].sum()/Fs2[m].sum(); k3_ = B3_[m].sum()/Fs3[m].sum(); f2 = Fd2[m].sum()*kd2_[i]/ks2_[i]; f3 = Fd3[m].sum()*kd3_[i]/ks3_[i]; T_ = Fdt_[m].sum()
        v.append((f2*k2_ + f3*k3_)/T_); e2_ = k2_*np.sqrt(1/B2_[m].sum() + 0.25/max(Sp2[m].sum(), 1) + 0.25/max(Sm2[m].sum(), 1)); e3_ = k3_*np.sqrt(1/B3_[m].sum() + 0.25/max(Sp3[m].sum(), 1) + 0.25/max(Sm3[m].sum(), 1)); e_.append(np.sqrt((f2*e2_)**2 + (f3*e3_)**2)/T_)
    return np.array(v), np.array(e_)
KCsimv, _ = k_two()   # with the simulation ratio: a sensitivity variant only, the simulation does not model the combinatorial background of the data
def k_data_only(edges=kcf):   # central value: event mixing in the data, classes n = 2 and n >= 3 weighted with their share of the like-sign formula in the data
    kd2_, e2_ = k_of_mix(X('kclass_data_c2_pm'), X('kclass_data_c2_pp'), X('kclass_data_c2_mm'), edges); kd3_, e3_ = k_of_mix(*pools('data', 'nge3'), edges); v, e_ = [], []
    for i, (a_, b_) in enumerate(zip(edges[:-1], edges[1:])):
        m = inw(a_, b_); f2 = Fd2[m].sum(); f3 = Fd3[m].sum(); T_ = Fdt_[m].sum(); v.append((f2*kd2_[i] + f3*kd3_[i])/T_); e_.append(np.sqrt((f2*e2_[i])**2 + (f3*e3_[i])**2)/T_)
    return np.array(v), np.array(e_)
KC, KCs = k_data_only()
KFc = k_final(kcf); kdav = 0.5*(KFc['kda'] + KFc['kdn']); KCsys = np.sqrt((KC - kdav)**2 + (0.5*(KFc['kda'] - KFc['kdn']))**2); KCtot = np.sqrt(KCs**2 + KCsys**2)
def agauss(M, a, mu, sR, c, sL=0.015):
    sg = np.where(M < mu, sL, sR); return c + a*np.exp(-0.5*((M - mu)/sg)**2)
def binavg_(p_): return np.array([agauss(np.linspace(a_, b_, 41), *p_).mean() for a_, b_ in zip(kcf[:-1], kcf[1:])])
def fit_k(y, err): return curve_fit(lambda _x, a, mu, sR, c: binavg_((a, mu, sR, c)), xcf, y, p0=(.18, .135, .07, 1.01), sigma=err, absolute_sigma=True, bounds=([0, .1, .005, .9], [1, .2, .3, 1.1]))
cen58 = 0.5*(lo + hi); mask58 = cen58 >= 0.05 - 1e-9
grid58 = np.array([np.linspace(l_, h_, 41) for l_, h_ in zip(lo, hi)])
kfun = lambda q: np.where(mask58, agauss(grid58, *q).mean(axis=1), 1.)
pfit, Cfit = fit_k(KC, KCtot); chi2fit = float((((KC - binavg_(pfit))/KCtot)**2).sum())
rng_ = np.random.default_rng(7); toys_ = []
for _ in range(300):
    try: toys_.append(kfun(fit_k(KC + rng_.normal(0, KCs), KCtot)[0]))
    except Exception: pass
toys_ = np.array(toys_); pu_, _ = fit_k(KC + KCsys, KCtot); pd_, _ = fit_k(KC - KCsys, KCtot); kup, kdn = kfun(pu_), kfun(pd_)
bstat = toys_.std(axis=0); bsys = 0.5*np.abs(kup - kdn); band58 = np.hypot(bstat, bsys); kfit58 = kfun(pfit)
fig, ax = plt.subplots(figsize=(8.2, 5.2))
ax.errorbar(kx, KF['kt'], xerr=kw, yerr=KF['ekt'], fmt='o', ms=6.5, color=INK, lw=1.2, label='k$_{true}$: tło prawdziwe / wzór', zorder=4)
for pool, lab, c, mk in (('all', 'k$_{mix}$, pula: wszystkie zdarzenia', C1, 's'), ('nge3', 'k$_{mix}$, pula: zdarzenia z n ≥ 3', C2, '^')):
    v_, e_ = k_of_mix(*pools('sim', pool)); ax.errorbar(kx, v_, xerr=kw, yerr=e_, fmt=mk + '-', ms=5.5, color=c, lw=0.9, elinewidth=1.1, label=lab)
ax.axhline(1, color=INK2, lw=1, ls='--'); ax.axvspan(0, 0.05, color=GRID, alpha=0.6, lw=0); ax.text(0.025, 1.26, 'nie akceptancja', fontsize=11.5, color=INK2, ha='center', rotation=90, va='top')
ax_mass(ax, 0.8); ax.set_ylim(0.75, 1.3); ax.set_xlabel(MLAB); ax.set_ylabel('k'); ax.legend(loc='upper right', fontsize=13)
save(fig, 'k_sim')
fig, ax = plt.subplots(figsize=(8.2, 5.2))
order_ = np.argsort(cen58); xx_ = cen58[(cen58 >= 0.05) & (cen58 <= 0.8)]; mm_ = np.isin(cen58, xx_)
ax.fill_between(cen58[mm_], (kfit58 - band58)[mm_], (kfit58 + band58)[mm_], color=RED, alpha=0.20, lw=0, label='pasmo niepewności dopasowania')
ax.plot(cen58[mm_], kfit58[mm_], color=RED, lw=2.4, label='k: gładkie dopasowanie')
ax.errorbar(xcf, KC, xerr=wcf, yerr=KCtot, fmt='o', ms=6.5, color=INK, lw=1.2, label='k z mieszania w danych: klasy n = 2 i n ≥ 3', zorder=5)
ax.plot(kx[kx > 0.05], KF['kda'][kx > 0.05], 's-', ms=5, color=C1, lw=0.9, label='k$_{mix}^{dane}$, pula: wszystkie zdarzenia')
ax.plot(kx[kx > 0.05], KF['kdn'][kx > 0.05], '^-', ms=5, color=C2, lw=0.9, label='k$_{mix}^{dane}$, pula: n ≥ 3')
ax.axhline(1, color=INK2, lw=1, ls='--'); ax_mass(ax, 0.8); ax.set_xlim(0.04, 0.8); ax.set_ylim(0.94, 1.16); ax.set_xlabel(MLAB); ax.set_ylabel('k'); ax.legend(loc='upper right', fontsize=12.5)
save(fig, 'k_data')
k58 = kfit58.copy()
with open(os.path.join(HERE, '../../plots/output/k_factor_data.txt'), 'w') as f_:   # columns: bin low, bin high, k, band; used by plots/mass_ee_signal_vs_sim_channels.C
    for l_, h_, k_, b_ in zip(lo, hi, k58, band58): f_.write(f'{l_:.4f} {h_:.4f} {k_:.5f} {b_:.5f}\n')

# ---- why the formula is low below 0.3 GeV/c2 (simulation, isBest >= 0)
fig, ax = plt.subplots(figsize=(7.2, 5.0))
r, e, m = ratio_series('all'); ax.errorbar(cx[m], r[m], xerr=cw[m], yerr=e[m], fmt='o', ms=6, color=INK, lw=1.2, label='wzór / tło prawdziwe', zorder=4)
for pool, lab, c, mk in (('all', '1 / k$_{mix}$, pula: wszystkie', C1, 's'), ('nge3', '1 / k$_{mix}$, pula: n ≥ 3', C2, '^')):
    v_, e_ = k_of_mix(*pools('sim', pool, 'all')); mm_ = kx >= 0.05; ax.plot(kx[mm_], 1/v_[mm_], mk + '-', ms=5.5, color=c, lw=1.0, label=lab)
ax.axhline(1, color=INK2, lw=1, ls='--'); ax_mass(ax, 0.8); ax.set_ylim(0.8, 1.3); ax.set_xlabel(MLAB); ax.set_ylabel('stosunek'); ax.legend(loc='upper right', fontsize=14)
save(fig, 'why_low')

# ---- why the choice is a coin toss: D of the pairs and who wins
de = x_['hD_edges']; De = np.linspace(de[0], de[1], int(de[2]) + 1)
fig, (a, b_) = plt.subplots(1, 2, figsize=(9.6, 6.6), gridspec_kw={'width_ratios': [1.25, 1], 'wspace': 0.32})
for k_, lab, c in (('gen', 'pary prawdziwe', C1), ('lep', 'tło z leptonów', RED), ('fake', 'tło z fałszywym śladem', INK2)):
    y = X(f'hD_{k_}_n3'); a.stairs(y/y.sum(), De, color=c, lw=2.2, label=f'{lab}', baseline=None)
a.set_xlim(-1.0, 1.0); a.set_ylim(0, 0.19); a.xaxis.set_major_formatter(FuncFormatter(lambda v, p: ('%g' % v).replace('.', ',').replace('-', '−'))); a.yaxis.set_major_formatter(comma); a.set_xlabel('D [ns]'); a.set_ylabel('ułamek par / 50 ps'); a.legend(loc='upper left', fontsize=13); a.set_title('Trzy leptony w zdarzeniu', loc='left', fontsize=15)
W_, T_ = X('choice_win'), X('choice_tot')
idx = [8 + 1, 8 + 2, 8 + 3]; labs = ['ten sam\ndetektor', 'prawdziwa\nw RPC', 'rywal\nw RPC']; cols = [INK2, C3, C2]
vals = [W_[i]/T_[i] for i in idx]; errs = [np.sqrt(v*(1 - v)/T_[i]) for v, i in zip(vals, idx)]
b_.bar(range(3), vals, yerr=errs, color=cols, width=0.62, capsize=4); b_.axhline(0.5, color=INK, lw=1.2, ls='--'); b_.set_xticks(range(3)); b_.set_xticklabels(labs, fontsize=12.5); b_.set_ylim(0, 0.85); b_.yaxis.set_major_formatter(comma); b_.set_ylabel('P(wybrana para prawdziwa)'); b_.grid(axis='x', visible=False); b_.set_title('Kto wygrywa', loc='left', fontsize=15)
for i, v in enumerate(vals): b_.text(i, v + 0.035, ('%.2f' % v).replace('.', ','), ha='center', fontsize=14, color=INK)
save(fig, 'choice')

# ---- experimental data (PT3, vertex cut, all from existing outputs)
dN, dNb, dPa, dPb, dMa, dMb = X('data_all_pm'), X('data_best_pm'), X('data_all_pp'), X('data_best_pp'), X('data_all_mm'), X('data_best_mm')
gm = lambda p_, m_: np.where((p_ > 0) & (m_ > 0), 2*np.sqrt(np.maximum(p_, 0)*np.maximum(m_, 0)), 0.)
gme = lambda p_, m_: np.where((p_ > 0) & (m_ > 0), gm(p_, m_)*0.5*np.sqrt(1/np.maximum(p_, 1) + 1/np.maximum(m_, 1)), 0.)
dCa, dCb, eCa, eCb = gm(dPa, dMa), gm(dPb, dMb), gme(dPa, dMa), gme(dPb, dMb)
dSa, dSb, eSa, eSb = dN - dCa, dNb - dCb, np.sqrt(dN + eCa**2), np.sqrt(dNb + eCb**2)
fig, axs = plt.subplots(2, 2, figsize=(11.6, 6.0), sharex=True, gridspec_kw={'height_ratios': [3.0, 1.15], 'hspace': 0.07, 'wspace': 0.05})
for iax, (n_, c_, s_, ec_, es_, t) in enumerate(((dN, dCa, dSa, eCa, eSa, 'isBest ≥ 0'), (dNb, dCb, dSb, eCb, eSb, 'isBest == 1'))):
    ax, axr = axs[0][iax], axs[1][iax]
    ok = n_ > 0; ax.errorbar(cen[ok], n_[ok]/wid[ok], yerr=np.sqrt(n_[ok])/wid[ok], fmt='o', ms=4.3, color='k', lw=1, label='wszystkie pary e⁺e⁻')
    ok = c_ > 0; ax.errorbar(cen[ok], c_[ok]/wid[ok], yerr=ec_[ok]/wid[ok], fmt='o', ms=4.3, color=RED, lw=1, label='tło ze wzoru 2√(N$_{++}$N$_{--}$)')
    ok = s_ > 0; ax.errorbar(cen[ok], s_[ok]/wid[ok], yerr=es_[ok]/wid[ok], fmt='o', ms=4.3, color='#0000D0', lw=1, label='sygnał = wszystkie − tło')
    ax_mass(ax, 1.0, True); ax.set_ylim(1e3, 2e10); ax.set_title(t, loc='left', fontsize=17)
    with np.errstate(divide='ignore', invalid='ignore'): rr_ = s_/c_
    ok = c_ > 0; axr.plot(cen[ok], np.maximum(rr_[ok], 1e-3), 'o', ms=4, color='#0000D0'); axr.axhline(1, color=INK2, lw=1, ls='--'); ax_mass(axr, 1.0, True); axr.set_ylim(0.03, 300); axr.set_xlabel(MLAB)
    if iax == 1: ax.tick_params(labelleft=False); axr.tick_params(labelleft=False); axr.set_xticks([.2, .4, .6, .8, 1.0])
axs[0][0].set_ylabel('pary / (GeV/c$^2$)'); axs[1][0].set_ylabel('sygnał / tło', fontsize=13)
h_, l_ = axs[0][1].get_legend_handles_labels(); fig.legend(h_, l_, loc='lower center', ncol=3, fontsize=13.5, bbox_to_anchor=(0.5, 0.955), handletextpad=0.5, columnspacing=1.6)
save(fig, 'data_spectra')
dc = np.array([0, .01, .02, .03, .05, .08, .11, .14, .18, .22, .26, .30, .36, .42, .50, .60, .70, .80, .95, 1.40])
def dsum(y): return np.array([y[inw(a, b)].sum() for a, b in zip(dc[:-1], dc[1:])])
dx = 0.5*(dc[:-1] + dc[1:]); dw = 0.5*(dc[1:] - dc[:-1])
fig, ax = plt.subplots(figsize=(8.0, 5.0))
series_ = ((dN, dNb, 'e⁺e⁻', 'k', 'o'), (dPa, dPb, 'e⁺e⁺', C2, '^'), (dMa, dMb, 'e⁻e⁻', C5, 'v'))
fr_ = {}
for a_, b_, lab, c, mk in series_:
    A0, B0 = dsum(a_), dsum(b_); ok = A0 > 30; f_ = B0[ok]/A0[ok]; fr_[lab] = (dx[ok], f_); ax.errorbar(dx[ok], f_, xerr=dw[ok], yerr=np.sqrt(np.maximum(f_*(1 - f_), 1/A0[ok])/A0[ok]), fmt=mk, ms=5.5, color=c, lw=1.1, label=lab)
A0, B0 = dsum(dCa), dsum(dCb); ok = A0 > 30; ax.plot(dx[ok], B0[ok]/A0[ok], 's--', ms=4.5, color=RED, lw=1.2, label='tło ze wzoru: isBest == 1 / isBest ≥ 0')
ax_mass(ax, 1.0); ax.set_ylim(0.55, 1.03); ax.set_xlabel(MLAB); ax.set_ylabel('isBest == 1 / isBest ≥ 0'); ax.legend(loc='lower left', fontsize=13)
save(fig, 'data_kept')
def emu_cb():
    tot = np.zeros(len(lo))
    for c in range(2, 6):
        g_ = gm(X(f'emu_pp_c{c}'), X(f'emu_mm_c{c}')); cr, wn = X(f'emu_cross_c{c}'), X(f'emu_win_c{c}')
        with np.errstate(divide='ignore', invalid='ignore'): eps = np.where(cr > 0, wn/cr, 1. if c == 2 else wn.sum()/max(1., cr.sum()))
        tot += eps*g_
    return tot
dCe = emu_cb()
dSk = dN - k58*dCa; dSe = dNb - k58*dCe
k58s = np.ones(len(lo))
for a_, b_, v_ in zip(kcf[:-1], kcf[1:], KCsimv): k58s[inw(a_, b_)] = v_   # only to show how much a simulation-based ratio would change the signal
# display grid of this figure: 5 MeV/c2 below 0.14, 20 MeV/c2 above (5 MeV histograms of the same data)
fe_ = np.linspace(0, 1.4, 281); de_ = np.concatenate([np.arange(0, 0.14 - 1e-9, 0.005), np.arange(0.14, 0.80 + 1e-9, 0.02)])
reb_ = lambda y: np.array([y[(fe_[:-1] >= a_ - 1e-9) & (fe_[1:] <= b_ + 1e-9)].sum() for a_, b_ in zip(de_[:-1], de_[1:])])
Fx = lambda k: reb_(X(k)); dc_ = 0.5*(de_[:-1] + de_[1:]); dw_ = de_[1:] - de_[:-1]
fN, fP, fM, fNb, fPb, fMb = Fx('fine_all_pm'), Fx('fine_all_pp'), Fx('fine_all_mm'), Fx('fine_best_pm'), Fx('fine_best_pp'), Fx('fine_best_mm')
fC, fCb, feC = gm(fP, fM), gm(fPb, fMb), gme(fP, fM)
fCe = np.zeros(len(dc_))
for c in range(2, 6):
    g_ = gm(Fx(f'fine_pp_c{c}'), Fx(f'fine_mm_c{c}')); cr, wn = Fx(f'fine_cross_c{c}'), Fx(f'fine_win_c{c}')
    with np.errstate(divide='ignore', invalid='ignore'): eps = np.where(cr > 0, wn/cr, 1. if c == 2 else wn.sum()/max(1., cr.sum()))
    fCe += eps*g_
fk = np.array([agauss(np.linspace(a_, b_, 41), *pfit).mean() if a_ >= 0.05 - 1e-9 else 1. for a_, b_ in zip(de_[:-1], de_[1:])])
fb = np.where(dc_ >= 0.05, np.interp(dc_, cen58, band58), 0.)
fSk, fSb, fSe, fbS, feS = fN - fk*fC, fNb - fCb, fNb - fk*fCe, fb*fC, np.sqrt(fN + (fk*feC)**2)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.6, 5.9), sharex=True, gridspec_kw={'height_ratios': [2.6, 1], 'hspace': 0.07})
for y_, lab, c, mk in ((fSk, 'isBest ≥ 0:  N$_{+-}$ − k · 2√(N$_{++}$N$_{--}$)', 'k', 'o'), (fSb, 'isBest == 1:  N$_{+-}$ − 2√(N$_{++}$N$_{--}$)', RED, 's'), (fSe, 'isBest == 1:  N$_{+-}$ − k · Σ$_n$ ε$_n$ 2√(N$_{++,n}$N$_{--,n}$)', C3, '^')):
    ok = y_ > 0; a1.errorbar(dc_[ok], y_[ok]/dw_[ok], yerr=feS[ok]/dw_[ok], fmt=mk, ms=3.8, color=c, lw=0.9, label=lab)
ok = fSk > 0; a1.fill_between(dc_[ok], np.maximum(fSk - fbS, 1)[ok]/dw_[ok], (fSk + fbS)[ok]/dw_[ok], color='k', alpha=0.18, lw=0, label='niepewność k')
ax_mass(a1, 0.8, True); a1.set_ylim(2e4, 3e10); a1.set_ylabel('sygnał / (GeV/c$^2$)'); a1.legend(loc='upper right', fontsize=12.5)
ok_ = fSk > 300
a2.fill_between(dc_[ok_], 1 - (fbS/fSk)[ok_], 1 + (fbS/fSk)[ok_], color='k', alpha=0.18, lw=0)
a2.plot(dc_[ok_], fSb[ok_]/fSk[ok_], 's', ms=3.8, color=RED); a2.plot(dc_[ok_], fSe[ok_]/fSk[ok_], '^', ms=3.8, color=C3); a2.axhline(1, color=INK2, lw=1, ls='--')
ax_mass(a2, 0.8); a2.set_ylim(0.0, 1.5); a2.set_xlabel(MLAB); a2.set_ylabel('do czarnych', fontsize=12.5)
save(fig, 'data_signal')

# ---- k and eps_n in other observables (window 0.14-0.30 GeV/c2)
PV = [('mm', 'masa brakująca e⁺e⁻ [GeV/c$^2$]', 38, 1.6, 3.5), ('thcms', 'θ pary w CMS [°]', 18, 0., 180.), ('cosh', '|cos θ| e⁺ w układzie pary', 20, 0., 1.), ('pt', 'p$_T$ pary [GeV/c]', 25, 0., 1.0)]
def kx_(w_, v_, pool):
    pm, pp, mm = X(f'mix_w{w_}_{v_}_{pool}_pm'), X(f'mix_w{w_}_{v_}_{pool}_pp'), X(f'mix_w{w_}_{v_}_{pool}_mm')
    g_ = 2*np.sqrt(np.maximum(pp, 0)*np.maximum(mm, 0))
    with np.errstate(divide='ignore', invalid='ignore'): k_ = pm/g_
    return np.where((pm > 3000) & (g_ > 0), k_, np.nan), np.nansum(pm)/np.nansum(g_)
def ex_(w_, v_, c_):
    a_, b_ = X(f'eps_w{w_}_{v_}_{c_}_all'), X(f'eps_w{w_}_{v_}_{c_}_best')
    with np.errstate(divide='ignore', invalid='ignore'): e_ = b_/a_
    return np.where(a_ > 300, e_, np.nan), b_.sum()/a_.sum()
fig, axs = plt.subplots(2, 4, figsize=(14.2, 6.0), sharex='col', gridspec_kw={'hspace': 0.09, 'wspace': 0.2})
for j, (v_, lab, nb_, lo_, hi_) in enumerate(PV):
    e_ = np.linspace(lo_, hi_, nb_ + 1); c_ = 0.5*(e_[:-1] + e_[1:]); a0, a1_ = axs[0][j], axs[1][j]
    for pool, plab, col, mk in (('all', 'pula: wszystkie zdarzenia', C1, 'o'), ('nge3', 'pula: n ≥ 3', C2, '^')):
        y_, m_ = kx_(0, v_, pool); a0.plot(c_, y_, mk + '-', ms=4.5, lw=1, color=col, label=plab); a0.axhline(m_, color=col, lw=0.9, ls=':')
    a0.axhline(1, color=INK2, lw=1, ls='--'); a0.set_ylim(0.75, 1.55); a0.yaxis.set_major_formatter(comma)
    for c__, plab, col, mk in (('n3', 'n = 3', C5, 's'), ('n4', 'n = 4', C3, 'v')):
        y_, m_ = ex_(0, v_, c__); a1_.plot(c_, y_, mk + '-', ms=4.5, lw=1, color=col, label=plab); a1_.axhline(m_, color=col, lw=0.9, ls=':')
    a1_.set_ylim(0.0, 0.8); a1_.yaxis.set_major_formatter(comma); a1_.xaxis.set_major_formatter(comma); a1_.set_xlabel(lab, fontsize=14); a1_.set_xlim(lo_, hi_)
axs[0][0].set_ylabel('k$_{mix}$'); axs[1][0].set_ylabel('isBest == 1 / wszystkie', fontsize=13); axs[0][0].legend(loc='upper right', fontsize=11.5); axs[1][0].legend(loc='upper right', fontsize=11.5, ncol=2)
save(fig, 'projections')

# ---- equations
EQ = {'cb': r'\mathrm{CB}_{\text{wz\'or}} = 2\sqrt{N_{++}\,N_{--}}',
      's': r'S_{\text{wz\'or}} = N_{+-} - \mathrm{CB}_{\text{wz\'or}}',
      'rdec': r'R_{\mathrm{isBest}==1} = R_{\mathrm{isBest}\geq 0}\cdot\frac{f_{\mathrm{jedn}}}{f_{\mathrm{CB}}}',
      'derive': r'\begin{aligned} N_{++} &= \tfrac12\langle N(N-1)\rangle\,\varepsilon_+^2, \qquad N_{--} = \tfrac12\langle N(N-1)\rangle\,\varepsilon_-^2\\ N_{+-}^{\mathrm{CB}} &= \langle N(N-1)\rangle\,\varepsilon_+\varepsilon_- = 2\sqrt{N_{++}N_{--}} \end{aligned}',
      'chi2': r'\chi^2 = \frac{D^2}{4}\left(\frac{1}{\sigma_1^2}+\frac{1}{\sigma_2^2}\right), \qquad D = (t_1 - t_2) - (t_1^{\mathrm{p}} - t_2^{\mathrm{p}})',
      'kmix': r'k_{\mathrm{mix}}(M) = \frac{N_{+-}^{\mathrm{mix}}(M)}{2\sqrt{N_{++}^{\mathrm{mix}}(M)\,N_{--}^{\mathrm{mix}}(M)}}',
      'kfinal': r'k(M) = \frac{1}{F}\sum_{n} F_{n}\,k_{\mathrm{mix},n}, \qquad F = 2\sqrt{N_{++}N_{--}}'}
TEX = r'''\documentclass{article}\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage{amsmath}\pagestyle{empty}\begin{document}\noindent$\displaystyle %s$\end{document}'''
from PIL import Image
sizes = {}
tmp = tempfile.mkdtemp()
for k, e_ in EQ.items():
    open(f'{tmp}/{k}.tex', 'w').write(TEX % e_)
    r_ = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory', tmp, f'{tmp}/{k}.tex'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace')
    if r_.returncode != 0: raise SystemExit(f'pdflatex failed for {k}:\n' + r_.stdout[-1500:])
    out_ = os.path.join(HERE, 'eq', f'{k}.png')
    subprocess.run(['mutool', 'draw', '-r', '500', '-o', out_, f'{tmp}/{k}.pdf', '1'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    im = Image.open(out_).convert('RGB'); arr = np.asarray(im); ink = np.where(arr.min(axis=2) < 250)
    pad = 6; box = (max(ink[1].min() - pad, 0), max(ink[0].min() - pad, 0), min(ink[1].max() + pad + 1, im.width), min(ink[0].max() + pad + 1, im.height))
    im = im.crop(box); g_ = np.asarray(im.convert('L'), dtype=np.uint8)
    rgba = np.zeros(g_.shape + (4,), dtype=np.uint8); rgba[..., :3] = (0x11, 0x18, 0x27); rgba[..., 3] = 255 - g_
    Image.fromarray(rgba, 'RGBA').save(out_); sizes[k] = [round(im.width/500, 3), round(im.height/500, 3)]
shutil.rmtree(tmp)
json.dump({'eq_sizes': sizes}, open(os.path.join(HERE, 'data/deck_data.json'), 'w'), indent=1)

# ---- numbers for the tables
pl = lambda v, nd=3: (f'{v:.{nd}f}').replace('.', ',').replace('-', '−')
th = lambda v: f'{int(round(v)):,}'.replace(',', ' ')
def win(a, b): return (lo >= a - 1e-9) & (hi <= b + 1e-9)
T = {}
rows = []
for a, b_ in [(0, .02), (.02, .05), (.05, .14), (.14, .30), (.30, .56), (.14, .80)]:
    m = win(a, b_); row = {'w': f'{pl(a,2)}–{pl(b_,2)}'}
    for v in ('all', 'best'):
        Nv, Sv, Bv, CBv, e2 = N(v)[m].sum(), Sg(v)[m].sum(), B(v)[m].sum(), CB(v)[m].sum(), (CBe(v)[m]**2).sum()
        R = CBv/Bv; eR = R*np.sqrt(e2/CBv**2 + 1/Bv); row[v] = dict(N=th(Nv), S=th(Sv), B=th(Bv), CB=th(CBv), R=pl(R), eR=pl(eR), Sr=pl((Nv - CBv)/Sv), raw=[Nv, Sv, Bv, CBv, R])
    rows.append(row)
T['windows'] = rows
# decomposition
dec = []
for a, b_ in [(.05, .14), (.14, .30), (.30, .56), (.14, .80)]:
    m = win(a, b_); RA = CB('all')[m].sum()/B('all')[m].sum(); RB = CB('best')[m].sum()/B('best')[m].sum(); fL = CB('best')[m].sum()/CB('all')[m].sum(); fB = B('best')[m].sum()/B('all')[m].sum()
    dec.append(dict(w=f'{pl(a,2)}–{pl(b_,2)}', RA=pl(RA), fL=pl(fL), fB=pl(fB), RB=pl(RB), check=pl(RA*fL/fB)))
T['dec'] = dec
# classes in 0.14-0.80
m = win(.14, .80); cl = []
BT = sum(A(f'hcB_all_all_oa0_c{c}')[m].sum() for c in range(2, 6))
for c in range(2, 6):
    ba = A(f'hcB_all_all_oa0_c{c}')[m].sum(); bb = A(f'hcB_all_best_oa0_c{c}')[m].sum()
    pa = A(f'hc_all_all_oa0_pp_c{c}')[m].sum(); pb = A(f'hc_all_best_oa0_pp_c{c}')[m].sum(); ma = A(f'hc_all_all_oa0_mm_c{c}')[m].sum(); mb = A(f'hc_all_best_oa0_mm_c{c}')[m].sum()
    cl.append(dict(c=c, share=pl(100*ba/BT, 1), keptB=pl(bb/ba), keptP=pl(pb/pa) if pa > 0 else '–', keptM=pl(mb/ma) if ma > 0 else '–', B=th(ba)))
T['classes'] = cl
# channels
TOT = {}
chans = {'pi0': [0], 'eta': [1], 'cpi0': [2], 'ceta': [3], 'coth': [4]}
for k, cats in chans.items(): TOT[k] = th(cat('all', 0, 'pm', cats).sum()); TOT[k+'_raw'] = float(cat('all', 0, 'pm', cats).sum())
T['channels'] = TOT
# like-sign composition in 0.14-0.80
def comp(p):
    out = {}; tot = sum(cat('all', 0, p, [c])[m].sum() for c in range(15))
    for name, cats in (('dd', [6]), ('dc', [7]), ('cc', [8]), ('sec', [9, 10]), ('fake', [11, 12, 13])): out[name] = pl(100*cat('all', 0, p, cats)[m].sum()/tot, 1)
    out['tot'] = th(tot); return out
T['ls'] = {'pp': comp('pp'), 'mm': comp('mm')}
tb = cat('all', 0, 'pm', [6,7,8,9,10,11,12,13])[m].sum(); T['cbcomp'] = {n: pl(100*cat('all',0,'pm',c)[m].sum()/tb,1) for n,c in (('dd',[6]),('dc',[7]),('cc',[8]),('sec',[9,10]),('fake',[11,12,13]))}
# isBest >= 0 against isBest == 1
m_ = win(.14, .80); T['compare'] = {}
for v in ('all', 'best'):
    Nv, Sv, Bv, CBv, e2 = N(v)[m_].sum(), Sg(v)[m_].sum(), B(v)[m_].sum(), CB(v)[m_].sum(), (CBe(v)[m_]**2).sum(); R = CBv/Bv
    T['compare'][v] = dict(N=th(Nv), S=th(Sv), B=th(Bv), CB=th(CBv), pp=th(Npp(v)[m_].sum()), mm=th(Nmm(v)[m_].sum()), R=pl(R), eR=pl(R*np.sqrt(e2/CBv**2 + 1/Bv)), Sr=pl((Nv - CBv)/Sv), fom=pl(Sv/np.sqrt(Nv + e2), 1))
leg = d['legs']; T['pid'] = {'nep': th(leg[0]), 'nem': th(leg[1]), 'in1': pl(100*(d['inside_ep'][0]+d['inside_em'][0])/(leg[0]+leg[1]),1), 'in3': pl(100*(d['inside_ep'][1]+d['inside_em'][1])/(leg[0]+leg[1]),1), 'in5': pl(100*(d['inside_ep'][2]+d['inside_em'][2])/(leg[0]+leg[1]),1)}
T['counts'] = dict(ev=th(d['n_events_pairs'][0]), pairs=th(d['n_events_pairs'][1]), pp=th(Npp('all').sum()), mm=th(Nmm('all').sum()))
# parents
t_ = d['genuine_tally']; g = [(int(t_[i]), int(t_[i+1]), int(t_[i+2]), t_[i+3]) for i in range(0, len(t_), 4)]
def share(par, gen): return sum(c for p, g_, gp, c in g if p == par and g_ in gen)
pt = sum(c for p, _, _, c in g if p == 7)
ct = sum(c for p, _, _, c in g if p == 1)
codes = [('π⁰ pierwotne (kod 707)', {707}), ('η (kod 1717)', {1717}), ('ρ⁺ → π⁺π⁰ (kod 420807)', {420807}), ('ω → π⁺π⁻π⁰ (kod 52080908)', {52080908}), ('ρ⁻ → π⁻π⁰ (kod 430907)', {430907})]
T['src'] = [dict(lab=lab, pi0=pl(share(7, c)/pt*100, 2), conv=pl(share(1, c)/ct*100, 2)) for lab, c in codes]
T['src'].append(dict(lab='pozostałe', pi0=pl(100 - sum(share(7, c) for _, c in codes)/pt*100, 2), conv=pl(100 - sum(share(1, c) for _, c in codes)/ct*100, 2)))
T['src_tot'] = dict(pi0=th(pt), conv=th(ct))
T['eta_direct'] = pl(share(17, {1717})/sum(c for p, _, _, c in g if p == 17)*100, 2)

# ---- classes: overview plot of the "stays" columns
cl_ = T['classes']; ns = [2, 3, 4, 5]; num_ = lambda v: float(v.replace(',', '.'))
fig, ax = plt.subplots(figsize=(5.4, 3.5))
for k_, lab, c, mk in (('keptB', 'tło e⁺e⁻', RED, 's'), ('keptP', 'e⁺e⁺', C2, '^'), ('keptM', 'e⁻e⁻', C5, 'v')):
    ax.plot(ns, [num_(r_[k_]) if r_[k_] != '–' else np.nan for r_ in cl_], mk + '-', ms=8, lw=2.2, color=c, label=lab)
ax.plot(ns, [1, 0.5, 0.25, 1/6], color=INK2, lw=1.2, ls=':', label='1 / liczba hipotez e⁺e⁻')
ax.set_xticks(ns); ax.set_xticklabels(['2', '3', '4', '5+']); ax.set_ylim(0, 1.1); ax.yaxis.set_major_formatter(comma); ax.set_xlabel('liczba leptonów w zdarzeniu, n'); ax.set_ylabel('zostaje'); ax.legend(loc='lower left', fontsize=11.5)
save(fig, 'classes')
# ---- simulation: per-source check of the formula (isBest >= 0)
why = []
kWh = np.array([.05, .14, .20, .30, .56]); ksmA, _ = k_of_mix(*pools('sim', 'all', 'all'), kWh); ksmN, _ = k_of_mix(*pools('sim', 'nge3', 'all'), kWh)
for i_, (a_, b2_) in enumerate(zip(kWh[:-1], kWh[1:])):
    mw = win(a_, b2_); Rt = CB('all')[mw].sum()/B('all')[mw].sum()
    pm_ = cat('all', 0, 'pm', [6, 7, 8])[mw].sum(); pp_ = cat('all', 0, 'pp', [6, 7, 8])[mw].sum(); mm_ = cat('all', 0, 'mm', [6, 7, 8])[mw].sum()
    sec_ = cat('all', 0, 'pm', [9])[mw].sum()/cat('all', 0, 'pm', [6, 7, 8, 9, 10, 11, 12, 13])[mw].sum()
    why.append(dict(w=f'{pl(a_,2)}–{pl(b2_,2)}', R=pl(Rt), Rlep=pl(2*np.sqrt(pp_*mm_)/pm_), ppmm=pl(pp_/mm_, 2), kA=pl(ksmA[i_]), kN=pl(ksmN[i_]), sec=pl(100*sec_, 1)))
T['why'] = why
# ---- numbers used in the text
T['ls_pp_mm'] = {f'{a_}-{b2_}': pl(cat('all', 0, 'pp', [6, 7, 8])[win(a_, b2_)].sum()/cat('all', 0, 'mm', [6, 7, 8])[win(a_, b2_)].sum(), 2) for a_, b2_ in ((.05, .14), (.14, .20), (.20, .30), (.30, .56))}
# ---- k factor table
KT = k_final(np.array([.05, .14, .20, .30, .56, .80])); rk = []
for i_, (a_, b2_) in enumerate(zip([.05, .14, .20, .30, .56], [.14, .20, .30, .56, .80])):
    rk.append(dict(w=f'{pl(a_,2)}–{pl(b2_,2)}', kt=pl(KT['kt'][i_]), ekt=pl(KT['ekt'][i_]), ksa=pl(KT['ksa'][i_]), ksn=pl(KT['ksn'][i_]), kda=pl(KT['kda'][i_]), kdn=pl(KT['kdn'][i_]), k1=pl(KT['k1'][i_]), k2=pl(KT['k2'][i_]), k=pl(KT['k'][i_]), err=pl(KT['tot'][i_])))
T['k'] = rk
e58 = band58.copy()
def s_err(a_, b_):
    m_ = win(a_, b_); st_ = (toys_*dCa)[:, m_].sum(axis=1).std(); sy_ = 0.5*abs(((kup - kdn)*dCa)[m_].sum()); return float(np.hypot(st_, sy_))
def k_w(a_, b_):
    m_ = win(a_, b_); C_ = dCa[m_]; return (k58[m_]*C_).sum()/C_.sum(), s_err(a_, b_)/C_.sum()
for r_, (a_, b2_) in zip(rk, zip([.05, .14, .20, .30, .56], [.14, .20, .30, .56, .80])): kw_ = k_w(a_, b2_); r_['kw'] = pl(kw_[0]); r_['ekw'] = pl(kw_[1])
kw_ = k_w(.14, .80); T['k_int'] = dict(k=pl(kw_[0]), err=pl(kw_[1]))
kw_ = k_w(.14, .30); T['k_int30'] = dict(k=pl(kw_[0]), err=pl(kw_[1]))
# ---- data table
dT = []
for a_, b2_ in ((.05, .14), (.14, .30), (.30, .56), (.56, .80), (.14, .80)):
    m_w = win(a_, b2_); Na, Nb, Ca, Cb, Pa, Pb, Ma, Mb = dN[m_w].sum(), dNb[m_w].sum(), dCa[m_w].sum(), dCb[m_w].sum(), dPa[m_w].sum(), dPb[m_w].sum(), dMa[m_w].sum(), dMb[m_w].sum()
    kCa = (k58*dCa)[m_w].sum(); Ce = (k58*dCe)[m_w].sum(); Skm = Na - (k58s*dCa)[m_w].sum()
    Sa = Na - Ca; Sk = Na - kCa; Sb = Nb - Cb; Se = Nb - Ce
    dT.append(dict(eS=th(s_err(a_, b2_)), SkN=float(Sk), SaN=float(Sa), SbN=float(Sb), SkmN=float(Skm), w=f'{pl(a_,2)}–{pl(b2_,2)}', Na=th(Na), Nb=th(Nb), Ca=th(Ca), Cb=th(Cb), Sa=th(Sa), Sk=th(Sk), Skm=th(Skm), Sb=th(Sb), Se=th(Se), Ce=th(Ce), rb=pl(Sb/Sk), re=pl(Se/Sk), rnoK=pl(Sa/Sk), eps=pl(Nb/Na), epp=pl(Pb/Pa), emm=pl(Mb/Ma), rC=pl(Cb/Ca), kav=pl(kCa/Ca)))
T['data'] = dT
T['kfit'] = dict(a=pl(pfit[0]), mu=pl(pfit[1]), sR=pl(pfit[2]), c=pl(pfit[3]), chi2=pl(chi2fit, 1), ndf=str(len(KC) - 4))
# trigger: PT2 against PT3 (two pools, as in the simple procedure); simulation PT2 = bit 12 = practically all events
we_ = np.array([.05, .14, .20, .30, .56, .80])
kmx_ = lambda p_: k_of_mix(X(p_ + '_pm'), X(p_ + '_pp'), X(p_ + '_mm'), we_)[0]
ka2_, ka3_ = 0.5*(kmx_('kpool_data2_all') + kmx_('kpool_data2_nge3')), 0.5*(kmx_('kpool_data_all') + kmx_('kpool_data_nge3'))
T['trig'] = [dict(w=f'{pl(a_,2)}–{pl(b2_,2)}', k2=pl(ka2_[i_]), k3=pl(ka3_[i_])) for i_, (a_, b2_) in enumerate(zip(we_[:-1], we_[1:]))]
# k(M, pT)
NM_, NP_ = 7, 6; Y_ = lambda k_: np.array(x_[k_], float).reshape(NM_, NP_); g2_ = lambda p_, m_: np.where((p_ > 0) & (m_ > 0), 2*np.sqrt(np.maximum(p_, 0)*np.maximum(m_, 0)), 0.)
B2d, Fs2d_, Fd2d_ = Y_('t2_sim_B'), g2_(Y_('t2_sim_pp'), Y_('t2_sim_mm')), g2_(Y_('t2_data_pp'), Y_('t2_data_mm'))
wm_ = [.05, .14, .16, .20, .26, .34, .56, .80]
T['cpt'] = []
for i_ in range(NM_):
    vals_ = []
    for pool in ('all', 'nge3'):
        ms_, gs_, md_, gd_ = Y_(f'k2_sim_{pool}_pm')[i_], g2_(Y_(f'k2_sim_{pool}_pp')[i_], Y_(f'k2_sim_{pool}_mm')[i_]), Y_(f'k2_data_{pool}_pm')[i_], g2_(Y_(f'k2_data_{pool}_pp')[i_], Y_(f'k2_data_{pool}_mm')[i_])
        with np.errstate(divide='ignore', invalid='ignore'): k2_ = md_/gd_
        ok_ = np.isfinite(k2_) & (Fd2d_[i_] > 0); eff_ = (k2_[ok_]*Fd2d_[i_][ok_]).sum()/Fd2d_[i_][ok_].sum(); vals_.append(eff_/(md_.sum()/gd_.sum()))
    T['cpt'].append(dict(w=f'{pl(wm_[i_],2)}–{pl(wm_[i_+1],2)}', all=pl(vals_[0]), n3=pl(vals_[1])))
tg_ = lambda m_, i_: x_[f'tag_{m_}_win'][i_]/x_[f'tag_{m_}_tot'][i_]
T['tag'] = dict(all_d=pl(tg_('data', 0)), same_d=pl(tg_('data', 1)), rpc_d=pl(tg_('data', 2)), tof_d=pl(tg_('data', 3)), all_s=pl(tg_('sim', 0)), same_s=pl(tg_('sim', 1)), rpc_s=pl(tg_('sim', 2)), tof_s=pl(tg_('sim', 3)), n=th(x_['tag_data_tot'][0]))
# eps_n tables: by window, shares of the formula, trigger comparison, like-sign check
def eps_w(c_, a_, b_, cross='emu_cross', win_='emu_win'):
    m_ = win(a_, b_); return X(f'{win_}_c{c_}')[m_].sum()/X(f'{cross}_c{c_}')[m_].sum()
EW = [(.05, .14), (.14, .30), (.30, .80)]
mS = win(.14, .80); shr = [gm(X(f'emu_pp_c{c}'), X(f'emu_mm_c{c}'))[mS].sum() for c in range(2, 6)]
T['eps'] = [dict(n='2' if c == 2 else ('3' if c == 3 else ('4' if c == 4 else '5 i więcej')), hyp='1' if c == 2 else ('2' if c == 3 else ('3 lub 4' if c == 4 else '6 i więcej')), e=[('1,000' if c == 2 else pl(eps_w(c, a_, b_))) for a_, b_ in EW], share=pl(100*shr[c - 2]/sum(shr), 0)) for c in range(2, 6)]
T['eps_tr'] = {f'n{c}': [dict(w=f'{pl(a_,2)}–{pl(b_,2)}', p3=pl(eps_w(c, a_, b_)), p2=pl(eps_w(c, a_, b_, 'emu2_cross', 'emu2_win'))) for a_, b_ in EW] for c in (3, 4)}
def cl_(p_, c_):
    m_ = win(.14, .80); d_ = X(f'cl_{p_}_best_c{c_}')[m_].sum()/X(f'cl_{p_}_all_c{c_}')[m_].sum(); e_ = X(f'cl_{p_}_crosswin_c{c_}')[m_].sum()/X(f'cl_{p_}_cross_c{c_}')[m_].sum(); return dict(d=pl(d_), e=pl(e_))
T['eps_ls'] = {p_: {f'n{c_}': cl_(p_, c_) for c_ in (3, 4)} for p_ in ('pp', 'mm')}
# signal per lepton class and the signal kept by isBest == 1
def cls_(a_, b_):
    m_ = win(a_, b_); N_ = np.array([X(f'emu_pm_c{c}')[m_].sum() for c in range(2, 6)]); S_ = np.array([(X(f'emu_pm_c{c}') - k58*gm(X(f'emu_pp_c{c}'), X(f'emu_mm_c{c}')))[m_].sum() for c in range(2, 6)])
    e_ = np.array([1.] + [eps_w(c, a_, b_) for c in (3, 4, 5)]); return N_, S_, float((S_*e_).sum()/S_.sum())
N_c, S_c, _ = cls_(.14, .80)
T['cls'] = [dict(n=n_, pairs=pl(100*N_c[i]/N_c.sum(), 1), sig=pl(100*S_c[i]/S_c.sum(), 1), sn=pl(100*S_c[i]/N_c[i], 1)) for i, n_ in enumerate(('2', '3', '4', '5 i więcej'))]
T['keep'] = [dict(w=r_['w'], meas=r_['re'], exp=pl(cls_(a_, b2_)[2])) for r_, (a_, b2_) in zip(T['data'], ((.05, .14), (.14, .30), (.30, .56), (.56, .80), (.14, .80)))]
T['weightCB'] = dict(naive=th(dCb[mS].sum()), w=th(dCe[mS].sum()), r=pl(dCe[mS].sum()/dCb[mS].sum(), 2))
# k and eps_n in other observables: spread across bins
def spread(y_):
    v_ = y_[np.isfinite(y_)]; return (pl(v_.min(), 2), pl(v_.max(), 2))
T['proj'] = {}
for w_ in (0, 1):
    T['proj'][str(w_)] = {}
    for v_, *_ in PV:
        ka_, ma_ = kx_(w_, v_, 'all'); kn_, mn_ = kx_(w_, v_, 'nge3'); e3_, m3_ = ex_(w_, v_, 'n3'); e4_, m4_ = ex_(w_, v_, 'n4')
        T['proj'][str(w_)][v_] = dict(k=spread(ka_), kn=spread(kn_), km=pl(ma_), e3=spread(e3_), e3m=pl(m3_), e4=spread(e4_), e4m=pl(m4_))
ch = lambda w_, d_: (x_['choice_win'][4*w_ + d_], x_['choice_tot'][4*w_ + d_])
T['choice'] = {k_: dict(p=pl(ch(2, i_)[0]/ch(2, i_)[1]), n=th(ch(2, i_)[1])) for i_, k_ in enumerate(('all', 'same', 'grpc', 'gtof'))}
T['choice_rpc'] = dict(sim=pl(x_['choice_rpc'][0]/x_['choice_rpc'][1]), simT=pl(x_['rpc_sim_PT3'][0]/x_['rpc_sim_PT3'][1]), data=pl(x_['rpc_data_PT3'][0]/x_['rpc_data_PT3'][1]), dataN=th(x_['rpc_data_PT3'][1]))
Dc_ = 0.5*(De[:-1] + De[1:])
def g2rms(k_):
    y_ = X(f'hD_{k_}_n3'); m0 = (y_*Dc_).sum()/y_.sum(); return float(np.sqrt((y_*(Dc_ - m0)**2).sum()/y_.sum()))
T['D'] = dict(gen=pl(g2rms('gen')), lep=pl(g2rms('lep')), fake=pl(g2rms('fake')))
json.dump(T, open(os.path.join(HERE, 'data/deck_tables.json'), 'w'), ensure_ascii=False, indent=1)
for k_ in ('cls', 'keep', 'eps', 'eps_tr', 'eps_ls', 'weightCB', 'proj', 'k', 'k_int', 'k_int30', 'kfit', 'trig', 'cpt', 'tag', 'why', 'data', 'choice', 'choice_rpc', 'D', 'ls_pp_mm'): print(k_, json.dumps(T[k_], ensure_ascii=False))

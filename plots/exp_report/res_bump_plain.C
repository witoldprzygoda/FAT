// res_bump_plain.C -- COPY of plots/exp_report/res_bump_cnt.C made by the short-report figure agent (integ6, 2026-11-16):
// internal catalogue labels in the figure text replaced by plain names (glossary_v2.md); outputs redirected:
// figures of the scope -> docs/latex/<dir>/<base>_plain.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/run/res
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ6/shortfig/mk_plain.py; the physics code is unchanged.
// res_bump_cnt.C -- COPY of plots/exp_report/res_bump.C made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/res_bump
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// res_bump.C  --  the bump below the eta mass (eta -> pi+pi-pi0, photon lost): data only
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22, 29 days), channel pi+ pi- e+ e-. EXPERIMENTAL
// DATA ONLY. Chapter "Results without simulation" (prefix res), task 2.
//
// What it does
//   1. edge step (method of plots/scan_bump.C, engine scan::FitEdgeStep): in the
//      M(pi+pi-) slices 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40 the spectrum
//      N - CB is fitted with pol2 [+ eta Gaussian] + h x step, step = smeared
//      kinematic edge M4max(M(pipi)) of eta -> pi+pi-pi0 with the photon of
//      pi0 -> e+e-gamma lost, position FIXED (M4 scale 0.99277, M(pipi) scale
//      0.9970), smearing 5.4 MeV; h = height of the drop in rows per 2.5 MeV;
//      E = sum of h over the four slices (a visibility measure, NOT a yield);
//      slice 0.40-0.41 reported, not in E; control slice 0.42-0.46 above the
//      kinematic limit m_eta - m_pi0 = 0.413 (edge set by hand at 0.5150, the
//      M(pipi) cut of the bump-nominal selections is dropped there);
//      model spread: half range of h over order 2, 3 x half range 35, 45, 55 MeV
//      x smearing 4.4, 5.4, 6.4 MeV (fits with p > 0.05), added linearly over slices;
//   2. edge position fitted free: the factor between the true limit and the
//      measured edge scanned per slice, weighted mean against the eta peak;
//   3. bump window 0.4550-0.5225: N, CB, N - CB (contains the bump AND the
//      physical background);
//   4. ECAL-tagged eta (frozen definition of plots/scan_bump.C; photon energy
//      rescaled so that M(e+e-gamma) = m_pi0): content T, S/B, free width,
//      model spread; pi0 in M(e+e-gamma) of the tagged eta; share of the tagged
//      eta with M4 inside the bump window.
// Selections (rows = stored combinations; cor_ok == 1; FAT: vertex z > -500 mm, PT3):
//   pre = BUMP (frozen)  = cor_ok && isBest >= 0 && vz > -160 mm
//   bump nominal         = pre + oa_lab < 40 deg + 0.30 < M(pipi) < 0.41 GeV/c2
//   bump nominal+vertex  = bump nominal + both pions in the vertex fit
//   variants: isBest >= 0 (main) / isBest == 1; without / with OA(e+e-) > 2 deg
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root (scan::kCompact)
// Outputs: docs/latex/tables_exp/res_bump.txt, res_bump.tex, res_tagged.tex
//          docs/latex/figs_exp/res_bump_slices.pdf/.png, res_bump_tagged.pdf/.png
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/res/res_bump.root
//
// Run (from /home/przygoda/HADES/FAT/FAT, about 1 min):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_res", kTRUE)' 'plots/exp_report/res_bump.C+'
// ============================================================================
#include "TDatime.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStopwatch.h"

#include "res_common_plain.h"

namespace resbp {

using namespace res;
using scan::StepCfg;
using scan::StepFit;

// ---- slices (as plots/scan_bump.C) ---------------------------------------------
const int kNSl = 6, kNSig = 4;
const double kSlLo[kNSl] = {0.30, 0.34, 0.37, 0.39, 0.40, 0.42};
const double kSlHi[kNSl] = {0.34, 0.37, 0.39, 0.40, 0.41, 0.46};
const char* const kSlName[kNSl] = {"0.30-0.34", "0.34-0.37", "0.37-0.39", "0.39-0.40", "0.40-0.41", "0.42-0.46"};
const double kCtrlEdge = 0.5150;
inline StepCfg SliceCfg(int i, double lo, double hi) {
    StepCfg c;
    c.mpLo = lo; c.mpHi = hi;
    if (i == 3 || i == 4) c.halfRange = 0.035;
    if (i == 5) { c.useEdge = true; c.edge = kCtrlEdge; }
    return c;
}

struct BSel {
    std::string key, tex;
    Sel f;        // full selection
    Sel fCtrl;    // the same without the M(pipi) cut (control slice)
    double mpLo, mpHi;   // M(pipi) cut (-1: none)
};
std::vector<BSel> Selections() {
    std::vector<BSel> v;
    v.push_back({"pre (= BUMP)", "pre ($=$ BUMP)", [](const Row& r) { return Pre(r); }, [](const Row& r) { return Pre(r); }, -1, -1});
    v.push_back({"bump nominal", "bump nominal", [](const Row& r) { return Pre(r) && BumpNomKin(r); },
                 [](const Row& r) { return Pre(r) && r.oa_lab < 40.f; }, 0.30, 0.41});
    v.push_back({"bump nominal+vertex", "bump nominal + vertex", [](const Row& r) { return Pre(r) && BumpNomKin(r) && PiVtx(r); },
                 [](const Row& r) { return Pre(r) && r.oa_lab < 40.f && PiVtx(r); }, 0.30, 0.41});
    return v;
}
inline void SliceBounds(const BSel& s, int i, double& lo, double& hi) {
    lo = kSlLo[i]; hi = kSlHi[i];
    if (i == 5) return;                       // control: the M(pipi) cut is dropped
    if (s.mpLo >= 0) lo = std::max(lo, s.mpLo);
    if (s.mpHi >= 0) hi = std::min(hi, s.mpHi);
}

struct EdgeRes {
    StepFit f[kNSl]; bool ok[kNSl];
    double half[kNSig]; int kept[kNSig];
    double E = 0, eE = 0, spread = 0, chi2 = 0; int ndf = 0;
    scan::WinSum W;
};
EdgeRes Measure(const std::vector<Spectra>& sl, const BSel& s) {
    EdgeRes R;
    double v = 0;
    for (int i = 0; i < kNSl; ++i) {
        R.ok[i] = false;
        double lo, hi; SliceBounds(s, i, lo, hi);
        if (hi - lo < 0.0045 || sl[i].nsel[0] < 50) continue;
        R.f[i] = scan::FitEdgeStep(sl[i], SliceCfg(i, lo, hi));
        R.ok[i] = R.f[i].ok;
        if (!R.ok[i] || i >= kNSig) continue;
        R.E += R.f[i].h; v += R.f[i].eh * R.f[i].eh; R.chi2 += R.f[i].chi2; R.ndf += R.f[i].ndf;
    }
    R.eE = std::sqrt(v);
    // model spread per slice, added linearly
    for (int i = 0; i < kNSig; ++i) {
        R.half[i] = 0; R.kept[i] = 0;
        if (!R.ok[i]) continue;
        double lo, hi; SliceBounds(s, i, lo, hi);
        double mn = R.f[i].h, mx = R.f[i].h;
        for (int ord = 2; ord <= 3; ++ord)
            for (int ih = 0; ih < 3; ++ih)
                for (int is = -1; is <= 1; ++is) {
                    StepCfg c = SliceCfg(i, lo, hi);
                    c.order = ord; c.halfRange = 0.035 + 0.010 * ih; c.s = scan::kEdgeSigma + 0.001 * is;
                    const StepFit f = scan::FitEdgeStep(sl[i], c);
                    if (!f.ok || !(TMath::Prob(f.chi2, f.ndf) > 0.05)) continue;
                    R.kept[i]++;
                    mn = std::min(mn, f.h); mx = std::max(mx, f.h);
                }
        R.half[i] = 0.5 * (mx - mn);
        R.spread += R.half[i];
    }
    return R;
}

// free edge position (copy of plots/scan_bump.C FreeEdge)
struct EdgeFree {
    bool ok = false;
    double scale = 0, eLo = 0, eHi = 0, chi2min = 0, chi2nom = 0; int ndf = 0;
    double edge = 0, pred = 0, h = 0, eh = 0;
};
EdgeFree FreeEdge(const Spectra& sp, int i, double lo, double hi) {
    EdgeFree r;
    const StepCfg c0 = SliceCfg(i, lo, hi);
    const StepFit nom = scan::FitEdgeStep(sp, c0);
    if (!nom.ok) return r;
    const double ec = scan::M4max(0.5 * (lo + hi) / scan::kScalePiPi);
    r.pred = scan::kScale * ec; r.chi2nom = nom.chi2; r.ndf = nom.ndf;
    std::vector<double> sc, cs, hh, eh;
    int ib = -1;
    for (double s = 0.9700; s < 1.01501; s += 0.00025) {
        StepCfg v = c0;
        v.scale = s; v.lo = nom.lo; v.hi = nom.hi; v.withEta = nom.withEta ? 1 : 0;
        const StepFit f = scan::FitEdgeStep(sp, v);
        if (!f.ok) continue;
        sc.push_back(s); cs.push_back(f.chi2); hh.push_back(f.h); eh.push_back(f.eh);
        if (ib < 0 || f.chi2 < cs[ib]) ib = (int)sc.size() - 1;
    }
    if (ib <= 0 || ib >= (int)sc.size() - 1) return r;
    r.scale = sc[ib]; r.chi2min = cs[ib]; r.h = hh[ib]; r.eh = eh[ib];
    int a = ib, b = ib;
    while (a > 0 && cs[a] < r.chi2min + 1.0) --a;
    while (b < (int)sc.size() - 1 && cs[b] < r.chi2min + 1.0) ++b;
    r.eLo = r.scale - sc[a]; r.eHi = sc[b] - r.scale;
    r.edge = r.scale * ec;
    r.ok = true;
    return r;
}

void Txt(double x, double y, const char* t, double size = 0.050, int col = kBlack, int align = 11) {
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size); l.SetTextColor(col); l.SetTextAlign(align);
    l.DrawLatex(x, y, t);
}

} // namespace resbp

void res_bump_plain() {
    using namespace resbp;
    TStopwatch swAll;
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE); gSystem->mkdir(kWork, kTRUE); gSystem->mkdir(kPlainScratch, kTRUE);
    Style();

    Data dE, dG;
    {
        Data d;
        if (!scan::Load(d)) return;
        Filter(d, dE, [](const Row& r) { return Pre(r) && r.m4_fb >= 0.30f && r.m4_fb < 0.66f; });
        Filter(d, dG, [](const Row& r) { return Pre(r) && HasGam(r) && r.m_eeg >= 0.05f && r.m_eeg < 0.25f; });
    }
    printf("rows: M4 sample %zu / %zu / %zu, photon sample %zu / %zu / %zu\n", dE.s[0].size(), dE.s[1].size(), dE.s[2].size(),
           dG.s[0].size(), dG.s[1].size(), dG.s[2].size());

    // ------------------------------------------------------------------------
    // 1. edge steps for 3 selections x 4 variants
    // ------------------------------------------------------------------------
    const std::vector<BSel> bs = Selections();
    const int nB = (int)bs.size(), nV = 4;
    std::vector<Sel> sels;
    for (int b = 0; b < nB; ++b)
        for (int v = 0; v < nV; ++v) {
            for (int i = 0; i < kNSl; ++i) {
                const float lo = kSlLo[i], hi = kSlHi[i];
                const Sel base = Variant(i == 5 ? bs[b].fCtrl : bs[b].f, v);
                sels.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && base(r); });
            }
            sels.push_back(Variant(bs[b].f, v));       // whole selection (bump window)
        }
    const int per = kNSl + 1;
    const std::vector<Spectra> sp = scan::FillMany(dE, sels);
    std::vector<EdgeRes> R(nB * nV);
    for (int b = 0; b < nB; ++b)
        for (int v = 0; v < nV; ++v) {
            const int k = b * nV + v;
            std::vector<Spectra> sl(sp.begin() + k * per, sp.begin() + k * per + kNSl);
            R[k] = Measure(sl, bs[b]);
            R[k].W = sp[k * per + kNSl].win[1];
        }
    // free edge (isBest >= 0, no OA cut) for pre and bump nominal + vertex
    EdgeFree fe[2][kNSig];
    double feMean[2] = {0, 0}, feErr[2] = {0, 0};
    for (int q = 0; q < 2; ++q) {
        const int b = q == 0 ? 0 : 2, k = b * nV;
        double sw = 0, sx = 0;
        for (int i = 0; i < kNSig; ++i) {
            double lo, hi; SliceBounds(bs[b], i, lo, hi);
            fe[q][i] = FreeEdge(sp[k * per + i], i, lo, hi);
            if (!fe[q][i].ok) continue;
            const double e = 0.5 * (fe[q][i].eLo + fe[q][i].eHi), w = e > 0 ? 1 / (e * e) : 0;
            sw += w; sx += w * fe[q][i].scale;
        }
        if (sw > 0) { feMean[q] = sx / sw; feErr[q] = 1 / std::sqrt(sw); }
    }

    // ------------------------------------------------------------------------
    // 2. ECAL-tagged eta
    // ------------------------------------------------------------------------
    std::vector<Sel> ts;
    for (int b = 0; b < nB; ++b) for (int v = 0; v < nV; ++v) { const Sel f = Variant(bs[b].f, v); ts.push_back([f](const Row& r) { return Tagged(r) && f(r); }); }
    ts.push_back([](const Row& r) { return Tagged(r) && r.m4 >= scan::kBumpLo && r.m4 < scan::kBumpHi; });       // M4 in the bump window
    ts.push_back([](const Row& r) { return Tagged(r) && r.m4 < scan::kBumpHi; });                                // M4 below 0.5225
    const std::vector<Spectra> st = scan::FillMany(dG, ts, TagOpt());
    std::vector<PeakFit> T(nB * nV); std::vector<double> Tsp(nB * nV); std::vector<int> Tkept(nB * nV); std::vector<scan::FreeFit> Tff(nB * nV);
    for (int k = 0; k < nB * nV; ++k) { T[k] = scan::FitPeakMain(st[k], TagCfg()); Tsp[k] = TagSpread(st[k], &Tkept[k]); Tff[k] = scan::FitPeakFree(st[k], TagCfg()); }
    const PeakFit Tbw = scan::FitPeakMain(st[nB * nV], TagCfg()), Tlo = scan::FitPeakMain(st[nB * nV + 1], TagCfg());
    // pi0 in M(e+e-gamma): eta content (m5c fit) in 5 MeV bins of M(e+e-gamma), pre, isBest >= 0
    const int nPi = 32; const double piLo = 0.060, piBw = 0.005;
    std::vector<double> piS(nPi), piE(nPi);
    double pi0Par[3], pi0Err[3], pi0Chi2; int pi0Ndf;
    {
        std::vector<Sel> sb;
        for (int i = 0; i < nPi; ++i) {
            const float a = piLo + i * piBw, b = a + piBw;
            sb.push_back([=](const Row& r) { return r.m_eeg >= a && r.m_eeg < b; });
        }
        scan::FillOpt o = TagOpt();
        o.x = [](const Row& r) { return M5c(r); };
        const std::vector<Spectra> s = scan::FillMany(dG, sb, o);
        TGraphErrors g;
        for (int i = 0; i < nPi; ++i) {
            const PeakFit f = scan::FitPeakMain(s[i], TagCfg());
            piS[i] = f.S; piE[i] = f.eS;
            g.SetPoint(i, piLo + (i + 0.5) * piBw, f.S); g.SetPointError(i, 0, f.eS);
        }
        TF1 fg("fg", "[0]*0.005/(sqrt(2*pi)*[2])*exp(-0.5*((x-[1])/[2])^2)", 0.095, 0.200);
        fg.SetParameters(3500, 0.146, 0.012);
        TFitResultPtr fr = g.Fit(&fg, "RQSN");
        for (int k = 0; k < 3; ++k) { pi0Par[k] = fr->Parameter(k); pi0Err[k] = fr->ParError(k); }
        pi0Chi2 = fr->Chi2(); pi0Ndf = fr->Ndf();
    }

    // ------------------------------------------------------------------------
    // tables
    // ------------------------------------------------------------------------
    OpenTxt("res_bump.txt");
    P("# res_bump.txt -- written by plots/exp_report/res_bump.C (%s)\n", TDatime().AsSQLString());
    P("# EXPERIMENTAL DATA, p+p 4.53 GeV (feb22, 29 days). Rows = stored combinations; cor_ok == 1; FAT: vertex z > -500 mm, PT3 (trigbit 8192).\n");
    P("# Selections: pre (= BUMP of the simulation study) = cor_ok && isBest >= 0 && vz > -160 mm; bump nominal = pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41;\n");
    P("#   bump nominal+vertex = bump nominal + both pions in the vertex fit. Variants isBest >= 0 / == 1, without / with OA(e+e-) > 2 deg.\n");
    P("# Edge step: N - CB in M(pipi) slices = pol2 [+ eta Gaussian 0.5439 / 5.4 MeV] + h x smeared step at the FIXED predicted edge\n");
    P("#   (M4 scale %.5f, M(pipi) scale %.4f, smearing 5.4 MeV; fit range edge -+ 45 MeV, -+ 35 MeV for 0.39-0.40 and 0.40-0.41).\n", scan::kScale, scan::kScalePiPi);
    P("#   h in rows per 2.5 MeV; E = sum of h over 0.30-0.34, 0.34-0.37, 0.37-0.39, 0.39-0.40 (a visibility measure, NOT a yield).\n");
    P("#   spread = half range of h over 18 variations (order 2, 3; half range 35, 45, 55 MeV; smearing 4.4, 5.4, 6.4 MeV; p > 0.05), added linearly.\n");
    P("#   control slice 0.42-0.46 (above m_eta - m_pi0): edge set by hand at %.4f; there the M(pipi) cut of the bump-nominal selections is dropped.\n", kCtrlEdge);
    P("# W = N - CB in the bump window 0.4550-0.5225 (bump + physical background; CB = 2 sqrt(N++ N--), k = 1).\n\n");
    P("A. EDGE STEPS\n");
    P("%-20s %-16s | %-37s | %-11s %-11s %-11s %-11s | %-11s | %-11s | %9s\n", "selection", "variant", "E +- stat +- spread (chi2/ndf)",
      kSlName[0], kSlName[1], kSlName[2], kSlName[3], kSlName[4], "ctrl " , "W = N-CB");
    for (int v = 0; v < nV; ++v)
        for (int b = 0; b < nB; ++b) {
            const EdgeRes& r = R[b * nV + v];
            P("%-20s %-16s | %5.0f +- %3.0f +- %3.0f (%5.1f/%-3d) %6s |", bs[b].key.c_str(), VariantName(v), r.E, r.eE, r.spread, r.chi2, r.ndf, "");
            for (int i = 0; i < kNSl; ++i) {
                if (i == kNSig || i == kNSig + 1) P(" |");
                if (r.ok[i]) P(" %4.0f+-%-4.0f", r.f[i].h, r.f[i].eh); else P(" %10s", "-");
            }
            P(" | %7.0f+-%-4.0f\n", r.W.sig(), r.W.sigErr());
        }
    P("\nA2. Per slice (isBest >= 0, no OA cut): edge, fit range, eta term, h, err, spread, background at the edge, chi2/ndf, p, dchi2 (step)\n");
    for (int b = 0; b < nB; ++b) {
        const EdgeRes& r = R[b * nV];
        P("  %s\n", bs[b].key.c_str());
        for (int i = 0; i < kNSl; ++i) {
            if (!r.ok[i]) { P("    %-10s  no fit\n", kSlName[i]); continue; }
            const StepFit& f = r.f[i];
            P("    %-10s  edge %.4f  range %.4f-%.4f eta %d | h %6.1f +- %5.1f spread %5.1f | bkg %6.0f h/bkg %.3f | %5.1f/%-3d p %.3f dchi2 %6.1f%s\n",
              kSlName[i], f.edge, f.lo, f.hi, (int)f.withEta, f.h, f.eh, i < kNSig ? r.half[i] : 0., f.bkgAtEdge, f.bkgAtEdge > 0 ? f.h / f.bkgAtEdge : 0.,
              f.chi2, f.ndf, TMath::Prob(f.chi2, f.ndf), f.dchi2, i == 4 ? "  (not in E)" : (i == 5 ? "  (control)" : ""));
        }
        P("    window: N %.0f, CB %.0f (CB/N %.3f), N - CB %.0f +- %.0f\n", r.W.n[0], r.W.cb, r.W.cbShare(), r.W.sig(), r.W.sigErr());
    }
    P("\n   ratios (isBest >= 0, no OA cut): E(bump nominal)/E(pre) = %.3f, E(bump nominal+vertex)/E(pre) = %.3f; isBest == 1 / >= 0 for pre: %.3f\n",
      R[nV].E / R[0].E, R[2 * nV].E / R[0].E, R[1].E / R[0].E);
    P("   OA > 2 pass fraction of E: pre %.3f, bump nominal %.3f, bump nominal+vertex %.3f\n", R[2].E / R[0].E, R[nV + 2].E / R[nV].E, R[2 * nV + 2].E / R[2 * nV].E);
    P("\nB. EDGE POSITION FITTED FREE (factor between the true limit and the measured edge; nominal %.5f from the eta peak +- 0.0009)\n", scan::kScale);
    for (int q = 0; q < 2; ++q) {
        P("  %s, isBest >= 0\n", q == 0 ? "pre" : "bump nominal+vertex");
        for (int i = 0; i < kNSig; ++i) {
            const EdgeFree& f = fe[q][i];
            if (!f.ok) { P("    %-10s no minimum inside 0.970-1.015\n", kSlName[i]); continue; }
            P("    %-10s factor %.4f -%.4f +%.4f | edge fit %.4f predicted %.4f diff %+.4f | chi2 min %.1f nominal %.1f ndf %d | h %.1f +- %.1f\n",
              kSlName[i], f.scale, f.eLo, f.eHi, f.edge, f.pred, f.edge - f.pred, f.chi2min, f.chi2nom, f.ndf, f.h, f.eh);
        }
        P("    weighted mean factor %.4f +- %.4f (eta peak %.4f +- 0.0009; difference %+.4f)\n", feMean[q], feErr[q], scan::kScale, feMean[q] - scan::kScale);
    }
    P("\nC. ECAL-TAGGED ETA (pi0 window %.3f <= M(e+e-gamma) < %.3f, m_ee < %.2f; m5c with the photon energy rescaled to M(e+e-gamma) = m_pi0;\n", kPi0Lo, kPi0Hi, kMeeMax);
    P("   P_main with the frozen eta shape %.4f / %.1f MeV + pol2 in %.2f-%.2f, CB in 20 MeV bins; spread over 18 variations; free = position and width free)\n",
      kTagMean, 1e3 * kTagSigma, kTagLo, kTagHi);
    P("%-20s %-16s | %7s %6s %5s | %6s %5s %5s %3s | %6s %5s %9s | %8s %6s %6s | %s\n", "selection", "variant", "rows", "N win", "CB", "T", "stat", "sprd", "fit",
      "S/B", "S/err", "chi2/ndf", "Y_T [nb]", "stat", "model", "free fit: position, width [MeV], S");
    for (int v = 0; v < nV; ++v)
        for (int b = 0; b < nB; ++b) {
            const int k = b * nV + v;
            P("%-20s %-16s | %7.0f %6.0f %5.0f | %6.0f %5.0f %5.0f %3d | %6.3f %5.1f %5.1f/%-3d | %8.4f %6.4f %6.4f |", bs[b].key.c_str(), VariantName(v),
              st[k].nsel[0], T[k].N, T[k].CB, T[k].S, T[k].eS, Tsp[k], Tkept[k], T[k].SoB, T[k].F, T[k].chi2, T[k].ndf, Yield(T[k].S), Yield(T[k].eS), Yield(Tsp[k]));
            if (Tff[k].ok) P(" %.4f +- %.4f, %.2f +- %.2f, %.0f +- %.0f\n", Tff[k].mean, Tff[k].eMean, 1e3 * Tff[k].sigma, 1e3 * Tff[k].eSigma, Tff[k].S, Tff[k].eS);
            else P(" n/a\n");
        }
    double fbw, efbw, flo, eflo;
    Ratio(Tbw.S, Tbw.eS, T[0].S, T[0].eS, fbw, efbw); Ratio(Tlo.S, Tlo.eS, T[0].S, T[0].eS, flo, eflo);
    P("   tagged eta with M4 in the bump window 0.4550-0.5225: %.0f +- %.0f, share %.3f; with M4 < 0.5225: %.0f +- %.0f, share %.3f (pre, isBest >= 0;\n",
      Tbw.S, Tbw.eS, fbw, Tlo.S, Tlo.eS, flo);
    P("   errors of the shares are not given: the samples are nested)\n");
    P("   efficiency of the cuts for the tagged eta (isBest >= 0): bump nominal %.3f, bump nominal+vertex %.3f; isBest == 1 / >= 0 (pre) %.3f\n",
      T[nV].S / T[0].S, T[2 * nV].S / T[0].S, T[1].S / T[0].S);
    P("   model spread of a slice = 0 means that no variation (and not the nominal fit) passed p > 0.05 there, see A2.\n");
    P("   pi0 of the tagged eta (eta content in 5 MeV bins of M(e+e-gamma), Gaussian fit 0.095-0.200): position %.4f +- %.4f, width %.2f +- %.2f MeV,\n",
      pi0Par[1], pi0Err[1], 1e3 * pi0Par[2], 1e3 * pi0Err[2]);
    P("   content %.0f +- %.0f, chi2/ndf %.1f/%d (m_pi0 = 0.1350: the ECAL photon energy is not calibrated; frozen window 0.119-0.173 = +-2 sigma)\n",
      pi0Par[0], pi0Err[0], pi0Chi2, pi0Ndf);
    for (int i = 0; i < nPi; ++i) P("     M(eeg) %.3f-%.3f: %7.1f +- %5.1f\n", piLo + i * piBw, piLo + (i + 1) * piBw, piS[i], piE[i]);
    P("\nWHY NO MODEL-FREE BUMP YIELD: the bump (eta -> pi+pi-pi0 with an e+e- pair from the pi0 and the photon lost) has no shape of its own in M4:\n");
    P("it is a broad distribution from about 0.40 up to the kinematic edge (0.49-0.52), lying on a multi-pion background whose shape under the\n");
    P("bump is not known (no side band below the bump; above the edge lies the eta peak; at the edge the drop h is only 0.2-0.5 of the level of\n");
    P("N - CB). A smooth polynomial absorbs most of the bump,\n");
    P("so N - CB - background in the bump window depends entirely on the background model. Only the sharp drop at the predicted edge is\n");
    P("model-independent: h measures the bump density just below the edge, per M(pipi) slice. Turning E (or W) into a number of eta -> pi+pi-pi0\n");
    P("decays needs the bump shape below the edge (Dalitz-plot matrix element, pi0 Dalitz / conversion mix, acceptance): this comes from the\n");
    P("simulation and is done in the next report. The ECAL-tagged eta T is a clean peak but counts only the events with a detected photon.\n");
    P("\nrun time %.0f s\n", swAll.RealTime());
    CloseTxt();

    // ---- LaTeX ----
    {
        Tex t("res_bump.tex");
        t("%% res_bump.tex -- written by plots/exp_report/res_bump.C (numbers: tables_exp/res_bump.txt)\n");
        t("\\begin{table}[H]\n\\centering\\scriptsize\n\\setlength{\\tabcolsep}{2pt}\n");
        t("\\caption{Edge step of the bump ($\\eta\\to\\pi^+\\pi^-\\pi^0$, $\\pi^0\\to e^+e^-\\gamma$ with the photon lost) in slices of $M(\\pi^+\\pi^-)$ [GeV/$c^2$], "
          "experimental data. $h$ = height of the drop of $N-\\mathrm{CB}$ at the predicted kinematic edge (position fixed, smearing 5.4~MeV/$c^2$) in rows per 2.5~MeV/$c^2$; "
          "$E=\\sum h$ over the first four slices $\\pm$ stat $\\pm$ model spread; $E$ is a visibility measure, not a yield. "
          "0.40--0.41: not in $E$; 0.42--0.46: control slice above the kinematic limit (edge set at 0.5150, $M_{\\pi\\pi}$ cut dropped). "
          "$W=N-\\mathrm{CB}$ in the bump window 0.4550--0.5225~GeV/$c^2$ (bump and physical background). "
          "Rows: all stored combinations (isBest $\\geq0$) or isBest $=1$, corrected momenta, PT3, vertex $z>-500$~mm, $v_z>-160$~mm; not normalised to luminosity. "
          "Edge position fitted free (factor measured/true $M_4$, isBest $\\geq0$): pre $%.4f\\pm%.4f$, bump nominal + vertex $%.4f\\pm%.4f$; "
          "from the $\\eta$ peak $0.9927\\pm0.0009$ ($0.54388/0.547862$). "
          "The statistical errors are fit errors; the calibration with independent halves of the data (Table~\\ref{tab:fix_stat}) finds them too small by a factor 1.2. "
          "Bump nominal + vertex, slice 0.39--0.40: the nominal fit has $p=0.005$ and no variation passes $p>0.05$, so its model spread is undefined (counted as 0) "
          "and the model error of that $E$ is a lower limit.}\n", feMean[0], feErr[0], feMean[1], feErr[1]);
        t("\\label{tab:res_bump}\n\\begin{tabular}{llcrrrrrrr}\n\\toprule\n");
        t("Selection & isBest & $E\\pm$stat$\\pm$model & 0.30--0.34 & 0.34--0.37 & 0.37--0.39 & 0.39--0.40 & 0.40--0.41 & 0.42--0.46 & $W$ \\\\\n\\midrule\n");
        for (int v = 0; v < nV; ++v) {
            if (v == 2) t("\\midrule\n\\multicolumn{10}{l}{\\emph{with OA$(e^+e^-)>2^\\circ$ (not a fixed cut)}} \\\\\n");
            for (int b = 0; b < nB; ++b) {
                const EdgeRes& r = R[b * nV + v];
                t("%s & %s & $%.0f\\pm%.0f\\pm%.0f$", bs[b].tex.c_str(), (v & 1) ? "$=1$" : "$\\geq0$", r.E, r.eE, r.spread);
                for (int i = 0; i < kNSl; ++i) { if (r.ok[i]) t(" & $%.0f\\pm%.0f$", r.f[i].h, r.f[i].eh); else t(" & --"); }
                t(" & $%.0f\\pm%.0f$ \\\\\n", r.W.sig(), r.W.sigErr());
            }
            if (v == 0 || v == 2) t("\\addlinespace\n");
        }
        t("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    }
    {
        Tex t("res_tagged.tex");
        t("%% res_tagged.tex -- written by plots/exp_report/res_bump.C (numbers: tables_exp/res_bump.txt)\n");
        t("\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n");
        t("\\caption{ECAL-tagged $\\eta$ (reference sample of the bump): rows with an ECAL photon and $0.119\\leq M(e^+e^-\\gamma)<0.173$~GeV/$c^2$, $m_{ee}<0.10$~GeV/$c^2$; "
          "$m_{5c}=M(\\pi^+\\pi^-e^+e^-\\gamma)$ with the photon energy rescaled so that $M(e^+e^-\\gamma)=m_{\\pi^0}$. "
          "$T$ = content of the frozen Gaussian (0.5463~GeV/$c^2$, 6.2~MeV/$c^2$) + pol2 fit of $N-\\mathrm{CB}$ in 0.49--0.62 $\\pm$ stat $\\pm$ model spread; "
          "$S/B$ in $\\pm1.9\\sigma$; width from the free fit. $Y_T=T\\times\\frac{64}{63}/L$ in visible rows per nb$^{-1}$, $L=%.1f$~nb$^{-1}$ (uncorrected; luminosity uncertainty $\\pm%.1f$\\%% not included). "
          "$\\pi^0$ of the tagged $\\eta$ in $M(e^+e^-\\gamma)$: $%.4f\\pm%.4f$~GeV/$c^2$, $\\sigma=%.1f\\pm%.1f$~MeV/$c^2$ (photon energy not calibrated).}\n",
          kLdata, 100 * kLrelTot, pi0Par[1], pi0Err[1], 1e3 * pi0Par[2], 1e3 * pi0Err[2]);
        t("\\label{tab:res_tagged}\n\\begin{tabular}{llcccc}\n\\toprule\n");
        t("Selection & isBest & $T\\pm$stat$\\pm$model & $S/B$ & width [MeV/$c^2$] & $Y_T$ [rows/nb$^{-1}$] \\\\\n\\midrule\n");
        for (int v = 0; v < nV; ++v) {
            if (v == 2) t("\\midrule\n\\multicolumn{6}{l}{\\emph{with OA$(e^+e^-)>2^\\circ$}} \\\\\n");
            for (int b = 0; b < nB; ++b) {
                const int k = b * nV + v;
                t("%s & %s & $%.0f\\pm%.0f\\pm%.0f$ & %.2f & ", bs[b].tex.c_str(), (v & 1) ? "$=1$" : "$\\geq0$", T[k].S, T[k].eS, Tsp[k], T[k].SoB);
                if (Tff[k].ok) t("$%.1f\\pm%.1f$", 1e3 * Tff[k].sigma, 1e3 * Tff[k].eSigma); else t("n/a");
                t(" & $%.3f\\pm%.3f$ \\\\\n", Yield(T[k].S), Yield(T[k].eS));
            }
        }
        t("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    }

    // ------------------------------------------------------------------------
    // figures
    // ------------------------------------------------------------------------
    TFile fr(Form("%s/res_bump.root", kWork), "RECREATE");
    {
        // [review fix 2026-11-14] 4 columns x 3 rows (larger text): row 1 pre, row 2 bump nominal + vertex, row 3 the two control
        // slices and a key; panel letters (a)-(d) pre, (e)-(h) bump nominal + vertex, (i), (j) control slices
        TCanvas* c = new TCanvas("c_sl", "", 2400, 2200);
        const int bsel[2] = {0, 2};
        const int isl[5] = {0, 1, 2, 3, 5};
        TGraphErrors* gK = nullptr; TGraph* gfK = nullptr; TGraph* g0K = nullptr;
        for (int row = 0; row < 2; ++row)
            for (int col = 0; col < 5; ++col) {
                const int b = bsel[row], i = isl[col], k = b * nV;
                const EdgeRes& r = R[k];
                const int pr = col < 4 ? row : 2, pc = col < 4 ? col : row;
                const char letter = col < 4 ? (char)('a' + row * 4 + col) : (char)('i' + row);
                c->cd();
                TPad* p = new TPad(Form("ps_%d_%d", row, col), "", pc / 4.0, 1.0 - (pr + 1) / 3.0, (pc + 1) / 4.0, 1.0 - pr / 3.0);
                p->SetLeftMargin(0.22); p->SetRightMargin(0.04); p->SetTopMargin(0.11); p->SetBottomMargin(0.16);
                p->Draw(); p->cd();
                const Spectra& s = sp[k * per + i];
                if (!r.ok[i]) continue;
                const StepFit& f = r.f[i];
                const double x1 = f.lo - 0.01, x2 = f.hi + 0.01;
                TGraphErrors* g = new TGraphErrors();
                double mn = 1e30, mx = -1e30;
                for (int q = 0; q < s.nbin; ++q) {
                    const double x = s.x(q);
                    if (x < x1 || x > x2) continue;
                    const double y = 1e-3 * s.sig[q], e = 1e-3 * std::sqrt(s.sigVar[q]);
                    const int n = g->GetN(); g->SetPoint(n, x, y); g->SetPointError(n, 0, e);
                    mn = std::min(mn, y - e); mx = std::max(mx, y + e);
                }
                TH1D* fh = new TH1D(Form("fs_%d_%d", row, col), "", 10, x1, x2); fh->SetDirectory(nullptr); fh->SetStats(0);
                fh->SetMinimum(mn - 0.05 * (mx - mn)); fh->SetMaximum(mx + 0.70 * (mx - mn));
                fh->GetXaxis()->SetTitle("M(#pi^{+}#pi^{#minus}e^{+}e^{#minus}) [GeV/c^{2}]");
                fh->GetYaxis()->SetTitle("10^{3} counts / 2.5 MeV/c^{2}");
                fh->GetXaxis()->SetNdivisions(505); fh->GetYaxis()->SetNdivisions(505);
                for (TAxis* a : {fh->GetXaxis(), fh->GetYaxis()}) { a->SetLabelSize(0.062); a->SetTitleSize(0.060); }
                fh->GetYaxis()->SetTitleOffset(1.55); fh->GetXaxis()->SetTitleOffset(1.10);
                fh->Draw("axis");
                g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetMarkerColor(kBlue + 1); g->SetLineColor(kBlue + 1);
                TGraph* gf = new TGraph(); TGraph* g0 = new TGraph();
                for (size_t q = 0; q < f.tx.size(); ++q) { gf->SetPoint((int)q, f.tx[q], 1e-3 * f.fy[q]); g0->SetPoint((int)q, f.tx[q], 1e-3 * f.fy0[q]); }
                gf->SetLineColor(kMagenta + 1); gf->SetLineWidth(3); g0->SetLineColor(kGray + 2); g0->SetLineStyle(2); g0->SetLineWidth(2);
                g0->Draw("l same"); g->Draw("pz same"); gf->Draw("l same");
                gK = g; gfK = gf; g0K = g0;
                TLine l; l.SetLineColor(kGreen + 2); l.SetLineStyle(7); l.SetLineWidth(2);
                l.DrawLine(f.edge, fh->GetMinimum(), f.edge, fh->GetMinimum() + 0.58 * (fh->GetMaximum() - fh->GetMinimum()));
                Txt(0.22, 0.915, Form("(%c) %s, M_{#pi#pi} %s", letter, row == 0 ? "base" : "+ #alpha_{lab}<40#circ, #pi vtx", kSlName[i]), 0.047);
                Txt(0.25, 0.80, Form("h = %.0f #pm %.0f", f.h, f.eh), 0.066, kBlue + 1);
                Txt(0.25, 0.725, Form("#chi^{2}/ndf = %.1f/%d", f.chi2, f.ndf), 0.056);
                Txt(0.25, 0.655, Form("#Delta#chi^{2} = %.1f", f.dchi2), 0.056);
                if (i == 5) Txt(0.58, 0.655, "control", 0.056, kRed + 1);
                if (row == 1 && i == 3 && TMath::Prob(f.chi2, f.ndf) < 0.05) Txt(0.58, 0.655, Form("p = %.3f", TMath::Prob(f.chi2, f.ndf)), 0.056, kRed + 1);
            }
        c->cd();
        TPad* pk = new TPad("ps_key", "", 0.5, 0.0, 1.0, 1.0 / 3.0);
        pk->SetMargin(0.02, 0.02, 0.02, 0.02); pk->Draw(); pk->cd();
        TLegend* lg = new TLegend(0.04, 0.50, 0.60, 0.94);
        lg->SetTextSize(0.060); lg->SetFillStyle(0);
        lg->AddEntry(gK, "N #minus CB (data, isBest #geq 0)", "pe");
        lg->AddEntry(gfK, "fit: pol2 [+ #eta Gaussian] + h #times step", "l");
        lg->AddEntry(g0K, "the same fit without its step", "l");
        TLine* le = new TLine(0, 0, 1, 1); le->SetLineColor(kGreen + 2); le->SetLineStyle(7); le->SetLineWidth(2);
        lg->AddEntry(le, "predicted kinematic edge (fixed)", "l");
        lg->Draw();
        Txt(0.05, 0.38, "(a)#minus(d): base selection; (e)#minus(h): + #alpha_{lab}<40#circ, #pi vtx fit;", 0.052);
        Txt(0.05, 0.30, "(i), (j): control slices M_{#pi#pi} 0.42#minus0.46 above the", 0.052);
        Txt(0.05, 0.22, "kinematic limit (edge set by hand, no step expected).", 0.052);
        Txt(0.05, 0.14, "h in counts per 2.5 MeV/c^{2}; edge smearing 5.4 MeV/c^{2}; CB with k = 1", 0.052);
        Save(c, "res_bump_slices");
    }
    {
        TCanvas* c = new TCanvas("c_tag", "", 2000, 800);
        const Spectra& s = st[0];
        const PeakFit& f = T[0];
        for (int col = 0; col < 2; ++col) {
            c->cd();
            TPad* p = new TPad(Form("pt_%d", col), "", col / 2.0, 0, (col + 1) / 2.0, 1);
            p->SetLeftMargin(0.15); p->SetRightMargin(0.03); p->SetTopMargin(0.10); p->SetBottomMargin(0.14);
            p->Draw(); p->cd();
            if (col == 0) {
                TH1D* hN = s.Hist("us", "tag_N"); TH1D* hC = s.Hist("cb", "tag_CB"); TH1D* hS = s.Hist("sig", "tag_sig");
                hN->SetMarkerStyle(20); hN->SetMarkerSize(0.8); hN->SetMarkerColor(kBlack); hN->SetLineColor(kBlack);
                hC->SetMarkerStyle(20); hC->SetMarkerSize(0.8); hC->SetMarkerColor(kRed + 1); hC->SetLineColor(kRed + 1);
                hS->SetMarkerStyle(24); hS->SetMarkerSize(0.8); hS->SetMarkerColor(kBlue + 1); hS->SetLineColor(kBlue + 1);
                // wide mass spectrum: log y (style rule; the CB is visible only on a log scale)
                // [verifier fix 2026-11-14: was linear with minimum 0]
                p->SetLogy();
                hN->GetXaxis()->SetRangeUser(0.42, 0.78);
                hN->SetMaximum(6.0 * hN->GetMaximum()); hN->SetMinimum(3);
                hN->GetXaxis()->SetTitle("m_{5c} = M(#pi^{+}#pi^{#minus}e^{+}e^{#minus}#gamma), E_{#gamma} rescaled [GeV/c^{2}]");
                hN->GetYaxis()->SetTitle("counts / 2.5 MeV/c^{2}");
                hN->GetXaxis()->SetNdivisions(505);
                hN->Draw("pe"); hS->Draw("pe same"); hC->Draw("pe same");
                TLegend* lg = new TLegend(0.18, 0.68, 0.50, 0.88); lg->SetTextSize(0.042);   // top left: empty on the log scale
                lg->AddEntry(hN, "N (e^{+}e^{#minus} counts)", "pe"); lg->AddEntry(hC, "CB", "pe"); lg->AddEntry(hS, "N #minus CB", "pe");
                lg->Draw();
                Txt(0.15, 0.93, "(a) ECAL-tagged counts, pre, isBest #geq 0", 0.046);
                hN->Write(); hC->Write(); hS->Write();
            } else {
                TGraphErrors* g = new TGraphErrors();
                double mn = 1e30, mx = -1e30;
                for (int q = 0; q < s.nbin; ++q) {
                    const double x = s.x(q);
                    if (x < 0.47 || x > 0.64) continue;
                    const int n = g->GetN(); g->SetPoint(n, x, s.sig[q]); g->SetPointError(n, 0, std::sqrt(s.sigVar[q]));
                    mn = std::min(mn, s.sig[q] - std::sqrt(s.sigVar[q])); mx = std::max(mx, s.sig[q] + std::sqrt(s.sigVar[q]));
                }
                TH1D* fh = new TH1D("ft", "", 10, 0.47, 0.64); fh->SetDirectory(nullptr); fh->SetStats(0);
                fh->SetMinimum(0); fh->SetMaximum(mx * 1.35);
                fh->GetXaxis()->SetTitle("m_{5c} [GeV/c^{2}]"); fh->GetYaxis()->SetTitle("N#minusCB [counts / 2.5 MeV/c^{2}]");
                fh->GetXaxis()->SetNdivisions(505);
                fh->Draw("axis");
                g->SetMarkerStyle(20); g->SetMarkerSize(0.9); g->SetMarkerColor(kBlue + 1); g->SetLineColor(kBlue + 1);
                TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
                for (int q = 0; q <= 200; ++q) {
                    const double x = kTagLo + q * (kTagHi - kTagLo) / 200;
                    gt->SetPoint(q, x, f.Total(x)); gb->SetPoint(q, x, f.Bkg(x));
                }
                gt->SetLineColor(kMagenta + 1); gt->SetLineWidth(3); gb->SetLineColor(kGray + 2); gb->SetLineStyle(2); gb->SetLineWidth(2);
                gb->Draw("l same"); g->Draw("pz same"); gt->Draw("l same");
                Txt(0.15, 0.93, "(b) N #minus CB: frozen Gaussian + pol2", 0.046);
                Txt(0.19, 0.83, Form("T = %.0f #pm %.0f #pm %.0f", f.S, f.eS, Tsp[0]), 0.046, kBlue + 1);
                Txt(0.19, 0.77, Form("S/B = %.2f, #chi^{2}/ndf = %.1f/%d", f.SoB, f.chi2, f.ndf), 0.042);
                if (Tff[0].ok) Txt(0.19, 0.71, Form("free: %.4f GeV/c^{2}, #sigma = %.1f #pm %.1f MeV/c^{2}", Tff[0].mean, 1e3 * Tff[0].sigma, 1e3 * Tff[0].eSigma), 0.042);
            }
        }
        Save(c, "res_bump_tagged");
    }
    fr.Close();
    printf("res_bump: done in %.0f s\n", swAll.RealTime());
}

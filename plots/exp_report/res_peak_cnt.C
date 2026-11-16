// res_peak_cnt.C -- COPY of plots/exp_report/res_peak.C made by figure agent A (integ4, 2026-11-15) with
// the figure strings 'rows' -> 'counts' (analyst's terminology) and the outputs of the copy redirected:
// figures of the scope -> docs/latex/figs_exp/<base>_cnt.{pdf,png}; everything else -> /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/res_peak
// Built by /tmp/claude-1002/-home-przygoda-HADES-FAT-FAT/846dd0d5-f6e8-4a0f-b2f4-a6bcc50f1087/scratchpad/final/integ4/figA/work/mk_cnt.py; the physics code is unchanged.
// ============================================================================
// res_peak.C  --  eta peak in M(pi+pi-e+e-): final numbers without simulation
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22, 29 days), channel pi+ pi- e+ e-. EXPERIMENTAL
// DATA ONLY. Chapter "Results without simulation" (prefix res), tasks 1 and 5.
//
// What it does
//   1. eta peak for the final candidate selections (pre, pre + vertex, peak
//      nominal, peak nominal + vertex, DIRECT, PEAK, X17 region), each for
//      isBest >= 0 (main) and isBest == 1 (comparison), with and without
//      OA(e+e-) > 2 deg:
//        P_main (chi2 fit of N - CB, frozen Gaussian 0.5439 GeV/c2 / 5.4 MeV +
//        pol2, 0.5225-0.6100, S = whole Gaussian), model spread (half range of
//        S over the 48 variations with p > 0.05), Poisson fit of N (control),
//        like-sign null tests, CB x 1.03 / x 0.97, free fit (position and
//        width), visible yield per luminosity Y = S x 64/63 / L (nb);
//   2. guard against sculpting: fake peak content that each selection makes
//      from a smooth background (pass fraction of mixed rows = pions of one
//      event with the lepton pair of another, as in plots/scan_peak.C);
//   3. mass scale: eta free fit / m_eta and K0s -> pi+pi- in M(pi+pi-)
//      (displaced pion crossing point, as in plots/scan_reference.C);
//   4. ECAL-tagged eta T (frozen definition of plots/scan_bump.C) and the ratio
//      peak content / tagged content (data only, uncorrected);
//   5. per-period stability (periods I, IIA, IIB, III) of the tagged yield per
//      nb^-1 and of the ratio peak / tagged (the peak yield per period itself is
//      in the samples chapter, plots/exp_report/smp_*; the numbers agree exactly).
//
// Selections (rows = stored combinations; cor_ok == 1; in FAT already vertex
// z > -500 mm and PT3 trigbit == 8192):
//   pre                 = cor_ok && isBest >= 0 && vz > -160 mm
//   pre + vertex        = pre + both pions in the vertex fit (pip_isUsedVertex && pim_isUsedVertex)
//   peak nominal        = pre + oa_lab < 55 deg + M(pipi) > 0.30 GeV/c2
//   peak nominal+vertex = peak nominal + both pions in the vertex fit
//   DIRECT              = pre + dz_pi < 40 mm           (frozen in the simulation study)
//   PEAK                = pre + dz_pi < 40 mm + MM4 < 2.85 GeV/c2 (frozen in the simulation study)
//   X17 region          = pre + 15.5 <= M(ee) < 24.5 MeV/c2 + both pions in the vertex fit
//                         (the frozen X17 selection without its M4 window; the peak is fitted)
//   dz_pi = max(|pip_z - vz|, |pim_z - vz|). Windows are the DATA windows.
// Luminosity: L = 6079.2 nb^-1 (runs present in the data), per day from the
// per-run table (res_common.h, section F).
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root (scan::kCompact)
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/res/lumi_coverage_perrun.csv
// Outputs: docs/latex/tables_exp/res_peak.txt, res_peak.tex, res_scale.tex, res_period.tex
//          docs/latex/figs_exp/res_peak_overview.pdf/.png, res_peak_fits.pdf/.png, res_tagged_period.pdf/.png
//          /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/res/res_peak.root
//
// Run (from /home/przygoda/HADES/FAT/FAT, about 2 min, 6 GB memory):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_res", kTRUE)' 'plots/exp_report/res_peak.C+'
// ============================================================================
#include "TBox.h"
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

#include "res_common_cnt.h"

namespace respk {

using namespace res;

// ---- mixed rows (copy of the method of plots/scan_peak.C) --------------------
// (pi+pi-) of one unlike-sign row with (e+e-) of a row of another event group,
// 20 fixed offsets in the row order; M4, oa_lab and MM4 recomputed from the
// corrected momenta; OA and M(ee) taken from the lepton row, the pion and vertex
// variables from the pion row. A mixed row has no eta.
const int kNOff = 20;
const int kOff[kNOff] = {1009, 2003, 3001, 4001, 5003, 6007, 7001, 8009, 9001, 10007,
                         11003, 12007, 13001, 14009, 15013, 16001, 17011, 18013, 19001, 20011};
struct Pool {
    const std::vector<Row>* rows = nullptr;
    std::vector<int> idx;
    std::vector<float> pp, ee;
    double pb = 0, eIni = 0;
};
void BuildPool(const Data& d, Pool& mp) {
    const double mPi = 1e3 * scan::kMpi, mE = 1e3 * scan::kMe, mP = 1e3 * scan::kMp, d2r = TMath::DegToRad();
    const double eb = scan::kTbeam + mP;
    mp.pb = std::sqrt(eb * eb - mP * mP); mp.eIni = eb + mP;
    mp.rows = &d.s[0];
    const std::vector<Row>& v = d.s[0];
    for (size_t i = 0; i < v.size(); ++i) if (v[i].cor_ok && Pre(v[i])) mp.idx.push_back((int)i);
    const size_t n = mp.idx.size();
    mp.pp.resize(4 * n); mp.ee.resize(4 * n);
    auto add = [&](float* o, double p, double th, double ph, double m) {
        const double st = std::sin(th * d2r);
        o[0] += (float)(p * st * std::cos(ph * d2r)); o[1] += (float)(p * st * std::sin(ph * d2r));
        o[2] += (float)(p * std::cos(th * d2r)); o[3] += (float)std::sqrt(p * p + m * m);
    };
    for (size_t k = 0; k < n; ++k) {
        const Row& r = v[mp.idx[k]];
        float* a = &mp.pp[4 * k]; float* b = &mp.ee[4 * k];
        a[0] = a[1] = a[2] = a[3] = 0; b[0] = b[1] = b[2] = b[3] = 0;
        add(a, r.pip_p, r.pip_theta, r.pip_phi, mPi); add(a, r.pim_p, r.pim_theta, r.pim_phi, mPi);
        add(b, r.ep_p, r.ep_theta, r.ep_phi, mE); add(b, r.em_p, r.em_theta, r.em_phi, mE);
    }
}
inline void Kin(const Pool& mp, size_t k, size_t j, double& m4, double& oal, double& mm) {
    const float* a = &mp.pp[4 * k]; const float* b = &mp.ee[4 * j];
    const double x = (double)a[0] + b[0], y = (double)a[1] + b[1], z = (double)a[2] + b[2], e = (double)a[3] + b[3];
    const double p2 = x * x + y * y + z * z, m2 = e * e - p2;
    m4 = 1e-3 * (m2 > 0 ? std::sqrt(m2) : -std::sqrt(-m2));
    const double pa = std::sqrt((double)a[0] * a[0] + (double)a[1] * a[1] + (double)a[2] * a[2]);
    const double pc = std::sqrt((double)b[0] * b[0] + (double)b[1] * b[1] + (double)b[2] * b[2]);
    const double cs = ((double)a[0] * b[0] + (double)a[1] * b[1] + (double)a[2] * b[2]) / (pa * pc);
    oal = std::acos(std::max(-1., std::min(1., cs))) * TMath::RadToDeg();
    const double q2 = (mp.eIni - e) * (mp.eIni - e) - x * x - y * y - (mp.pb - z) * (mp.pb - z);
    mm = 1e-3 * (q2 > 0 ? std::sqrt(q2) : -std::sqrt(-q2));
}
// pass counts of the selections on the mixed rows in the fit range; halves = offset parity
struct MixSpec { std::vector<std::vector<double>> pas[2]; std::vector<std::vector<double>> tot[2]; };
MixSpec MixedFill(const Pool& mp, const std::vector<Sel>& sels, const std::vector<Sel>& bases, long& nmix) {
    MixSpec ms;
    for (int h = 0; h < 2; ++h) {
        ms.pas[h].assign(sels.size(), std::vector<double>(scan::kNbin, 0.));
        ms.tot[h].assign(bases.size(), std::vector<double>(scan::kNbin, 0.));
    }
    const std::vector<Row>& v = *mp.rows;
    const size_t n = mp.idx.size();
    nmix = 0;
    double m4, oal, mm;
    for (size_t k = 0; k < n; ++k) {
        const Row& a = v[mp.idx[k]];
        for (int o = 0; o < kNOff; ++o) {
            const size_t j = (k + kOff[o]) % n;
            const Row& b = v[mp.idx[j]];
            if (a.grp == b.grp) continue;
            Kin(mp, k, j, m4, oal, mm);
            if (!(m4 >= scan::kFitLo && m4 < scan::kFitHi)) continue;
            const int ib = (int)((m4 - scan::kMlo) / scan::kBinW);
            Row m = a;
            m.oa = b.oa; m.mee = b.mee; m.m4 = (float)m4; m.oa_lab = (float)oal; m.mm4 = (float)mm;
            const int h = o & 1;
            nmix++;
            for (size_t s = 0; s < bases.size(); ++s) if (bases[s](m)) ms.tot[h][s][ib] += 1;
            for (size_t s = 0; s < sels.size(); ++s) if (sels[s](m)) ms.pas[h][s][ib] += 1;
        }
    }
    return ms;
}
// Gaussian (frozen shape) + pol2 fitted by linear least squares to y with variances w in the fit range
double GaussCoef(const std::vector<double>& y, const std::vector<double>& w, std::vector<double>* coef = nullptr) {
    std::vector<double> F, yy, var;
    const double x0 = 0.5 * (scan::kFitLo + scan::kFitHi), dx = 0.5 * (scan::kFitHi - scan::kFitLo);
    int npts = 0;
    for (int i = 0; i < scan::kNbin; ++i) {
        const double x1 = scan::kMlo + i * scan::kBinW, xc = x1 + 0.5 * scan::kBinW;
        if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
        const double t = (xc - x0) / dx;
        F.push_back(scan::GaussFrac(x1, x1 + scan::kBinW, scan::kPeakMean, scan::kPeakSigma));
        F.push_back(1.); F.push_back(t); F.push_back(t * t);
        yy.push_back(y[i]); var.push_back(std::max(w[i], 1e-9));
        npts++;
    }
    const scan::LinRes r = scan::LinFit(npts, 4, yy.data(), var.data(), F.data());
    if (!r.ok) return 0;
    if (coef) {
        coef->assign(scan::kNbin, 0.);
        int q = 0;
        for (int i = 0; i < scan::kNbin; ++i) {
            const double xc = scan::kMlo + (i + 0.5) * scan::kBinW;
            if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
            double a = 0;
            for (int l = 0; l < 4; ++l) a += r.cov[0][l] * F[q * 4 + l] / var[q];
            (*coef)[i] = a; q++;
        }
    }
    return r.par[0];
}
// y = smooth background of the base (pol2 of the Poisson fit of N) x pass fraction of the mixed rows
double FakeOne(const PeakFit& smooth, const std::vector<double>& pas, const std::vector<double>& tot, double* ebin = nullptr) {
    std::vector<double> y(scan::kNbin, 0.), w(scan::kNbin, 1.), ef(scan::kNbin, 0.), coef;
    for (int i = 0; i < scan::kNbin; ++i) {
        const double xc = scan::kMlo + (i + 0.5) * scan::kBinW;
        if (xc < scan::kFitLo || xc >= scan::kFitHi) continue;
        const double b = smooth.Bkg(xc);
        const double f = tot[i] > 0 ? pas[i] / tot[i] : 0;
        y[i] = b * f; w[i] = std::max(y[i], 1.);
        ef[i] = tot[i] > 0 ? b * std::sqrt(std::max(f * (1 - f), 0.) / tot[i]) : 0;
    }
    const double s = GaussCoef(y, w, &coef);
    if (ebin) {
        double e2 = 0;
        for (int i = 0; i < scan::kNbin; ++i) e2 += coef[i] * coef[i] * ef[i] * ef[i];
        *ebin = std::sqrt(e2);
    }
    return s;
}
std::vector<double> AddV(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> c(a.size());
    for (size_t i = 0; i < a.size(); ++i) c[i] = a[i] + b[i];
    return c;
}

// ---- one measurement ---------------------------------------------------------
struct Meas {
    std::string sel; int var = 0;
    Spectra sp;
    PeakFit pm, pp, lsP, lsM;
    scan::Spread sd;
    scan::FreeFit ff;
    double dCBup = 0, dCBdn = 0;       // S(CB x 1.03) - S, S(CB x 0.97) - S
    double fake = 0, eFake = 0; bool hasFake = false;
};

// drawing helpers
void Txt(double x, double y, const char* t, double size = 0.050, int col = kBlack, int align = 11) {
    TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size); l.SetTextColor(col); l.SetTextAlign(align);
    l.DrawLatex(x, y, t);
}
TH1D* Frame(const char* name, double x0, double x1, double y0, double y1, const char* xt, const char* yt) {
    TH1D* f = new TH1D(name, "", 10, x0, x1); f->SetDirectory(nullptr);
    f->SetMinimum(y0); f->SetMaximum(y1); f->SetStats(0);
    f->GetXaxis()->SetTitle(xt); f->GetYaxis()->SetTitle(yt);
    f->GetXaxis()->SetNdivisions(505);
    return f;
}

} // namespace respk

void res_peak_cnt() {
    using namespace respk;
    TStopwatch swAll;
    gSystem->mkdir(kFigOut, kTRUE); gSystem->mkdir(kTabOut, kTRUE); gSystem->mkdir(kWork, kTRUE); gSystem->mkdir(kCntScratch, kTRUE);
    Style();

    Data d;
    if (!scan::Load(d)) return;

    std::map<int, double> Lday;
    if (!LumiPerDay(Lday)) return;
    double Lper[kNPer] = {0, 0, 0, 0}, Lsum = 0;
    for (const auto& e : Lday) { const int p = Period(e.first); if (p >= 0) Lper[p] += e.second; Lsum += e.second; }

    // ------------------------------------------------------------------------
    // 1. peak fits for all selections and variants
    // ------------------------------------------------------------------------
    std::vector<NSel> sl = {SelPre(), SelPreVtx(), SelPeakNom(), SelPeakNomVtx(), SelDIRECT(), SelPEAK(), SelX17()};
    const int nS = (int)sl.size(), nV = 4;
    std::vector<Sel> sels;
    for (int s = 0; s < nS; ++s) for (int v = 0; v < nV; ++v) sels.push_back(Variant(sl[s].f, v));
    // per period: pre and peak nominal + vertex (isBest >= 0); the rows of all M4 per period
    const int iPer0 = (int)sels.size();
    for (int p = 0; p < kNPer; ++p) {
        sels.push_back([p](const Row& r) { return Pre(r) && Period(r.file_id) == p; });
        sels.push_back([p](const Row& r) { return Pre(r) && PeakNomKin(r) && PiVtx(r) && Period(r.file_id) == p; });
    }
    TStopwatch sw;
    const std::vector<Spectra> sp = scan::FillMany(d, sels);
    printf("FillMany of %zu selections: %.1f s\n", sels.size(), sw.RealTime());

    std::vector<Meas> M(nS * nV);
    for (int s = 0; s < nS; ++s)
        for (int v = 0; v < nV; ++v) {
            Meas& m = M[s * nV + v];
            m.sel = sl[s].key; m.var = v; m.sp = sp[s * nV + v];
            m.pm = scan::FitPeakMain(m.sp);
            m.pp = scan::FitPeakPois(m.sp);
            m.lsP = scan::FitPeakLikeSign(m.sp, 1);
            m.lsM = scan::FitPeakLikeSign(m.sp, 2);
            m.sd = scan::PeakModelSpread(m.sp);
            m.ff = scan::FitPeakFree(m.sp);
            m.dCBup = scan::FitPeakMain(ScaleCB(m.sp, 1.03)).S - m.pm.S;
            m.dCBdn = scan::FitPeakMain(ScaleCB(m.sp, 0.97)).S - m.pm.S;
        }

    // ------------------------------------------------------------------------
    // 2. guard against sculpting (isBest >= 0; variants without and with OA > 2)
    // ------------------------------------------------------------------------
    long nmix = 0; double mixTime = 0;
    {
        TStopwatch sm;
        Pool pool; BuildPool(d, pool);
        std::vector<Sel> ms, bases = {[](const Row& r) { return Pre(r); }, [](const Row& r) { return Pre(r) && r.oa > 2.f; }};
        for (int s = 0; s < nS; ++s) { ms.push_back(Variant(sl[s].f, 0)); ms.push_back(Variant(sl[s].f, 2)); }
        const MixSpec mx = MixedFill(pool, ms, bases, nmix);
        const PeakFit smooth[2] = {scan::FitPeakPois(M[0].sp), scan::FitPeakPois(M[2].sp)};   // pre, pre + OA > 2
        for (int s = 0; s < nS; ++s)
            for (int q = 0; q < 2; ++q) {
                const int is = 2 * s + q;
                double eb = 0;
                const double f = FakeOne(smooth[q], AddV(mx.pas[0][is], mx.pas[1][is]), AddV(mx.tot[0][q], mx.tot[1][q]), &eb);
                const double a = FakeOne(smooth[q], mx.pas[0][is], mx.tot[0][q]), b = FakeOne(smooth[q], mx.pas[1][is], mx.tot[1][q]);
                Meas& m = M[s * nV + (q ? 2 : 0)];
                m.fake = f; m.eFake = std::max(eb, 0.5 * std::fabs(a - b)); m.hasFake = true;
            }
        mixTime = sm.RealTime();
        printf("mixed rows in the fit range: %ld (%.1f s)\n", nmix, mixTime);
    }

    // ------------------------------------------------------------------------
    // 3. K0s -> pi+pi- (displaced pion crossing point), mass scale of M(pi+pi-)
    // ------------------------------------------------------------------------
    struct K0 { double pos = 0, ePos = 0, w = 0, eW = 0, n = 0, eN = 0, chi2 = 0; int ndf = 0; };
    K0 k0[3];
    TH1D* hK[3];
    {
        const char* kn[3] = {"all rows (no vz cut)", "pre", "pre, isBest == 1"};
        for (int i = 0; i < 3; ++i) { hK[i] = new TH1D(Form("hK0s_%d", i), kn[i], 100, 0.42, 0.57); hK[i]->SetDirectory(nullptr); }
        for (const auto& r : d.us()) {
            if (!r.cor_ok || !(r.dz_pipi > 20.f && r.doca_pipi < 10.f)) continue;
            hK[0]->Fill(r.mpipi);
            if (Pre(r)) { hK[1]->Fill(r.mpipi); if (r.isBest == 1) hK[2]->Fill(r.mpipi); }
        }
        for (int i = 0; i < 3; ++i) {
            TF1 f(Form("fK%d", i), "gausn(0)*0.0015+pol2(3)", 0.455, 0.54);
            f.SetParameters(0.03 * hK[i]->Integral(), 0.496, 0.006, hK[i]->GetBinContent(hK[i]->FindBin(0.47)), 0, 0);
            f.SetParLimits(1, 0.485, 0.505); f.SetParLimits(2, 0.003, 0.012);
            TFitResultPtr rr = hK[i]->Fit(&f, "RQSN0L");
            k0[i].pos = rr->Parameter(1); k0[i].ePos = rr->ParError(1); k0[i].w = rr->Parameter(2); k0[i].eW = rr->ParError(2);
            k0[i].n = rr->Parameter(0); k0[i].eN = rr->ParError(0); k0[i].chi2 = rr->Chi2(); k0[i].ndf = rr->Ndf();
        }
    }

    // ------------------------------------------------------------------------
    // 4. ECAL-tagged eta and 5. per period
    // ------------------------------------------------------------------------
    std::vector<Sel> ts = {[](const Row& r) { return Tagged(r) && Pre(r); },
                           [](const Row& r) { return Tagged(r) && Pre(r) && r.isBest == 1; },
                           [](const Row& r) { return Tagged(r) && Pre(r) && r.oa > 2.f; },
                           [](const Row& r) { return Tagged(r) && Pre(r) && r.isBest == 1 && r.oa > 2.f; }};
    for (int p = 0; p < kNPer; ++p) ts.push_back([p](const Row& r) { return Tagged(r) && Pre(r) && Period(r.file_id) == p; });
    const std::vector<Spectra> st = scan::FillMany(d, ts, TagOpt());
    PeakFit tf[4]; double tSpread[4]; int tKept[4]; scan::FreeFit tFree[4];
    for (int v = 0; v < 4; ++v) { tf[v] = scan::FitPeakMain(st[v], TagCfg()); tSpread[v] = TagSpread(st[v], &tKept[v]); tFree[v] = scan::FitPeakFree(st[v], TagCfg()); }
    PeakFit sPer[2][kNPer], tPer[kNPer];
    double rowsPer[kNPer];
    for (int p = 0; p < kNPer; ++p) {
        sPer[0][p] = scan::FitPeakMain(sp[iPer0 + 2 * p]);
        sPer[1][p] = scan::FitPeakMain(sp[iPer0 + 2 * p + 1]);
        tPer[p] = scan::FitPeakMain(st[4 + p], TagCfg());
        rowsPer[p] = sp[iPer0 + 2 * p].nsel[0];
    }

    // ------------------------------------------------------------------------
    // tables
    // ------------------------------------------------------------------------
    OpenTxt("res_peak.txt");
    P("# res_peak.txt -- written by plots/exp_report/res_peak.C (%s)\n", TDatime().AsSQLString());
    P("# EXPERIMENTAL DATA, p+p 4.53 GeV (feb22, 29 days). Rows = stored combinations; cor_ok == 1; FAT: vertex z > -500 mm, PT3 (trigbit 8192).\n");
    P("# Selections:\n");
    for (int s = 0; s < nS; ++s) P("#   %-20s = %s\n", sl[s].key.c_str(), sl[s].def.c_str());
    P("#   variants: isBest >= 0 (all rows, main) / isBest == 1 (comparison); without / with OA(e+e-) > 2 deg (OA > 2 is NOT a fixed cut)\n");
    P("# P_main: chi2 fit of N - CB, Gaussian 0.5439 GeV/c2 / 5.4 MeV (frozen) + pol2 in 0.5225-0.6100, 2.5 MeV bins; S = whole Gaussian.\n");
    P("# CB = 2 sqrt(N++ N--) in 10 MeV bins, interpolated to 2.5 MeV; k = 1. N, CB = rows in the DATA peak window 0.53364-0.55416.\n");
    P("# spread = half of (max - min) of S over the 48 model variations with p > 0.05 (orders 1-4, ranges 0.520-0.590 ... 0.505-0.650, width +-0.5 MeV).\n");
    P("# S/B = S inside the window / (polynomial + CB in the window); F = S/stat (compares selections, NOT a significance).\n");
    P("# Y = S x 64/63 / L, L = %.1f nb^-1 (runs present in the data; all runs %.1f): visible rows per nb^-1 = nb, not corrected for acceptance or efficiency.\n", kLdata, kLall);
    P("# Luminosity uncertainty (separate, not in the errors of Y): +-%.1f %% total (HADES syst %.2f %%, Argonne fit %.2f %%, normalisation %.1f %%);\n",
      100 * kLrelTot, 100 * kLrelSyst, 100 * kLrelFit, 100 * kLrelNorm);
    P("#   one-sided: if the PT3 live fraction relative to PT2 (%.3f) applies, L is lower by up to %.1f %% and Y higher by up to %.1f %%.\n",
      kLlive, 100 * (1 - kLlive), 100 * (1 / kLlive - 1));
    P("# dS(CB+-3%%) = S(CB x 1.03) - S and S(CB x 0.97) - S. S_pois = Poisson fit of unsubtracted N (control). S(++), S(--) = like-sign null tests.\n");
    P("# S_fake = fake Gaussian content made by the selection from a smooth background (mixed rows, relative to pre or pre + OA > 2; isBest >= 0 only).\n");
    P("# Luminosity per period [nb^-1]: I %.2f, IIA %.2f, IIB %.2f, III %.2f; sum %.2f.\n\n", Lper[0], Lper[1], Lper[2], Lper[3], Lsum);

    P("A. ETA PEAK, P_main\n");
    P("%-20s %-16s | %7s %6s %6s | %6s %5s %6s %3s | %6s %5s %5s %9s | %6s %6s | %8s %6s %6s | %6s %4s | %5s %4s %5s %4s | %6s %5s\n",
      "selection", "variant", "N", "CB", "CB/N", "S", "stat", "spread", "fit", "S/B", "F", "F_tot", "chi2/ndf", "dS+3%", "dS-3%",
      "Y [nb]", "stat", "model", "S_pois", "err", "S(++)", "err", "S(--)", "err", "S_fake", "err");
    for (int v = 0; v < nV; ++v) {
        for (int s = 0; s < nS; ++s) {
            const Meas& m = M[s * nV + v];
            const double ftot = m.pm.S / std::sqrt(m.pm.eS * m.pm.eS + m.sd.half * m.sd.half);
            P("%-20s %-16s | %7.0f %6.0f %6.3f | %6.0f %5.0f %6.0f %3d | %6.4f %5.2f %5.2f %5.1f/%-3d | %+6.0f %+6.0f | %8.4f %6.4f %6.4f | %6.0f %4.0f | %5.0f %4.0f %5.0f %4.0f |",
              m.sel.c_str(), VariantName(v), m.pm.N, m.pm.CB, m.pm.N > 0 ? m.pm.CB / m.pm.N : 0., m.pm.S, m.pm.eS, m.sd.half, m.sd.npass,
              m.pm.SoB, m.pm.F, ftot, m.pm.chi2, m.pm.ndf, m.dCBup, m.dCBdn, Yield(m.pm.S), Yield(m.pm.eS), Yield(m.sd.half),
              m.pp.S, m.pp.eS, m.lsP.S, m.lsP.eS, m.lsM.S, m.lsM.eS);
            if (m.hasFake) P(" %+6.0f %5.0f\n", m.fake, m.eFake); else P(" %6s\n", "-");
        }
        P("\n");
    }
    P("   isBest >= 0 / isBest == 1 ratio of S (no OA cut): ");
    for (int s = 0; s < nS; ++s) P("%s %.3f; ", sl[s].key.c_str(), M[s * nV + 1].pm.S != 0 ? M[s * nV].pm.S / M[s * nV + 1].pm.S : 0.);
    P("\n   OA > 2 pass fraction of S (isBest >= 0): ");
    for (int s = 0; s < nS; ++s) P("%s %.3f; ", sl[s].key.c_str(), M[s * nV].pm.S != 0 ? M[s * nV + 2].pm.S / M[s * nV].pm.S : 0.);
    P("\n   mixed rows used for S_fake: %ld in the fit range (%.0f s); threshold used in scan_peak.C: |S_fake| > 0.25 stat.\n", nmix, mixTime);

    P("\nB. FREE FIT (P_main, position and width free, frozen range and pol2; errors from delta chi2 = 1, eS includes the shape)\n");
    P("%-20s %-16s | %8s %7s | %6s %5s | %6s %5s | %8s %7s | %9s\n", "selection", "variant", "position", "err", "width", "err", "S", "err",
      "scale", "err", "chi2/ndf");
    for (int v = 0; v < nV; ++v)
        for (int s = 0; s < nS; ++s) {
            const Meas& m = M[s * nV + v];
            if (!m.ff.ok) { P("%-20s %-16s | no minimum (free fit does not converge)\n", m.sel.c_str(), VariantName(v)); continue; }
            P("%-20s %-16s | %8.5f %7.5f | %6.2f %5.2f | %6.0f %5.0f | %8.5f %7.5f | %5.1f/%-3d\n", m.sel.c_str(), VariantName(v),
              m.ff.mean, m.ff.eMean, 1e3 * m.ff.sigma, 1e3 * m.ff.eSigma, m.ff.S, m.ff.eS, m.ff.mean / scan::kMeta, m.ff.eMean / scan::kMeta,
              m.ff.chi2, m.ff.ndf);
        }
    P("   frozen: position 0.5439 (free fit on oa_lab < 50, no vz cut: 0.54388 +- 0.00049, width 5.44 +- 0.60 MeV), M4 scale 0.99274 +- 0.00090.\n");
    P("\nC. K0s -> pi+pi- in M(pi+pi-) (pion crossing point displaced: dz_pipi > 20 mm, doca_pipi < 10 mm; gausn + pol2, 0.455-0.540, likelihood)\n");
    for (int i = 0; i < 3; ++i)
        P("   %-22s: position %.5f +- %.5f, width %.2f +- %.2f MeV, content %.0f +- %.0f, chi2/ndf %.1f/%d, scale M(pipi) %.5f +- %.5f\n",
          hK[i]->GetTitle(), k0[i].pos, k0[i].ePos, 1e3 * k0[i].w, 1e3 * k0[i].eW, k0[i].n, k0[i].eN, k0[i].chi2, k0[i].ndf,
          k0[i].pos / scan::kMK0, k0[i].ePos / scan::kMK0);
    P("   isBest >= 0 rows: one pion pair can enter several times (errors slightly too small). Frozen M(pipi) scale 0.9970 (reference 0.99691 +- 0.00036).\n");

    P("\nD. ECAL-TAGGED ETA (pi0 window 0.119 <= M(e+e-gamma) < 0.173, m_ee < 0.10; m5c, frozen 0.5463 / 6.2 MeV + pol2, 0.49-0.62; CB in 20 MeV bins)\n");
    P("%-24s | %6s %6s %6s | %6s %5s %6s %3s | %6s %5s %9s | %8s %6s %6s | %s\n", "selection", "rows", "N win", "CB", "T", "stat", "spread", "fit",
      "S/B", "S/err", "chi2/ndf", "Y [nb]", "stat", "model", "free fit: position, width [MeV], S");
    const char* tn[4] = {"pre, isBest>=0", "pre, isBest==1", "pre + OA>2, isBest>=0", "pre + OA>2, isBest==1"};
    for (int v = 0; v < 4; ++v) {
        P("%-24s | %6.0f %6.0f %6.0f | %6.0f %5.0f %6.0f %3d | %6.3f %5.1f %5.1f/%-3d | %8.4f %6.4f %6.4f |", tn[v], st[v].nsel[0], tf[v].N, tf[v].CB,
          tf[v].S, tf[v].eS, tSpread[v], tKept[v], tf[v].SoB, tf[v].F, tf[v].chi2, tf[v].ndf, Yield(tf[v].S), Yield(tf[v].eS), Yield(tSpread[v]));
        if (tFree[v].ok) P(" %.4f +- %.4f, %.2f +- %.2f, %.0f +- %.0f\n", tFree[v].mean, tFree[v].eMean, 1e3 * tFree[v].sigma, 1e3 * tFree[v].eSigma, tFree[v].S, tFree[v].eS);
        else P(" n/a\n");
    }
    double rST[4], erST[4];
    // overlap of the two samples: unlike-sign pre rows in the P_main range that are also tagged rows in the m5c fit range
    double nFit = 0, nTagFit = 0, nBoth = 0, nBothW = 0, nWin = 0;
    for (const auto& r : d.us()) {
        if (!r.cor_ok || !Pre(r)) continue;
        const bool a = r.m4 >= scan::kFitLo && r.m4 < scan::kFitHi;
        bool t = false;
        if (Tagged(r)) { const double m5 = M5c(r); t = m5 >= kTagLo && m5 < kTagHi; }
        if (a) nFit++;
        if (t) nTagFit++;
        if (a && t) nBoth++;
        if (r.m4 >= scan::kPeakLo && r.m4 < scan::kPeakHi) { nWin++; if (t) nBothW++; }
    }
    P("\n   overlap: pre rows in the P_main range %.0f, tagged rows in the m5c fit range %.0f, in both %.0f (%.2f %% / %.2f %%);\n", nFit, nTagFit, nBoth,
      100 * nBoth / nFit, 100 * nBoth / nTagFit);
    P("            rows in the peak window %.0f, of them tagged in the m5c range %.0f (%.2f %%)\n", nWin, nBothW, 100 * nBothW / nWin);
    P("\n   RATIO peak content / tagged content (same selection; uncorrected; the two fits are taken as independent, see the overlap above):\n");
    for (int v = 0; v < 4; ++v) {
        Ratio(M[v].pm.S, M[v].pm.eS, tf[v].S, tf[v].eS, rST[v], erST[v]);
        const double mod = std::sqrt(std::pow(M[v].sd.half / M[v].pm.S, 2) + std::pow(tSpread[v] / tf[v].S, 2)) * rST[v];
        P("   %-24s S = %.0f +- %.0f +- %.0f, T = %.0f +- %.0f +- %.0f, S/T = %.3f +- %.3f (stat) +- %.3f (model)\n", tn[v], M[v].pm.S, M[v].pm.eS,
          M[v].sd.half, tf[v].S, tf[v].eS, tSpread[v], rST[v], erST[v], mod);
    }

    P("\nE. PER PERIOD (isBest >= 0; P_main with the frozen shape per period, stat errors only)\n");
    P("%-6s %9s | %8s %7s | %6s %5s %8s %6s | %6s %5s %8s %6s | %6s %5s %8s %6s | %6s %5s\n", "period", "L [nb^-1]", "rows", "rows/nb",
      "S pre", "stat", "Y [nb]", "stat", "S nv", "stat", "Y [nb]", "stat", "T", "stat", "Y_T [nb]", "stat", "S/T", "stat");
    double yP[3][kNPer], eyP[3][kNPer], rP[kNPer], erP[kNPer];
    for (int p = 0; p < kNPer; ++p) {
        const PeakFit* f3[3] = {&sPer[0][p], &sPer[1][p], &tPer[p]};
        for (int q = 0; q < 3; ++q) { yP[q][p] = f3[q]->S * kPT3 / Lper[p]; eyP[q][p] = f3[q]->eS * kPT3 / Lper[p]; }
        Ratio(sPer[0][p].S, sPer[0][p].eS, tPer[p].S, tPer[p].eS, rP[p], erP[p]);
        P("%-6s %9.2f | %8.0f %7.2f | %6.0f %5.0f %8.4f %6.4f | %6.0f %5.0f %8.4f %6.4f | %6.0f %5.0f %8.4f %6.4f | %6.3f %5.3f\n", kPerName[p], Lper[p],
          rowsPer[p], rowsPer[p] / Lper[p], sPer[0][p].S, sPer[0][p].eS, yP[0][p], eyP[0][p], sPer[1][p].S, sPer[1][p].eS, yP[1][p], eyP[1][p],
          tPer[p].S, tPer[p].eS, yP[2][p], eyP[2][p], rP[p], erP[p]);
    }
    const char* qn[4] = {"Y peak (pre)", "Y peak (peak nominal + vertex)", "Y tagged (pre)", "S/T (pre)"};
    double wm[4], ewm[4], c2[4];
    for (int q = 0; q < 4; ++q) {
        double sw = 0, sx = 0;
        for (int p = 0; p < kNPer; ++p) {
            const double y = q < 3 ? yP[q][p] : rP[p], e = q < 3 ? eyP[q][p] : erP[p];
            if (e > 0) { sw += 1 / (e * e); sx += y / (e * e); }
        }
        wm[q] = sx / sw; ewm[q] = 1 / std::sqrt(sw); c2[q] = 0;
        for (int p = 0; p < kNPer; ++p) {
            const double y = q < 3 ? yP[q][p] : rP[p], e = q < 3 ? eyP[q][p] : erP[p];
            c2[q] += (y - wm[q]) * (y - wm[q]) / (e * e);
        }
        P("   %-32s weighted mean %.4f +- %.4f, chi2/ndf = %.2f/%d (p = %.3f)\n", qn[q], wm[q], ewm[q], c2[q], kNPer - 1, TMath::Prob(c2[q], kNPer - 1));
    }
    P("   rows/nb = all selected unlike-sign rows of pre (any M4) per nb^-1 (normalisation check). Sector 2 had a trigger deficit on days 034-044 (period I),\n");
    P("   not in L. The model spread of S is common to all periods and is not included here.\n");
    P("\nrun time %.0f s\n", swAll.RealTime());
    CloseTxt();

    // ---- LaTeX: main peak table ----
    {
        Tex t("res_peak.tex");
        t("%% res_peak.tex -- written by plots/exp_report/res_peak.C (numbers: tables_exp/res_peak.txt)\n");
        t("\\begin{table}[H]\n\\centering\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n");
        t("\\caption{The $\\eta$ peak in $M(\\pi^+\\pi^-e^+e^-)$ for the final candidate selections, experimental data. "
          "$P_\\mathrm{main}$: $\\chi^2$ fit of $N-\\mathrm{CB}$ with the frozen Gaussian (0.5439~GeV/$c^2$, $\\sigma=5.4$~MeV/$c^2$) + pol2 in 0.5225--0.6100~GeV/$c^2$; "
          "$S$ = whole Gaussian $\\pm$ stat $\\pm$ model spread (old envelope: half range over the 48 variations with $p>0.05$; the revised model uncertainty, the rms over an edge-safe grid of variations, is in Table~\\ref{tab:fix_spread}). "
          "stat: fit error, per-row Poisson (no event-group term; checked with independent halves, Table~\\ref{tab:fix_stat}). "
          "$N$, CB: rows in the data peak window 0.53364--0.55416~GeV/$c^2$, CB $=2\\sqrt{N_{++}N_{--}}$, $k=1$. "
          "$\\Delta S_\\mathrm{CB}$: change of $S$ for CB $\\times 1.03$. "
          "$Y=S\\times\\frac{64}{63}/L$ with $L=%.1f$~nb$^{-1}$: visible rows per nb$^{-1}$ (rows, not events; not a cross section), not corrected for acceptance or efficiency; "
          "luminosity uncertainty $\\pm%.1f$\\%% (and one-sided $+%.1f$\\%% on $Y$ if the PT3 live fraction applies) not included. "
          "Selections: pre $=$ cor\\_ok, $v_z>-160$~mm; pre + vertex $=$ pre + both pions in the event-vertex fit; peak nominal $=$ pre + $\\alpha_\\mathrm{lab}<55^\\circ$ + $M_{\\pi\\pi}>0.30$~GeV/$c^2$; "
          "DIRECT $=$ pre + $\\Delta z_\\pi<40$~mm; PEAK $=$ DIRECT + $MM_4<2.85$~GeV/$c^2$; X17 region $=$ pre + $15.5\\leq M_{ee}<24.5$~MeV/$c^2$ + both pions in the vertex fit. "
          "Rows: all stored combinations (isBest $\\geq0$) or isBest $=1$, corrected momenta, PT3, vertex $z>-500$~mm. "
          "$F=S/\\sigma_\\mathrm{stat}$ compares selections, it is not a significance.}\n", kLdata, 100 * kLrelTot, 100 * (1 / kLlive - 1));
        t("\\label{tab:res_peak}\n");
        t("\\begin{tabular}{llrrcrrrc}\n\\toprule\n");
        t("Selection & isBest & $N$ & CB/$N$ & $S\\pm$stat$\\pm$model & $S/B$ & $F$ & $\\Delta S_\\mathrm{CB}$ & $Y$ [rows/nb$^{-1}$] \\\\\n\\midrule\n");
        for (int v = 0; v < nV; ++v) {
            if (v == 2) t("\\midrule\n\\multicolumn{9}{l}{\\emph{with OA$(e^+e^-)>2^\\circ$ (not a fixed cut)}} \\\\\n");
            for (int s = 0; s < nS; ++s) {
                const Meas& m = M[s * nV + v];
                t("%s & %s & %.0f & %.3f & $%.0f\\pm%.0f\\pm%.0f$ & %.3f & %.1f & %+.0f & $%.3f\\pm%.3f\\pm%.3f$ \\\\\n", sl[s].tex.c_str(),
                  (v & 1) ? "$=1$" : "$\\geq0$", m.pm.N, m.pm.N > 0 ? m.pm.CB / m.pm.N : 0., m.pm.S, m.pm.eS, m.sd.half, m.pm.SoB, m.pm.F,
                  std::round(m.dCBup) + 0.0, Yield(m.pm.S), Yield(m.pm.eS), Yield(m.sd.half));
            }
            if (v == 0 || v == 2) t("\\addlinespace\n");
        }
        t("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    }
    // ---- LaTeX: mass scale, free fits, sculpting guard ----
    {
        Tex t("res_scale.tex");
        t("%% res_scale.tex -- written by plots/exp_report/res_peak.C (numbers: tables_exp/res_peak.txt)\n");
        t("\\begin{table}[H]\n\\centering\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n");
        t("\\caption{Free fit of the $\\eta$ peak (position and width free, $P_\\mathrm{main}$, 0.5225--0.6100~GeV/$c^2$, pol2; $S$ error includes the shape), "
          "the $M_4$ mass scale $=$ position$/m_\\eta$ ($m_\\eta=0.547862$~GeV/$c^2$), and the fake peak content $S_\\mathrm{fake}$ that each selection makes "
          "from a smooth background (pass fraction of mixed rows: pions of one event, lepton pair of another; relative to pre). "
          "Data, isBest $\\geq0$. Last rows: K$^0_S\\to\\pi^+\\pi^-$ in $M(\\pi^+\\pi^-)$ with a displaced pion crossing point "
          "($\\Delta z_{\\pi\\pi}>20$~mm, DOCA $<10$~mm); there the scale column is the $M(\\pi^+\\pi^-)$ scale $=$ position$/m_{K^0}$.}\n");
        t("\\label{tab:res_scale}\n\\begin{tabular}{lccccc}\n\\toprule\n");
        t("Selection & position [GeV/$c^2$] & width [MeV/$c^2$] & $S$ (free) & $M_4$ scale & $S_\\mathrm{fake}$ \\\\\n\\midrule\n");
        for (int s = 0; s < nS; ++s) {
            const Meas& m = M[s * nV];
            // S_fake of pre is 0 by construction (pre is the reference of the pass fractions) [verifier fix: was printed as -0 +- 0]
            const std::string fk = s == 0 ? std::string("0 (ref.)") : std::string(Form("$%+.0f\\pm%.0f$", m.fake, m.eFake));
            if (m.ff.ok)
                t("%s & $%.5f\\pm%.5f$ & $%.2f\\pm%.2f$ & $%.0f\\pm%.0f$ & $%.4f\\pm%.4f$ & %s \\\\\n", sl[s].tex.c_str(), m.ff.mean, m.ff.eMean,
                  1e3 * m.ff.sigma, 1e3 * m.ff.eSigma, m.ff.S, m.ff.eS, m.ff.mean / scan::kMeta, m.ff.eMean / scan::kMeta, fk.c_str());
            else t("%s & n/a & n/a & n/a & n/a & $%+.0f\\pm%.0f$ \\\\\n", sl[s].tex.c_str(), m.fake, m.eFake);
        }
        t("\\midrule\n");
        const char* kt[3] = {"K$^0_S$, all rows (no $v_z$ cut)", "K$^0_S$, pre", "K$^0_S$, pre, isBest $=1$"};
        for (int i = 0; i < 3; ++i)
            t("%s & $%.5f\\pm%.5f$ & $%.2f\\pm%.2f$ & $%.0f\\pm%.0f$ & $%.5f\\pm%.5f$ & -- \\\\\n", kt[i], k0[i].pos, k0[i].ePos,
              1e3 * k0[i].w, 1e3 * k0[i].eW, k0[i].n, k0[i].eN, k0[i].pos / scan::kMK0, k0[i].ePos / scan::kMK0);
        t("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    }
    // ---- LaTeX: tagged, ratio, periods ----
    {
        Tex t("res_period.tex");
        t("%% res_period.tex -- written by plots/exp_report/res_peak.C (numbers: tables_exp/res_peak.txt)\n");
        t("\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{4pt}\n");
        t("\\caption{ECAL-tagged $\\eta$ content $T$ (rows with a photon, $0.119\\leq M(e^+e^-\\gamma)<0.173$~GeV/$c^2$, $m_{ee}<0.10$; frozen $m_{5c}$ peak "
          "0.5463~GeV/$c^2$, 6.2~MeV/$c^2$, pol2) and the ratio of the $\\eta$-peak content $S$ ($P_\\mathrm{main}$) to $T$ per period of the beam time, "
          "selection pre, isBest $\\geq0$; $Y_T=T\\times\\frac{64}{63}/L_\\mathrm{period}$ in visible rows per nb$^{-1}$; statistical errors only (the model spread is common). "
          "$S$ per period is identical to Table~\\ref{tab:smp_peak_periods}, where the peak yield per period is given. Last rows: weighted mean and $\\chi^2$/ndf of a constant. "
          "Sector 2 had a trigger deficit on days 034--044 (period I) that is not in $L$; its measured row loss, 3\\,\\%%, is too small to explain the low $Y_T$ of period I.}\n");
        t("\\label{tab:res_period}\n\\begin{tabular}{lrrccc}\n\\toprule\n");
        t("Period & $L$ [nb$^{-1}$] & $S$ & $T$ & $Y_T$ [rows/nb$^{-1}$] & $S/T$ \\\\\n\\midrule\n");
        for (int p = 0; p < kNPer; ++p)
            t("%s & %.1f & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ & $%.4f\\pm%.4f$ & $%.2f\\pm%.2f$ \\\\\n", kPerName[p], Lper[p],
              sPer[0][p].S, sPer[0][p].eS, tPer[p].S, tPer[p].eS, yP[2][p], eyP[2][p], rP[p], erP[p]);
        t("\\midrule\n");
        t("all & %.1f & $%.0f\\pm%.0f$ & $%.0f\\pm%.0f$ & $%.4f\\pm%.4f$ & $%.2f\\pm%.2f$ \\\\\n", Lsum, M[0].pm.S, M[0].pm.eS,
          tf[0].S, tf[0].eS, Yield(tf[0].S), Yield(tf[0].eS), rST[0], erST[0]);
        t("weighted mean & & & & $%.4f\\pm%.4f$ & $%.2f\\pm%.2f$ \\\\\n", wm[2], ewm[2], wm[3], ewm[3]);
        t("$\\chi^2$/ndf ($p$) & & & & %.1f/%d (%.2f) & %.1f/%d (%.2f) \\\\\n",
          c2[2], kNPer - 1, TMath::Prob(c2[2], kNPer - 1), c2[3], kNPer - 1, TMath::Prob(c2[3], kNPer - 1));
        t("\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    }

    // ------------------------------------------------------------------------
    // figures
    // ------------------------------------------------------------------------
    TFile fr(Form("%s/res_peak.root", kWork), "RECREATE");
    // ---- overview: pre, M4 0.35-0.95, log y ----
    {
        TCanvas* c = new TCanvas("c_ov", "", 1500, 1000);
        c->SetLogy(); c->SetLeftMargin(0.12); c->SetRightMargin(0.03); c->SetTopMargin(0.07);
        const Spectra& s0 = M[0].sp;
        TH1D* hN = s0.Hist("us", "ov_N"); TH1D* hC = s0.Hist("cb", "ov_CB"); TH1D* hS = s0.Hist("sig", "ov_sig");
        TH1D* f = Frame("ov_fr", 0.35, 0.95, 30, 3e5, "M(#pi^{+}#pi^{#minus}e^{+}e^{#minus}) [GeV/c^{2}]", "counts / 2.5 MeV/c^{2}");
        f->GetYaxis()->SetTitleOffset(1.0);
        f->Draw("axis");
        for (TH1D* h : {hN, hC, hS}) h->SetMarkerSize(0.7);
        hN->SetMarkerStyle(20); hN->SetMarkerColor(kBlack); hN->SetLineColor(kBlack);
        hC->SetMarkerStyle(20); hC->SetMarkerColor(kRed + 1); hC->SetLineColor(kRed + 1);
        hS->SetMarkerStyle(20); hS->SetMarkerColor(kBlue + 1); hS->SetLineColor(kBlue + 1);
        hN->Draw("pe same"); hS->Draw("pe same"); hC->Draw("pe same");
        TLine l; l.SetLineStyle(2); l.SetLineWidth(2);
        l.SetLineColor(kGreen + 2); l.DrawLine(scan::kBumpLo, 30, scan::kBumpLo, 1.2e5); l.DrawLine(scan::kBumpHi, 30, scan::kBumpHi, 1.2e5);
        l.SetLineColor(kMagenta + 1); l.DrawLine(scan::kPeakLo, 30, scan::kPeakLo, 1.2e5); l.DrawLine(scan::kPeakHi, 30, scan::kPeakHi, 1.2e5);
        TLegend* lg = new TLegend(0.55, 0.16, 0.95, 0.40);
        lg->SetTextSize(0.040);
        lg->AddEntry(hN, "N (e^{+}e^{#minus} counts)", "pe");
        lg->AddEntry(hC, "CB = 2#sqrt{N_{++}N_{#minus#minus}}", "pe");
        lg->AddEntry(hS, "N #minus CB", "pe");
        lg->Draw();
        Txt(0.12, 0.945, "data, selection pre (v_{z}>#minus160 mm), all stored combinations (isBest #geq 0), PT3", 0.040);
        Txt(0.16, 0.86, "green: bump window 0.4550#minus0.5225; magenta: peak window 0.53364#minus0.55416 GeV/c^{2}", 0.034);
        Save(c, "res_peak_overview");
        hN->Write(); hC->Write(); hS->Write();
    }
    // ---- peak fits of the five final selections ----
    // [review fix 2026-11-14] 3 + 2 layout (larger text), Y in visible rows per nb^-1, only the statistical error printed
    // (the model uncertainty quoted in the report is the edge-safe rms of fix_checks.C, not the old envelope)
    {
        const int idx[5] = {0, 3, 4, 5, 6};
        const char* lab[5] = {"pre", "peak nominal + vertex", "DIRECT", "PEAK", "X17 band"};
        TCanvas* c = new TCanvas("c_fits", "", 2400, 2300);
        // common y range of the lower rows
        double ymin = 0, ymax = 0;
        for (int k = 0; k < 5; ++k) {
            const Meas& m = M[idx[k] * nV];
            for (int i = 0; i < m.sp.nbin; ++i) {
                const double x = m.sp.x(i);
                if (x < 0.515 || x > 0.615) continue;
                const double r = m.sp.sig[i] - m.pm.Bkg(x), e = std::sqrt(m.sp.sigVar[i]);
                ymin = std::min(ymin, r - e); ymax = std::max(ymax, r + e);
            }
        }
        ymax *= 1.35; ymin *= 1.15;
        TGraphErrors* gLeg = nullptr; TGraph* gtLeg = nullptr; TGraph* gbLeg = nullptr;
        for (int k = 0; k < 5; ++k) {
            const Meas& m = M[idx[k] * nV];
            const int col = k % 3, blk = k / 3;
            const double yb = 0.5 - 0.5 * blk, ym = yb + 0.5 * 0.47, yt = yb + 0.5;   // top pad 53 %, bottom 47 % of a block
            for (int row = 0; row < 2; ++row) {
                c->cd();
                TPad* p = new TPad(Form("pf_%d_%d", k, row), "", col / 3.0, row == 0 ? ym : yb, (col + 1) / 3.0, row == 0 ? yt : ym);
                p->SetLeftMargin(0.215); p->SetRightMargin(0.03); p->SetTopMargin(row == 0 ? 0.12 : 0.03); p->SetBottomMargin(row == 0 ? 0.03 : 0.21);
                p->Draw(); p->cd();
                TGraphErrors* g = new TGraphErrors();
                double lo = 1e9, hi = -1e9;
                for (int i = 0; i < m.sp.nbin; ++i) {
                    const double x = m.sp.x(i);
                    if (x < 0.515 || x > 0.615) continue;
                    const double sc = row == 0 ? 1e-3 : 1.0;
                    const double y = sc * (row == 0 ? m.sp.sig[i] : m.sp.sig[i] - m.pm.Bkg(x)), e = sc * std::sqrt(m.sp.sigVar[i]);
                    const int n = g->GetN(); g->SetPoint(n, x, y); g->SetPointError(n, 0, e);
                    lo = std::min(lo, y - e); hi = std::max(hi, y + e);
                }
                const double y0 = row == 0 ? lo - 0.05 * (hi - lo) : ymin, y1 = row == 0 ? hi + 0.62 * (hi - lo) : ymax;
                TH1D* f = Frame(Form("fr_%d_%d", k, row), 0.515, 0.615, y0, y1, "M(#pi^{+}#pi^{#minus}e^{+}e^{#minus}) [GeV/c^{2}]",
                                row == 0 ? "10^{3} counts / 2.5 MeV/c^{2}" : "counts / 2.5 MeV/c^{2}");
                const double ts = row == 0 ? 0.068 : 0.076;
                f->GetXaxis()->SetLabelSize(row == 0 ? 0 : ts); f->GetXaxis()->SetTitleSize(row == 0 ? 0 : ts);
                f->GetYaxis()->SetLabelSize(ts); f->GetYaxis()->SetTitleSize(ts); f->GetYaxis()->SetTitleOffset(row == 0 ? 1.45 : 1.38);
                f->GetXaxis()->SetTitleOffset(1.10);
                f->GetYaxis()->SetMaxDigits(5); f->GetYaxis()->SetNdivisions(505); f->GetXaxis()->SetNdivisions(505);
                f->Draw("axis");
                g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetMarkerColor(kBlue + 1); g->SetLineColor(kBlue + 1);
                TGraph* gt = new TGraph(); TGraph* gb = new TGraph();
                for (int q = 0; q <= 200; ++q) {
                    const double x = scan::kFitLo + q * (scan::kFitHi - scan::kFitLo) / 200;
                    gt->SetPoint(q, x, row == 0 ? 1e-3 * m.pm.Total(x) : m.pm.Total(x) - m.pm.Bkg(x));
                    gb->SetPoint(q, x, row == 0 ? 1e-3 * m.pm.Bkg(x) : 0.);
                }
                gt->SetLineColor(kMagenta + 1); gt->SetLineWidth(3);
                gb->SetLineColor(kGray + 2); gb->SetLineStyle(2); gb->SetLineWidth(2);
                gb->Draw("l same"); g->Draw("pz same"); gt->Draw("l same");
                if (row == 0) {
                    gLeg = g; gtLeg = gt; gbLeg = gb;
                    Txt(0.215, 0.905, Form("(%c) %s: N #minus CB", 'a' + k, lab[k]), 0.072);
                    Txt(0.24, 0.79, Form("S = %.0f #pm %.0f (stat)", m.pm.S, m.pm.eS), 0.068, kBlue + 1);
                    Txt(0.24, 0.69, Form("#chi^{2}/ndf = %.1f/%d, S/B = %.3f", m.pm.chi2, m.pm.ndf, m.pm.SoB), 0.062);
                } else {
                    TLine l; l.SetLineStyle(3); l.SetLineColor(kGray + 1);
                    l.DrawLine(0.515, 0, 0.615, 0);
                    Txt(0.24, 0.88, Form("Y = %.3f #pm %.3f counts/nb^{#minus1}", Yield(m.pm.S), Yield(m.pm.eS)), 0.072, kBlue + 1);
                    Txt(0.60, 0.79, "N #minus CB #minus pol2", 0.066);
                }
            }
        }
        // key in the free slot
        c->cd();
        TPad* pk = new TPad("pf_key", "", 2 / 3.0, 0.0, 1.0, 0.5);
        pk->SetMargin(0.02, 0.02, 0.02, 0.02); pk->Draw(); pk->cd();
        TLegend* lg = new TLegend(0.06, 0.66, 0.98, 0.94);
        lg->SetTextSize(0.050); lg->SetFillStyle(0);
        lg->AddEntry(gLeg, "data, isBest #geq 0", "pe");
        lg->AddEntry(gtLeg, "P_{main} fit", "l");
        lg->AddEntry(gbLeg, "pol2 (top row)", "l");
        lg->Draw();
        const char* key[8] = {"P_{main}: Gaussian at 0.5439 GeV/c^{2},", "#sigma = 5.4 MeV/c^{2} (frozen), + pol2", "in 0.5225#minus0.6100 GeV/c^{2}",
                              "CB = 2#sqrt{N_{++}N_{#minus#minus}}, k = 1", "Y = S #times 64/63 / L,", "L = 6079.2 nb^{#minus1}",
                              "(visible counts per nb^{#minus1})", "errors: statistical only"};
        for (int q = 0; q < 8; ++q) Txt(0.06, 0.57 - 0.07 * q, key[q], 0.048);
        Save(c, "res_peak_fits");
    }
    // ---- per period: tagged eta and the ratio S/T (the peak per period is in the samples chapter, smp_peak_periods) ----
    {
        TCanvas* c = new TCanvas("c_per", "", 1500, 750);
        const char* ytit[2] = {"Y_{T} = T #times 64/63 / L [counts/nb^{#minus1}]", "S / T"};
        const char* ptit[2] = {"(a) ECAL-tagged #eta, pre", "(b) #eta peak / ECAL-tagged #eta, pre"};
        for (int q = 0; q < 2; ++q) {
            const int iq = q == 0 ? 2 : 3;
            double y[kNPer], e[kNPer];
            for (int pp = 0; pp < kNPer; ++pp) { y[pp] = q == 0 ? yP[2][pp] : rP[pp]; e[pp] = q == 0 ? eyP[2][pp] : erP[pp]; }
            c->cd();
            TPad* p = new TPad(Form("pp_%d", q), "", q / 2.0, 0, (q + 1) / 2.0, 1);
            p->SetLeftMargin(0.20); p->SetRightMargin(0.03); p->SetTopMargin(0.16); p->SetBottomMargin(0.14);
            p->Draw(); p->cd();
            double lo = 1e9, hi = -1e9;
            for (int pp = 0; pp < kNPer; ++pp) { lo = std::min(lo, y[pp] - e[pp]); hi = std::max(hi, y[pp] + e[pp]); }
            TH1D* f = new TH1D(Form("fper_%d", q), "", kNPer, 0, kNPer); f->SetDirectory(nullptr); f->SetStats(0);
            for (int pp = 0; pp < kNPer; ++pp) f->GetXaxis()->SetBinLabel(pp + 1, kPerName[pp]);
            f->SetMinimum(std::max(0., lo - 0.3 * (hi - lo))); f->SetMaximum(hi + 0.45 * (hi - lo));
            f->GetYaxis()->SetTitle(ytit[q]); f->GetXaxis()->SetTitle("period");
            f->GetXaxis()->SetLabelSize(0.07); f->GetYaxis()->SetTitleOffset(1.55);
            f->Draw("axis");
            TGraphErrors* g = new TGraphErrors();
            for (int pp = 0; pp < kNPer; ++pp) { g->SetPoint(pp, pp + 0.5, y[pp]); g->SetPointError(pp, 0, e[pp]); }
            g->SetMarkerStyle(20); g->SetMarkerSize(1.4); g->SetMarkerColor(q == 0 ? kGreen + 2 : kBlue + 1); g->SetLineColor(g->GetMarkerColor());
            g->SetLineWidth(2);
            TBox* band = new TBox(0, wm[iq] - ewm[iq], kNPer, wm[iq] + ewm[iq]);
            band->SetFillColorAlpha(kGray, 0.5); band->Draw();
            TLine l; l.SetLineStyle(2); l.SetLineColor(kGray + 2); l.SetLineWidth(2); l.DrawLine(0, wm[iq], kNPer, wm[iq]);
            g->Draw("pz same");
            Txt(0.20, 0.93, ptit[q], 0.050);
            Txt(0.20, 0.86, Form(q == 0 ? "mean %.3f #pm %.3f counts/nb^{#minus1}, #chi^{2}/ndf = %.1f/%d" : "mean %.2f #pm %.2f, #chi^{2}/ndf = %.1f/%d",
                                 wm[iq], ewm[iq], c2[iq], kNPer - 1), q == 0 ? 0.040 : 0.045);
        }
        Save(c, "res_tagged_period");
    }
    fr.Close();
    printf("res_peak: done in %.0f s\n", swAll.RealTime());
}

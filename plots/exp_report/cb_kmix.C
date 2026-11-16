// ============================================================================
// cb_kmix.C  --  event-mixing k factor of the combinatorial background in
//                M(pi+pi-e+e-), experimental data
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// The like-sign estimate CB = 2 sqrt(N++ N--) equals the unlike-sign pairs of
// independent leptons only if e+ and e- have the same kinematic distributions
// (shape) or, more precisely, if the unlike-sign pair acceptance equals the
// geometric mean of the like-sign ones. The correction k (CB = k 2 sqrt(N++ N--))
// is measured, as in the e+e- study (plots/k_mix_pools.C, branch pp45_epem_sim),
// from pairs of leptons taken from DIFFERENT events:
//   k(M) = N+-_mix(M) / (2 sqrt(N++_mix(M) N--_mix(M)))
// For the four-particle channel the mixed combination keeps the pion pair and
// ONE lepton of the current event E and takes the other lepton from an earlier
// event B of the same class:
//   N+-_mix: (pi+pi- e+)_E e-_B and (pi+pi- e-)_E e+_B
//   N++_mix: (pi+pi- e+)_E e+_B ;  N--_mix: (pi+pi- e-)_E e-_B
// The same selection (cuts on the mixed combination) is applied as for the rows.
//
// Event pool: every event (adjacent rows of one group) of the e+e- tree with
// selection pre (cor_ok == 1, isBest >= 0, vz > -160 mm; FAT: vertex z > -500 mm,
// PT3); tracks = distinct (p, theta, phi) of the rows with cor_ok == 1; pion
// pairs = the (pi+, pi-) pairs that occur in the rows; a lepton that is the same
// track as one of the pions (equal theta and phi) is not combined with that pair.
// Variants (all with a rolling buffer, B earlier than E):
//   A     main: 5 classes of the event vertex z (20 mm: -160..-140 ... -80..),
//         the last 10 events of the class; all selections, M(pi+pi-) slices of
//         the bump; statistical error from 10 sub-samples (event index mod 10)
//   Anoz  one class for all z (last 10 events)
//   Ange3 only events with >= 3 identified leptons in E and B (their leptons
//         are mostly unpaired, as in the "nge3" pool of the e+e- study)
//   B     pions from E, BOTH leptons from two different earlier events of the
//         class (last 5 events, unordered pairs): no lepton is correlated with
//         the pions; tests whether the pion-lepton correlation matters
//   C     lepton pairs only (E lepton x B lepton) -> k_ee(M(e+e-)), transferred
//         to M4 by averaging k_ee over the M(e+e-) of the like-sign rows in each
//         10 MeV bin of M4 (k from the lepton pair alone; "do the pions matter?")
// Then k_A(M4) per selection (and per bump slice) is applied to the data CB and
// the eta peak S (P_main), the bump edge-step sum E and the bump-window N - CB
// (W) are recomputed.
//
// Selections (cb_common.h; rows = all stored combinations, isBest >= 0)
//   pre, pre + OA>2, pre + vertex (both pions in the event-vertex fit),
//   DIRECT (pre + dz_pi < 40 mm), PEAK (DIRECT + MM4 < 2.85 GeV/c2),
//   peak nominal (pre + oa_lab < 55 + M(pipi) > 0.30), + vertex, + OA>2,
//   bump nominal (pre + oa_lab < 40 + 0.30 < M(pipi) < 0.41), + vertex, + OA>2.
//   No k for isBest == 1 (the plain like-sign formula is biased there for
//   another reason, see cb_checks.C, section B).
//
// Inputs : /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/scan_v2/scan_compact.root
//          (scan::kCompact) through ../ScanEngine.h and cb_common.h
// Outputs: /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/cb/cb_kmix.root
//            (k_<variant>_<selection> in 10 MeV bins, k_A_<sel>_sl<i>, k_ee, mixed spectra)
//          /home/przygoda/HADES/FAT/FAT/docs/latex/tables_exp/cb_kmix.txt, cb_kmix.tex
//          /home/przygoda/HADES/FAT/FAT/docs/latex/figs_exp/cb_kmix.pdf (+ .png)
//
// Run (from /home/przygoda/HADES/FAT/FAT) BEFORE cb_checks.C, about 1 min, 3 GB of memory:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_cb", kTRUE)' 'plots/exp_report/cb_kmix.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStopwatch.h"

#include "cb_common.h"

namespace kmx {

using namespace cbc;

const double kMe = 0.51099895, kMpi = 139.56995, kMp = 938.27231;   // MeV
const double kD2R = 1.74532925199432955e-02, kR2D = 57.2957795130823229;

struct V4 { double x = 0, y = 0, z = 0, e = 0; };
inline V4 Mk(double p, double th, double ph, double m) {
    const double t = th * kD2R, f = ph * kD2R;
    V4 v; v.x = p * std::sin(t) * std::cos(f); v.y = p * std::sin(t) * std::sin(f); v.z = p * std::cos(t);
    v.e = std::sqrt(p * p + m * m); return v;
}
inline V4 operator+(const V4& a, const V4& b) { V4 v; v.x = a.x + b.x; v.y = a.y + b.y; v.z = a.z + b.z; v.e = a.e + b.e; return v; }
inline V4 operator-(const V4& a, const V4& b) { V4 v; v.x = a.x - b.x; v.y = a.y - b.y; v.z = a.z - b.z; v.e = a.e - b.e; return v; }
inline double Mass(const V4& v) { const double m2 = v.e * v.e - v.x * v.x - v.y * v.y - v.z * v.z; return m2 > 0 ? std::sqrt(m2) : 0.0; }
// signed mass as TLorentzVector::M() (negative for m2 < 0), used for the missing mass
inline double MassS(const V4& v) { const double m2 = v.e * v.e - v.x * v.x - v.y * v.y - v.z * v.z; return m2 >= 0 ? std::sqrt(m2) : -std::sqrt(-m2); }
inline double Angle(const V4& a, const V4& b) {
    const double na = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z), nb = std::sqrt(b.x * b.x + b.y * b.y + b.z * b.z);
    if (!(na > 0 && nb > 0)) return 0;
    double c = (a.x * b.x + a.y * b.y + a.z * b.z) / (na * nb);
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) * kR2D;
}
V4 Initial() {
    const double eb = scan::kTbeam + kMp, pb = std::sqrt(eb * eb - kMp * kMp);
    V4 v; v.z = pb; v.e = eb + kMp; return v;
}

struct Lep { V4 v; float th, ph, p; };
struct Pair { V4 sum; double mpipi; bool vtx; float dz; float thp, php, thm, phm; };
struct Ev { std::vector<Pair> pr; std::vector<Lep> lp[2]; int zc = 0; int nlep = 0; long idx = 0; };

// ---- mixed-combination selections (same cuts as cb_common.h, on the mixed combination)
const int kNSelM = 11;
const char* const kSelKey[kNSelM] = {"pre", "preoa", "prevtx", "direct", "peaksim", "peaknom", "peakvtx", "bumpnom", "bumpvtx", "peakoa", "bumpoa"};
inline unsigned SelBits(double oa, double oaLab, double mpipi, bool vtx, float dz, double mm4) {
    unsigned b = 1u;                                            // pre
    const bool o2 = oa > 2.0, dzc = dz < 40.f, pk = oaLab < 55.0 && mpipi > 0.30, bm = oaLab < 40.0 && mpipi > 0.30 && mpipi < 0.41;
    if (o2) b |= 1u << 1;
    if (vtx) b |= 1u << 2;
    if (dzc) b |= 1u << 3;
    if (dzc && mm4 < 2.85) b |= 1u << 4;
    if (pk) b |= 1u << 5;
    if (pk && vtx) b |= 1u << 6;
    if (bm) b |= 1u << 7;
    if (bm && vtx) b |= 1u << 8;
    if (pk && o2) b |= 1u << 9;
    if (bm && o2) b |= 1u << 10;
    return b;
}
const int kSlSel[3] = {0, 7, 8};          // selections with M(pipi)-slice spectra: pre, bump nominal, bump nominal + vertex

// ---- accumulators: [sub][combo 0 +-, 1 ++, 2 --][bin] ---------------------------------
const int kNsub = 10;
struct Acc {
    int nb = 0, nsub = 1; double lo = 0, hi = 1;
    std::vector<double> h;
    void Init(int n, double a, double b, int ns) { nb = n; lo = a; hi = b; nsub = ns; h.assign((size_t)ns * 3 * n, 0.); }
    inline void Fill(int sub, int c, double x) {
        if (!(x >= lo && x < hi)) return;
        const int b = (int)((x - lo) / (hi - lo) * nb);
        h[((size_t)(sub % nsub) * 3 + c) * nb + b] += 1;
    }
    double At(int sub, int c, int b) const { return h[((size_t)sub * 3 + c) * nb + b]; }
    double Sum(int c, int b) const { double s = 0; for (int q = 0; q < nsub; ++q) s += At(q, c, b); return s; }
    // k in wide bins of reb fine bins; sub < 0: all sub-samples
    std::vector<double> K(int reb, int sub = -1) const {
        const int nw = nb / reb; std::vector<double> k(nw, 1.0);
        for (int j = 0; j < nw; ++j) {
            double s[3] = {0, 0, 0};
            for (int c = 0; c < 3; ++c)
                for (int i = j * reb; i < (j + 1) * reb; ++i) s[c] += sub < 0 ? Sum(c, i) : At(sub, c, i);
            k[j] = (s[1] > 0 && s[2] > 0) ? s[0] / (2.0 * std::sqrt(s[1] * s[2])) : 1.0;
        }
        return k;
    }
    double KRange(double a, double b, int sub = -1) const {
        double s[3] = {0, 0, 0};
        for (int i = 0; i < nb; ++i) {
            const double x = lo + (i + 0.5) * (hi - lo) / nb;
            if (x < a || x >= b) continue;
            for (int c = 0; c < 3; ++c) s[c] += sub < 0 ? Sum(c, i) : At(sub, c, i);
        }
        return (s[1] > 0 && s[2] > 0) ? s[0] / (2.0 * std::sqrt(s[1] * s[2])) : 1.0;
    }
    double Total(int c) const { double s = 0; for (int i = 0; i < nb; ++i) s += Sum(c, i); return s; }
};
// error of a k-derived number from the spread over the sub-samples
template <class F> void SubErr(const Acc& a, F f, double& val, double& err) {
    val = f(-1);
    if (a.nsub < 2) { err = 0; return; }
    double s = 0, s2 = 0;
    for (int q = 0; q < a.nsub; ++q) { const double v = f(q); s += v; s2 += v * v; }
    const double m = s / a.nsub, var = std::max(0., s2 / a.nsub - m * m);
    err = std::sqrt(var / (a.nsub - 1));
}

} // namespace kmx

// ============================================================================
void cb_kmix(int depth = 10, int depthB = 5) {
    using namespace kmx;
    TStopwatch swAll;
    gSystem->mkdir(kWork, kTRUE);
    OutFile() = fopen(Form("%s/cb_kmix.txt", kTabOut), "w");
    P("cb_kmix.txt -- written by plots/exp_report/cb_kmix.C (experimental data only)\n");
    P("input %s\n", scan::kCompact);
    P("k(M4) = N+-_mix / (2 sqrt(N++_mix N--_mix)); mixed combination = (pi+ pi- l1) of event E + l2 of an earlier event B.\n");
    P("pool depth %d (variants A, Anoz, Ange3), %d (variant B); vertex classes 20 mm from -160 mm.\n\n", depth, depthB);

    Data d;
    if (!scan::Load(d)) return;
    const V4 ini = Initial();

    // ---- accumulators ----
    const int nb = scan::kNbin;
    std::vector<Acc> A(kNSelM), Anoz(kNSelM), Ange3(1), Bv(1);
    for (int s = 0; s < kNSelM; ++s) { A[s].Init(nb, scan::kMlo, scan::kMhi, kNsub); Anoz[s].Init(nb, scan::kMlo, scan::kMhi, 1); }
    Ange3[0].Init(nb, scan::kMlo, scan::kMhi, kNsub); Bv[0].Init(nb, scan::kMlo, scan::kMhi, kNsub);
    std::vector<std::vector<Acc>> Asl(3, std::vector<Acc>(kNSig));
    for (int q = 0; q < 3; ++q) for (int i = 0; i < kNSig; ++i) Asl[q][i].Init(nb, scan::kMlo, scan::kMhi, kNsub);
    Acc aMee, aOa, aC;            // variant A, pre: M(ee) and OA of the mixed pair; variant C: lepton pairs only
    aMee.Init(160, 0.0, 0.80, kNsub); aOa.Init(60, 0.0, 120.0, kNsub); aC.Init(160, 0.0, 0.80, kNsub);

    // ---- validation of the kinematics on the rows themselves ----
    double maxd[5] = {0, 0, 0, 0, 0};
    long nVal = 0;
    for (const Row& r : d.s[0]) {
        if (!r.cor_ok || !Pre(r)) continue;
        const V4 pp = Mk(r.pip_p, r.pip_theta, r.pip_phi, kMpi), pm = Mk(r.pim_p, r.pim_theta, r.pim_phi, kMpi);
        const V4 ep = Mk(r.ep_p, r.ep_theta, r.ep_phi, kMe), em = Mk(r.em_p, r.em_theta, r.em_phi, kMe);
        const V4 v4 = pp + pm + ep + em, vee = ep + em, vpp = pp + pm;
        const double c[5] = {Mass(v4) / 1000., Mass(vee) / 1000., Angle(ep, em), Angle(vpp, vee), MassS(ini - v4) / 1000.};
        const double s[5] = {r.m4, r.mee, r.oa, r.oa_lab, r.mm4};
        for (int k = 0; k < 5; ++k) maxd[k] = std::max(maxd[k], std::fabs(c[k] - s[k]));
        ++nVal;
    }
    P("validation on %ld unlike-sign rows (pre): max |recomputed - stored|: M4 %.2g, M(ee) %.2g GeV/c2, OA %.2g, oa_lab %.2g deg, MM4 %.2g GeV/c2\n\n",
      nVal, maxd[0], maxd[1], maxd[2], maxd[3], maxd[4]);

    // ---- one mixed combination ----
    long nComb = 0;
    auto fillCombo = [&](int sub, int c, const Pair& pr, const Lep& a, const Lep& b, std::vector<Acc>* acc, bool slices, bool pairHists,
                         unsigned mask) {
        const V4 vee = a.v + b.v, v4 = pr.sum + vee;
        const double m4 = Mass(v4) / 1000.;
        const double oa = Angle(a.v, b.v);
        if (pairHists) { aMee.Fill(sub, c, Mass(vee) / 1000.); aOa.Fill(sub, c, oa); }
        if (!(m4 >= scan::kMlo && m4 < scan::kMhi)) return;
        ++nComb;
        const double oaLab = Angle(pr.sum, vee);
        const double mm4 = (mask & ((1u << 4))) ? MassS(ini - v4) / 1000. : 0.;
        const unsigned bits = SelBits(oa, oaLab, pr.mpipi, pr.vtx, pr.dz, mm4) & mask;
        for (int s = 0; s < (int)acc->size(); ++s) if (bits & (1u << s)) (*acc)[s].Fill(sub, c, m4);
        if (slices) {
            const int sl = SliceOf(pr.mpipi);
            if (sl >= 0) for (int q = 0; q < 3; ++q) if (bits & (1u << kSlSel[q])) Asl[q][sl].Fill(sub, c, m4);
        }
    };
    auto overlap = [](const Pair& pr, const Lep& l) {
        return (l.th == pr.thp && l.ph == pr.php) || (l.th == pr.thm && l.ph == pr.phm);
    };
    // E (pions + lepton) x lepton of B
    auto mixPair = [&](const Ev& E, const Ev& B, std::vector<Acc>* acc, bool slices, bool pairHists, unsigned mask) {
        for (const Pair& pr : E.pr)
            for (int qa = 0; qa < 2; ++qa)                     // charge of E's lepton: 0 = e+, 1 = e-
                for (const Lep& a : E.lp[qa]) {
                    if (overlap(pr, a)) continue;
                    for (int qb = 0; qb < 2; ++qb) {
                        const int c = (qa != qb) ? 0 : (qa == 0 ? 1 : 2);
                        for (const Lep& b : B.lp[qb]) fillCombo((int)(E.idx % kNsub), c, pr, a, b, acc, slices, pairHists, mask);
                    }
                }
    };
    // variant B: pions of E, leptons of X and Y (X earlier than Y or not: unordered pair)
    auto mixThree = [&](const Ev& E, const Ev& X, const Ev& Y) {
        for (const Pair& pr : E.pr)
            for (int qa = 0; qa < 2; ++qa)
                for (const Lep& a : X.lp[qa])
                    for (int qb = 0; qb < 2; ++qb) {
                        const int c = (qa != qb) ? 0 : (qa == 0 ? 1 : 2);
                        for (const Lep& b : Y.lp[qb]) {
                            const V4 v4 = pr.sum + a.v + b.v;
                            const double m4 = Mass(v4) / 1000.;
                            Bv[0].Fill((int)(E.idx % kNsub), c, m4);
                        }
                    }
    };
    // variant C: lepton pairs only
    auto mixLep = [&](const Ev& E, const Ev& B) {
        for (int qa = 0; qa < 2; ++qa)
            for (const Lep& a : E.lp[qa])
                for (int qb = 0; qb < 2; ++qb) {
                    const int c = (qa != qb) ? 0 : (qa == 0 ? 1 : 2);
                    for (const Lep& b : B.lp[qb]) aC.Fill((int)(E.idx % kNsub), c, Mass(a.v + b.v) / 1000.);
                }
    };

    // ---- stream over the event groups of the e+e- tree ----
    const int nZ = 5;
    std::vector<std::deque<Ev>> bufZ(nZ);
    std::deque<Ev> bufAll, bufN3;
    long nEv = 0, nEvN3 = 0, nPairs = 0, nLep = 0;
    const std::vector<Row>& us = d.s[0];
    size_t i0 = 0;
    TStopwatch sw;
    while (i0 < us.size()) {
        size_t i1 = i0;
        while (i1 < us.size() && us[i1].grp == us[i0].grp) ++i1;
        Ev E;
        bool any = false;
        std::vector<std::pair<int, int>> seenPair;
        std::vector<Lep> pips, pims;
        for (size_t i = i0; i < i1; ++i) {
            const Row& r = us[i];
            if (!r.cor_ok || !Pre(r)) continue;
            any = true;
            E.zc = std::max(0, std::min(nZ - 1, (int)((r.vz + 160.f) / 20.f)));
            E.nlep = r.n_lep1 + r.n_lep2;
            auto findOrAdd = [](std::vector<Lep>& v, float p, float th, float ph, double m) {
                for (size_t k = 0; k < v.size(); ++k) if (v[k].p == p && v[k].th == th && v[k].ph == ph) return (int)k;
                Lep l; l.p = p; l.th = th; l.ph = ph; l.v = Mk(p, th, ph, m); v.push_back(l); return (int)v.size() - 1;
            };
            const int ip = findOrAdd(pips, r.pip_p, r.pip_theta, r.pip_phi, kMpi);
            const int im = findOrAdd(pims, r.pim_p, r.pim_theta, r.pim_phi, kMpi);
            findOrAdd(E.lp[0], r.ep_p, r.ep_theta, r.ep_phi, kMe);
            findOrAdd(E.lp[1], r.em_p, r.em_theta, r.em_phi, kMe);
            const std::pair<int, int> key(ip, im);
            if (std::find(seenPair.begin(), seenPair.end(), key) == seenPair.end()) {
                seenPair.push_back(key);
                Pair pr; pr.sum = pips[ip].v + pims[im].v; pr.mpipi = Mass(pr.sum) / 1000.;
                pr.vtx = r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; pr.dz = r.dz_pi;
                pr.thp = r.pip_theta; pr.php = r.pip_phi; pr.thm = r.pim_theta; pr.phm = r.pim_phi;
                E.pr.push_back(pr);
            }
        }
        i0 = i1;
        if (!any) continue;
        E.idx = nEv++;
        nPairs += E.pr.size(); nLep += E.lp[0].size() + E.lp[1].size();
        // variant A (vertex classes), with slices and pair histograms
        for (const Ev& B : bufZ[E.zc]) mixPair(E, B, &A, true, true, 0x7FFu);
        // variant Anoz (pre, peak nominal, bump nominal)
        for (const Ev& B : bufAll) mixPair(E, B, &Anoz, false, false, (1u << 0) | (1u << 5) | (1u << 7));
        // variant Ange3
        if (E.nlep >= 3) {
            ++nEvN3;
            for (const Ev& B : bufN3) mixPair(E, B, &Ange3, false, false, 1u);
        }
        // variant B: unordered pairs among the last depthB events of the class
        {
            const std::deque<Ev>& bz = bufZ[E.zc];
            const int n = (int)bz.size(), first = std::max(0, n - depthB);
            for (int x = first; x < n; ++x) for (int y = x + 1; y < n; ++y) mixThree(E, bz[x], bz[y]);
        }
        // variant C
        for (const Ev& B : bufZ[E.zc]) mixLep(E, B);
        // buffers
        bufZ[E.zc].push_back(E); if ((int)bufZ[E.zc].size() > depth) bufZ[E.zc].pop_front();
        bufAll.push_back(E); if ((int)bufAll.size() > depth) bufAll.pop_front();
        if (E.nlep >= 3) { bufN3.push_back(E); if ((int)bufN3.size() > depth) bufN3.pop_front(); }
        if (nEv % 500000 == 0) printf("cb_kmix: %ld events, %.0f s\n", nEv, sw.RealTime()), sw.Continue();
    }
    P("events %ld (>= 3 leptons: %ld); pion pairs per event %.3f, leptons per event %.3f; mixed combinations in 0.35-0.95 (A, Anoz, Ange3): %ld\n",
      nEv, nEvN3, (double)nPairs / nEv, (double)nLep / nEv, nComb);
    P("integrated k (0.35-0.95): A pre %.4f, Anoz pre %.4f, Ange3 %.4f, B %.4f; mixed counts A pre: +- %.0f, ++ %.0f, -- %.0f; ratio --/++ %.3f\n\n",
      A[0].KRange(0.35, 0.95), Anoz[0].KRange(0.35, 0.95), Ange3[0].KRange(0.35, 0.95), Bv[0].KRange(0.35, 0.95), A[0].Total(0), A[0].Total(1), A[0].Total(2),
      A[0].Total(2) / A[0].Total(1));

    // ======================================================================
    // data spectra and the CB-weighted k per region
    // ======================================================================
    const std::vector<SD> sels = Selections();
    std::vector<int> sIdx(kNSelM);
    std::vector<Sel> sv;
    for (int s = 0; s < kNSelM; ++s) { sIdx[s] = Find(sels, kSelKey[s]); sv.push_back(sels[sIdx[s]].sel); }
    scan::FillOpt fo;
    fo.windows.push_back(std::make_pair(0.35, 0.455));
    fo.windows.push_back(std::make_pair(0.60, 0.70));
    fo.windows.push_back(std::make_pair(0.70, 0.95));
    fo.windows.push_back(std::make_pair(0.35, 0.95));
    const std::vector<Spectra> sp = scan::FillMany(d, sv, fo);
    const int reb = scan::kCBreb, nw = nb / reb;
    // CB-weighted mean of k over [a, b): sum_j k_j CB_j / sum_j CB_j (data CB of the same selection)
    auto kAvg = [&](const std::vector<double>& k, const Spectra& s, double a, double b) {
        double num = 0, den = 0;
        for (int j = 0; j < nw; ++j) {
            const double x = scan::kMlo + (j + 0.5) * reb * scan::kBinW;
            if (x < a || x >= b) continue;
            num += k[j] * s.cbW[j]; den += s.cbW[j];
        }
        return den > 0 ? num / den : 1.0;
    };
    const double rA[6] = {0.35, scan::kBumpLo, 0.5225, 0.56, 0.60, 0.70};   // 0.35-0.455, bump, peak region 0.5225-0.56, 0.56-0.60, 0.60-0.70, 0.70-0.95
    const double rB[6] = {0.455, scan::kBumpHi, 0.56, 0.60, 0.70, 0.95};
    const char* rN[6] = {"0.35-0.455", "bump 0.455-0.5225", "0.5225-0.56", "0.56-0.60", "0.60-0.70", "0.70-0.95"};
    P("A. k_A (vertex classes, depth %d), CB-weighted mean per M4 region [GeV/c2] (10 MeV bins, data CB of the selection),\n", depth);
    P("   error = spread over 10 sub-samples / sqrt(10)\n");
    P("%-24s", "selection");
    for (int q = 0; q < 6; ++q) P(" | %-17s", rN[q]);
    P(" | %-15s\n", "0.35-0.95");
    std::vector<std::vector<double>> kA(kNSelM);
    for (int s = 0; s < kNSelM; ++s) {
        kA[s] = A[s].K(reb);
        P("%-24s", sels[sIdx[s]].name.c_str());
        for (int q = 0; q <= 6; ++q) {
            const double a = q < 6 ? rA[q] : 0.35, b = q < 6 ? rB[q] : 0.95;
            double v, e;
            SubErr(A[s], [&](int sub) { return kAvg(A[s].K(reb, sub), sp[s], a, b); }, v, e);
            P(" | %7.4f +- %6.4f", v, e);
        }
        P("\n");
    }
    // ---- variants for pre ----
    P("\nB. VARIANTS for pre (CB-weighted mean per region):\n");
    {
        const std::vector<double> kNoz = Anoz[0].K(reb), kN3 = Ange3[0].K(reb), kB = Bv[0].K(reb);
        // variant C: k_ee per M(ee) bin (5 MeV) -> average over the like-sign rows of each 10 MeV M4 bin
        const std::vector<double> kee = aC.K(1);
        std::vector<double> sumK[2], cnt[2];
        for (int q = 0; q < 2; ++q) { sumK[q].assign(nw, 0.); cnt[q].assign(nw, 0.); }
        for (int t = 1; t <= 2; ++t)
            for (const Row& r : d.s[t]) {
                if (!r.cor_ok || !Pre(r)) continue;
                if (!(r.m4 >= scan::kMlo && r.m4 < scan::kMhi)) continue;
                const int j = (int)((r.m4 - scan::kMlo) / (reb * scan::kBinW));
                int ib = (int)(r.mee / 0.005); ib = std::max(0, std::min(159, ib));
                sumK[t - 1][j] += kee[ib]; cnt[t - 1][j] += 1;
            }
        std::vector<double> kC(nw, 1.0);
        for (int j = 0; j < nw; ++j) {
            const double a = cnt[0][j] > 0 ? sumK[0][j] / cnt[0][j] : 1.0, b = cnt[1][j] > 0 ? sumK[1][j] / cnt[1][j] : 1.0;
            kC[j] = 0.5 * (a + b);
        }
        const std::vector<double>* kv[5] = {&kA[0], &kNoz, &kN3, &kB, &kC};
        const char* vn[5] = {"A (z classes)", "Anoz (one pool)", "Ange3 (>=3 leptons)", "B (both leptons foreign)", "C (lepton pair only)"};
        for (int v = 0; v < 5; ++v) {
            P("%-26s", vn[v]);
            for (int q = 0; q <= 6; ++q) {
                const double a = q < 6 ? rA[q] : 0.35, b = q < 6 ? rB[q] : 0.95;
                P(" | %7.4f", kAvg(*kv[v], sp[0], a, b));
            }
            if (v == 2 || v == 3) {
                const Acc& ac = v == 2 ? Ange3[0] : Bv[0];
                double val, e; SubErr(ac, [&](int sub) { return kAvg(ac.K(reb, sub), sp[0], 0.35, 0.95); }, val, e);
                P("   (0.35-0.95 error %.4f)", e);
            }
            P("\n");
        }
        // k_ee in M(ee) classes (variant C) and k_A in OA / M(ee) of the mixed pair (pre)
        P("\nC. k_ee(M(ee)) of lepton pairs from different events (variant C, as the e+e- study) and k_A against M(ee), OA of the mixed pair (pre):\n");
        const double eLo[] = {0.0, 0.01, 0.02, 0.05, 0.10, 0.14, 0.20, 0.30, 0.50};
        const double eHi[] = {0.01, 0.02, 0.05, 0.10, 0.14, 0.20, 0.30, 0.50, 0.80};
        for (int q = 0; q < 9; ++q) {
            double v1, e1, v2, e2;
            SubErr(aC, [&](int sub) { return aC.KRange(eLo[q], eHi[q], sub); }, v1, e1);
            SubErr(aMee, [&](int sub) { return aMee.KRange(eLo[q], eHi[q], sub); }, v2, e2);
            P("   M(ee) %4.2f-%4.2f: k_ee %.4f +- %.4f | k_A(4-body, pre) %.4f +- %.4f\n", eLo[q], eHi[q], v1, e1, v2, e2);
        }
        {
            double v1, e1; SubErr(aC, [&](int sub) { return aC.KRange(0.14, 0.80, sub); }, v1, e1);
            P("   k_ee(0.14-0.80) = %.4f +- %.4f (e+e- study, data, PT3, classes: 1.032 +- 0.007)\n", v1, e1);
        }
        const double oLo[] = {0, 2, 4, 6, 10, 15, 20, 30, 50, 80}, oHi[] = {2, 4, 6, 10, 15, 20, 30, 50, 80, 120};
        for (int q = 0; q < 10; ++q) {
            double v, e; SubErr(aOa, [&](int sub) { return aOa.KRange(oLo[q], oHi[q], sub); }, v, e);
            P("   OA %3.0f-%3.0f deg: k_A %.4f +- %.4f\n", oLo[q], oHi[q], v, e);
        }

        // ---- shape test: real like-sign rows against mixed like-sign (pre), normalised in 0.35-0.95 ----
        P("\nD. SHAPE TEST (pre): real like-sign rows / mixed like-sign combinations (variant A), both normalised to unit sum in 0.35-0.95;\n");
        P("   CB = 2 sqrt(N++ N--) / mixed +- (normalised). 50 MeV regions:\n");
        const Spectra& s0 = sp[0];
        double tr[3] = {0, 0, 0}, tm[3] = {0, 0, 0};
        for (int i = 0; i < nb; ++i) { tr[1] += s0.n[1][i]; tr[2] += s0.n[2][i]; tr[0] += s0.cb[i]; for (int c = 0; c < 3; ++c) tm[c] += A[0].Sum(c, i); }
        for (int q = 0; q < 12; ++q) {
            const double a = 0.35 + 0.05 * q, b = a + 0.05;
            double r[3] = {0, 0, 0}, m[3] = {0, 0, 0};
            for (int i = 0; i < nb; ++i) {
                const double x = scan::kMlo + (i + 0.5) * scan::kBinW;
                if (x < a || x >= b) continue;
                r[1] += s0.n[1][i]; r[2] += s0.n[2][i]; r[0] += s0.cb[i];
                for (int c = 0; c < 3; ++c) m[c] += A[0].Sum(c, i);
            }
            P("   %4.2f-%4.2f: ++ %.3f +- %.3f | -- %.3f +- %.3f | CB/mixed+- %.3f\n", a, b, (r[1] / tr[1]) / (m[1] / tm[1]), std::sqrt(r[1]) / tr[1] / (m[1] / tm[1]),
              (r[2] / tr[2]) / (m[2] / tm[2]), std::sqrt(r[2]) / tr[2] / (m[2] / tm[2]), (r[0] / tr[0]) / (m[0] / tm[0]));
        }

        // ======================================================================
        // E. effect of k_A on S, E, W
        // ======================================================================
        P("\nE. EFFECT OF k_A(M4) (per 10 MeV bin, the selection's own k) ON THE ETA PEAK S, THE BUMP EDGE-STEP SUM E AND W = N - CB (bump window)\n");
        P("   selection | S stat -> S with k (relative change) | W -> W with k (relative change) | CB/N peak window: k = 1 -> k | CB/N bump window\n");
        TFile fo2(kKmixFile, "RECREATE");
        for (int s = 0; s < kNSelM; ++s) {
            const PeakFit a = scan::FitPeakMain(sp[s]);
            const Spectra m = Modify(sp[s], 1.0, 0, &kA[s]);
            const PeakFit b = scan::FitPeakMain(m);
            P("%-24s | %7.0f %5.0f -> %7.0f (%+6.3f) | W %8.0f -> %8.0f (%+6.3f) | %.4f -> %.4f | %.4f -> %.4f\n", sels[sIdx[s]].name.c_str(), a.S, a.eS, b.S,
              (b.S - a.S) / a.S, sp[s].win[1].sig(), m.win[1].sig(), (m.win[1].sig() - sp[s].win[1].sig()) / sp[s].win[1].sig(), sp[s].win[0].cbShare(),
              m.win[0].cbShare(), sp[s].win[1].cbShare(), m.win[1].cbShare());
            // histograms of k with sub-sample errors
            TH1D* h = new TH1D(Form("k_A_%s", kSelKey[s]), ";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];k", nw, scan::kMlo, scan::kMhi);
            for (int j = 0; j < nw; ++j) {
                double v, e; SubErr(A[s], [&](int sub) { return A[s].K(reb, sub)[j]; }, v, e);
                h->SetBinContent(j + 1, v); h->SetBinError(j + 1, e);
            }
            h->Write();
            for (int c = 0; c < 3; ++c) {
                const char* cn[3] = {"pm", "pp", "mm"};
                TH1D* hm = new TH1D(Form("mix_A_%s_%s", kSelKey[s], cn[c]), "", nb, scan::kMlo, scan::kMhi);
                for (int i = 0; i < nb; ++i) hm->SetBinContent(i + 1, A[s].Sum(c, i));
                hm->Write();
            }
        }
        // bump E with the slice k
        P("\n   bump edge-step sum E with k_A of the selection AND slice (rows per 2.5 MeV):\n");
        for (int q = 0; q < 3; ++q) {
            const int s = kSlSel[q];
            std::vector<Sel> ss;
            for (int i = 0; i < kNSig; ++i) {
                const Sel base = sv[s]; const float lo = kSlLo[i], hi = kSlHi[i];
                ss.push_back([=](const Row& r) { return r.mpipi >= lo && r.mpipi < hi && base(r); });
            }
            const std::vector<Spectra> sl = scan::FillMany(d, ss);
            std::vector<std::vector<double>> ksl(kNSig);
            for (int i = 0; i < kNSig; ++i) {
                ksl[i] = Asl[q][i].K(reb);
                TH1D* h = new TH1D(Form("k_A_%s_sl%d", kSelKey[s], i), ";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];k", nw, scan::kMlo, scan::kMhi);
                for (int j = 0; j < nw; ++j) {
                    double v, e; SubErr(Asl[q][i], [&](int sub) { return Asl[q][i].K(reb, sub)[j]; }, v, e);
                    h->SetBinContent(j + 1, v); h->SetBinError(j + 1, e);
                }
                h->Write();
            }
            const EVal e0 = EdgeSum(sl), e1 = EdgeSum(ModifyAll(sl, 1.0, 0, &ksl));
            std::vector<std::vector<double>> kall(kNSig, kA[s]);
            const EVal e2 = EdgeSum(ModifyAll(sl, 1.0, 0, &kall));
            P("   %-24s | E %7.1f +- %5.1f -> %7.1f (slice k, %+6.3f), %7.1f (selection k, %+6.3f)\n", sels[sIdx[s]].name.c_str(), e0.E, e0.eE, e1.E,
              (e1.E - e0.E) / e0.E, e2.E, (e2.E - e0.E) / e0.E);
            P("      slice k in the edge region (CB-weighted 0.47-0.53):");
            for (int i = 0; i < kNSig; ++i) P(" %.4f", kAvg(ksl[i], sl[i], 0.47, 0.53));
            P("\n");
        }
        // variant histograms
        // bin errors: spread over the sub-samples where they exist (Ange3, B), else 0
        auto writeK = [&](const char* name, const std::vector<double>& k, const Acc* ac) {
            TH1D* h = new TH1D(name, ";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];k", nw, scan::kMlo, scan::kMhi);
            for (int j = 0; j < nw; ++j) {
                h->SetBinContent(j + 1, k[j]);
                double v = k[j], e = 0;
                if (ac) SubErr(*ac, [&](int sub) { return ac->K(reb, sub)[j]; }, v, e);
                h->SetBinError(j + 1, e);
            }
            h->Write();
            return h;
        };
        TH1D* hNoz = writeK("k_Anoz_pre", kNoz, nullptr);
        TH1D* hN3 = writeK("k_Ange3_pre", kN3, &Ange3[0]);
        TH1D* hB = writeK("k_B_pre", kB, &Bv[0]);
        TH1D* hC = writeK("k_C_pre", kC, nullptr);
        for (int s : {5, 7}) writeK(Form("k_Anoz_%s", kSelKey[s]), Anoz[s].K(reb), nullptr);
        TH1D* hee = new TH1D("k_ee", ";M(e^{+}e^{-}) [GeV/c^{2}];k_{ee}", 32, 0.0, 0.80);
        {
            for (int j = 0; j < 32; ++j) {
                double v, e; SubErr(aC, [&](int sub) { return aC.K(5, sub)[j]; }, v, e);
                hee->SetBinContent(j + 1, v); hee->SetBinError(j + 1, e);
            }
        }
        hee->Write();
        TH1D* hoa = new TH1D("k_A_pre_oa", ";OA(e^{+}e^{-}) [deg];k", 30, 0.0, 120.);
        for (int j = 0; j < 30; ++j) { double v, e; SubErr(aOa, [&](int sub) { return aOa.K(2, sub)[j]; }, v, e); hoa->SetBinContent(j + 1, v); hoa->SetBinError(j + 1, e); }
        hoa->Write();
        TH1D* hmee = new TH1D("k_A_pre_mee", ";M(e^{+}e^{-}) [GeV/c^{2}];k", 32, 0.0, 0.80);
        for (int j = 0; j < 32; ++j) { double v, e; SubErr(aMee, [&](int sub) { return aMee.K(5, sub)[j]; }, v, e); hmee->SetBinContent(j + 1, v); hmee->SetBinError(j + 1, e); }
        hmee->Write();

        // ======================================================================
        // figure
        // ======================================================================
        Style();
        TCanvas* c = new TCanvas("cb_kmix", "", 1500, 520);
        c->Divide(3, 1, 0.002, 0.002);
        auto frame = [](const char* n, double a, double b, double lo, double hi, const char* xt, const char* yt) {
            TH1D* f = new TH1D(n, "", 10, a, b); f->SetMinimum(lo); f->SetMaximum(hi);
            f->GetXaxis()->SetTitle(xt); f->GetYaxis()->SetTitle(yt); f->GetYaxis()->SetTitleOffset(1.35); f->GetXaxis()->SetNdivisions(505);
            f->Draw("AXIS"); return f;
        };
        auto band = [](double a, double b, double y0, double y1) {
            for (double y : {y0, y1}) { TLine* l = new TLine(a, y, b, y); l->SetLineStyle(2); l->SetLineColor(kGray + 1); l->Draw(); }
            TLine* l = new TLine(a, 1, b, 1); l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->Draw();
        };
        const char* xm = "M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}]";
        // (a) variants for pre
        c->cd(1); gPad->SetLeftMargin(0.16);
        frame("frk1", 0.35, 0.95, 0.85, 1.15, xm, "k = N_{+-}^{mix} / 2#sqrt{N_{++}^{mix}N_{--}^{mix}}");
        band(0.35, 0.95, 0.97, 1.03);
        TH1D* hA = (TH1D*)gDirectory->Get("k_A_pre");
        TLegend* l1 = new TLegend(0.40, 0.17, 0.97, 0.45); l1->SetNColumns(1); l1->SetBorderSize(0); l1->SetFillStyle(0); l1->SetTextSize(0.038);
        {
            struct It { TH1D* h; int col, mk; const char* n; };
            const It it[5] = {{hA, kBlack, 20, "A: (#pi#pi l)_{E} + l_{B}"}, {hNoz, kGray + 1, 24, "A, one z pool"}, {hN3, kOrange + 7, 26, "A, #geq 3 leptons"},
                              {hB, kGreen + 2, 25, "B: #pi#pi_{E} + l_{X} + l_{Y}"}, {hC, kMagenta + 1, 27, "C: k_{ee} via like-sign M_{ee}"}};
            for (int q = 0; q < 5; ++q) {
                TGraphErrors* g = new TGraphErrors();
                for (int j = 0; j < nw; ++j) {
                    g->SetPoint(j, it[q].h->GetBinCenter(j + 1) + (q - 2) * 0.0012, it[q].h->GetBinContent(j + 1));
                    g->SetPointError(j, 0, it[q].h->GetBinError(j + 1));
                }
                Deco(g, it[q].col, it[q].mk, 0.6); g->SetLineWidth(1); g->Draw("P same");
                l1->AddEntry(g, it[q].n, "p");
            }
        }
        l1->Draw();
        Txt(0.17, 0.935, "(a) data, pre: event-mixing k, variants", 0.046);
        // (b) selections, variant A
        c->cd(2); gPad->SetLeftMargin(0.16);
        frame("frk2", 0.35, 0.95, 0.85, 1.15, xm, "k (variant A)");
        band(0.35, 0.95, 0.97, 1.03);
        TLegend* l2 = new TLegend(0.40, 0.17, 0.97, 0.47); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.038);
        {
            const int ss[6] = {0, 1, 3, 6, 7, 8};
            const int col[6] = {kBlack, kGray + 1, kGreen + 2, kAzure + 7, kOrange + 7, kRed + 1};
            const int mk[6] = {20, 24, 22, 25, 26, 32};
            for (int q = 0; q < 6; ++q) {
                TH1D* h = (TH1D*)gDirectory->Get(Form("k_A_%s", kSelKey[ss[q]]));
                TGraphErrors* g = new TGraphErrors();
                for (int j = 0; j < nw; ++j) { g->SetPoint(j, h->GetBinCenter(j + 1) + (q - 2.5) * 0.001, h->GetBinContent(j + 1)); g->SetPointError(j, 0, h->GetBinError(j + 1)); }
                Deco(g, col[q], mk[q], 0.6); g->SetLineWidth(1); g->Draw("P same");
                l2->AddEntry(g, sels[sIdx[ss[q]]].name.c_str(), "p");
            }
        }
        l2->Draw();
        Txt(0.17, 0.935, "(b) data: k for the selections", 0.046);
        // (c) k against M(ee) and OA
        c->cd(3); gPad->SetLeftMargin(0.16);
        frame("frk3", 0.0, 0.80, 0.85, 1.25, "M(e^{+}e^{-}) [GeV/c^{2}]", "k_{ee}");
        band(0.0, 0.80, 0.97, 1.03);
        {
            TGraphErrors* g1 = new TGraphErrors();
            for (int j = 0; j < 32; ++j) { g1->SetPoint(j, hee->GetBinCenter(j + 1), hee->GetBinContent(j + 1)); g1->SetPointError(j, 0, hee->GetBinError(j + 1)); }
            Deco(g1, kMagenta + 1, 20, 0.7); g1->SetLineWidth(1);
            g1->Draw("P same");
            // e+e- channel, data, PT3, classes n = 2 and >= 3 (memory note project_isbest_choice_and_k, 2026-11-09)
            const double eLo[4] = {0.05, 0.14, 0.20, 0.30}, eHi[4] = {0.14, 0.20, 0.30, 0.56};
            const double eK[4] = {1.029, 1.069, 1.011, 1.000}, eE[4] = {0.011, 0.007, 0.006, 0.008};
            TGraphErrors* g2 = new TGraphErrors();
            for (int q = 0; q < 4; ++q) { g2->SetPoint(q, 0.5 * (eLo[q] + eHi[q]), eK[q]); g2->SetPointError(q, 0.5 * (eHi[q] - eLo[q]), eE[q]); }
            Deco(g2, kRed, 21, 0.9); g2->SetLineWidth(2);
            g2->Draw("P same");
            TLegend* l3 = new TLegend(0.19, 0.70, 0.97, 0.90); l3->SetBorderSize(0); l3->SetFillStyle(0); l3->SetTextSize(0.038);
            l3->AddEntry(g1, "this channel, mixed lepton pairs (C)", "pe");
            l3->AddEntry(g2, "e^{+}e^{-} study (data, k_{mix} with classes)", "pe");
            l3->Draw();
        }
        Txt(0.17, 0.935, "(c) data: k_{ee} of the lepton pair alone", 0.046);
        Save(c, "cb_kmix");
        c->Write();
        fo2.Close();
    }

    // ---- LaTeX table: CB-weighted k per region and selection, with the effect on S, W ----
    {
        FILE* t = fopen(Form("%s/cb_kmix.tex", kTabOut), "w");
        fprintf(t, "%% written by plots/exp_report/cb_kmix.C\n");
        fprintf(t, "\\begin{table}[H]\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}\n");
        fprintf(t, "\\caption{Event-mixing correction $k=N_{+-}^\\text{mix}/(2\\sqrt{N_{++}^\\text{mix}N_{--}^\\text{mix}})$ of the like-sign CB, experimental data, "
                   "isBest $\\geq 0$. Mixed combination: pion pair and one lepton of an event, the other lepton from one of the last %d events of the same "
                   "vertex-$z$ class (20~mm); the selection is applied to the mixed combination. Mean of $k$ weighted with the data CB of the selection in "
                   "regions of $M(\\pi^+\\pi^-e^+e^-)$ [GeV/$c^2$]; errors from 10 sub-samples. Last columns: relative change of the $\\eta$-peak content $S$ "
                   "($P_\\text{main}$) and of $N-\\text{CB}$ in the bump window ($W$) when the CB is multiplied by $k(M)$. Not normalised to luminosity.}\n", depth);
        fprintf(t, "\\label{tab:cb_kmix}\n\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{lrrrrrrrr}\n\\toprule\n");   // [review fix] 0.56-0.60 column added (quoted in the text)
        fprintf(t, "selection & 0.35--0.455 & bump & 0.5225--0.56 & 0.56--0.60 & 0.60--0.70 & 0.70--0.95 & $\\Delta S/S$ & $\\Delta W/W$\\\\\n\\midrule\n");
        for (int s = 0; s < kNSelM; ++s) {
            fprintf(t, "%s", TexName(sels[sIdx[s]].name).c_str());
            for (int q : {0, 1, 2, 3, 4, 5}) {
                double v, e; SubErr(A[s], [&](int sub) { return kAvg(A[s].K(reb, sub), sp[s], rA[q], rB[q]); }, v, e);
                fprintf(t, " & $%.3f\\pm%.3f$", v, e);
            }
            const PeakFit a = scan::FitPeakMain(sp[s]);
            const Spectra m = Modify(sp[s], 1.0, 0, &kA[s]);
            const PeakFit b = scan::FitPeakMain(m);
            fprintf(t, " & %+.2f\\%% & %+.2f\\%%\\\\\n", 100. * (b.S - a.S) / a.S, 100. * (m.win[1].sig() - sp[s].win[1].sig()) / sp[s].win[1].sig());
        }
        fprintf(t, "\\bottomrule\n\\end{tabular}}\n\\end{table}\n");
        fclose(t);
    }
    P("\nrun time %.1f s\n", swAll.RealTime());
    fclose(OutFile()); OutFile() = nullptr;
}

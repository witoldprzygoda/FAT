// ============================================================================
// wp2_sector.C  --  WP2: sector-2 trigger deficit of days 034-044 in the pi+pi-e+e- rows (data)
// ============================================================================
// What it does (plain words):
//   1. Every track gets a HADES sector from its azimuth (phi of the compact is in (-180, 180]; it is
//      brought to [0, 360) first): sector = (floor(phi/60) + 5) mod 6
//      (phi 180-240 deg = sector 2). The boundaries are checked on the gaps (frames) of the
//      measured phi distribution.
//   2. Rows of the data compact (us, pp, mm) are counted per day, per M4 window and per
//      "sector use" category (row has a track in sector s; no track in sector 2; one or two
//      tracks in sector 2; the sector-2 track is a pion or a lepton), for the selections
//      S0 (nominal), pre, none and S0 with isBest == 1. Errors: event-group variance.
//   3. Rates = rows per nb^-1 with the luminosity of the analysed runs of each day
//      (wp2_lumi_perday.txt, L6_an). The deficit of sector 2 is measured as
//          eps_single = R_2(early)/R_2(late)                (rows with a sector-2 track)
//          eps_double = eps_single / [R_no2(early)/R_no2(late)]  (relative to rows without one)
//      early = days 034-044 (the deficit days), late = 045-067; variant late = 045-046 only
//      (same period I). Also per day and per 10-deg phi bin (the deficit must sit at 180-240).
//   4. Like-sign rows give CB = 2 sqrt(N++ N--) per category and window, so the same ratios
//      are also formed for N - CB in the bump (B) and peak (P) windows.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_sector.root (histograms for Fig. N1),
//         docs/latex/tables_all/wp2_sector.txt (all numbers), wp2_sector_eps.txt (machine-readable summary)
// Run:    root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)' \
//              'plots/all_report/wp2_sector.C+O'
// ============================================================================
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include "TSystem.h"

#include "defs.h"

namespace wp2s {

const char* const kOutDir = "/hdd3/przygoda/HADES/FAT_out/all_report/wp2";
const char* const kTabDir = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all";
const char* const kLumiDay = "/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp2_lumi_perday.txt";

// the compact stores phi = TVector3::Phi() in deg, (-180, 180] (FAT PParticle::phi()); HADES phi = that mod 360
inline float Phi360(float phi) { return phi < 0.f ? phi + 360.f : phi; }
inline int Sec(float phi) {
    if (!(phi == phi)) return -1;
    phi = Phi360(phi);
    int k = (int)std::floor(phi / 60.);
    return ((k + 5) % 6 + 6) % 6;
}

const int kND = 29;
const int kDays[kND] = {34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
                        59, 60, 61, 62, 63, 64, 65, 66, 67};
inline int DayIdx(int d) { for (int i = 0; i < kND; ++i) if (kDays[i] == d) return i; return -1; }
inline bool Early(int d) { return d >= 34 && d <= 45; }   // phi histograms: deficit span 034-045 (day 045 has deficit blocks)

// selections
enum { kS0 = 0, kPre, kSnone, kS0best, kNS };
const char* kSName[kNS] = {"S0", "pre", "none", "S0_best"};
// windows
enum { kWall = 0, kWB, kWP, kWCRU, kWlo, kWhi, kNW };
const char* kWName[kNW] = {"all M4", "B 0.455-0.5225", "P 0.5336-0.5542", "CR-U 0.565-0.70", "M4 < 0.455", "M4 >= 0.70"};
inline void Wins(float m4, int* w, int& n) {
    n = 0;
    w[n++] = kWall;
    if (m4 >= alldefs::kBlo && m4 < alldefs::kBhi) w[n++] = kWB;
    if (m4 >= alldefs::kPlo && m4 < alldefs::kPhi) w[n++] = kWP;
    if (m4 >= 0.565f && m4 < 0.700f) w[n++] = kWCRU;
    if (m4 < 0.455f) w[n++] = kWlo;
    if (m4 >= 0.700f) w[n++] = kWhi;
}
// categories
// A_k = no sector-2 track and k distinct sectors used by the four tracks; B_k = >= 1 sector-2 track and k distinct OTHER
// sectors; A34 = A_3 + A_4; Imm = rows with >= 3 distinct non-2 sectors (A_3, A_4, B_3): they keep >= 3 fired sectors
// without sector 2 (immune if PT3 is a ">= 3 sectors" condition); NImm = the other rows
enum { kCtot = 0, kCs0, kCs1, kCs2, kCs3, kCs4, kCs5, kCno2, kC2one, kC2two, kC2pi, kC2lep, kC2both,
       kCx0, kCx1, kCx2unused, kCx3, kCx4, kCx5,
       kA1, kA2, kA3, kA4, kA34, kB0, kB1, kB2, kB3, kImm, kNImm, kNC };
const char* kCName[kNC] = {"total", "uses s0", "uses s1", "uses s2", "uses s3", "uses s4", "uses s5", "no s2 track",
                           "one s2 track", ">=2 s2 tracks", "s2: pion(s) only", "s2: lepton(s) only", "s2: pion and lepton",
                           "s0 not s2", "s1 not s2", "-", "s3 not s2", "s4 not s2", "s5 not s2",
                           "A1 no s2, 1 sector", "A2 no s2, 2 sectors", "A3 no s2, 3 sectors", "A4 no s2, 4 sectors",
                           "A34 no s2, >=3 sectors", "B0 s2 only", "B1 s2 + 1 other", "B2 s2 + 2 others", "B3 s2 + 3 others",
                           "Imm >=3 non-2 sectors", "NImm <3 non-2 sectors"};

const int kNCell = kNS * 3 * kND * kNW * kNC;
inline int Cell(int s, int smp, int d, int w, int c) { return (((s * 3 + smp) * kND + d) * kNW + w) * kNC + c; }

struct Acc {
    std::vector<double> sum, var, buf;
    std::vector<int> touched;
    Acc() : sum(kNCell, 0.), var(kNCell, 0.), buf(kNCell, 0.) {}
    void Add(int i) { if (buf[i] == 0.) touched.push_back(i); buf[i] += 1.; sum[i] += 1.; }
    void Flush() { for (int i : touched) { var[i] += buf[i] * buf[i]; buf[i] = 0.; } touched.clear(); }
};

std::map<int, double> ReadLumi() {
    std::map<int, double> L;
    std::ifstream in(kLumiDay);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.compare(0, 3, "day") == 0) continue;
        std::stringstream ss(line);
        std::string tok;
        std::vector<std::string> t;
        while (std::getline(ss, tok, ',')) t.push_back(tok);
        if (t.size() < 6) continue;
        L[std::stoi(t[0])] = std::stod(t[5]);   // L6_an
    }
    return L;
}

}  // namespace wp2s

void wp2_sector() {
    using namespace wp2s;
    TStopwatch sw;
    gSystem->mkdir(kOutDir, kTRUE);
    std::map<int, double> Lday = ReadLumi();
    if ((int)Lday.size() != kND) { printf("wp2_sector: luminosity table has %zu days\n", Lday.size()); return; }
    double Lsum = 0;
    for (auto& x : Lday) Lsum += x.second;
    printf("wp2_sector: L(analysed runs) = %.3f nb^-1\n", Lsum);

    scan::Data D;
    if (!scan::Load(D)) return;

    Acc A;
    // phi histograms (S0, us): all tracks, early / late; per track type; day x phi map
    TH1D* hPhiE = new TH1D("hPhiE", "S0 e^{+}e^{-} rows, all four tracks, days 034-045;#phi [deg];tracks / 0.5 deg", 720, 0, 360);
    TH1D* hPhiL = new TH1D("hPhiL", "S0 e^{+}e^{-} rows, all four tracks, days 046-067;#phi [deg];tracks / 0.5 deg", 720, 0, 360);
    TH1D* hPhiT[4];
    const char* tn[4] = {"pip", "pim", "ep", "em"};
    for (int k = 0; k < 4; ++k) hPhiT[k] = new TH1D(Form("hPhi_%s", tn[k]), Form("S0, %s;#phi [deg];tracks / 0.5 deg", tn[k]), 720, 0, 360);
    TH2D* hDayPhi = new TH2D("hDayPhi", "S0 e^{+}e^{-} rows;day index;#phi [deg]", kND, 0, kND, 36, 0, 360);
    TH2D* hDayPhiLS = new TH2D("hDayPhiLS", "S0 like-sign rows;day index;#phi [deg]", kND, 0, kND, 36, 0, 360);
    long nBadPhi = 0, nPhiOut = 0;

    for (int smp = 0; smp < 3; ++smp) {
        const std::vector<scan::Row>& v = D.s[smp];
        int cur = -1;
        for (size_t i = 0; i < v.size(); ++i) {
            const scan::Row& r = v[i];
            if (r.grp != cur) { A.Flush(); cur = r.grp; }
            const int di = DayIdx(r.file_id);
            if (di < 0) continue;
            bool pass[kNS];
            const bool corOk = r.cor_ok == 1;
            pass[kS0] = alldefs::PassSel(alldefs::kSelS0, r, corOk, r.dz_pi, r.mm4);
            pass[kPre] = alldefs::PassSel(alldefs::kSelPre, r, corOk, r.dz_pi, r.mm4);
            pass[kSnone] = alldefs::PassSel(alldefs::kSelNone, r, corOk, r.dz_pi, r.mm4);
            pass[kS0best] = pass[kS0] && r.isBest == 1;
            if (!pass[kSnone] && !pass[kPre]) continue;
            const float phi[4] = {Phi360(r.pip_phi), Phi360(r.pim_phi), Phi360(r.ep_phi), Phi360(r.em_phi)};
            int sec[4];
            bool bad = false;
            for (int k = 0; k < 4; ++k) {
                sec[k] = Sec(phi[k]);
                if (sec[k] < 0) bad = true;
                if (phi[k] < 0.f || phi[k] >= 360.f) nPhiOut++;
            }
            if (bad) { nBadPhi++; continue; }
            int use = 0, n2 = 0, n2pi = 0, n2lep = 0;
            for (int k = 0; k < 4; ++k) {
                use |= 1 << sec[k];
                if (sec[k] == 2) { ++n2; if (k < 2) ++n2pi; else ++n2lep; }
            }
            int cats[kNC], nc = 0;
            cats[nc++] = kCtot;
            for (int s = 0; s < 6; ++s) if (use & (1 << s)) cats[nc++] = kCs0 + s;
            int kOther = 0;
            for (int q = 0; q < 6; ++q) if (q != 2 && (use & (1 << q))) ++kOther;
            if (n2 == 0) { cats[nc++] = kA1 + (kOther - 1); if (kOther >= 3) cats[nc++] = kA34; }
            else cats[nc++] = kB0 + kOther;
            cats[nc++] = kOther >= 3 ? kImm : kNImm;
            if (n2 == 0) {
                cats[nc++] = kCno2;
                for (int s = 0; s < 6; ++s) if (s != 2 && (use & (1 << s))) cats[nc++] = kCx0 + s;
            } else {
                cats[nc++] = n2 == 1 ? kC2one : kC2two;
                cats[nc++] = n2lep == 0 ? kC2pi : (n2pi == 0 ? kC2lep : kC2both);
            }
            int w[kNW], nw = 0;
            Wins(r.m4, w, nw);
            for (int s = 0; s < kNS; ++s) {
                if (!pass[s]) continue;
                for (int a = 0; a < nw; ++a)
                    for (int c = 0; c < nc; ++c) A.Add(Cell(s, smp, di, w[a], cats[c]));
            }
            if (pass[kS0]) {
                for (int k = 0; k < 4; ++k) {
                    if (smp == 0) {
                        (Early(r.file_id) ? hPhiE : hPhiL)->Fill(phi[k]);
                        hPhiT[k]->Fill(phi[k]);
                        hDayPhi->Fill(di + 0.5, phi[k]);
                    } else {
                        hDayPhiLS->Fill(di + 0.5, phi[k]);
                    }
                }
            }
        }
        A.Flush();
    }
    printf("wp2_sector: filled in %.1f s; rows with NaN phi skipped %ld; phi outside [0,360): %ld\n", sw.RealTime(), nBadPhi, nPhiOut);
    sw.Continue();

    // ------------------------------------------------------------------ helpers
    auto Sum = [&](int s, int smp, int w, int c, bool (*dsel)(int)) {
        double n = 0, v = 0, L = 0;
        for (int d = 0; d < kND; ++d) {
            if (!dsel(kDays[d])) continue;
            n += A.sum[Cell(s, smp, d, w, c)];
            v += A.var[Cell(s, smp, d, w, c)];
            L += Lday[kDays[d]];
        }
        return std::array<double, 3>{n, v, L};
    };
    // N - CB of a group of days: CB = 2 sqrt(N++ N--) of the summed like-sign counts (window-level estimate)
    struct Val { double n, v, L; };
    auto Sig = [&](int s, int w, int c, bool (*dsel)(int), bool sub) {
        auto u = Sum(s, 0, w, c, dsel);
        Val x{u[0], u[1], u[2]};
        if (sub) {
            auto p = Sum(s, 1, w, c, dsel), m = Sum(s, 2, w, c, dsel);
            const double cb = (p[0] > 0 && m[0] > 0) ? 2. * std::sqrt(p[0] * m[0]) : 0.;
            // d CB/d Npp = sqrt(Nmm/Npp), d CB/d Nmm = sqrt(Npp/Nmm)
            const double vcb = (p[0] > 0 && m[0] > 0) ? (m[0] / p[0]) * p[1] + (p[0] / m[0]) * m[1] : 0.;
            x.n -= cb;
            x.v += vcb;
        }
        return x;
    };
    auto IsEarly = [](int d) { return d >= 34 && d <= 44; };
    auto IsLate = [](int d) { return d >= 45; };
    auto IsLateI = [](int d) { return d == 45 || d == 46; };
    auto IsLateII = [](int d) { return d >= 47 && d <= 53; };
    auto IsLateIII = [](int d) { return d >= 59; };
    // ratio of rates (a/La)/(b/Lb) with uncorrelated errors (different days)
    auto RR = [](const Val& a, const Val& b, double& r, double& e) {
        r = (a.n / a.L) / (b.n / b.L);
        e = r * std::sqrt(a.v / (a.n * a.n) + b.v / (b.n * b.n));
    };

    FILE* tab = fopen(Form("%s/wp2_sector.txt", kTabDir), "w");
    FILE* sum = fopen(Form("%s/wp2_sector_eps.txt", kTabDir), "w");
    fprintf(tab, "# wp2_sector.C: sector-2 trigger deficit (days 034-044) in the pi+pi-e+e- rows; data scan_v2 compact\n");
    fprintf(tab, "# sector = (floor(phi/60) + 5) mod 6; rates = rows / L6_an [nb^-1] (wp2_lumi_perday.txt); errors = event-group variance\n");
    fprintf(tab, "# early = 034-044 (L %.3f), late = 045-067 (L %.3f), lateI = 045-046 (L %.3f)\n",
            Sum(0, 0, 0, 0, IsEarly)[2], Sum(0, 0, 0, 0, IsLate)[2], Sum(0, 0, 0, 0, IsLateI)[2]);
    fprintf(sum, "# wp2_sector.C summary: key,value,error (eps = efficiency of rows with a sector-2 track in days 034-044)\n");

    // ---------------------------------------------------------------- phi gaps
    fprintf(tab, "\n[phi gaps] S0 e+e- rows, all tracks, 0.5 deg bins: content around the sector boundaries k x 60 deg\n");
    fprintf(tab, "# boundary  min position within +-6 deg  content(+-1 deg)/mean content(+-10..20 deg)\n");
    {
        TH1D* h = (TH1D*)hPhiE->Clone("hPhiAll");
        h->Add(hPhiL);
        for (int k = 0; k < 6; ++k) {
            const double b = 60. * k;
            double cmin = 1e30, xmin = 0, cin = 0, cout = 0;
            int nin = 0, nout = 0;
            for (int j = 1; j <= 720; ++j) {
                double x = h->GetBinCenter(j);
                double dx = x - b;
                if (dx > 180) dx -= 360;
                if (dx < -180) dx += 360;
                const double c = h->GetBinContent(j);
                if (std::fabs(dx) <= 6 && c < cmin) { cmin = c; xmin = dx; }
                if (std::fabs(dx) <= 1) { cin += c; ++nin; }
                if (std::fabs(dx) >= 10 && std::fabs(dx) <= 20) { cout += c; ++nout; }
            }
            fprintf(tab, "  %5.0f deg  %+5.2f deg  %.3f   (sectors %d | %d)\n", b, xmin, (cin / nin) / (cout / nout), Sec(b - 1e-3), Sec(b + 1e-3));
        }
        delete h;
    }

    // ---------------------------------------------------------------- per sector early/late (all M4, S0)
    fprintf(tab, "\n[sector use] early/late rate ratios, all M4 (N = unlike-sign rows; N-CB with the like-sign estimate)\n");
    fprintf(tab, "# sel  category            N_early   N_late   R_early   R_late   ratio_early/late (stat)   ratio/ratio(no s2)\n");
    for (int s = 0; s < kNS; ++s) {
        Val no2e = Sig(s, kWall, kCno2, IsEarly, false), no2l = Sig(s, kWall, kCno2, IsLate, false);
        double rno2, eno2;
        RR(no2e, no2l, rno2, eno2);
        for (int c = 0; c < kNC; ++c) {
            if (c == kCx2unused) continue;
            Val e = Sig(s, kWall, c, IsEarly, false), l = Sig(s, kWall, c, IsLate, false);
            if (e.n <= 0 || l.n <= 0) continue;
            double r, er;
            RR(e, l, r, er);
            fprintf(tab, "  %-7s %-20s %9.0f %9.0f %8.4f %8.4f   %.4f +- %.4f   %.4f\n", kSName[s], kCName[c], e.n, l.n, e.n / e.L,
                    l.n / l.L, r, er, r / rno2);
        }
    }

    // ---------------------------------------------------------------- eps: single, double, by window, N and N-CB, lateI
    fprintf(tab, "\n[eps] efficiency of rows with a sector-2 track in days 034-044\n");
    fprintf(tab, "# eps_single = R_2(early)/R_2(late); eps_double = eps_single / [R_no2(early)/R_no2(late)]\n");
    fprintf(tab, "# sel     window              quantity  late-ref  eps_single (stat)    ref ratio no-s2 (stat)   eps_double (stat)   f_2(late) = share with s2\n");
    struct LateDef { const char* tag; bool (*f)(int); };
    LateDef lates[4] = {{"045-067", IsLate}, {"045-046", IsLateI}, {"047-053", IsLateII}, {"059-067", IsLateIII}};
    for (int s = 0; s < kNS; ++s) {
        for (int w = 0; w < kNW; ++w) {
            for (int sub = 0; sub < 2; ++sub) {
                for (int li = 0; li < 4; ++li) {
                    if (li > 0 && !(w == kWall && sub == 0)) continue;
                    Val e2 = Sig(s, w, kCs2, IsEarly, sub), l2 = Sig(s, w, kCs2, lates[li].f, sub);
                    Val e0 = Sig(s, w, kCno2, IsEarly, sub), l0 = Sig(s, w, kCno2, lates[li].f, sub);
                    Val lt = Sig(s, w, kCtot, lates[li].f, sub);
                    if (e2.n <= 0 || l2.n <= 0 || e0.n <= 0 || l0.n <= 0) continue;
                    double r2, e2r, r0, e0r;
                    RR(e2, l2, r2, e2r);
                    RR(e0, l0, r0, e0r);
                    const double dbl = r2 / r0, edbl = dbl * std::sqrt(std::pow(e2r / r2, 2) + std::pow(e0r / r0, 2));
                    const double f2 = l2.n / lt.n;
                    fprintf(tab, "  %-7s %-18s %-8s %-8s %.4f +- %.4f     %.4f +- %.4f     %.4f +- %.4f    %.4f\n", kSName[s], kWName[w],
                            sub ? "N-CB" : "N", lates[li].tag, r2, e2r, r0, e0r, dbl, edbl, f2);
                    fprintf(sum, "%s|%d|%s|%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", kSName[s], w, sub ? "NmCB" : "N", lates[li].tag, r2, e2r,
                            r0, e0r, dbl, edbl, f2);
                }
            }
        }
    }

    // ---------------------------------------------------------------- per category of the sector-2 tracks (S0, all M4)
    fprintf(tab, "\n[eps by kind] S0, all M4, N: eps_double for rows by the kind of their sector-2 track(s)\n");
    {
        Val e0 = Sig(kS0, kWall, kCno2, IsEarly, false), l0 = Sig(kS0, kWall, kCno2, IsLate, false);
        double r0, e0r;
        RR(e0, l0, r0, e0r);
        const int cs[5] = {kC2one, kC2two, kC2pi, kC2lep, kC2both};
        for (int c : cs) {
            Val e = Sig(kS0, kWall, c, IsEarly, false), l = Sig(kS0, kWall, c, IsLate, false);
            double r, er;
            RR(e, l, r, er);
            fprintf(tab, "  %-22s N_early %8.0f  eps_double %.4f +- %.4f   share of s2 rows (late) %.4f\n", kCName[c], e.n, r / r0,
                    (r / r0) * std::sqrt(std::pow(er / r, 2) + std::pow(e0r / r0, 2)), l.n / Sig(kS0, kWall, kCs2, IsLate, false).n);
            fprintf(sum, "kind|%s,%.6f,%.6f\n", kCName[c], r / r0, (r / r0) * std::sqrt(std::pow(er / r, 2) + std::pow(e0r / r0, 2)));
        }
    }

    // ---------------------------------------------------------------- trigger model: eps per sector pattern, reference A34
    fprintf(tab, "\n[trigger model] eps_c = (n_c/n_ref)(early)/(n_c/n_ref)(late), early 034-045 and 034-044, late 046-067 / 047-067, reference A34\n");
    fprintf(tab, "# (no sector-2 track, >= 3 distinct sectors); share = n_c/n_tot in the late days; also the reference rate ratio per nb^-1\n");
    {
        auto E45 = [](int d) { return d >= 34 && d <= 45; };
        auto L46 = [](int d) { return d >= 46; };
        auto E44 = [](int d) { return d >= 34 && d <= 44; };
        auto L47 = [](int d) { return d >= 47; };
        struct Sp { const char* tag; bool (*e)(int); bool (*l)(int); };
        Sp sps[2] = {{"034-045 vs 046-067", E45, L46}, {"034-044 vs 047-067", E44, L47}};
        for (int s = 0; s < kNS; ++s) {
            if (s == kSnone) continue;
            for (int w : {kWall, kWB, kWP}) {
                for (int sub = 0; sub < 2; ++sub) {
                    for (int q = 0; q < 2; ++q) {
                        if (q == 1 && !(s == kS0 && w == kWall && sub == 0)) continue;
                        Val re = Sig(s, w, kA34, sps[q].e, sub), rl = Sig(s, w, kA34, sps[q].l, sub);
                        Val tl = Sig(s, w, kCtot, sps[q].l, sub);
                        double rr0, err0;
                        RR(re, rl, rr0, err0);
                        fprintf(tab, "  %-7s %-16s %-4s %s: reference A34 rate ratio early/late %.4f +- %.4f\n", kSName[s], kWName[w],
                                sub ? "N-CB" : "N", sps[q].tag, rr0, err0);
                        const int cc[13] = {kA1, kA2, kA3, kA4, kB0, kB1, kB2, kB3, kImm, kNImm, kCs2, kC2one, kC2two};
                        for (int c : cc) {
                            Val e = Sig(s, w, c, sps[q].e, sub), l = Sig(s, w, c, sps[q].l, sub);
                            if (e.n <= 0 || l.n <= 0) continue;
                            double r, er;
                            RR(e, l, r, er);
                            const double ep = r / rr0, eep = ep * std::sqrt(std::pow(er / r, 2) + std::pow(err0 / rr0, 2));
                            fprintf(tab, "      %-24s share(late) %.4f   eps %.4f +- %.4f\n", kCName[c], l.n / tl.n, ep, eep);
                            fprintf(sum, "tm|%s|%d|%s|%d|%s,%.6f,%.6f,%.6f\n", kSName[s], w, sub ? "NmCB" : "N", q, kCName[c], ep, eep, l.n / tl.n);
                        }
                    }
                }
            }
        }
    }

    // ---------------------------------------------------------------- per day (S0, all M4, N)
    fprintf(tab, "\n[per day] S0, all M4, N (unlike-sign rows): rates per nb^-1 and eps_day = [R_2(day)/R_2(late)] / [R_no2(day)/R_no2(late)]\n");
    fprintf(tab, "# day  L6_an   N_tot   R_tot   R_s0 R_s1 R_s2 R_s3 R_s4 R_s5   R_no2   eps_day (stat)\n");
    TH1D* hR[8];
    const char* hrn[8] = {"tot", "s0", "s1", "s2", "s3", "s4", "s5", "no2"};
    for (int k = 0; k < 8; ++k) hR[k] = new TH1D(Form("hRate_%s", hrn[k]), Form("S0 rows per nb^{-1}, %s;day index;rows / nb^{-1}", hrn[k]), kND, 0, kND);
    TH1D* hEps = new TH1D("hEpsDay", "eps_day (double ratio, S0, N);day index;#varepsilon", kND, 0, kND);
    TH1D* hEpsLS = new TH1D("hEpsDayLS", "eps_day (double ratio, S0, like-sign N++ + N--);day index;#varepsilon", kND, 0, kND);
    {
        Val l2 = Sig(kS0, kWall, kCs2, IsLate, false), l0 = Sig(kS0, kWall, kCno2, IsLate, false);
        auto lsSum = [&](int c, bool (*f)(int)) {
            auto p = Sum(kS0, 1, kWall, c, f), m = Sum(kS0, 2, kWall, c, f);
            return Val{p[0] + m[0], p[1] + m[1], p[2]};
        };
        Val l2s = lsSum(kCs2, IsLate), l0s = lsSum(kCno2, IsLate);
        for (int d = 0; d < kND; ++d) {
            const double L = Lday[kDays[d]];
            auto c = [&](int smp, int cat) { return A.sum[Cell(kS0, smp, d, kWall, cat)]; };
            auto vv = [&](int smp, int cat) { return A.var[Cell(kS0, smp, d, kWall, cat)]; };
            const int cats8[8] = {kCtot, kCs0, kCs1, kCs2, kCs3, kCs4, kCs5, kCno2};
            for (int k = 0; k < 8; ++k) {
                hR[k]->SetBinContent(d + 1, c(0, cats8[k]) / L);
                hR[k]->SetBinError(d + 1, std::sqrt(vv(0, cats8[k])) / L);
                hR[k]->GetXaxis()->SetBinLabel(d + 1, Form("%03d", kDays[d]));
            }
            Val d2{c(0, kCs2), vv(0, kCs2), L}, d0{c(0, kCno2), vv(0, kCno2), L};
            double r2, e2r, r0, e0r;
            RR(d2, l2, r2, e2r);
            RR(d0, l0, r0, e0r);
            const double eps = r2 / r0, eeps = eps * std::sqrt(std::pow(e2r / r2, 2) + std::pow(e0r / r0, 2));
            hEps->SetBinContent(d + 1, eps);
            hEps->SetBinError(d + 1, eeps);
            hEps->GetXaxis()->SetBinLabel(d + 1, Form("%03d", kDays[d]));
            Val d2s{c(1, kCs2) + c(2, kCs2), vv(1, kCs2) + vv(2, kCs2), L}, d0s{c(1, kCno2) + c(2, kCno2), vv(1, kCno2) + vv(2, kCno2), L};
            double q2, eq2, q0, eq0;
            RR(d2s, l2s, q2, eq2);
            RR(d0s, l0s, q0, eq0);
            hEpsLS->SetBinContent(d + 1, q2 / q0);
            hEpsLS->SetBinError(d + 1, (q2 / q0) * std::sqrt(std::pow(eq2 / q2, 2) + std::pow(eq0 / q0, 2)));
            hEpsLS->GetXaxis()->SetBinLabel(d + 1, Form("%03d", kDays[d]));
            fprintf(tab, "  %03d %8.3f %8.0f %7.2f  %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f  %7.2f   %.4f +- %.4f\n", kDays[d], L, c(0, kCtot),
                    c(0, kCtot) / L, c(0, kCs0) / L, c(0, kCs1) / L, c(0, kCs2) / L, c(0, kCs3) / L, c(0, kCs4) / L, c(0, kCs5) / L,
                    c(0, kCno2) / L, eps, eeps);
            fprintf(sum, "day|%03d,%.6f,%.6f\n", kDays[d], eps, eeps);
        }
    }

    // ---------------------------------------------------------------- phi profile of the deficit (S0 us, 10-deg bins)
    fprintf(tab, "\n[phi profile] S0 e+e- tracks: (034-045 / 046-067 track rate per nb^-1) per 10-deg phi bin, normalised to the mean of the bins outside 180-240\n");
    TH1D* hPhiRatio = new TH1D("hPhiRatio", "early/late track rate per 10 deg;#phi [deg];(early/late) / mean(outside sector 2)", 36, 0, 360);
    {
        auto E45 = [](int d) { return d >= 34 && d <= 45; };
        auto L46 = [](int d) { return d >= 46; };
        const double Le = Sum(0, 0, 0, 0, E45)[2], Ll = Sum(0, 0, 0, 0, L46)[2];
        double ref = 0;
        int nref = 0;
        std::vector<double> rr(36), er(36);
        for (int b = 0; b < 36; ++b) {
            double ne = 0, nl = 0;
            for (int j = 1; j <= 20; ++j) { ne += hPhiE->GetBinContent(20 * b + j); nl += hPhiL->GetBinContent(20 * b + j); }
            rr[b] = (ne / Le) / (nl / Ll);
            er[b] = rr[b] * std::sqrt(1. / ne + 1. / nl);   // track counts, Poisson (tracks of one row are correlated: indicative)
            if (b < 18 || b >= 24) { ref += rr[b]; ++nref; }
        }
        ref /= nref;
        for (int b = 0; b < 36; ++b) {
            hPhiRatio->SetBinContent(b + 1, rr[b] / ref);
            hPhiRatio->SetBinError(b + 1, er[b] / ref);
            fprintf(tab, "  %3d-%3d deg (sector %d): %.4f +- %.4f\n", 10 * b, 10 * b + 10, Sec(10. * b + 5.), rr[b] / ref, er[b] / ref);
        }
        fprintf(tab, "# reference = mean early/late ratio of the 30 bins outside 180-240 deg: %.4f\n", ref);
    }

    fprintf(tab, "\n# wp2_sector.C finished in %.1f s\n", sw.RealTime());
    fclose(tab);
    fclose(sum);

    TFile* fo = TFile::Open(Form("%s/wp2_sector.root", kOutDir), "RECREATE");
    hPhiE->Write(); hPhiL->Write();
    for (int k = 0; k < 4; ++k) hPhiT[k]->Write();
    hDayPhi->Write(); hDayPhiLS->Write();
    for (int k = 0; k < 8; ++k) hR[k]->Write();
    hEps->Write(); hEpsLS->Write(); hPhiRatio->Write();
    // the raw accumulators for later use: per (sel, sample, day, window, category)
    TH1D* hA = new TH1D("hAcc", "cells (sel,sample,day,window,category); see wp2_sector.C Cell()", kNCell, 0, kNCell);
    for (int i = 0; i < kNCell; ++i) { hA->SetBinContent(i + 1, A.sum[i]); hA->SetBinError(i + 1, std::sqrt(A.var[i])); }
    hA->Write();
    fo->Close();
    printf("wp2_sector: done in %.1f s; tables %s/wp2_sector.txt, wp2_sector_eps.txt\n", sw.RealTime(), kTabDir);
}

// ============================================================================
// wp2_win.C  --  WP2: target-window (empty-target) contribution f_win to the pi+pi-e+e- rows
// ============================================================================
// Plain words: on days 034, 047 and 065 the beam also ran on the EMPTY target cell (3195 A field; 7 + 16 + 104 runs,
// 70 + 160 + 1040 DSTs, /lustre/hades/dst/feb22/gen4/4500MeV/3195A/emptytarget/). Everything recorded then comes from
// material (cell windows, foils, residual gas, beam line), not from the LH2. The material in front of the cell
// (vertex z = -200..-170 mm, outside the LH2 cell -150..-100 mm) is hit by the same beam per beam particle in both
// kinds of runs, so the number of rows with a vertex there measures the recorded beam flux of each sample (same
// trigger, same dead time as the yields). Then, day by day,
//     b     = U_full / U_empty                       (flux ratio, full/empty, of the same day)
//     f_win = b x Y_empty / Y_full                   (share of the full-target yield Y that comes from material)
// for every selection, M4 window and quantity (N, like-sign, N - CB with CB = 2 sqrt(N++ N--) of the window).
// The three days are combined with the luminosity of the period each one represents (034 -> I, 047 -> IIA + IIB,
// 065 -> III), as in the luminosity analysis (Task E, window_fraction_final).
// Normalisers U (vertex windows in front of the cell):
//   0 ch     channel rows (us + pp + mm, base selection = isBest >= 0, vz > -500, PT3, cor_ok), vz in [-200,-170)  NOMINAL
//   1 chw    the same, vz in [-250,-170)
//   2 chfar  the same, vz in [-500,-250)
//   3 eenear EpEm_ID rows (PT3, isBest == 1; wp2_upstream.root), vz in [-200,-170)
//   4 eefar  EpEm_ID rows, vz in [-700,-250) (the region of the luminosity analysis' normalisers U1 + U2)
// The simulated LH2 rows put 1e-5 of their vertices into [-200,-170) (wp2_simvz.root): the nominal window is LH2-free.
// Peak signal: S_P = (N - CB)_P - (N - CB)_SB x w_P / w_SB (side bands 0.5250-0.5336 and 0.5542-0.5628, linear).
// Inputs: data compact scan_v2 (full) and wp2/empty/scan2/scan_compact_empty.root (empty), wp2_upstream.root.
// Output: /hdd3/.../wp2/wp2_win.root (vz and M4 spectra), docs/latex/tables_all/wp2_win.txt (per day and combined)
//         and wp2_win_summary.txt
// Run: root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp2", kTRUE)' \
//           'plots/all_report/wp2_win.C+O'
// ============================================================================
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"

#include "defs.h"

namespace wp2w {
const char* const kEmpty = "/hdd3/przygoda/HADES/FAT_out/all_report/wp2/empty/scan2/scan_compact_empty.root";
const char* const kUpFile = "/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_upstream.root";
const char* const kSimVz = "/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_simvz.root";
const int kDays[3] = {34, 47, 65};
// luminosity of the period each empty-target day represents (analysed runs, wp2_lumi_perperiod.txt)
const double kLper[3] = {2809.453, 372.610 + 1078.022, 1819.110};
const char* kPerName[3] = {"I", "IIA+IIB", "III"};
inline int DI(int d) { return d == 34 ? 0 : d == 47 ? 1 : d == 65 ? 2 : -1; }
enum { sS0 = 0, sPre, sNone, sS0b1, sVz150, sVz170, sDz40, sVtx, kNS };   // sVtx = S0 without the vz > -160 cut, spectra only
const char* kSN[kNS] = {"S0", "pre", "none", "S0b1", "vz150", "vz170", "dz40", "S0 without vz cut"};
enum { wAll = 0, wB, wP, wCRU, wLo, wHi, wB0, wSB, wBP, wFit, kNW };
const char* kWN[kNW] = {"all M4", "B", "P", "CR-U", "M4 < 0.455", "M4 >= 0.70", "0.455-0.565", "P side bands", "B+P", "0.40-0.70"};
enum { uCh = 0, uChW, uChFar, uEeNear, uEeFar, kNU };
const char* kUN[kNU] = {"ch[-200,-170)", "ch[-250,-170)", "ch[-500,-250)", "ee[-200,-170)", "ee[-700,-250)"};
const char* kUT[kNU] = {"ch", "chw", "chfar", "eenear", "eefar"};
struct Acc {
    double n[2][3][3][kNS][kNW] = {}, v[2][3][3][kNS][kNW] = {};   // [target][sample][day][sel][win]
    double u[2][3][3][3] = {}, uv[2][3][3][3] = {};               // channel normalisers [target][sample][day][ch, chw, chfar]
};
struct Q { double y = 0, v = 0; };
}  // namespace wp2w

void wp2_win() {
    using namespace wp2w;
    Acc* A = new Acc();
    TFile* fo = TFile::Open("/hdd3/przygoda/HADES/FAT_out/all_report/wp2/wp2_win.root", "RECREATE");
    TH1D* hvz[2][3][3];   // [target][sample][day]: base + pion vertex, all vz
    TH1D* hm4[2][3][3];   // [target][sample][day]: S0, M4
    TH1D* hvzB[2][3][3];  // [target][sample][day]: base selection, all vz (normaliser region)
    for (int t = 0; t < 2; ++t)
        for (int s = 0; s < 3; ++s)
            for (int d = 0; d < 3; ++d) {
                hvz[t][s][d] = new TH1D(Form("vz_%s_%s_%03d", t ? "empty" : "full", scan::SampleName(s), kDays[d]),
                                        ";v_{z} [mm];rows / 2 mm", 300, -500, 100);
                hvzB[t][s][d] = new TH1D(Form("vzbase_%s_%s_%03d", t ? "empty" : "full", scan::SampleName(s), kDays[d]),
                                         ";v_{z} [mm];rows / 2 mm", 300, -500, 100);
                hm4[t][s][d] = new TH1D(Form("m4_%s_%s_%03d", t ? "empty" : "full", scan::SampleName(s), kDays[d]),
                                        ";M(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];rows / 10 mm", 80, 0.30, 1.10);
                hvz[t][s][d]->Sumw2(); hvzB[t][s][d]->Sumw2(); hm4[t][s][d]->Sumw2();
            }
    for (int t = 0; t < 2; ++t) {
        scan::Data D;
        if (t == 0) { if (!scan::Load(D)) return; }
        else { if (!scan::Load(D, kEmpty)) return; }
        for (int s = 0; s < 3; ++s) {
            const std::vector<scan::Row>& v = D.s[s];
            int cur = -1, curDay = -1;
            double g[kNS][kNW] = {}, gu[3] = {};
            auto flush = [&]() {
                if (curDay >= 0) {
                    for (int a = 0; a < kNS; ++a) for (int w = 0; w < kNW; ++w) A->v[t][s][curDay][a][w] += g[a][w] * g[a][w];
                    for (int q = 0; q < 3; ++q) A->uv[t][s][curDay][q] += gu[q] * gu[q];
                }
                for (int a = 0; a < kNS; ++a) for (int w = 0; w < kNW; ++w) g[a][w] = 0;
                for (int q = 0; q < 3; ++q) gu[q] = 0;
            };
            for (size_t i = 0; i < v.size(); ++i) {
                const scan::Row& r = v[i];
                if (r.grp != cur) { flush(); cur = r.grp; curDay = DI(r.file_id); }
                const int d = DI(r.file_id);
                if (d < 0) continue;
                const bool corOk = r.cor_ok == 1;
                if (!alldefs::PassBase(r, corOk)) continue;
                // channel normalisers
                const bool inU[3] = {r.vz >= -200.f && r.vz < -170.f, r.vz >= -250.f && r.vz < -170.f, r.vz >= -500.f && r.vz < -250.f};
                for (int q = 0; q < 3; ++q) if (inU[q]) { A->u[t][s][d][q] += 1; gu[q] += 1; }
                hvzB[t][s][d]->Fill(r.vz);
                bool pass[kNS];
                pass[sS0] = alldefs::PassSel(alldefs::kSelS0, r, corOk, r.dz_pi, r.mm4);
                pass[sPre] = alldefs::PassSel(alldefs::kSelPre, r, corOk, r.dz_pi, r.mm4);
                pass[sNone] = true;
                pass[sS0b1] = pass[sS0] && r.isBest == 1;
                pass[sVz150] = alldefs::PassSel(alldefs::kSelVz150, r, corOk, r.dz_pi, r.mm4);
                pass[sVz170] = alldefs::PassSel(alldefs::kSelVz170, r, corOk, r.dz_pi, r.mm4);
                pass[sDz40] = pass[sS0] && r.dz_pi < 40.f;   // S0 + DIRECT (dz_pi < 40 mm)
                pass[sVtx] = alldefs::PiVtx(r);
                if (pass[sVtx]) hvz[t][s][d]->Fill(r.vz);
                if (pass[sS0]) hm4[t][s][d]->Fill(r.m4);
                int w[kNW], nw = 0;
                w[nw++] = wAll;
                if (r.m4 >= alldefs::kBlo && r.m4 < alldefs::kBhi) w[nw++] = wB;
                if (r.m4 >= alldefs::kPlo && r.m4 < alldefs::kPhi) w[nw++] = wP;
                if (r.m4 >= 0.565f && r.m4 < 0.700f) w[nw++] = wCRU;
                if (r.m4 < 0.455f) w[nw++] = wLo;
                if (r.m4 >= 0.700f) w[nw++] = wHi;
                if (r.m4 >= 0.455f && r.m4 < 0.565f) w[nw++] = wB0;
                if ((r.m4 >= 0.5250f && r.m4 < 0.53364f) || (r.m4 >= 0.55416f && r.m4 < 0.5628f)) w[nw++] = wSB;
                if ((r.m4 >= alldefs::kBlo && r.m4 < alldefs::kBhi) || (r.m4 >= alldefs::kPlo && r.m4 < alldefs::kPhi)) w[nw++] = wBP;
                if (r.m4 >= 0.400f && r.m4 < 0.700f) w[nw++] = wFit;
                for (int a = 0; a < kNS; ++a) {
                    if (!pass[a] || a == sVtx) continue;
                    for (int k = 0; k < nw; ++k) { A->n[t][s][d][a][w[k]] += 1; g[a][w[k]] += 1; }
                }
            }
            flush();
        }
    }
    // normalisers U[target][day][kNU] and their variances
    double U[2][3][kNU] = {}, Uv[2][3][kNU] = {};
    for (int t = 0; t < 2; ++t)
        for (int d = 0; d < 3; ++d)
            for (int q = 0; q < 3; ++q)
                for (int s = 0; s < 3; ++s) { U[t][d][q] += A->u[t][s][d][q]; Uv[t][d][q] += A->uv[t][s][d][q]; }
    {
        TFile* fu = TFile::Open(kUpFile, "READ");
        for (int t = 0; t < 2; ++t)
            for (int d = 0; d < 3; ++d) {
                TH1D* h = fu ? (TH1D*)fu->Get(Form("h_EpEm_ID_%s_%03d", t ? "empty" : "full", kDays[d])) : nullptr;
                if (!h) { printf("wp2_win: EpEm_ID normaliser missing (t %d, day %d)\n", t, kDays[d]); continue; }
                const double lo[2] = {-200, -700}, hi[2] = {-170, -250};
                for (int q = 0; q < 2; ++q) {
                    const int b1 = h->FindBin(lo[q] + 1e-3), b2 = h->FindBin(hi[q] - 1e-3);
                    U[t][d][uEeNear + q] = h->Integral(b1, b2);
                    Uv[t][d][uEeNear + q] = U[t][d][uEeNear + q];   // EpEm_ID isBest == 1: one row per event
                }
            }
        if (fu) fu->Close();
    }
    // LH2 leak into the nominal normaliser window, from the simulation
    double simLeak = -1, simIn = -1;
    {
        TFile* fs = TFile::Open(kSimVz, "READ");
        TH1D* h = fs ? (TH1D*)fs->Get("simvz_T12") : nullptr;
        if (h) {
            simLeak = h->Integral(h->FindBin(-199.999), h->FindBin(-170.001));
            simIn = h->Integral(h->FindBin(-159.999), h->GetNbinsX());
        }
        if (fs) fs->Close();
    }
    // quantity per target/day/sel/win: 0 N, 1 LS, 2 N-CB, 3 peak signal (only for wP)
    auto get = [&](int t, int d, int a, int w, int qn) -> Q {
        Q q;
        auto one = [&](int ww, double& y, double& vv) {
            const double nu = A->n[t][0][d][a][ww], np = A->n[t][1][d][a][ww], nm = A->n[t][2][d][a][ww];
            const double vu = A->v[t][0][d][a][ww], vp = A->v[t][1][d][a][ww], vm = A->v[t][2][d][a][ww];
            if (qn == 0) { y = nu; vv = vu; return; }
            if (qn == 1) { y = np + nm; vv = vp + vm; return; }
            const double cb = (np > 0 && nm > 0) ? 2 * std::sqrt(np * nm) : 0;
            y = nu - cb;
            vv = vu + ((np > 0 && nm > 0) ? (nm / np) * vp + (np / nm) * vm : 0);
        };
        if (qn <= 2) { one(w, q.y, q.v); return q; }
        double yp, vp, ys, vs;
        one(wP, yp, vp);
        one(wSB, ys, vs);
        const double k = (0.55416 - 0.53364) / ((0.53364 - 0.5250) + (0.5628 - 0.55416));
        q.y = yp - k * ys;
        q.v = vp + k * k * vs;
        return q;
    };
    const char* kQN[4] = {"N", "LS", "N-CB", "S_P"};
    FILE* o = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp2_win.txt", "w");
    FILE* sm = fopen("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp2_win_summary.txt", "w");
    fprintf(o, "# wp2_win.C: target-window (empty-target) share f_win = b Y_E / Y_F, b = U_F / U_E (same day); isBest >= 0, PT3, cor_ok\n");
    fprintf(o, "# empty-target runs: day 034 (70 DSTs), 047 (160), 065 (1040), 3195 A; full = all analysed runs of the same day\n");
    fprintf(o, "# normalisers U (rows with a vertex in front of the cell): ");
    for (int q = 0; q < kNU; ++q) fprintf(o, "%d %s%s", q, kUN[q], q + 1 < kNU ? "; " : "\n");
    fprintf(o, "#   ch = channel rows us+pp+mm (base selection), ee = EpEm_ID PT3 isBest==1 rows; NOMINAL = ch[-200,-170)\n");
    fprintf(o, "# errors: event-group variance (Y, U_ch) and Poisson (U_ee), full and empty samples\n");
    fprintf(o, "# combination over days: f = sum_d (L_p(d)/L) f_d, L_p = %.3f (I, day 034), %.3f (IIA+IIB, day 047), %.3f (III, day 065)\n",
            kLper[0], kLper[1], kLper[2]);
    fprintf(o, "# simulated LH2 (T12, S0 without the vz cut): share in [-200,-170) / share in vz > -160 = %.2e\n\n",
            simIn > 0 ? simLeak / simIn : -1.);
    fprintf(o, "[normalisers] day: U_full U_empty  b = U_F/U_E +- stat   for each normaliser\n");
    double b[3][kNU], eb[3][kNU];
    for (int d = 0; d < 3; ++d) {
        fprintf(o, "  %03d", kDays[d]);
        for (int q = 0; q < kNU; ++q) {
            b[d][q] = U[1][d][q] > 0 ? U[0][d][q] / U[1][d][q] : 0;
            eb[d][q] = b[d][q] * std::sqrt((U[0][d][q] > 0 ? Uv[0][d][q] / (U[0][d][q] * U[0][d][q]) : 0) +
                                           (U[1][d][q] > 0 ? Uv[1][d][q] / (U[1][d][q] * U[1][d][q]) : 0));
            fprintf(o, " | %s %.0f %.0f %.3f +- %.3f", kUT[q], U[0][d][q], U[1][d][q], b[d][q], eb[d][q]);
        }
        fprintf(o, "\n");
    }
    fprintf(o, "  ratios b(q)/b(ch): ");
    for (int d = 0; d < 3; ++d) {
        fprintf(o, "%03d:", kDays[d]);
        for (int q = 1; q < kNU; ++q) fprintf(o, " %s %.3f", kUT[q], b[d][q] / b[d][0]);
        fprintf(o, d < 2 ? "; " : "\n");
    }
    fprintf(o, "\n[f_win] sel | window | quantity | day: Y_full Y_empty | f [%%] +- stat for the normalisers ch chw chfar eenear eefar\n");
    fprintf(sm, "# wp2_win.C summary: f_win [%%] combined over the three days (L-weighted), nominal normaliser ch[-200,-170)\n");
    fprintf(sm, "# key = sel|window|quantity ; f_ch, stat ; f_eenear ; f_eefar ; f_chw ; f_chfar ; f per day (ch) 034 047 065 ; Y_full Y_empty (sum)\n");
    const double Lsum = kLper[0] + kLper[1] + kLper[2];
    for (int a = 0; a < sVtx; ++a)
        for (int w = 0; w < kNW; ++w)
            for (int qn = 0; qn < 4; ++qn) {
                if (qn == 3 && w != wP) continue;
                double fc[kNU] = {}, efc2[kNU] = {}, fd0[3] = {}, YFs = 0, YEs = 0;
                for (int d = 0; d < 3; ++d) {
                    const Q F = get(0, d, a, w, qn), E = get(1, d, a, w, qn);
                    YFs += F.y; YEs += E.y;
                    fprintf(o, "  %-5s %-13s %-5s %03d: %9.0f %8.1f |", kSN[a], kWN[w], kQN[qn], kDays[d], F.y, E.y);
                    for (int q = 0; q < kNU; ++q) {
                        const double f = (F.y != 0) ? b[d][q] * E.y / F.y : 0;
                        const double ef = (F.y != 0) ? std::sqrt(std::pow(b[d][q] * std::sqrt(E.v) / F.y, 2) +
                                                                 std::pow(f * std::sqrt(F.v) / F.y, 2) +
                                                                 std::pow(f * eb[d][q] / (b[d][q] > 0 ? b[d][q] : 1), 2))
                                                     : 0;
                        fprintf(o, " %7.3f +- %6.3f", 100 * f, 100 * ef);
                        fc[q] += kLper[d] / Lsum * f;
                        efc2[q] += std::pow(kLper[d] / Lsum * ef, 2);
                        if (q == 0) fd0[d] = f;
                    }
                    fprintf(o, "\n");
                }
                fprintf(o, "  %-5s %-13s %-5s all: %9.0f %8.1f |", kSN[a], kWN[w], kQN[qn], YFs, YEs);
                for (int q = 0; q < kNU; ++q) fprintf(o, " %7.3f +- %6.3f", 100 * fc[q], 100 * std::sqrt(efc2[q]));
                fprintf(o, "\n");
                fprintf(sm, "%s|%s|%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.0f,%.1f\n", kSN[a], kWN[w], kQN[qn], 100 * fc[0],
                        100 * std::sqrt(efc2[0]), 100 * fc[uEeNear], 100 * fc[uEeFar], 100 * fc[uChW], 100 * fc[uChFar], 100 * fd0[0],
                        100 * fd0[1], 100 * fd0[2], YFs, YEs);
            }
    fclose(o);
    fclose(sm);
    fo->cd();
    for (int t = 0; t < 2; ++t) for (int s = 0; s < 3; ++s) for (int d = 0; d < 3; ++d) { hvz[t][s][d]->Write(); hvzB[t][s][d]->Write(); hm4[t][s][d]->Write(); }
    fo->Close();
    printf("wp2_win: done; tables docs/latex/tables_all/wp2_win.txt, wp2_win_summary.txt\n");
}

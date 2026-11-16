// ============================================================================
// vc2recon_bedge.C  --  independent verification of package c2recon: the per-slice edge estimator (bedge) with
//                       the verifier's own QED-weighted C2 templates
// ============================================================================
// Written 2026-11-15 (prefix vc2recon, verifier). New file. bedge_ana.C is #included unchanged (its fit code and the
// nominal configuration); bedge_ana() is never called and nothing is written into the bedge directory. The C2
// templates are built HERE from the verifier's row list vc2_bedge_rows.txt (vc2recon_fit.py: bedge bins of the nominal
// configuration and the weight multipliers of every variant), never from c2recon files.
//   1. nominal reproduction (E_data, E_sim, linear C2 / C3 parts per class);
//   2. the own isotropic C2 template against sim_S.bin [nom][C2] bin by bin;
//   3. for every weight variant v: SS[nom][C2] -> SS - T_iso + T_v, full re-evaluation: linear C2 step per class with the
//      data covariance (-> effective factor), E_sim with the WP6 own weights (-> f directly), the bedge formula;
//   4. per-bin linear coefficients of the class fits (unit responses with the data covariance) for the row bootstrap.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/vc2_bedge.txt
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon/aclic", kTRUE)' \
//        'plots/all_report/vc2recon_bedge.C+O'
// ============================================================================
#include <fstream>
#include <sstream>
#include "bedge_ana.C"

namespace vc2b {
const char* const kOutV = "/hdd3/przygoda/HADES/FAT_out/all_report/vc2recon";
FILE* gO = nullptr;
#define VO(...) do { printf(__VA_ARGS__); if (gO) fprintf(gO, __VA_ARGS__); fflush(stdout); } while (0)
}

void vc2recon_bedge() {
    using namespace bedgea;
    using namespace vc2b;
    gO = fopen(Form("%s/vc2_bedge.txt", kOutV), "w");
    if (!ReadArr(Form("%s/data_N.bin", kOut), DN, SizeDN()) || !ReadArr(Form("%s/sim_S.bin", kOut), SS, SizeSS()) ||
        !ReadArr(Form("%s/sim_hq.bin", kOut), HQ, SizeHq())) return;
    InitSets();
    for (int k = 0; k < 3; ++k) CF[2][k] = 0.5 * (CF[0][k] + CF[1][k]);
    Var nom; nom.tag = "nominal";
    Res R0;
    Evaluate(nom, R0);
    double q0[kNQ]; Quant(MakeInp(R0), q0);
    for (int k = 0; k < 3; ++k)
        VO("NOM k %d Ed %.6f Es %.6f linT12 %.6f linC2 %.6f linC3 %.6f\n", k, R0.c[k].Ed, R0.c[k].Es, R0.c[k].lin[0], R0.c[k].lin[1], R0.c[k].lin[2]);
    VO("NOMQ f_all %.6f fc_rec %.6f (CF E1 %.4f E2 %.4f)\n", q0[0], q0[5], CF[2][1], CF[2][2]);

    // ---- the verifier's rows: e p ib W m_1 .. m_nv
    std::ifstream in(Form("%s/vc2_bedge_rows.txt", kOutV));
    std::string line; std::vector<std::string> vn;
    struct Rw { int e, p, ib; double W; std::vector<double> m; };
    std::vector<Rw> rows;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        if (line[0] == '#') { std::string t; ss >> t; while (ss >> t) vn.push_back(t); continue; }
        Rw r; ss >> r.e >> r.p >> r.ib >> r.W;
        r.m.resize(vn.size());
        for (size_t j = 0; j < vn.size(); ++j) ss >> r.m[j];
        rows.push_back(r);
    }
    const int nv = (int)vn.size();
    VO("ROWS %zu variants %d\n", rows.size(), nv);
    const size_t nb = (size_t)kNCls * kNPP * kNM4;
    auto Ix = [](int e, int p, int i) { return ((size_t)e * kNPP + p) * kNM4 + i; };
    std::vector<std::vector<double>> T(nv + 1, std::vector<double>(nb, 0.));
    for (const Rw& r : rows) {
        if (r.e < 0 || r.p < 0 || r.ib < 0) continue;
        T[0][Ix(r.e, r.p, r.ib)] += r.W;
        for (int j = 0; j < nv; ++j) T[j + 1][Ix(r.e, r.p, r.ib)] += r.W * r.m[j];
    }
    // ---- check against sim_S.bin
    double maxd = 0, maxdFit = 0, sumS = 0, sumT = 0;
    for (int e = 0; e < kNCls; ++e)
        for (int p = 0; p < kNPP; ++p)
            for (int i = 0; i < kNM4; ++i) {
                const double s = SS[IdxSS(sNom, gC2, e, p, i)], t = T[0][Ix(e, p, i)];
                maxd = std::max(maxd, std::fabs(s - t));
                if (p >= 6 && p < 26 && i >= 20 && i < 108) { maxdFit = std::max(maxdFit, std::fabs(s - t)); sumS += s; sumT += t; }
            }
    VO("CHECK own iso C2 template vs sim_S.bin[nom][C2]: max |diff| all bins %.3g, in M(pipi) 0.30-0.40 x M4 0.40-0.62: %.3g (sums %.6f / %.6f)\n",
       maxd, maxdFit, sumS, sumT);

    // ---- variants
    std::vector<double> orig(nb);
    for (int e = 0; e < kNCls; ++e) for (int p = 0; p < kNPP; ++p) for (int i = 0; i < kNM4; ++i) orig[Ix(e, p, i)] = SS[IdxSS(sNom, gC2, e, p, i)];
    for (int v = 0; v <= nv; ++v) {
        for (int e = 0; e < kNCls; ++e) for (int p = 0; p < kNPP; ++p) for (int i = 0; i < kNM4; ++i)
            SS[IdxSS(sNom, gC2, e, p, i)] = orig[Ix(e, p, i)] - T[0][Ix(e, p, i)] + T[v][Ix(e, p, i)];
        Res Rv; Evaluate(nom, Rv);
        const double l1 = Rv.c[1].lin[1], l2 = Rv.c[2].lin[1];
        const double b1 = R0.c[1].lin[1], b2 = R0.c[2].lin[1];
        const double red = ((b1 - l1) + (b2 - l2)) / (R0.c[1].lin[0] + R0.c[2].lin[0]);
        const double fForm = R0.c[0].Ed / (R0.c[0].Es * (1 - red));
        const double fDir = Rv.c[0].Ed / Rv.c[0].Es;
        VO("VAR %s c_E1 %.6f c_E2 %.6f c_cls %.6f c_allfit %.6f linC2 %.6f %.6f %.6f Es %.6f %.6f %.6f f_formula %.6f f_direct %.6f\n",
           v == 0 ? "iso" : vn[v - 1].c_str(), l1 / b1, l2 / b2, (l1 + l2) / (b1 + b2), Rv.c[0].lin[1] / R0.c[0].lin[1], Rv.c[0].lin[1], l1, l2,
           Rv.c[0].Es, Rv.c[1].Es, Rv.c[2].Es, fForm, fDir);
    }
    for (int e = 0; e < kNCls; ++e) for (int p = 0; p < kNPP; ++p) for (int i = 0; i < kNM4; ++i) SS[IdxSS(sNom, gC2, e, p, i)] = orig[Ix(e, p, i)];
    // the applied class factors, and the direct own-weight E_sim with them
    {
        for (int e = 0; e < kNCls; ++e) for (int p = 0; p < kNPP; ++p) for (int i = 0; i < kNM4; ++i)
            SS[IdxSS(sNom, gC2, e, p, i)] = orig[Ix(e, p, i)] * (e == 0 ? CF[2][1] : CF[2][2]);
        Res Ra; Evaluate(nom, Ra);
        VO("APPLIED c_cls %.6f c_allfit %.6f f_direct %.6f\n", (Ra.c[1].lin[1] + Ra.c[2].lin[1]) / (R0.c[1].lin[1] + R0.c[2].lin[1]),
           Ra.c[0].lin[1] / R0.c[0].lin[1], Ra.c[0].Ed / Ra.c[0].Es);
        for (int e = 0; e < kNCls; ++e) for (int p = 0; p < kNPP; ++p) for (int i = 0; i < kNM4; ++i) SS[IdxSS(sNom, gC2, e, p, i)] = orig[Ix(e, p, i)];
    }
    // ---- unit responses of the class fits (data covariance) for every (class fit k, slice, M4 bin)
    const std::vector<ESlice>& set = gSets[0];
    for (int k = 0; k < 3; ++k)
        for (size_t s = 0; s < set.size(); ++s) {
            std::vector<double> y(kNM4, 0.);
            for (int i = 0; i < kNM4; ++i) {
                std::fill(y.begin(), y.end(), 0.); y[i] = 1.0;
                const scan::StepFit fl = FitStep(SimSpec(y.data(), kNM4, scan::kMlo, kBwD), R0.c[k].cS[s], nullptr, &R0.c[k].Vd[s]);
                if (fl.ok && fl.h != 0) fprintf(gO, "COEF %d %zu %d %d %d %.10g\n", k, s, PB(set[s].lo), PB(set[s].hi), i, fl.h);
            }
        }
    VO("done\n");
    fclose(gO);
}

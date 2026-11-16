// ============================================================================
// wp8_syst.C  --  WP8: refit of the WP7 chain with the efficiency-model variants of S12
// ============================================================================
// Written by WP8 (prefix wp8), 2026-11-14. Reuses wp7_syst.C (Run, Cfg, templates, estimators; included, not changed,
// not executed) so that every quantity of Table S1 is recomputed exactly as in WP7 for:
//   ang_toy  : class C2 (bump pi0 Dalitz, HGeant isotropic gamma* decay) x the QED decay-angle factor of the toy
//              (wp8_toy_angle.txt, p_e > 0.06 GeV/c): E1 0.8886, E2 0.8675, F5g (>= 20 MeV) 0.8864, all-M(ee) windows 0.8800
//   ang_clos : the same with the factors of the full-simulation lep/std closure of pi0 Dalitz pairs (PLUTO with the QED
//              angle vs HGeant isotropic; closure_ratio_bins.txt of the simulation study, divided by the PLUTO/HGeant
//              generator ratio of the q spectrum): E1 0.869, E2 0.896, F5g 0.909, all 0.891
//   pty3x3, ptyPT, ptyY : (pT4, y4) reweighting of the eta production (wp8_pty.txt): bump classes C2, C3 x k_B, peak classes
//              C4, C5 and G_eta x k_P (3 x 3 cells as specified; pT-only and y-only as information)
// The simulation templates are the WP0 S0__nom F10/F5g templates scaled per component and M(ee) class; nothing else changes.
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp8/{wp8_syst_raw.txt, wp8_syst.log}
// Run (from /home/przygoda/HADES/FAT/FAT; about 1 min, < 3 GB):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp8", kTRUE)' 'plots/all_report/wp8_syst.C+O'
// ============================================================================
#include "wp7_syst.C"

namespace wp8y {
using namespace wp7;
// scale component t of a template: main histograms by class (h even = E1, odd = E2), F5g by f5, aligned windows per cell,
// non-aligned all-cell windows (P, side bands, tot) by fAll
void ScaleComp(Tmpl& T, int t, double fE1, double fE2, double f5, double fAll) {
    const int nb = T.nb;
    for (int h = 0; h < 10; ++h) {
        const double f = (h % 2 == 0) ? fE1 : fE2;
        for (int x = 0; x < nb; ++x) { T.s[t][(size_t)h * nb + x] *= f; T.vs[t][(size_t)h * nb + x] *= f * f; }
        for (int w = 0; w < kNW; ++w) T.ws[t][(size_t)h * kNW + w] *= f;
    }
    for (int j = 0; j < 5; ++j)
        for (int x = 0; x < nb; ++x) { T.s5[t][(size_t)j * nb + x] *= f5; T.vs5[t][(size_t)j * nb + x] *= f5 * f5; }
    for (int w = 0; w < kNW; ++w) {
        if (kWins[w].aligned && T.ws[t][w] == T.ws[t][w]) {
            double s = 0; for (int h = 0; h < 10; ++h) s += T.ws[t][(size_t)h * kNW + w];
            T.wa[t][w] = s;
        } else T.wa[t][w] *= fAll;
    }
}
std::string KeyOf(const Cfg& c) {
    return Form("%d|%s|%s|%s|%.6f|%d", c.src, c.simSel.c_str(), c.simVar.c_str(), c.L.tag.c_str(), c.simScale, (int)c.fact);
}
bool ReadPty(const char* path, std::map<std::string, std::pair<double, double>>& kk) {   // name -> (kB, kP)
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream is(line);
        std::string a; is >> a;
        if (a == "pty") {
            std::string k; double v; std::map<std::string, double> m;
            while (is >> k >> v) m[k] = v;
            kk["pty3x3"] = {m["kB"], m["kP"]};
        } else if (a == "pty1d") {
            std::string nm; is >> nm;
            std::string k; double v; std::map<std::string, double> m;
            while (is >> k >> v) m[k] = v;
            const double kP = 1. / (1. + m["rel"]), kB = 1. / (1. + m["aB_rel"]);
            kk[std::string("pty") + nm] = {kB, kP};
        }
    }
    return true;
}
}  // namespace wp8y

void wp8_syst() {
    using namespace alldefs;
    using namespace wp7;
    using namespace wp8y;
    TStopwatch total;
    const char* out = "/hdd3/przygoda/HADES/FAT_out/all_report/wp8";
    gSystem->mkdir(out, kTRUE);
    gLog = fopen(Form("%s/wp8_syst.log", out), "w");
    gRaw = fopen(Form("%s/wp8_syst_raw.txt", out), "w");
    gRepF = nullptr;
    fprintf(gRaw, "# wp8_syst.C (WP8, 2026-11-14): WP7 chain (wp7_syst.C Run) with the S12 efficiency-model variants; S0, isBest >= 0\n");
    fprintf(gRaw, "# cfg <tag> <family> chi2 | text ; q <tag> <quantity> value stat N CB A F mc nrep\n");
    if (!LoadRows(gRows)) { P("rows missing\n"); return; }
    gFH = TFile::Open(kHistFile);
    Cfg nom; nom.tag = "nom"; nom.fam = "nom"; nom.mc = false; nom.text = "nominal (= WP7 nominal)";
    Tmpl& T0 = GetTmpl(nom);
    Res rNom = Run(nom); PrintRes(rNom); WriteRes(rNom);
    P("CHECK vs WP7 nominal: a_eta %.8f (0.38228397) BR_G %.6e (2.6546e-4) Y_bump %.3f (16389.273) Y_peak %.3f (5087.9747) DR %.6f (1.3974303) sigma_eta %.6f (0.43744039)\n",
      rNom.q["a_eta"].v, rNom.q["BR_G"].v, rNom.q["Y_bump"].v, rNom.q["Y_peak"].v, rNom.q["DR"].v, rNom.q["sigma_eta"].v);

    struct V { std::string tag, text; double f[5][4]; };   // per signal component (C2 C3 C4 C5 G_eta): E1, E2, F5g, all
    std::vector<V> vars;
    {
        V v; v.tag = "ang_toy"; v.text = "C2 x QED decay-angle factor (toy): E1 0.8886, E2 0.8675, F5g 0.8864, all 0.8800";
        for (int t = 0; t < 5; ++t) for (int k = 0; k < 4; ++k) v.f[t][k] = 1.;
        v.f[0][0] = 0.8886; v.f[0][1] = 0.8675; v.f[0][2] = 0.8864; v.f[0][3] = 0.8800;
        vars.push_back(v);
    }
    {
        V v; v.tag = "ang_clos"; v.text = "C2 x QED decay-angle factor (lep/std closure, q-spectrum corrected): E1 0.869, E2 0.896, F5g 0.909, all 0.891";
        for (int t = 0; t < 5; ++t) for (int k = 0; k < 4; ++k) v.f[t][k] = 1.;
        v.f[0][0] = 0.869; v.f[0][1] = 0.896; v.f[0][2] = 0.909; v.f[0][3] = 0.891;
        vars.push_back(v);
    }
    std::map<std::string, std::pair<double, double>> kk;
    if (!ReadPty("/home/przygoda/HADES/FAT/FAT/docs/latex/tables_all/wp8_pty.txt", kk)) P("!! wp8_pty.txt missing\n");
    for (const char* nm : {"pty3x3", "ptypT", "ptyy"}) {
        if (!kk.count(nm)) continue;
        V v; v.tag = nm;
        const double kB = kk[nm].first, kP = kk[nm].second;
        v.text = Form("(pT4, y4) reweighting %s: C2, C3 x k_B = %.4f; C4, C5, G_eta x k_P = %.4f", nm, kB, kP);
        for (int t = 0; t < 5; ++t) for (int k = 0; k < 4; ++k) v.f[t][k] = (t <= 1) ? kB : kP;
        vars.push_back(v);
    }
    double sc = 1.0;
    for (const V& v : vars) {
        sc += 1e-6;
        Cfg c = nom; c.tag = v.tag; c.fam = "S12"; c.text = v.text; c.simScale = sc; c.mc = false;
        Tmpl Tm = T0;
        for (int t = 0; t < 5; ++t) ScaleComp(Tm, t, v.f[t][0], v.f[t][1], v.f[t][2], v.f[t][3]);
        gTmpl[KeyOf(c)] = Tm;
        Res r = Run(c); PrintRes(r); WriteRes(r);
    }
    fclose(gRaw); gRaw = nullptr; fclose(gLog); gLog = nullptr;
    printf("wp8_syst done, %.0f s\n", total.RealTime());
}

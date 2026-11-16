// ============================================================================
// wp4_toys.C  --  WP4 stress closure: parametric toys and the estimators (prefix wp4)
// ============================================================================
// SPEC report_all, WP4. For every scenario (a)-(h) (+ informational (i)) and every
// toy: truth = half B of the simulation (non-eta background smoothed exp(pol3) per
// (slice, M(ee) class), eta classes raw) scaled to L_data, plus a CB-like component
// (data CB, spread rule); toy N = Poisson(truth) on the refined cells; toy like-sign
// counts = Poisson(data N++, N--) per (2.5 MeV bin, slice, M(ee) class); the toy CB is
// rebuilt from them with the WP0 rule. The estimators F-main, F-free, F-gold, B4 and
// the peak side bands (wp4_lib.h) are run with fit templates from half A.
//   cfg 0 "BA": truth half B, templates half A (SPEC nominal)
//   cfg 1 "AA": truth half A, templates half A (no template/truth MC mismatch: method bias)
//   cfg 2 "AB": truth half A, templates half B (swap: the MC mismatch changes sign)
//   cfg 3 "BB": truth half B, templates half B
//   cfg 4 "RF": truth full sample (A + B), templates = bootstrap replica r = toy mod 100 of the full sample
//               (hists_v1 _rep, WP0); sigma_MC per quantity = RMS of the 100 replica fits of the Asimov data;
//               toy -1 = Asimov with the full-sample templates (method at full MC statistics)
//   cfg 5 "RB": as RF, but with the double bootstrap of wp4_fill_boot.C: templates of toy t = parent replica
//               p = t mod P; the C children of p give the per-toy sigma_MC (RMS) and the bootstrap bias correction
//               est = 2 est_p - mean_c est_pc (the procedure the data analysis can apply with its own replicas)
// toy index -1 = Asimov (expected counts, like-sign = data).
// Output: /hdd3/przygoda/HADES/FAT_out/all_report/wp4/toys/toys_<cfg>.root (trees toys, truth)
//
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_wp4", kTRUE)' \
//        'plots/all_report/wp4_toys.C+O(0, 500)'
// ============================================================================
#include <cstdio>
#include <string>

#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"

#include "wp4_lib.h"

void wp4_toys(int cfg = 0, int ntoy = 500, int scenMask = -1, const char* tag = "") {
    using namespace wp4;
    TStopwatch sw;
    const char* cfgTag[6] = {"BA", "AA", "AB", "BB", "RF", "RB"};
    const int truthHalf[6] = {1, 0, 0, 1, -1, -1}, tmplHalf[6] = {0, 0, 1, 1, -1, -1};
    Grid G; G.Init();
    gSystem->mkdir(Form("%s/toys", kOutDir), kTRUE);
    FILE* log = fopen(Form("%s/toys/toys_%s%s.log", kOutDir, cfgTag[cfg], tag), "w");
    TFile* fh = TFile::Open(kHistFile);
    DataIn Din;
    if (!LoadData(fh, G, Din, log)) return;
    // ---- check the CB rebuild against hists_v1 (data S0)
    {
        Obs O; O.n.assign(kNh15 * (G.nR + 2), 0.); O.ls = Din.ls; Derive(G, O);
        double worst = 0;
        for (int k = 0; k < 10; ++k) {
            TH1* h = allio::Data(fh, "S0", "all", "CB", allio::F10(k / 2, k % 2));
            for (int x = 1; x <= G.nF; ++x) worst = std::max(worst, std::fabs(O.cb10[k][x] - h->GetBinContent(x)) / std::max(1e-9, h->GetBinContent(x)));
            delete h;
        }
        for (int j = 0; j < kNj; ++j) {
            TH1* h = allio::Data(fh, "S0", "all", "CB", allio::F5g(j));
            for (int x = 1; x <= G.nF; ++x) if (h->GetBinContent(x) > 0) worst = std::max(worst, std::fabs(O.cb5[j][x] - h->GetBinContent(x)) / h->GetBinContent(x));
            delete h;
        }
        TH1* w = allio::Data(fh, "S0", "all", "CB", "win");
        const double cw[4] = {O.CBB, O.CBP, O.CBSB, O.CBU};
        const double fw[4] = {allio::Win(w, kWB), allio::Win(w, kWP), allio::Win(w, kWSBlo) + allio::Win(w, kWSBhi), allio::Win(w, kWCRU)};
        for (int i = 0; i < 4; ++i) worst = std::max(worst, std::fabs(cw[i] - fw[i]) / fw[i]);
        fprintf(log, "CB rebuild vs hists_v1 (F10, F5g, windows B P SB CR-U): worst relative difference %.3g\n", worst);
        fprintf(log, "data CB windows: B %.2f P %.2f SB %.2f CRU %.2f; replica sigma B %.2f P %.2f SB %.2f CRU %.2f\n", cw[0], cw[1], cw[2], cw[3],
                std::sqrt(Din.covW[0][0]), std::sqrt(Din.covW[1][1]), std::sqrt(Din.covW[2][2] + Din.covW[3][3] + 2 * Din.covW[2][3]), std::sqrt(Din.covW[4][4]));
        printf("wp4_toys: CB rebuild worst relative difference %.3g\n", worst);
        delete w;
    }
    // ---- templates
    TFile* ft = TFile::Open(kTmplFile);
    Tmpl* Tfit = new Tmpl(); Tmpl* Ttr = new Tmpl(); Tmpl* Ttrf = new Tmpl();
    const bool rf = (cfg == 4), rb = (cfg == 5);
    if (!rf && !rb) {
        if (!LoadTmpl(ft, kMnom, tmplHalf[cfg], G, *Tfit) || !LoadTmpl(ft, kMnom, truthHalf[cfg], G, *Ttr) ||
            !LoadTmpl(ft, kMf, truthHalf[cfg], G, *Ttrf)) return;
    } else {
        if (!LoadTmplFull(ft, kMnom, G, *Tfit) || !LoadTmplFull(ft, kMnom, G, *Ttr) || !LoadTmplFull(ft, kMf, G, *Ttrf)) return;
    }
    FitTmpl* F = new FitTmpl();
    MakeFitTmpl(*Tfit, G, *F);
    const int nRep = 100;
    std::vector<FitTmpl*> Frep;
    std::vector<Tmpl*> Trep;
    if (rf && !LoadReplicas(fh, G, *Tfit, nRep, Frep, Trep)) return;
    int bP = 0, bC = 0;
    std::vector<FitTmpl*> Fb;
    std::vector<Tmpl*> Tb;
    if (rb) {
        if (!LoadBoot(Form("%s/wp4_boot.root", kOutDir), G, *Tfit, bP, bC, Fb, Tb)) { printf("wp4_toys: wp4_boot.root missing\n"); return; }
        fprintf(log, "double bootstrap: %d parents x %d children\n", bP, bC);
    }
    if (rf) {   // check: replica mean of the window sums against the full sample
        double m[kNSig] = {0};
        for (int r = 0; r < nRep; ++r) for (int t = 0; t < kNSig; ++t) m[t] += Frep[r]->W[t][0] / nRep;
        fprintf(log, "replicas: mean C2/C3/C4/C5/Geta in B over full sample: %.4f %.4f %.4f %.4f %.4f\n", m[0] / F->W[0][0], m[1] / F->W[1][0],
                m[2] / F->W[2][0], m[3] / F->W[3][0], m[4] / F->W[4][0]);
    }
    TruthSmooth* TS = new TruthSmooth(); TruthSmooth* TSf = new TruthSmooth();
    MakeTruthSmooth(*Ttr, G, *TS); MakeTruthSmooth(*Ttrf, G, *TSf);
    fprintf(log, "fit smoothing (half %s): F10 histogram: GD order/N_eff, GC order/N_eff\n", tmplHalf[cfg] == 0 ? "A" : tmplHalf[cfg] == 1 ? "B" : "A+B");
    for (int k = 0; k < 10; ++k)
        fprintf(log, "  S%dE%d GD %d/%.0f GC %d/%.0f\n", k / 2 + 1, k % 2 + 1, F->B.GD10[k].order, F->B.GD10[k].neff, F->B.GC10[k].order, F->B.GC10[k].neff);
    for (int j = 0; j < kNj; ++j)
        fprintf(log, "  F5g S%d GD %d%s/%.0f GC %d%s/%.0f\n", j + 1, F->B.GD5[j].order, F->B.GD5[j].scaled ? "s" : "", F->B.GD5[j].neff,
                F->B.GC5[j].order, F->B.GC5[j].scaled ? "s" : "", F->B.GC5[j].neff);
    fprintf(log, "truth smoothing (half %s): F10 histogram GD GC GDom GDnom GCom GCnom order/N_eff (s = parent shape scaled); F5g GD GC\n",
            truthHalf[cfg] == 0 ? "A" : truthHalf[cfg] == 1 ? "B" : "A+B");
    auto so = [](const Smooth& x) { return std::string(Form("%d%s/%.0f", x.order, x.scaled ? "s" : "", x.neff)); };
    for (int k = 0; k < 10; ++k)
        fprintf(log, "  S%dE%d  %s %s %s %s %s %s\n", k / 2 + 1, k % 2 + 1, so(TS->B.GD10[k]).c_str(), so(TS->B.GC10[k]).c_str(),
                so(TS->GDom[k]).c_str(), so(TS->GDnom[k]).c_str(), so(TS->GCom[k]).c_str(), so(TS->GCnom[k]).c_str());
    for (int j = 0; j < kNj; ++j) fprintf(log, "  F5g S%d  %s %s\n", j + 1, so(TS->B.GD5[j]).c_str(), so(TS->B.GC5[j]).c_str());
    // ---- output
    TFile* fo = TFile::Open(Form("%s/toys/toys_%s%s.root", kOutDir, cfgTag[cfg], tag), "RECREATE");
    TTree* tt = new TTree("toys", Form("WP4 toys, cfg %s", cfgTag[cfg]));
    int scen, toy, nitFM;
    double est[kNQ], err[kNQ], errd[kNQ], errs[kNQ];
    int ok[kNQ], rep;
    double chi2[4]; int ndf[4];
    double fmTh[13], saTh[9], corrAB, corrRK;
    tt->Branch("scen", &scen, "scen/I"); tt->Branch("toy", &toy, "toy/I"); tt->Branch("nitFM", &nitFM, "nitFM/I");
    tt->Branch("est", est, Form("est[%d]/D", kNQ)); tt->Branch("err", err, Form("err[%d]/D", kNQ));
    tt->Branch("errd", errd, Form("errd[%d]/D", kNQ)); tt->Branch("errs", errs, Form("errs[%d]/D", kNQ));
    tt->Branch("ok", ok, Form("ok[%d]/I", kNQ)); tt->Branch("rep", &rep, "rep/I");
    double mct[kNQ], est0[kNQ], cmean[kNQ];
    if (rb) {
        tt->Branch("mct", mct, Form("mct[%d]/D", kNQ));       // per-toy sigma_MC = RMS over the children
        tt->Branch("est0", est0, Form("est0[%d]/D", kNQ));    // uncorrected estimate (parent templates)
        tt->Branch("cmean", cmean, Form("cmean[%d]/D", kNQ)); // mean over the children
    }
    tt->Branch("chi2", chi2, "chi2[4]/D"); tt->Branch("ndf", ndf, "ndf[4]/I");
    tt->Branch("fmTh", fmTh, "fmTh[13]/D"); tt->Branch("saTh", saTh, "saTh[9]/D");
    tt->Branch("corrAB", &corrAB, "corrAB/D"); tt->Branch("corrRK", &corrRK, "corrRK/D");
    TTree* tu = new TTree("truth", "true values per scenario");
    double val[kNQ], cbTot, bkgTot, etaTot, mc[kNQ], mcMean[kNQ];
    tu->Branch("scen", &scen, "scen/I"); tu->Branch("val", val, Form("val[%d]/D", kNQ));
    tu->Branch("mc", mc, Form("mc[%d]/D", kNQ)); tu->Branch("mcMean", mcMean, Form("mcMean[%d]/D", kNQ));
    TTree* tr = new TTree("repasi", "Asimov fits with the bootstrap-replica templates (cfg RF)");
    tr->Branch("scen", &scen, "scen/I"); tr->Branch("rep", &rep, "rep/I"); tr->Branch("est", est, Form("est[%d]/D", kNQ));
    tr->Branch("ok", ok, Form("ok[%d]/I", kNQ));
    tu->Branch("cbTot", &cbTot, "cbTot/D"); tu->Branch("bkgTot", &bkgTot, "bkgTot/D"); tu->Branch("etaTot", &etaTot, "etaTot/D");

    for (int s = 0; s < kNScen; ++s) {
        if (!(scenMask & (1 << s))) continue;
        const Scen& S = kScen[s];
        Truth U;
        MakeTruth(S, S.mv == kMf ? *Ttrf : *Ttr, S.mv == kMf ? *TSf : *TS, Din, G, U);
        scen = s;
        for (int q = 0; q < kNQ; ++q) val[q] = U.val[q];
        cbTot = U.cbTot; bkgTot = U.bkgTot; etaTot = U.etaTot;
        for (int q = 0; q < kNQ; ++q) { mc[q] = 0; mcMean[q] = 0; }
        if (rf) {   // sigma_MC: the Asimov data fitted with each replica's templates
            Obs Oa;
            MakeToy(U, Din, G, nullptr, Oa);
            std::vector<double> sx(kNQ, 0.), sxx(kNQ, 0.);
            std::vector<int> nx(kNQ, 0);
            for (int r = 0; r < nRep; ++r) {
                Result R;
                RunEstimators(Oa, Din, *Frep[r], G, R);
                rep = r;
                for (int q = 0; q < kNQ; ++q) {
                    est[q] = R.est[q]; ok[q] = R.ok[q] ? 1 : 0;
                    if (R.ok[q]) { sx[q] += R.est[q]; sxx[q] += R.est[q] * R.est[q]; ++nx[q]; }
                }
                tr->Fill();
            }
            for (int q = 0; q < kNQ; ++q)
                if (nx[q] > 1) {
                    const double mm = sx[q] / nx[q];
                    mcMean[q] = mm - U.val[q];
                    mc[q] = std::sqrt(std::max(0., sxx[q] / nx[q] - mm * mm));
                }
            fprintf(log, "  sigma_MC (RMS of %d replica Asimov fits): a %.4f Ybump %.1f Ypeak %.1f rho %.4f BR_G %.3g B4 %.1f PSB %.1f; "
                         "mean shift a %+.4f Ybump %+.1f Ypeak %+.1f\n", nRep, mc[qFMa], mc[qFMYb], mc[qFMYp], mc[qFMrho], mc[qFGBRG], mc[qB4Yb],
                    mc[qPSBYp], mcMean[qFMa], mcMean[qFMYb], mcMean[qFMYp]);
        }
        tu->Fill();
        double dataTot = 0;
        for (int k = 0; k < 10; ++k) for (int x = 1; x <= G.nF; ++x) dataTot += Din.N10[k][x];
        fprintf(log, "scenario %s (%s): truth in 0.40-0.70: CB %.0f, background %.0f, eta %.0f, total %.0f (data S0 N %.0f); "
                     "true Y_bump %.1f, Y_peak %.1f\n", S.tag, S.text, U.cbTot, U.bkgTot, U.etaTot, U.cbTot + U.bkgTot + U.etaTot, dataTot,
                U.val[qFMYb], U.val[qFMYp]);
        int nfail = 0;
        for (int it = -1; it < ntoy; ++it) {
            TRandom3 rnd(1234567u + 100000u * (unsigned)cfg + 1000u * (unsigned)s + (unsigned)(it + 1) * 7919u);
            Obs O;
            MakeToy(U, Din, G, it < 0 ? nullptr : &rnd, O);
            Result R;
            rep = (rf && it >= 0) ? it % nRep : (rb && it >= 0) ? it % bP : -1;
            RunEstimators(O, Din, rep < 0 ? *F : rf ? *Frep[rep] : *Fb[rep], G, R);
            toy = it; nitFM = R.nitFM;
            for (int q = 0; q < kNQ; ++q) { est[q] = R.est[q]; err[q] = R.err[q]; errd[q] = R.errd[q]; errs[q] = R.errs[q]; ok[q] = R.ok[q] ? 1 : 0; }
            if (rb) {
                for (int q = 0; q < kNQ; ++q) { est0[q] = R.est[q]; mct[q] = 0; cmean[q] = R.est[q]; }
                if (rep >= 0) {
                    std::vector<double> sx(kNQ, 0.), sxx(kNQ, 0.);
                    std::vector<int> nx(kNQ, 0);
                    for (int c = 0; c < bC; ++c) {
                        Result Rc;
                        RunEstimators(O, Din, *Fb[bP + rep * bC + c], G, Rc);
                        for (int q = 0; q < kNQ; ++q) if (Rc.ok[q]) { sx[q] += Rc.est[q]; sxx[q] += Rc.est[q] * Rc.est[q]; ++nx[q]; }
                    }
                    for (int q = 0; q < kNQ; ++q) {
                        if (nx[q] < 2) { ok[q] = 0; continue; }
                        const double m = sx[q] / nx[q];
                        cmean[q] = m;
                        mct[q] = std::sqrt(std::max(0., sxx[q] / nx[q] - m * m));
                        est[q] = 2. * est0[q] - m;
                    }
                }
            }
            chi2[0] = R.chi2FM; chi2[1] = R.chi2FF; chi2[2] = R.chi2FG; chi2[3] = R.chi2SA;
            ndf[0] = R.ndfFM; ndf[1] = R.ndfFF; ndf[2] = R.ndfFG; ndf[3] = R.ndfSA;
            for (int k = 0; k < 13; ++k) fmTh[k] = R.fmTh[k];
            for (int k = 0; k < 9; ++k) saTh[k] = R.saTh[k];
            corrAB = R.rhoAbE; corrRK = R.rhoAbB;
            if (!R.ok[qFMa]) ++nfail;
            tt->Fill();
            if (it < 0) {
                fprintf(log, "  Asimov: FM a %.4f (true %.4f) kappa %.4f rho %.4f g %.4f Ybump %.1f (%.1f) Ypeak %.1f (%.1f) chi2 %.1f/%d; "
                             "B4 %.1f PSB %.1f; FF a2 %.4f; FG BR_G %.4g; SA chi2 %.1f/%d\n",
                        R.est[qFMa], U.val[qFMa], R.est[qFMkap], R.est[qFMrho], R.est[qFMg], R.est[qFMYb], U.val[qFMYb], R.est[qFMYp], U.val[qFMYp],
                        R.chi2FM, R.ndfFM, R.est[qB4Yb], R.est[qPSBYp], R.est[qFFa2], R.est[qFGBRG], R.chi2SA, R.ndfSA);
                fflush(log);
            }
        }
        printf("wp4_toys %s: scenario %s done (%d toys, %d F-main failures), %.1f s\n", cfgTag[cfg], S.tag, ntoy, nfail, sw.RealTime());
        fprintf(log, "  toys %d, F-main failures %d\n", ntoy, nfail);
        sw.Continue();
    }
    fo->cd(); tt->Write(); tu->Write(); if (rf) tr->Write();
    fo->Close();
    fclose(log);
    printf("wp4_toys %s: finished in %.1f s\n", cfgTag[cfg], sw.RealTime());
}

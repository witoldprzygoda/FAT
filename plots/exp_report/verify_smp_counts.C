// ============================================================================
// verify_smp_counts.C  --  independent re-count of the main numbers of the
//                          chapter "samples, luminosity and cut flow" (smp)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// Verification macro (written by the verifier of the smp chapter). It does NOT
// use scan_compact.root, ScanEngine.h or smp_common.h: it reads scan_nt of the
// three FAT outputs directly and applies the selections with its own code.
//   pre      = cor_ok == 1 && isBest >= 0 && eVertReco_z > -160 mm
//   vertex   = pip_isUsedVertex == 1 && pim_isUsedVertex == 1
//   dz_pi    = max(|pip_z - vz|, |pim_z - vz|) (float arithmetic)
//   peak nominal = pre && oa_lab_cor < 55 && m_pippim_cor > 0.30
//   bump nominal = pre && oa_lab_cor < 40 && 0.30 < m_pippim_cor < 0.41
//   DIRECT = pre && dz_pi < 40 ; PEAK = DIRECT && mm_pippimepem_cor < 2.85
//   X17 = pre && vertex && 15.5 <= 1000 m_ee_cor < 24.5 && M4 window
//         (sim 0.5333-0.5587 and data 0.53364-0.55416)
// Event group = adjacent rows with equal (file_id, runnumber, event).
// Outputs (directory /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/verify_smp):
//   counts.txt  : rows / events / isBest==1 rows per selection (also OA > 2),
//                 rows in the data windows; per-day rows after pre with the
//                 event-group variance; sector-2 shares; runs present;
//                 like-sign vz modes
//   hists.root  : M4 (2.5 MeV bins, 0.35-0.95) of us/pp/mm for pre and
//                 peak nominal + vertex per period (I, IIA, IIB, III, all),
//                 isBest >= 0
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_smp_verify", kTRUE)' 'plots/exp_report/verify_smp_counts.C+'
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TTree.h"

namespace vsc {
const char* kProd = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530";
const char* kOut  = "/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/verify_smp";
int Period(int d) { if (d >= 34 && d <= 46) return 0; if (d <= 48 && d >= 47) return 1; if (d >= 49 && d <= 53) return 2; if (d >= 59 && d <= 67) return 3; return -1; }
}

void verify_smp_counts()
{
    using namespace vsc;
    gSystem->mkdir(kOut, kTRUE);
    FILE* o = fopen(Form("%s/counts.txt", kOut), "w");
    TFile* fh = TFile::Open(Form("%s/hists.root", kOut), "RECREATE");
    const char* fn[3] = {"output_pippimepem_exp.root", "output_pippimepep_exp.root", "output_pippimemem_exp.root"};
    const char* hn[3] = {"us", "pp", "mm"};
    const int NSEL = 13;
    const char* sname[NSEL] = {"fat", "cor", "pre", "pre+vtx", "pre+vtx+dz", "DIRECT", "PEAK", "pnom", "pnomv", "bnom", "bnomv", "x17sim", "x17data"};
    const double pkLo = 0.5439 - 1.9 * 0.0054, pkHi = 0.5439 + 1.9 * 0.0054;   // 0.53364 - 0.55416
    const double buLo = 0.4550, buHi = 0.5225;
    // M4 histograms: [k][sel 0 pre / 1 pnomv][group 0..3 periods, 4 all]
    TH1D* h[3][2][5];
    for (int k = 0; k < 3; ++k) for (int s = 0; s < 2; ++s) for (int g = 0; g < 5; ++g) {
        h[k][s][g] = new TH1D(Form("m4_%s_%s_%d", hn[k], s == 0 ? "pre" : "pnomv", g), "", 240, 0.35, 0.95);
    }
    std::set<int> runsUS;
    for (int k = 0; k < 3; ++k) {
        TFile* f = TFile::Open(Form("%s/%s", kProd, fn[k]), "READ");
        TTree* t = (TTree*)f->Get("scan_nt");
        t->SetBranchStatus("*", 0);
        const char* br[] = {"file_id", "runnumber", "event", "isBest", "cor_ok", "eVertReco_z", "trigbit", "m_pippimepem_cor", "m_ee_cor",
                            "m_pippim_cor", "mm_pippimepem_cor", "oa", "oa_pippim_epem_lab_cor", "pip_isUsedVertex", "pim_isUsedVertex",
                            "pip_z", "pim_z", "ep_sector", "em_sector", "pip_sector", "pim_sector"};
        const int NB = sizeof(br) / sizeof(br[0]);
        std::vector<Float_t> v(NB, 0.f);
        for (int b = 0; b < NB; ++b) { t->SetBranchStatus(br[b], 1); t->SetBranchAddress(br[b], &v[b]); }
        const Long64_t n = t->GetEntries();
        // counters [variant 0 all OA, 1 OA>2][sel]
        double rows[2][NSEL] = {}, ev[2][NSEL] = {}, best[2][NSEL] = {}, wpk[2][NSEL] = {}, wbu[2][NSEL] = {};
        bool gpass[2][NSEL] = {};
        // per day after pre (variant 0): rows, sum k^2, isBest==1 rows
        std::map<int, double> dRows, dVar, dBest; int gk = 0; int gday = 0;
        // sectors after pre: [0 = days 034-044, 1 = others][0 leptons 1 pions][sector]
        double sec[2][2][6] = {}; double rowS2[2][2] = {};   // [group][0 rows, 1 with any sector-2 track]
        long long badTrig = 0, badVz = 0, nanM4cor1 = 0;
        // like-sign / unlike-sign vz histogram (1 mm), all FAT rows
        TH1D hvz(Form("vz_%s", hn[k]), "", 400, -260, 140);
        float lf = -1, lr = -1, le = -1; bool have = false;
        auto flush = [&]() {
            for (int a = 0; a < 2; ++a) for (int s = 0; s < NSEL; ++s) { if (gpass[a][s]) ev[a][s] += 1; gpass[a][s] = false; }
            if (gk > 0) { dVar[gday] += (double)gk * gk; }
            gk = 0;
        };
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            const float fid = v[0], run = v[1], evn = v[2], isB = v[3], cor = v[4], vz = v[5], trig = v[6];
            const float m4 = v[7], mee = v[8], mpp = v[9], mm4 = v[10], oa = v[11], oal = v[12];
            const float pvu = v[13], mvu = v[14], pz = v[15], mz = v[16];
            if (!have || fid != lf || run != lr || evn != le) { if (have) flush(); lf = fid; lr = run; le = evn; have = true; gday = (int)std::lround(fid); }
            if (trig != 8192.f) badTrig++;
            if (!(vz > -500.f)) badVz++;
            if (k == 0) runsUS.insert((int)std::lround(run));
            hvz.Fill(vz);
            const bool b1 = (isB == 1.f);
            const bool pre = (cor == 1.f) && (isB >= 0.f) && (vz > -160.f);
            if (cor == 1.f && !(m4 == m4)) nanM4cor1++;
            const bool vtx = (pvu == 1.f) && (mvu == 1.f);
            const float dzpi = std::max(std::fabs(pz - vz), std::fabs(mz - vz));
            const bool dz = dzpi < 40.f;
            const bool mmc = mm4 < 2.85f;
            const bool pkin = oal < 55.f && mpp > 0.30f;
            const bool bkin = oal < 40.f && mpp > 0.30f && mpp < 0.41f;
            const float meeM = mee * 1000.f;
            const bool x17ee = meeM >= 15.5f && meeM < 24.5f;
            const double m = m4;
            const bool wsim = m >= 0.5333 && m < 0.5587, wdat = m >= pkLo && m < pkHi;
            bool p[NSEL];
            p[0] = true; p[1] = cor == 1.f; p[2] = pre; p[3] = pre && vtx; p[4] = pre && vtx && dz; p[5] = pre && dz; p[6] = pre && dz && mmc;
            p[7] = pre && pkin; p[8] = pre && pkin && vtx; p[9] = pre && bkin; p[10] = pre && bkin && vtx;
            p[11] = pre && vtx && x17ee && wsim; p[12] = pre && vtx && x17ee && wdat;
            for (int a = 0; a < 2; ++a) {
                if (a == 1 && !(oa > 2.f)) continue;
                for (int s = 0; s < NSEL; ++s) if (p[s]) {
                    rows[a][s] += 1; if (b1) best[a][s] += 1; gpass[a][s] = true;
                    if (m >= pkLo && m < pkHi) wpk[a][s] += 1;
                    if (m >= buLo && m < buHi) wbu[a][s] += 1;
                }
            }
            if (pre) {
                const int d = (int)std::lround(fid);
                dRows[d] += 1; if (b1) dBest[d] += 1; gk++;
                const int g = (d <= 44) ? 0 : 1;
                bool s2 = false;
                for (int q = 0; q < 4; ++q) {
                    const int sc = (int)std::lround(v[17 + q]);
                    if (sc < 0 || sc > 5) continue;
                    sec[g][q < 2 ? 0 : 1][sc] += 1;
                    if (sc == 2) s2 = true;
                }
                rowS2[g][0] += 1; if (s2) rowS2[g][1] += 1;
                const int per = Period(d);
                if (m == m) {
                    h[k][0][per]->Fill(m); h[k][0][4]->Fill(m);
                    if (pkin && vtx) { h[k][1][per]->Fill(m); h[k][1][4]->Fill(m); }
                }
            }
        }
        if (have) flush();
        fprintf(o, "\n== %s (%s): %lld rows; rows with trigbit != 8192: %lld; rows with vz <= -500: %lld; cor_ok==1 rows with NaN M4: %lld\n",
                hn[k], fn[k], n, badTrig, badVz, nanM4cor1);
        for (int a = 0; a < 2; ++a) {
            fprintf(o, "-- variant %s\n", a == 0 ? "isBest>=0 (all OA)" : "OA > 2");
            for (int s = 0; s < NSEL; ++s)
                fprintf(o, "SEL %s %d %-12s rows %9.0f events %9.0f isBest1 %9.0f Npk %7.0f Nbu %7.0f\n", hn[k], a, sname[s], rows[a][s], ev[a][s], best[a][s], wpk[a][s], wbu[a][s]);
        }
        for (auto& q : dRows) fprintf(o, "DAY %s %d %.0f %.0f %.0f\n", hn[k], q.first, q.second, dVar[q.first], dBest[q.first]);
        for (int g = 0; g < 2; ++g) {
            for (int tt = 0; tt < 2; ++tt) {
                double tot = 0; for (int q = 0; q < 6; ++q) tot += sec[g][tt][q];
                fprintf(o, "SECSHARE %s %s %s", hn[k], g == 0 ? "034-044" : "045-067", tt == 0 ? "leptons" : "pions");
                for (int q = 0; q < 6; ++q) fprintf(o, " %.4f", sec[g][tt][q] / tot);
                fprintf(o, "\n");
            }
            fprintf(o, "SECROW %s %s rows %.0f with_sector2 %.0f frac %.4f\n", hn[k], g == 0 ? "034-044" : "045-067", rowS2[g][0], rowS2[g][1], rowS2[g][1] / rowS2[g][0]);
        }
        // vz mode (3 mm running sum) and fraction in -170..-150
        int bm = 1; double vm = 0;
        for (int b = 2; b < hvz.GetNbinsX(); ++b) { const double s3 = hvz.GetBinContent(b - 1) + hvz.GetBinContent(b) + hvz.GetBinContent(b + 1); if (s3 > vm) { vm = s3; bm = b; } }
        const double inN = hvz.Integral(hvz.FindBin(-169.5), hvz.FindBin(-150.5)), all = hvz.Integral(0, hvz.GetNbinsX() + 1);
        const double above160 = hvz.Integral(hvz.FindBin(-159.5), hvz.GetNbinsX() + 1);
        fprintf(o, "VZ %s mode %.1f mm; fraction in [-170,-150) %.4f; fraction vz > -160 (FAT rows) %.4f\n", hn[k], hvz.GetBinCenter(bm), inN / all, above160 / all);
        f->Close(); delete f;
    }
    fprintf(o, "\nRUNS_US %zu\n", runsUS.size());
    FILE* rr = fopen(Form("%s/runs_us.txt", kOut), "w");
    for (int r : runsUS) fprintf(rr, "%d\n", r);
    fclose(rr);
    fclose(o);
    fh->cd();
    for (int k = 0; k < 3; ++k) for (int s = 0; s < 2; ++s) for (int g = 0; g < 5; ++g) h[k][s][g]->Write();
    fh->Close();
    printf("verify_smp_counts done\n");
}

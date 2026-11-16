// ============================================================================
// verify_smp_x17cb.C  --  CB/N of the X17 region with and without the M4 window
//                         inside the selection (verification of the smp chapter)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// The smp cut flow quotes CB/N for the X17 region whose selection already contains
// the M4 window; the like-sign 10 MeV bins at the window edges are then only partly
// filled and the interpolated CB is distorted. Here the CB is formed from the
// selection WITHOUT the M4 cut (pre + both pions in the vertex fit + 15.5 <= M(ee)
// < 24.5 MeV/c2) and integrated over the data peak window 0.53364-0.55416 and the
// frozen simulation window 0.5333-0.5587 (ScanEngine FillMany, k = 1).
// Input: scan_v2/scan_compact.root. Output: printout only.
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_smp_verify", kTRUE)' 'plots/exp_report/verify_smp_x17cb.C+'
// ============================================================================
#include "../ScanEngine.h"

void verify_smp_x17cb()
{
    scan::Data d;
    if (!scan::Load(d)) return;
    auto pre = [](const scan::Row& r) { return r.cor_ok == 1 && r.isBest >= 0 && r.vz > -160.f; };
    auto vtx = [](const scan::Row& r) { return r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; };
    auto ee  = [](const scan::Row& r) { const float m = r.mee * 1000.f; return m >= 15.5f && m < 24.5f; };
    std::vector<scan::Sel> s;
    s.push_back([=](const scan::Row& r) { return pre(r) && vtx(r) && ee(r); });                                   // no M4 cut
    s.push_back([=](const scan::Row& r) { return pre(r) && vtx(r) && ee(r) && r.m4 >= 0.5333f && r.m4 < 0.5587f; }); // as in smp_flow (sim window)
    scan::FillOpt fo; fo.windows.push_back(std::make_pair(0.5333, 0.5587));
    std::vector<scan::Spectra> sp = scan::FillMany(d, s, fo);
    const char* nm[2] = {"X17 ee + vertex, no M4 cut", "X17 with the sim M4 window in the selection"};
    for (int i = 0; i < 2; ++i)
        printf("%-46s: data window N %.0f CB %.1f +- %.1f CB/N %.4f ; sim window N %.0f CB %.1f +- %.1f CB/N %.4f\n", nm[i],
               sp[i].win[0].n[0], sp[i].win[0].cb, std::sqrt(sp[i].win[0].cbVar), sp[i].win[0].cbShare(),
               sp[i].win[2].n[0], sp[i].win[2].cb, std::sqrt(sp[i].win[2].cbVar), sp[i].win[2].cbShare());
}

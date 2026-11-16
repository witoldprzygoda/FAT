// ============================================================================
// verify_smp_cbscale.C  --  sensitivity of the eta-peak content per period to a
//                           +-3 % change of the CB (verification of the smp chapter)
// ============================================================================
// HADES p+p, T = 4.53 GeV (feb22), channel pi+ pi- e+ e-. EXPERIMENTAL DATA ONLY.
// P_main (ScanEngine FitPeakMain: chi2 fit of N - CB, frozen Gaussian 0.5439 GeV/c2 /
// 5.4 MeV + pol2 in 0.5225-0.6100) with CB = k x 2 sqrt(N++ N--), k = 0.97, 1.00, 1.03,
// for pre and peak nominal + vertex (isBest >= 0), per period I, IIA, IIB, III and all
// days. s = S x 64/63 / L with L = START-free luminosity of the analysed runs per period
// (numbers of smp_lumi.txt). Input: scan_v2/scan_compact.root. Output: printout only.
// Run (from /home/przygoda/HADES/FAT/FAT):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/report_exp/aclic_smp_verify", kTRUE)' 'plots/exp_report/verify_smp_cbscale.C+'
// ============================================================================
#include "../ScanEngine.h"

void verify_smp_cbscale()
{
    scan::Data d;
    if (!scan::Load(d)) return;
    auto per = [](int day) { return (day <= 46) ? 0 : (day <= 48 ? 1 : (day <= 53 ? 2 : 3)); };
    auto pre = [](const scan::Row& r) { return r.cor_ok == 1 && r.isBest >= 0 && r.vz > -160.f; };
    auto pnv = [=](const scan::Row& r) { return pre(r) && r.oa_lab < 55.f && r.mpipi > 0.30f && r.pip_isUsedVertex == 1 && r.pim_isUsedVertex == 1; };
    const double L[5] = {2809.453, 372.610, 1078.022, 1819.110, 6079.195};
    const char* gn[5] = {"I", "IIA", "IIB", "III", "all"};
    std::vector<scan::Sel> s;
    for (int b = 0; b < 2; ++b)
        for (int g = 0; g < 5; ++g)
            s.push_back([=](const scan::Row& r) { return (g == 4 || per(r.file_id) == g) && (b == 0 ? pre(r) : pnv(r)); });
    const std::vector<scan::Spectra> sp = scan::FillMany(d, s);
    for (int b = 0; b < 2; ++b)
        for (int g = 0; g < 5; ++g) {
            double S[3];
            const double kf[3] = {0.97, 1.00, 1.03};
            for (int q = 0; q < 3; ++q) {
                scan::Spectra x = sp[b * 5 + g];
                for (auto& v : x.cbW) v *= kf[q];
                for (auto& v : x.cbWVar) v *= kf[q] * kf[q];
                for (int i = 0; i < x.nbin; ++i) {
                    x.cb[i] *= kf[q]; x.cbVar[i] *= kf[q] * kf[q];
                    x.sig[i] = x.n[0][i] - x.cb[i]; x.sigVar[i] = x.n[0][i] + x.cbVar[i];
                }
                S[q] = scan::FitPeakMain(x).S;
            }
            printf("CBSCALE %-6s %-4s S(k=0.97) %.1f S(k=1) %.1f S(k=1.03) %.1f  dS/S(+3%%) %+.4f  s(k=1) %.4f\n", b == 0 ? "pre" : "pnomv", gn[g],
                   S[0], S[1], S[2], S[2] / S[1] - 1, S[1] * 64. / 63. / L[g]);
        }
}

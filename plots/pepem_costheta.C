// pepem_costheta.C — cos(theta_CMS) of pe+e- system
// Usage: root -l -b -q plots/pepem_costheta.C

#include "PlotUtils.h"

void pepem_costheta() {

    PlotUtils pu("output_pepem.root", "output_pepep.root", "output_pemem.root");

    // cos(theta_CMS) of pe+e- (OA > 9, MM proton window)
    TH1D *all1, *cb1, *sig1;
    std::tie(all1, cb1, sig1) = pu.drawSignal("dilepton_nt", "pepem_costheta_cms",
                                                25, -1.0, 1.0,
                                                "oa>9 && m_ee>0.14 && mm_pepem_mass>=0.88 && mm_pepem_mass<=1.02",
                                                ";cos#theta_{CMS}(pe^{+}e^{-});Counts");

    auto* c1 = pu.drawTriple(all1, cb1, sig1,
                              "cos#theta_{CMS}(pe^{+}e^{-}) (OA > 9#circ, M_{ee} > 0.14, MM proton)",
                              "c_pepem_costheta");
    c1->SetLogy(0);
    all1->SetMaximum(100);
    all1->SetMinimum(0);
    c1->Update();
    pu.save(c1, "pepem_costheta_cms");

    double i_all = all1->Integral();
    double i_cb  = cb1->Integral();
    double i_sig = sig1->Integral();

    std::cout << "\n=== cos(theta_CMS) pe+e- (OA>9, MM proton) ===\n";
    std::cout << "  Integral:  all = " << i_all
              << "  CB = " << i_cb
              << "  sig = " << i_sig << "\n";

    std::cout << "\nDone. Check plots/output/\n";
}

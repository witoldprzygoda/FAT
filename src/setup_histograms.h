/**
 * @file setup_histograms.h
 * @brief Histogram definitions
 *
 * This file defines all histograms for the analysis.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_HISTOGRAMS_H
#define SETUP_HISTOGRAMS_H

#include "manager.h"
#include <iostream>

/**
 * @brief Define histograms for the analysis
 * @param mgr Manager object for histogram creation
 */
inline void setupHistograms(Manager& mgr) {
    std::cout << "Setting up histograms...\n";

    // Lepton momentum (RECONSTRUCTED + CORRECTED)
    mgr.create1D("ep_p", "e^{+} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "leptons");
    mgr.create1D("em_p", "e^{-} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "leptons");
    mgr.create1D("ep_p_cor", "e^{+} momentum (corr);p [MeV/c];Counts",
                 100, 0, 2000, "leptons");
    mgr.create1D("em_p_cor", "e^{-} momentum (corr);p [MeV/c];Counts",
                 100, 0, 2000, "leptons");

    // Momentum correction: delta_p vs p_reconstructed
    mgr.create2D("ep_dp_vs_p", "e^{+}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");
    mgr.create2D("em_dp_vs_p", "e^{-}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");

    // Dilepton invariant mass (REC + COR)
    mgr.create1D("mass_ee", "e^{+}e^{-} invariant mass;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");
    mgr.create1D("mass_ee_cor", "e^{+}e^{-} invariant mass (corr);M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");

    // e+e-gamma invariant mass — pi0 Dalitz candidate (only when ecal_mult == 1)
    mgr.create1D("mass_epemg", "e^{+}e^{-}#gamma invariant mass;M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "epemg");
    mgr.create1D("mass_epemg_cor", "e^{+}e^{-}#gamma invariant mass (corr);M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "epemg");

    // gamma+gamma invariant mass (only when ecal_mult == 2 with both passing quality)
    // ECAL clusters have no momentum correction — single histogram only.
    mgr.create1D("mass_gg", "#gamma#gamma invariant mass;M_{#gamma#gamma} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "ecal_only");

    // 4-gamma invariant mass under double-pi0 narrow constraint
    // (ecal_mult == 4, all gammas pass quality, both gg pairings in narrow pi0 window).
    // Filled across all 3 pairings, like a combinatorial sum.
    mgr.create1D("mass_gggg_pi0pi0",
                 "4#gamma invariant mass under (#pi^{0}, #pi^{0});M_{#gamma#gamma#gamma#gamma} [GeV/c^{2}];Counts",
                 200, 0.2, 2.2, "ecal_only");

    // Opening angle and before/after OA cut mass spectra (REC + COR)
    mgr.create1D("opening_angle", "e^{+}e^{-} opening angle;#theta_{open} [deg];Counts",
                 180, 0, 180, "dilepton");
    mgr.create1D("mass_ee_before_oa", "M_{ee} before OA cut;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");
    mgr.create1D("mass_ee_after_oa", "M_{ee} after OA cut;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");
    mgr.create1D("mass_ee_before_oa_cor", "M_{ee} before OA cut (corr);M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");
    mgr.create1D("mass_ee_after_oa_cor", "M_{ee} after OA cut (corr);M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");

    // CMS frame observables (after boost to beam-target CMS)
    mgr.create1D("rapidity_cms", "Dilepton rapidity in CMS;y_{CMS};Counts",
                 100, -2.0, 2.0, "cms");
    mgr.create1D("pt_cms", "Dilepton transverse momentum;p_{T} [MeV/c];Counts",
                 100, 0, 1500, "cms");
    mgr.create1D("theta_cms", "Dilepton polar angle in CMS;#theta_{CMS} [deg];Counts",
                 90, 0, 180, "cms");
    mgr.create2D("rapidity_vs_mass", "y_{CMS} vs M_{ee};M_{e^{+}e^{-}} [GeV/c^{2}];y_{CMS}",
                 100, 0.0, 1.0, 100, -2.0, 2.0, "cms");

    std::cout << "  Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

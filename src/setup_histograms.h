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

    // Lepton momentum
    mgr.create1D("ep_p", "e^{+} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "leptons");
    mgr.create1D("em_p", "e^{-} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "leptons");

    // Momentum correction: delta_p vs p_reconstructed
    mgr.create2D("ep_dp_vs_p", "e^{+}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");
    mgr.create2D("em_dp_vs_p", "e^{-}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");

    // Dilepton invariant mass
    mgr.create1D("mass_ee", "e^{+}e^{-} invariant mass;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");

    // Opening angle and before/after OA cut mass spectra
    mgr.create1D("opening_angle", "e^{+}e^{-} opening angle;#theta_{open} [deg];Counts",
                 180, 0, 180, "dilepton");
    mgr.create1D("mass_ee_before_oa", "M_{ee} before OA cut;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");
    mgr.create1D("mass_ee_after_oa", "M_{ee} after OA cut;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
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

    // pe+e- compound histograms (FWD proton + dilepton)
    mgr.create1D("mm_pepem", "MM(pe^{+}e^{-});MM(pe^{+}e^{-}) [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "pepem");
    mgr.create1D("m_ee_mm_pepem", "M_{ee} (OA>9, MM proton);M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 80, 0.0, 0.8, "pepem");
    mgr.create1D("pepem_inv_mass", "M(pe^{+}e^{-}) (OA>9, M_{ee}>0.14, MM proton);M(pe^{+}e^{-}) [GeV/c^{2}];Counts",
                 50, 0.8, 1.8, "pepem");
    mgr.create1D("pepem_cms_costheta", "cos#theta_{CMS}(pe^{+}e^{-}) (OA>9, M_{ee}>0.14, MM proton);cos#theta_{CMS};Counts",
                 25, -1.0, 1.0, "pepem");

    std::cout << "  Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

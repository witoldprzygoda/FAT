/**
 * @file setup_histograms.h
 * @brief Histogram definitions for the pi+pi- hadronic analysis
 *
 * Channels:
 *   - pippim       = pi+ pi-
 *   - gg           = gamma + gamma (from ECAL, ecal_mult == 2)
 *   - pippimgg     = pippim + gg          (for eta -> pi+ pi- pi0 -> pi+ pi- gamma gamma)
 *
 * All compound observables exist in two flavors: RECONSTRUCTED (default) and
 * CORRECTED (suffix "_cor"), drawn from KinematicType::CORRECTED of the inputs.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_HISTOGRAMS_H
#define SETUP_HISTOGRAMS_H

#include "manager.h"
#include <iostream>

inline void setupHistograms(Manager& mgr) {
    std::cout << "Setting up histograms...\n";

    // ========================================================================
    // Pion momentum (RECONSTRUCTED + CORRECTED)
    // ========================================================================
    mgr.create1D("pip_p", "#pi^{+} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "pions");
    mgr.create1D("pim_p", "#pi^{-} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "pions");

    mgr.create1D("pip_p_cor", "#pi^{+} momentum (corrected);p [MeV/c];Counts",
                 100, 0, 2000, "pions_cor");
    mgr.create1D("pim_p_cor", "#pi^{-} momentum (corrected);p [MeV/c];Counts",
                 100, 0, 2000, "pions_cor");

    // Momentum corrections
    mgr.create2D("pip_dp_vs_p", "#pi^{+}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");
    mgr.create2D("pim_dp_vs_p", "#pi^{-}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");

    // ========================================================================
    // Compound invariant masses
    // ========================================================================
    mgr.create1D("mass_pippim", "#pi^{+}#pi^{-} invariant mass;M_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 2.0, "compound");
    mgr.create1D("mass_gg", "#gamma#gamma invariant mass;M_{#gamma#gamma} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "compound");
    mgr.create1D("mass_pippimgg", "#pi^{+}#pi^{-}#gamma#gamma invariant mass;M_{#pi^{+}#pi^{-}#gamma#gamma} [GeV/c^{2}];Counts",
                 200, 0.2, 1.4, "compound");

    mgr.create1D("mass_pippim_cor", "#pi^{+}#pi^{-} invariant mass (cor);M_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 2.0, "compound_cor");
    mgr.create1D("mass_gg_cor", "#gamma#gamma invariant mass (cor);M_{#gamma#gamma} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "compound_cor");
    mgr.create1D("mass_pippimgg_cor", "#pi^{+}#pi^{-}#gamma#gamma invariant mass (cor);M_{#pi^{+}#pi^{-}#gamma#gamma} [GeV/c^{2}];Counts",
                 200, 0.2, 1.4, "compound_cor");

    // ========================================================================
    // Missing masses: MM(X) = beam + target - X
    // ========================================================================
    mgr.create1D("mm_pippim", "MM(#pi^{+}#pi^{-});MM(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses");
    mgr.create1D("mm_pippimgg", "MM(#pi^{+}#pi^{-}#gamma#gamma);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses");

    mgr.create1D("mm_pippim_cor", "MM(#pi^{+}#pi^{-}) (cor);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses_cor");
    mgr.create1D("mm_pippimgg_cor", "MM(#pi^{+}#pi^{-}#gamma#gamma) (cor);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses_cor");

    // ========================================================================
    // Opening angles (LAB)
    // ========================================================================
    mgr.create1D("oa_pip_pim", "Opening angle #pi^{+}#pi^{-};#theta_{open} [deg];Counts",
                 180, 0, 180, "angles");
    mgr.create1D("oa_gg", "Opening angle #gamma#gamma;#theta_{open} [deg];Counts",
                 180, 0, 180, "angles");

    // ========================================================================
    // CMS frame observables (after boost to beam-target CMS) — pippim
    // ========================================================================
    mgr.create1D("rapidity_cms", "Pippim rapidity in CMS;y_{CMS};Counts",
                 100, -2.0, 2.0, "cms");
    mgr.create1D("pt_cms", "Pippim transverse momentum;p_{T} [MeV/c];Counts",
                 100, 0, 1500, "cms");
    mgr.create1D("theta_cms", "Pippim polar angle in CMS;#theta_{CMS} [deg];Counts",
                 90, 0, 180, "cms");

    mgr.create1D("rapidity_cms_cor", "Pippim rapidity in CMS (cor);y_{CMS};Counts",
                 100, -2.0, 2.0, "cms_cor");
    mgr.create1D("pt_cms_cor", "Pippim transverse momentum (cor);p_{T} [MeV/c];Counts",
                 100, 0, 1500, "cms_cor");
    mgr.create1D("theta_cms_cor", "Pippim polar angle in CMS (cor);#theta_{CMS} [deg];Counts",
                 90, 0, 180, "cms_cor");

    std::cout << "  Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

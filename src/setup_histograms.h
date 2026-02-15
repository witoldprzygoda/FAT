/**
 * @file setup_histograms.h
 * @brief Histogram definitions
 *
 * This file defines all histograms for the pp analysis.
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

    // ppi0 invariant mass (proton + pi0 - delta system)
    mgr.create1D("ppi0_inv_mass", "M(p#pi^{0}) invariant mass;M(p#pi^{0}) [GeV/c^{2}];Counts",
                 1000, 1.0, 2.0, "pp");

    // pp invariant mass (proton + proton)
    mgr.create1D("pp_inv_mass", "M(pp) invariant mass;M(pp) [GeV/c^{2}];Counts",
                 1000, 1.8, 2.4, "pp");

    // pi0 missing mass (beam + target - proton1 - proton2)
    mgr.create1D("pi0_miss_mass", "MM(pp) #pi^{0} missing mass;MM(pp) [GeV/c^{2}];Counts",
                 1000, 0.0, 1.0, "pp");

    // pi0 missing mass squared (includes negative values due to resolution)
    mgr.create1D("pi0_miss_mass2", "MM^{2}(pp) #pi^{0} missing mass squared;MM^{2}(pp) [GeV^{2}/c^{4}];Counts",
                 1000, -0.5, 0.4, "pp");

    // ========================================================================
    // 2D Debug Histograms
    // ========================================================================

    // dphi vs tantan (anti-elastic check)
    mgr.create2D("dphi_vs_tantan",
                 "#Delta#phi vs tan#theta_{1}#upointtan#theta_{2};"
                 "#Delta#phi [deg];tan#theta_{1}#upointtan#theta_{2}",
                 1000, 160.0, 200.0, 1000, 0.0, 1.0, "debug");

    // Missing mass p1 vs missing mass p2
    mgr.create2D("mm_p1_vs_mm_p2",
                 "MM(p_{1}) vs MM(p_{2});"
                 "MM(p_{1}) [GeV/c^{2}];MM(p_{2}) [GeV/c^{2}]",
                 1000, 1.0, 1.6, 1000, 1.0, 1.6, "debug");

    // M(ppi0) vs cos_theta_cms(ppi0) - Dalitz-like
    mgr.create2D("mppi0_vs_costh_cms",
                 "M(p#pi^{0}) vs cos#theta_{CMS}(p#pi^{0});"
                 "M(p#pi^{0}) [GeV/c^{2}];cos#theta_{CMS}",
                 1000, 0.8, 1.8, 1000, -1.0, 1.0, "debug");

    // ========================================================================
    // PWA (Partial Wave Analysis) Histograms
    // ========================================================================

    // Group A: cos(theta) of single particles in CMS frame
    mgr.create1D("pwa_pi0_costh", "#pi^{0} cos#theta_{CMS};cos#theta_{CMS};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_p_costh", "Proton cos#theta_{CMS};cos#theta_{CMS};Counts",
                 40, -1.0, 1.0, "pwa");

    // Group B: Momenta in LAB frame
    mgr.create1D("pwa_pi0_p", "#pi^{0} momentum;p [GeV/c];Counts",
                 50, 0.0, 1.0, "pwa");
    mgr.create1D("pwa_p_p", "Proton momentum;p [GeV/c];Counts",
                 50, 0.0, 2.0, "pwa");

    // Group C: Invariant masses of compound systems
    mgr.create1D("pwa_ppi0_m", "M(p#pi^{0});M(p#pi^{0}) [GeV/c^{2}];Counts",
                 60, 1.0, 2.0, "pwa");
    mgr.create1D("pwa_pp_m", "M(pp);M(pp) [GeV/c^{2}];Counts",
                 60, 1.8, 2.4, "pwa");

    // Group D: Helicity distributions
    mgr.create1D("pwa_pi0_helicity", "#pi^{0} helicity angle;cos#theta_{helicity};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_p_helicity", "Proton helicity angle;cos#theta_{helicity};Counts",
                 40, -1.0, 1.0, "pwa");

    // Group E: Gottfried-Jackson distributions
    mgr.create1D("pwa_pi0_gj", "#pi^{0} GJ angle;cos#theta_{GJ};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_p_gj", "Proton GJ angle;cos#theta_{GJ};Counts",
                 40, -1.0, 1.0, "pwa");

    // ========================================================================
    // Additional histograms (D+)
    // ========================================================================

    // Delta (p + pi0) system
    mgr.create1D("mass_deltaP", "M(p#pi^{0}) delta mass;M(p#pi^{0}) [GeV/c^{2}];Counts",
                 100, 1.0, 2.0, "pwa");
    mgr.create1D("cos_theta_deltaP", "cos#theta(p#pi^{0}) delta;cos#theta_{CMS};Counts",
                 40, -1.0, 1.0, "pwa");

    // pi0
    mgr.create1D("cos_theta_pi0", "cos#theta(#pi^{0});cos#theta_{CMS};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("mass_pi0", "M(#pi^{0}) missing mass;M(#pi^{0}) [GeV/c^{2}];Counts",
                 100, 0.0, 0.5, "pwa");

    std::cout << "  Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

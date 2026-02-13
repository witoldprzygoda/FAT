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

    // ppip compound invariant mass (proton + pi+)
    mgr.create1D("ppip_inv_mass", "M(p#pi^{+}) invariant mass;M(p#pi^{+}) [GeV/c^{2}];Counts",
                 100, 1.0, 2.5, "ppip");

    // ppip missing mass (beam + target - proton - pi+)
    mgr.create1D("ppip_miss_mass", "MM(p#pi^{+}) missing mass;MM(p#pi^{+}) [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "ppip");

    // ========================================================================
    // PWA (Partial Wave Analysis) Histograms
    // ========================================================================

    // Group A: cos(theta) of single particles in CMS frame
    mgr.create1D("pwa_pip_costh", "#pi^{+} cos#theta_{CMS};cos#theta_{CMS};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_p_costh", "Proton cos#theta_{CMS};cos#theta_{CMS};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_n_costh", "Neutron cos#theta_{CMS};cos#theta_{CMS};Counts",
                 40, -1.0, 1.0, "pwa");

    // Group B: Momenta in LAB frame
    mgr.create1D("pwa_pip_p", "#pi^{+} momentum;p [GeV/c];Counts",
                 50, 0.0, 1.0, "pwa");
    mgr.create1D("pwa_p_p", "Proton momentum;p [GeV/c];Counts",
                 50, 0.0, 2.0, "pwa");
    mgr.create1D("pwa_n_p", "Neutron momentum;p [GeV/c];Counts",
                 50, 0.0, 2.0, "pwa");

    // Group C: Invariant masses of compound systems
    mgr.create1D("pwa_ppip_m", "M(p#pi^{+});M(p#pi^{+}) [GeV/c^{2}];Counts",
                 60, 1.0, 1.6, "pwa");
    mgr.create1D("pwa_npip_m", "M(n#pi^{+});M(n#pi^{+}) [GeV/c^{2}];Counts",
                 60, 0.9, 1.5, "pwa");
    mgr.create1D("pwa_pn_m", "M(pn);M(pn) [GeV/c^{2}];Counts",
                 60, 1.7, 2.3, "pwa");

    // Group D: Helicity distributions
    mgr.create1D("pwa_pip_helicity", "#pi^{+} helicity angle;cos#theta_{helicity};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_pipn_helicity", "#pi^{+}n helicity angle;cos#theta_{helicity};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_n_helicity", "Neutron helicity angle;cos#theta_{helicity};Counts",
                 40, -1.0, 1.0, "pwa");

    // Group E: Gottfried-Jackson distributions
    mgr.create1D("pwa_pip_gj", "#pi^{+} GJ angle;cos#theta_{GJ};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_pipn_gj", "#pi^{+}n GJ angle;cos#theta_{GJ};Counts",
                 40, -1.0, 1.0, "pwa");
    mgr.create1D("pwa_n_gj", "Neutron GJ angle;cos#theta_{GJ};Counts",
                 40, -1.0, 1.0, "pwa");

    std::cout << "  Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

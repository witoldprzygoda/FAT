/**
 * @file setup_histograms.h
 * @brief Histogram definitions - Tutorial
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

    // ========================================================================
    // STEP 2: First histograms - lepton momentum
    // ========================================================================
    // create1D(name, title, nbins, xmin, xmax, folder)
    //
    // Title format: "Display title;X-axis label;Y-axis label"
    // Folder organizes histograms in the output ROOT file

    // Positron momentum
    mgr.create1D("ep_p", "e^{+} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "leptons");

    // Electron momentum
    mgr.create1D("em_p", "e^{-} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "leptons");

    // ========================================================================
    // STEP 2b: 2D histograms - momentum correction
    // ========================================================================
    // create2D(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, folder)
    //
    // Show (p_corrected - p_reconstructed) vs p_reconstructed
    // This visualizes how much the energy-loss correction changes the momentum

    // Positron: momentum correction vs momentum
    mgr.create2D("ep_dp_vs_p", "e^{+}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");

    // Electron: momentum correction vs momentum
    mgr.create2D("em_dp_vs_p", "e^{-}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");

    // ========================================================================
    // STEP 3: Dilepton invariant mass
    // ========================================================================
    // The dilepton (e+e-) invariant mass is the key observable.
    // M = sqrt((E_ep + E_em)² - (p_ep + p_em)²)

    mgr.create1D("mass_ee", "e^{+}e^{-} invariant mass;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");

    // ========================================================================
    // STEP 5a: Opening angle histogram
    // ========================================================================
    // The opening angle between e+ and e- momentum vectors.
    // Filled BEFORE the opening angle cut, to see full distribution.

    mgr.create1D("opening_angle", "e^{+}e^{-} opening angle;#theta_{open} [deg];Counts",
                 180, 0, 180, "dilepton");

    // ========================================================================
    // STEP 5b: Before/after histograms for opening angle cut
    // ========================================================================
    // Compare mass distribution before and after the opening angle cut.
    // This shows what the cut removes from the mass spectrum.

    mgr.create1D("mass_ee_before_oa", "M_{ee} before OA cut;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");

    mgr.create1D("mass_ee_after_oa", "M_{ee} after OA cut;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");

    // ========================================================================
    // STEP 7: CMS (center of mass) histograms
    // ========================================================================
    // Variables after boosting dilepton to the beam-target CMS frame.
    // In CMS frame: total momentum is zero, rapidity is centered around 0.
    //
    // Key observables in CMS:
    //   y_cms    - rapidity (velocity-like variable, Lorentz additive)
    //   pt       - transverse momentum (invariant under z-boost)
    //   theta_cms - polar angle in CMS

    mgr.create1D("rapidity_cms", "Dilepton rapidity in CMS;y_{CMS};Counts",
                 100, -2.0, 2.0, "cms");

    mgr.create1D("pt_cms", "Dilepton transverse momentum;p_{T} [MeV/c];Counts",
                 100, 0, 1500, "cms");

    mgr.create1D("theta_cms", "Dilepton polar angle in CMS;#theta_{CMS} [deg];Counts",
                 90, 0, 180, "cms");

    // 2D: rapidity vs mass (useful for physics interpretation)
    mgr.create2D("rapidity_vs_mass", "y_{CMS} vs M_{ee};M_{e^{+}e^{-}} [GeV/c^{2}];y_{CMS}",
                 100, 0.0, 1.0, 100, -2.0, 2.0, "cms");

    std::cout << "  Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

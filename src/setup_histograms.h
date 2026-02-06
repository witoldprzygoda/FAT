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

    std::cout << "  Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

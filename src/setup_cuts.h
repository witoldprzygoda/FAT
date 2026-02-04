/**
 * @file setup_cuts.h
 * @brief Cut definitions for e+e- dilepton analysis with ECAL photons
 *
 * This file contains all cut definitions.
 * Edit this file to add/modify cuts for your analysis.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_CUTS_H
#define SETUP_CUTS_H

#include "cut_manager.h"
#include <iostream>

/**
 * @brief Define all cuts for the e+e- dilepton analysis
 * @param cuts CutManager object for cut definition
 *
 * EDIT THIS FUNCTION to customize cuts for your analysis.
 */
inline void setupCuts(CutManager& cuts) {
    std::cout << "Setting up cuts...\n";

    // ========================================================================
    // Range Cuts (min <= value <= max)
    // ========================================================================

    // Dilepton invariant mass cut
    // Example: select pi0 Dalitz region (M < 0.135 GeV)
    // Example: select eta Dalitz region (0.135 < M < 0.548 GeV)
    // Example: select omega/rho region (0.6 < M < 0.9 GeV)
    cuts.defineRangeCut("dilepton_mass", 0.0, 1.5, "Dilepton mass window [GeV]");

    // Opening angle cut (to reject conversion pairs at small angles)
    // Typical cut: > 9 degrees to reject gamma conversions
    cuts.defineRangeCut("opening_angle", 9.0, 180.0, "Opening angle cut [deg]");

    // ========================================================================
    // ECAL Photon Quality Cuts
    // ========================================================================

    // Photon minimum energy cut [MeV]
    // Typical: > 50-100 MeV to reject noise
    cuts.defineRangeCut("photon_energy", 50.0, 10000.0, "Photon energy window [MeV]");

    // Photon theta angular acceptance [deg]
    // ECAL covers roughly 12-45 degrees
    cuts.defineRangeCut("photon_theta", 12.0, 45.0, "Photon theta acceptance [deg]");

    // Example: Momentum cuts (uncomment if needed)
    // cuts.defineRangeCut("ep_momentum", 50.0, 2000.0, "e+ momentum range [MeV/c]");
    // cuts.defineRangeCut("em_momentum", 50.0, 2000.0, "e- momentum range [MeV/c]");

    // Example: Angular cuts (uncomment if needed)
    // cuts.defineRangeCut("ep_theta", 18.0, 85.0, "e+ theta range [deg]");
    // cuts.defineRangeCut("em_theta", 18.0, 85.0, "e- theta range [deg]");

    // ========================================================================
    // Trigger Cuts (uncomment if needed)
    // ========================================================================

    // cuts.defineTriggerCut("physics", 4, false, "PT3 trigger");
    // cuts.defineTriggerCut("lepton", 0x10, false, "Lepton trigger");

    // ========================================================================
    // Graphical Cuts (uncomment and edit path if needed)
    // ========================================================================

    // PID cuts for electron/positron identification
    // try {
    //     cuts.loadGraphicalCut("ep_pid", "cuts/ep_pid.root", "cutg_ep");
    //     cuts.loadGraphicalCut("em_pid", "cuts/em_pid.root", "cutg_em");
    // } catch (const std::exception& e) {
    //     std::cerr << "Warning: Could not load lepton PID cuts: " << e.what() << "\n";
    // }

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

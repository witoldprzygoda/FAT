/**
 * @file setup_cuts.h
 * @brief Cut definitions
 *
 * This file defines all cuts for the analysis.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_CUTS_H
#define SETUP_CUTS_H

#include "cut_manager.h"
#include <iostream>

/**
 * @brief Define cuts for the analysis
 * @param cuts CutManager object for cut definition
 */
inline void setupCuts(CutManager& cuts) {
    std::cout << "Setting up cuts...\n";

    // Event-level cuts (applied before particle creation)
    cuts.defineValueCut("isBest", 1, "Best candidate selection");
    cuts.defineMinCut("vertex_z", -500, "Vertex Z quality [mm]");

    // Opening angle cut (rejects close e+e- pairs)
    cuts.defineMinCut("opening_angle", 9.0, "Opening angle > 9 deg");

    // Dilepton mass cut (above pi0)
    cuts.defineMinCut("m_ee", 0.14, "M_ee > 0.14 GeV/c^2");

    // Forward detector time cut
    cuts.defineMaxCut("fwd_time", 40.0, "FWD time < 40 ns");

    // Missing mass of pe+e- (proton mass window)
    cuts.defineRangeCut("mm_pepem", 0.88, 1.02, "MM(pe+e-) proton window");

    // ECAL quality cuts (AND logic; passCutSet values must match this order)
    cuts.defineCutSet("ecal_quality", "ECAL particle quality")
        .addValueCut("ecal_pid", 1.0, "PID == 1")
        .addRangeCut("ecal_beta", 0.8, 1.2, "0.8 < beta < 1.2")
        .addMinCut("cluster_energy", 100.0, "Cluster energy > 100 MeV");

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

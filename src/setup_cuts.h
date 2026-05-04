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

    // Event-level cuts (applied before particle creation).
    // Simulation: trigger_PT3 and start_detector cuts are NOT applied because
    // SMASH-generated events have no trigger word and no LGAD response.
    cuts.defineValueCut("isBest", 1, "Best candidate selection");
    cuts.defineMinCut("vertex_z", -500, "Vertex Z quality [mm]");

    // Opening angle cuts (reject close e+e- pairs).
    // opening_angle_4 is the ACTIVE cut applied in data.
    // opening_angle_9 is defined but not applied (kept for future studies).
    cuts.defineMinCut("opening_angle_4", 4.0, "Opening angle > 4 deg (active)");
    cuts.defineMinCut("opening_angle_9", 9.0, "Opening angle > 9 deg (for future use)");

    // ECAL quality cuts (AND logic; passCutSet values must match this order)
    cuts.defineCutSet("ecal_quality", "ECAL particle quality")
        .addValueCut("ecal_pid", 1.0, "PID == 1")
        .addRangeCut("ecal_beta", 0.8, 1.2, "0.8 < beta < 1.2")
        .addMinCut("cluster_energy", 100.0, "Cluster energy > 100 MeV");

    // pi0 invariant-mass windows (applicable to either M(gg) or M(e+e-gamma)
    // for Dalitz studies). Wide and narrow versions defined separately.
    cuts.defineRangeCut("pi0_mass_window", 0.10, 0.18,
                        "M in [0.10, 0.18] GeV/c^2 (wide pi0 window)");
    cuts.defineRangeCut("pi0_mass_window_narrow", 0.125, 0.145,
                        "M in [0.125, 0.145] GeV/c^2 (narrow pi0 window)");

    // eta invariant-mass window — applicable on M(e+e-gamma) for the mult==1
    // Dalitz study and the mult==3 rotational analysis.
    cuts.defineRangeCut("eta_mass_window", 0.5, 0.6,
                        "M(e+e-gamma) in [0.5, 0.6] GeV/c^2 (eta window)");

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

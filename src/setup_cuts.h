/**
 * @file setup_cuts.h
 * @brief Cut definitions - Tutorial
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

    // ========================================================================
    // STEP 4a: Event-level cuts
    // ========================================================================
    // These cuts are applied at the BEGINNING of processEvent(),
    // BEFORE creating any particles. If an event fails, we skip it entirely.
    //
    // Cut types:
    //   defineValueCut(name, target, desc)  - value == target
    //   defineMinCut(name, min, desc)       - value > min
    //   defineMaxCut(name, max, desc)       - value < max
    //   defineRangeCut(name, min, max, desc) - min <= value <= max

    // Best candidate selection (isBest must equal 1)
    cuts.defineValueCut("isBest", 1, "Best candidate selection");

    // Vertex Z quality (must be > -500 mm)
    cuts.defineMinCut("vertex_z", -500, "Vertex Z quality [mm]");

    // ========================================================================
    // STEP 5a: Opening angle cut
    // ========================================================================
    // Applied AFTER creating particles but BEFORE combining them.
    // Rejects close e+e- pairs (e.g., from conversions or Dalitz decays).

    cuts.defineMinCut("opening_angle", 9.0, "Opening angle > 9 deg");

    // ========================================================================
    // STEP 8b: ECAL quality cuts
    // ========================================================================
    // CutSet for ECAL particle quality selection.
    // All cuts must pass (AND logic).
    //
    // Order matters! Values passed to passCutSet() must match this order:
    //   1. ecal_pid (exact match = 1)
    //   2. ecal_beta (range 0.8 - 1.2)
    //   3. cluster_energy (min > 100 MeV)

    cuts.defineCutSet("ecal_quality", "ECAL particle quality")
        .addValueCut("ecal_pid", 1.0, "PID == 1")
        .addRangeCut("ecal_beta", 0.8, 1.2, "0.8 < beta < 1.2")
        .addMinCut("cluster_energy", 100.0, "Cluster energy > 100 MeV");

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

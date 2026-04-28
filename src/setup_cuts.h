/**
 * @file setup_cuts.h
 * @brief Cut definitions for the pi+pi- hadronic analysis
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_CUTS_H
#define SETUP_CUTS_H

#include "cut_manager.h"
#include <iostream>

inline void setupCuts(CutManager& cuts) {
    std::cout << "Setting up cuts...\n";

    // Event-level cuts (applied first, before any compound construction)
    cuts.defineValueCut("isBest", 1, "Best candidate selection");
    cuts.defineMinCut("vertex_z", -500, "Vertex Z quality [mm]");

    // Trigger: PT3 only (trigbit == 8192)
    cuts.defineValueCut("trigger_PT3", 8192, "Trigger PT3 (trigbit == 8192)");

    // Start detector (LGAD) — defined as a cut set so future iteration/timing
    // conditions can be appended without changing the call site.
    // Only the first iteration (start_iteration == 3) is required for now.
    cuts.defineCutSet("start_detector", "Start detector (LGAD) response")
        .addValueCut("start_iteration", 3.0, "start_iteration == 3");

    // ECAL gamma quality (AND logic; passCutSet values must match this order)
    cuts.defineCutSet("ecal_quality", "ECAL particle quality")
        .addValueCut("ecal_pid", 1.0, "PID == 1")
        .addRangeCut("ecal_beta", 0.8, 1.2, "0.8 < beta < 1.2")
        .addMinCut("cluster_energy", 100.0, "Cluster energy > 100 MeV");

    // pi0 invariant-mass window for the gamma+gamma hypothesis (wide — for ntuple gating)
    cuts.defineRangeCut("pi0_mass_window", 0.1, 0.18,
                        "M(gg) in [0.10, 0.18] GeV/c^2 (wide pi0 window)");
    // Narrow pi0 window — used for histogram filling now that statistics is huge
    cuts.defineRangeCut("pi0_mass_window_narrow", 0.125, 0.145,
                        "M(gg) in [0.125, 0.145] GeV/c^2 (narrow pi0 window)");

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

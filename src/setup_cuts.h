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

    // Opening angle cuts (reject close e+e- pairs)
    // opening_angle_4 is the ACTIVE cut applied in data.
    // opening_angle_9 is defined but not applied (kept for future studies).
    cuts.defineMinCut("opening_angle_4", 4.0, "Opening angle > 4 deg (active)");
    cuts.defineMinCut("opening_angle_9", 9.0, "Opening angle > 9 deg (for future use)");

    // ECAL quality cuts (AND logic; passCutSet values must match this order)
    cuts.defineCutSet("ecal_quality", "ECAL particle quality")
        .addValueCut("ecal_pid", 1.0, "PID == 1")
        .addRangeCut("ecal_beta", 0.8, 1.2, "0.8 < beta < 1.2")
        .addMinCut("cluster_energy", 100.0, "Cluster energy > 100 MeV");

    // 4-body pi+pi-e+e- selection chain (AND logic; order matches passCutSet values)
    cuts.defineCutSet("pippimepem_selection", "pi+pi-e+e- 4-body selection")
        .addMaxCut("oa_pippim_epem_lab",  50.0, "OA_LAB((pi+pi-),(e+e-)) < 50 deg")
        .addMaxCut("m_pippim",            0.420, "M(pi+pi-) < 0.420 GeV/c^2")
        .addMinCut("oa_pippim_epem_rest", 140.0, "OA((pi+pi-),(e+e-)) in pippimepem rest frame > 140 deg");

    // MM(pi+pi-e+e-) slice windows [GeV/c^2] — 6 adjacent bins for scan studies.
    // Used in tandem with the selection chain (see main.cc).
    cuts.defineRangeCut("mm_slice_20_22", 2.0, 2.2, "MM(pi+pi-e+e-) in [2.0, 2.2] GeV/c^2");
    cuts.defineRangeCut("mm_slice_22_24", 2.2, 2.4, "MM(pi+pi-e+e-) in [2.2, 2.4] GeV/c^2");
    cuts.defineRangeCut("mm_slice_24_26", 2.4, 2.6, "MM(pi+pi-e+e-) in [2.4, 2.6] GeV/c^2");
    cuts.defineRangeCut("mm_slice_26_28", 2.6, 2.8, "MM(pi+pi-e+e-) in [2.6, 2.8] GeV/c^2");
    cuts.defineRangeCut("mm_slice_28_30", 2.8, 3.0, "MM(pi+pi-e+e-) in [2.8, 3.0] GeV/c^2");

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

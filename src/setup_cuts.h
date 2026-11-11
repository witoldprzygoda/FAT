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
    // isBestAll: all stored combinations (isBest >= 0). The flag takes the values
    // 0 and 1 and MinCut is a strict ">", hence the threshold -0.5.
    cuts.defineMinCut("isBestAll", -0.5, "isBest >= 0 (all stored combinations)");
    // isBest == 1 (best combination only): kept defined, but NOT applied in main.cc.
    cuts.defineValueCut("isBest", 1, "Best candidate selection (defined, not applied)");
    cuts.defineMinCut("vertex_z", -500, "Vertex Z quality [mm]");

    // Trigger: PT3 only (trigbit == 8192). PT2 (trigbit == 4096) is defined
    // for the trigger-bias correction counters; it is NOT applied as a cut.
    cuts.defineValueCut("trigger_PT3", 8192, "Trigger PT3 (trigbit == 8192)");
    cuts.defineValueCut("trigger_PT2", 4096, "Trigger PT2 (trigbit == 4096, downscale 64)");

    // No start-detector (LGAD) requirement: the betanew productions use the
    // relative-time reconstruction (beta_new); start_iteration is only stored.

    // Opening angle window of the e+e- pair. Not applied as an event cut —
    // stored as the oa_pass flag in the ntuples.
    cuts.defineRangeCut("opening_angle_2_6", 2.0, 6.0, "2 deg <= OA(e+e-) <= 6 deg (oa_pass)");

    // The same opening angle cut from below and from above, as two standalone
    // cuts (flags pass_oa_min, pass_oa_max in scan_nt). Both are strict, while
    // the window above is inclusive at 2 and 6 deg.
    cuts.defineMinCut("oa_min", 2.0, "OA(e+e-) > 2 deg");
    cuts.defineMaxCut("oa_max", 6.0, "OA(e+e-) < 6 deg");

    // ECAL quality cuts (AND logic; passCutSet values must match this order)
    cuts.defineCutSet("ecal_quality", "ECAL particle quality")
        .addValueCut("ecal_pid", 1.0, "PID == 1")
        .addRangeCut("ecal_beta", 0.8, 1.2, "0.8 < beta < 1.2")
        .addMinCut("cluster_energy", 100.0, "Cluster energy > 100 MeV");

    // 4-body pi+pi-e+e- selection: three standalone cuts (the former cut set
    // "pippimepem_selection", same thresholds). main.cc evaluates each of them
    // on every row and stores the decisions as flags; sel_pass is their AND.
    // The same cut objects are used for RECONSTRUCTED and CORRECTED values
    // (no "_cor" twins), so their cut-flow counters mix REC and COR.
    cuts.defineMaxCut("oa_pippim_epem_lab", 50.0, "OA_LAB((pi+pi-),(e+e-)) < 50 deg");
    cuts.defineMaxCut("m_pippim", 0.420, "M(pi+pi-) < 0.420 GeV/c^2");
    cuts.defineMinCut("oa_pippim_epem_eta_rest", 140.0,
                      "OA((pi+pi-),(e+e-)) in m_eta-constrained rest frame > 140 deg");

    // MM(pi+pi-e+e-) slice windows [GeV/c^2] — 6 adjacent bins for scan studies.
    // Used in tandem with the selection chain (see main.cc).
    cuts.defineRangeCut("mm_slice_20_22", 2.0, 2.2, "MM(pi+pi-e+e-) in [2.0, 2.2] GeV/c^2");
    cuts.defineRangeCut("mm_slice_22_24", 2.2, 2.4, "MM(pi+pi-e+e-) in [2.2, 2.4] GeV/c^2");
    cuts.defineRangeCut("mm_slice_24_26", 2.4, 2.6, "MM(pi+pi-e+e-) in [2.4, 2.6] GeV/c^2");
    cuts.defineRangeCut("mm_slice_26_28", 2.6, 2.8, "MM(pi+pi-e+e-) in [2.6, 2.8] GeV/c^2");
    cuts.defineRangeCut("mm_slice_28_30", 2.8, 3.0, "MM(pi+pi-e+e-) in [2.8, 3.0] GeV/c^2");

    // pi0 invariant-mass windows for the e+e-gamma Dalitz hypothesis.
    // Two variants:
    //   pi0_mass_window         — wide [0.10, 0.18] GeV/c^2 (defined for fallback /
    //                             systematic studies; NOT applied to pippimepemg_pass).
    //   pi0_mass_window_narrow  — narrow [0.125, 0.145] GeV/c^2 (±10 MeV around m_pi0).
    //                             This is the ACTIVE cut driving pippimepemg_pass —
    //                             matches the experimental analysis convention for
    //                             selecting the eta -> pi+pi- pi0 (Dalitz) candidate.
    cuts.defineRangeCut("pi0_mass_window", 0.1, 0.18,
                        "M(e+e-gamma) in [0.10, 0.18] GeV/c^2 (pi0 Dalitz, wide — fallback)");
    cuts.defineRangeCut("pi0_mass_window_narrow", 0.125, 0.145,
                        "M(e+e-gamma) in [0.125, 0.145] GeV/c^2 (pi0 Dalitz, narrow — ACTIVE)");

    // eta invariant-mass window — used for f1(1285) -> pi+pi- eta search.
    // The eta candidate is identified either via eta -> e+e-gamma Dalitz (mult==1)
    // or eta -> gamma gamma direct (mult==2), and the cut is applied on the
    // corresponding mass — M(e+e-gamma) or M(gamma gamma).
    cuts.defineRangeCut("eta_mass_window", 0.5, 0.6,
                        "M(eta candidate) in [0.50, 0.60] GeV/c^2");

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

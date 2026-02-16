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

    // Graphical cut: anti-elastic (reject events inside the cut)
    cuts.loadGraphicalCut("anti_elastic", "CUT_dphi_tantan_fitted.root", "cutg",
                         "Anti-elastic cut (dphi vs tantan)");

    // FWD proton beta cut
    cuts.defineMinCut("proton_fwd_cut", 0.9, "FWD beta > 0.9");

    // Missing mass proton cuts
    cuts.defineMinCut("mm_p1", 1.05, "Missing mass proton1 > 1.05 GeV");
    cuts.defineMinCut("mm_p2", 1.05, "Missing mass proton2 > 1.05 GeV");

    // Pi0 mass window cut
    cuts.defineRangeCut("pion0_mass_cut", 0.0, 0.23, "Pi0 mass window [GeV/c^2]");

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

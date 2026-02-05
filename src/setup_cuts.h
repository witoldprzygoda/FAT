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

    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

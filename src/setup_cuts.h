/**
 * @file setup_cuts.h
 * @brief Cut definitions - Tutorial Starting Point
 *
 * This file defines all cuts for the analysis.
 * Currently empty - we will add cuts step by step.
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
 *
 * TUTORIAL: We will add cuts here step by step.
 */
inline void setupCuts(CutManager& cuts) {
    std::cout << "Setting up cuts...\n";

    // ========================================================================
    // STEP 0: No cuts yet
    // ========================================================================
    // We will add cuts as we develop the analysis.

    std::cout << "  (No cuts defined yet - Step 0)\n";
}

#endif // SETUP_CUTS_H

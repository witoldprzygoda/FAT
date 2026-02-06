/**
 * @file setup_ntuples.h
 * @brief Ntuple definitions - Tutorial Starting Point
 *
 * This file defines output ntuples for the analysis.
 * Currently empty - we will add ntuples step by step.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_NTUPLES_H
#define SETUP_NTUPLES_H

#include "manager.h"
#include "analysis_config.h"
#include <iostream>

/**
 * @brief Setup output ntuples for the analysis
 * @param manager Reference to the Manager
 * @param config Reference to AnalysisConfig
 *
 * TUTORIAL: We will add ntuples here step by step.
 */
inline void setupNtuples(Manager& manager, const AnalysisConfig& config) {
    std::cout << "Setting up ntuples...\n";

    // ========================================================================
    // STEP 6: Output ntuple with dilepton variables
    // ========================================================================
    // DynamicHNtuple allows adding variables at any time via operator[].
    // Variables are automatically discovered and stored.
    //
    // Naming convention:
    //   ep_*  - positron variables
    //   em_*  - electron variables
    //   oa    - opening angle
    //   m_ee  - dilepton invariant mass

    manager.createDynamicNtuple("dilepton_nt", "Dilepton event data");

    std::cout << "  Created output ntuple 'dilepton_nt'\n";
}

#endif // SETUP_NTUPLES_H

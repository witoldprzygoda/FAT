/**
 * @file setup_ntuples.h
 * @brief Output ntuple definitions
 *
 * This file defines output ntuples for the analysis.
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
 * To add a new ntuple: manager.createDynamicNtuple("name", "description")
 */
inline void setupNtuples(Manager& manager, const AnalysisConfig& config) {
    std::cout << "Setting up ntuples...\n";

    // ppip compound ntuple (proton + pi+)
    manager.createDynamicNtuple("ppip_nt", "p + pi+ compound data");
    std::cout << "  Created output ntuple 'ppip_nt'\n";
}

#endif // SETUP_NTUPLES_H

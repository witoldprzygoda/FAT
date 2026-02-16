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

    // Forward Tracker ntuple (fwd_mult, fwd variables indexed _1 to _3)
    manager.createDynamicNtuple("fwdet_nt", "Forward Tracker data");
    std::cout << "  Created output ntuple 'fwdet_nt'\n";

    // pp compound ntuple (HADES proton + FT proton, one entry per FWD candidate)
    manager.createDynamicNtuple("pp_nt", "p(HADES) + p(FT) compound data");
    std::cout << "  Created output ntuple 'pp_nt'\n";
}

#endif // SETUP_NTUPLES_H

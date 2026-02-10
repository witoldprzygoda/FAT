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

    // Dilepton ntuple (ep_*, em_*, oa, m_ee, CMS variables, cut flags)
    manager.createDynamicNtuple("dilepton_nt", "Dilepton event data");
    std::cout << "  Created output ntuple 'dilepton_nt'\n";

    // ECAL ntuple (ecal_mult, cluster/ecal variables indexed _1 to _5)
    manager.createDynamicNtuple("ecal_nt", "ECAL detector data");
    std::cout << "  Created output ntuple 'ecal_nt'\n";

    // Forward Tracker ntuple (fwd_mult, fwd variables indexed _1 to _3)
    manager.createDynamicNtuple("fwdet_nt", "Forward Tracker data");
    std::cout << "  Created output ntuple 'fwdet_nt'\n";

    // e+e-gamma compound ntuple (one entry per passing gamma, for pi0 Dalitz)
    manager.createDynamicNtuple("epemg_nt", "e+e-gamma compound data");
    std::cout << "  Created output ntuple 'epemg_nt'\n";
}

#endif // SETUP_NTUPLES_H

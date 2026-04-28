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

    // e+e-gamma compound ntuple (mult==1, for pi0 Dalitz)
    manager.createDynamicNtuple("epemg_nt", "e+e-gamma compound data");
    std::cout << "  Created output ntuple 'epemg_nt'\n";

    // e+e-gg compound ntuple (mult==2: epem + gg, for omega -> e+e- pi0 search)
    manager.createDynamicNtuple("epemgg_nt", "e+e-gg compound data (mult==2)");
    std::cout << "  Created output ntuple 'epemgg_nt'\n";

    // e+e-ggg compound ntuple (mult==3: 3 rotational entries per event;
    // each row stores m_epemg / m_gg / m_epemggg of one rotation plus cut flags).
    manager.createDynamicNtuple("epemggg_nt", "e+e-ggg combinatorial data (mult==3, 3 rotations)");
    std::cout << "  Created output ntuple 'epemggg_nt'\n";
}

#endif // SETUP_NTUPLES_H

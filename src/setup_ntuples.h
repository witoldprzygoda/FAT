/**
 * @file setup_ntuples.h
 * @brief Output ntuple definitions for the pi+pi- hadronic analysis
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_NTUPLES_H
#define SETUP_NTUPLES_H

#include "manager.h"
#include "analysis_config.h"
#include <iostream>

inline void setupNtuples(Manager& manager, const AnalysisConfig& config) {
    std::cout << "Setting up ntuples...\n";

    // Pi+Pi- ntuple (RECONSTRUCTED-derived compound observables)
    manager.createDynamicNtuple("pippim_nt", "Pi+Pi- event data (RECONSTRUCTED)");
    std::cout << "  Created output ntuple 'pippim_nt'\n";

    // Mirror ntuple: same field names but compound observables computed from
    // CORRECTED kinematics (energy-loss corrected pion momenta).
    manager.createDynamicNtuple("pippim_nt_cor", "Pi+Pi- event data (CORRECTED)");
    std::cout << "  Created output ntuple 'pippim_nt_cor'\n";

    // ECAL ntuple (ecal_mult, cluster/ecal variables indexed _1 to _5)
    manager.createDynamicNtuple("ecal_nt", "ECAL detector data");
    std::cout << "  Created output ntuple 'ecal_nt'\n";

    // Forward Tracker ntuple (fwd_mult, fwd variables indexed _1 to _3)
    manager.createDynamicNtuple("fwdet_nt", "Forward Tracker data");
    std::cout << "  Created output ntuple 'fwdet_nt'\n";
}

#endif // SETUP_NTUPLES_H

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
#include <string>
#include <vector>

/**
 * @brief Setup output ntuples for the analysis
 * @param manager Reference to the Manager
 * @param config Reference to AnalysisConfig
 *
 * To add a new ntuple: manager.createDynamicNtuple("name", "description")
 */
inline void setupNtuples(Manager& manager, const AnalysisConfig& config) {
    std::cout << "Setting up ntuples...\n";

    // Pi+Pi-e+e- ntuple (ep_*, em_*, pip_*, pim_*, oa, m_ee, CMS variables, cut flags)
    // Default uses RECONSTRUCTED-derived compound observables.
    manager.createDynamicNtuple("pippimepem_nt", "Pi+Pi-e+e- event data (RECONSTRUCTED)");
    std::cout << "  Created output ntuple 'pippimepem_nt'\n";

    // Mirror ntuple: same field names but compound observables computed from
    // CORRECTED kinematics (energy-loss corrected momenta of the constituents).
    manager.createDynamicNtuple("pippimepem_nt_cor", "Pi+Pi-e+e- event data (CORRECTED)");
    std::cout << "  Created output ntuple 'pippimepem_nt_cor'\n";

    // ECAL ntuple (ecal_mult, cluster/ecal variables indexed _1 to _5)
    manager.createDynamicNtuple("ecal_nt", "ECAL detector data");
    std::cout << "  Created output ntuple 'ecal_nt'\n";

    // Forward Tracker ntuple (fwd_mult, fwd variables indexed _1 to _3)
    manager.createDynamicNtuple("fwdet_nt", "Forward Tracker data");
    std::cout << "  Created output ntuple 'fwdet_nt'\n";

    // e+e-gamma compound ntuple (one entry per passing gamma, for pi0 Dalitz)
    manager.createDynamicNtuple("epemg_nt", "e+e-gamma compound data");
    std::cout << "  Created output ntuple 'epemg_nt'\n";

    // Cut-scan ntuple: one row per combination after isBestAll + vertex + trigger,
    // no other pre-selection. Keys, track quantities, cut variables as values
    // (REC and _cor), cut flags, three best ECAL photons, two forward candidates.
    manager.createDynamicNtuple("scan_nt", "Pi+Pi-e+e- cut-scan data (one row per combination)");
    std::cout << "  Created output ntuple 'scan_nt'\n";
}

/**
 * @brief Names of all ntuples created in setupNtuples()
 *
 * main() uses the list to check the written output file (rows in the file ==
 * number of fills) and to remove intermediate files after a failure.
 * Keep it in step with setupNtuples(): an ntuple missing here is not checked.
 */
inline std::vector<std::string> outputNtupleNames() {
    return {"pippimepem_nt", "pippimepem_nt_cor", "ecal_nt", "fwdet_nt", "epemg_nt", "scan_nt"};
}

#endif // SETUP_NTUPLES_H

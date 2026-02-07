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

    // ========================================================================
    // STEP 8: ECAL ntuple
    // ========================================================================
    // Separate ntuple for ECAL detector variables.
    // Variables are indexed 1, 2, 3 for each ECAL hit.
    //
    // Naming convention:
    //   ecal_mult       - number of ECAL hits (0-3)
    //   ecal_beta_N     - velocity (beta = v/c)
    //   ecal_pid_N      - particle ID
    //   ecal_energy_N   - energy from calorimeter [MeV]
    //   ecal_theta_N    - polar angle [deg]
    //   ecal_phi_N      - azimuthal angle [deg]
    //   ecal_chi2_N     - fit quality

    manager.createDynamicNtuple("ecal_nt", "ECAL detector data");

    std::cout << "  Created output ntuple 'ecal_nt'\n";
}

#endif // SETUP_NTUPLES_H

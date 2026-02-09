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
    // Variables are indexed 1, 2, 3, 4, 5 for each ECAL hit.
    //
    // Naming convention:
    //   ecal_mult         - number of ECAL hits (0-5)
    //   cluster_energy_N  - cluster energy [MeV] (used for kinematics)
    //   cluster_theta_N   - cluster polar angle [deg] (used for kinematics)
    //   cluster_phi_N     - cluster azimuthal angle [deg] (used for kinematics)
    //   ecal_beta_N       - velocity (beta = v/c)
    //   ecal_pid_N        - particle ID
    //   ecal_energy_N     - energy from primary reconstruction [MeV]
    //   ecal_theta_N      - polar angle (primary reco) [deg]
    //   ecal_phi_N        - azimuthal angle (primary reco) [deg]
    //   ecal_chi2_N       - fit quality

    manager.createDynamicNtuple("ecal_nt", "ECAL detector data");

    std::cout << "  Created output ntuple 'ecal_nt'\n";

    // ========================================================================
    // STEP 9: Forward Tracker ntuple
    // ========================================================================
    // Separate ntuple for Forward Tracker detector variables.
    // Variables are indexed 1, 2, 3 for each FT hit.
    //
    // Naming convention:
    //   fwd_mult       - number of FT hits (0-3)
    //   fwd_p_N        - momentum [MeV/c]
    //   fwd_theta_N    - polar angle [deg]
    //   fwd_phi_N      - azimuthal angle [deg]
    //   fwd_beta_N     - velocity (beta = v/c)
    //   fwd_mass_N     - reconstructed mass [MeV/c^2]
    //   fwd_mass2_N    - mass squared [MeV^2/c^4]
    //   fwd_chi2_N     - fit chi-squared
    //   fwd_ndf_N      - degrees of freedom
    //   fwd_chi2ndf_N  - chi2/ndf
    //   fwd_r_N        - radial position [mm]
    //   fwd_z_N        - Z position [mm]

    manager.createDynamicNtuple("fwdet_nt", "Forward Tracker data");

    std::cout << "  Created output ntuple 'fwdet_nt'\n";
}

#endif // SETUP_NTUPLES_H

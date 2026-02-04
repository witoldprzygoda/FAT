/**
 * @file setup_ntuples.h
 * @brief Ntuple definitions for the e+e- dilepton analysis with ECAL photons
 *
 * This file defines all output ntuples for the analysis.
 * Similar to setup_histograms.h, it keeps ntuple definitions
 * separate from the main analysis logic.
 *
 * Usage:
 *   1. Define your ntuples in setupNtuples()
 *   2. Access them via manager.getDynamicNtuple("name")
 *   3. Fill them in processEvent with: nt["var"] = value; nt.fill();
 *
 * Key features of DynamicHNtuple:
 *   - Add variables at ANY time (no prebooking needed)
 *   - Uses TTree internally for dynamic schema
 *   - Converts to flat TNtuple at finalization
 *   - Missing values filled with configurable sentinel (default: -1)
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_NTUPLES_H
#define SETUP_NTUPLES_H

#include "manager.h"
#include "analysis_config.h"

/**
 * @brief Setup all output ntuples for e+e- dilepton analysis with ECAL
 *
 * Define your analysis ntuples here. Each ntuple:
 * - Uses intermediate TTree (supports dynamic variable addition)
 * - Converts to TNtuple at closeFile() for easy plotting
 * - Variables are alphabetically ordered in final ntuple
 *
 * @param manager Reference to the Manager
 * @param config Reference to AnalysisConfig (for missing_value and keep_intermediate options)
 *
 * Ntuple contents:
 *
 * nt_particles:
 *   - ep_p, ep_theta, ep_phi, ep_mass      (positron)
 *   - em_p, em_theta, em_phi, em_mass      (electron)
 *   - ee_p, ee_theta, ee_phi, ee_mass      (dilepton)
 *   - oa_epem                               (opening angle)
 *   - n_gamma, n_gamma_good                 (photon multiplicity)
 *   - g1_energy, g1_theta, g1_phi, g1_good  (gamma1)
 *   - g2_energy, g2_theta, g2_phi, g2_good  (gamma2)
 *   - g3_energy, g3_theta, g3_phi, g3_good  (gamma3)
 *   - weight
 *
 * nt_compound:
 *   - m_ee, ee_pt, ee_rapidity              (dilepton)
 *   - cos_th_ee_cms, cos_th_ep_cms, cos_th_em_cms  (CMS angles)
 *   - oa_epem                               (opening angle)
 *   - m_eeg1, m_eeg2, m_eeg3                (e+e-gamma invariant masses)
 *   - m_gg, oa_gg, pi0_energy, pi0_theta    (pi0 reconstruction from gamma-gamma)
 *   - ep_helicity, ep_gj, em_helicity       (PWA variables)
 *   - weight
 */
inline void setupNtuples(Manager& manager, const AnalysisConfig& config) {

    // ========================================================================
    // Ntuple 1: Basic particle observables (e+, e-, e+e-, gammas)
    // ========================================================================
    // Contains: momenta, angles, masses of individual leptons, dilepton,
    // and ECAL photon candidates

    manager.createDynamicNtuple(
        "nt_particles",
        "Lepton and photon observables",
        config.getMissingValue(),
        config.getKeepIntermediateTree()
    );

    // ========================================================================
    // Ntuple 2: Compound observables (dilepton kinematics, e+e-gamma, PWA)
    // ========================================================================
    // Contains: dilepton mass, pt, rapidity, CMS angles, e+e-gamma masses,
    // PWA variables

    manager.createDynamicNtuple(
        "nt_compound",
        "Dilepton and e+e-gamma observables",
        config.getMissingValue(),
        config.getKeepIntermediateTree()
    );

    // ========================================================================
    // Additional ntuples (examples - uncomment to use)
    // ========================================================================

    // Control distributions ntuple
    // manager.createDynamicNtuple(
    //     "nt_control",
    //     "Control distributions",
    //     config.getMissingValue(),
    //     config.getKeepIntermediateTree()
    // );

    // Efficiency ntuple (for acceptance studies)
    // manager.createDynamicNtuple(
    //     "nt_efficiency",
    //     "Efficiency variables",
    //     config.getMissingValue(),
    //     config.getKeepIntermediateTree()
    // );
}

#endif // SETUP_NTUPLES_H

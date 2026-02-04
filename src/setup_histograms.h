/**
 * @file setup_histograms.h
 * @brief Histogram definitions for e+e- dilepton analysis with ECAL photons
 *
 * This file contains all histogram definitions.
 * Edit this file to add/modify histograms for your analysis.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_HISTOGRAMS_H
#define SETUP_HISTOGRAMS_H

#include "manager.h"
#include <iostream>

/**
 * @brief Define all histograms for the e+e- dilepton analysis with ECAL
 * @param mgr Manager object for histogram creation
 *
 * EDIT THIS FUNCTION to customize histograms for your analysis.
 */
inline void setupHistograms(Manager& mgr) {
    std::cout << "Setting up histogram system...\n";

    // ========================================================================
    // Quality Control Histograms
    // ========================================================================

    // Positron mass check (should be ~0.511 MeV)
    mgr.create1D("mass_ep", "e^{+} mass;M [GeV/c^{2}];Counts",
                 200, 0.0, 0.01, "quality");

    // Electron mass check (should be ~0.511 MeV)
    mgr.create1D("mass_em", "e^{-} mass;M [GeV/c^{2}];Counts",
                 200, 0.0, 0.01, "quality");

    // Dilepton invariant mass (before cuts)
    mgr.create1D("mass_ee", "e^{+}e^{-} invariant mass;M [GeV/c^{2}];Counts",
                 500, 0.0, 1.5, "quality");

    // Dilepton invariant mass (after cuts)
    mgr.create1D("mass_ee_cut", "e^{+}e^{-} invariant mass (after cuts);M [GeV/c^{2}];Counts",
                 500, 0.0, 1.5, "quality");

    // Event vertex
    mgr.create1D("eVertX", "Event Vertex X;X [mm];Counts", 200, -100, 100, "quality");
    mgr.create1D("eVertY", "Event Vertex Y;Y [mm];Counts", 200, -100, 100, "quality");
    mgr.create1D("eVertZ", "Event Vertex Z;Z [mm];Counts", 500, -400, 100, "quality");

    // ========================================================================
    // Dilepton Mass Distributions
    // ========================================================================

    // Main dilepton mass histogram (fine binning for spectroscopy)
    mgr.create1D("mass_dilepton", "Dilepton invariant mass;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 500, 0.0, 1.5, "dilepton");

    // ========================================================================
    // LAB Frame Kinematics
    // ========================================================================

    // Momenta
    mgr.create1D("ep_p_lab", "e^{+} momentum (LAB);p [MeV/c];Counts",
                 100, 0, 2000, "lab/momentum");
    mgr.create1D("em_p_lab", "e^{-} momentum (LAB);p [MeV/c];Counts",
                 100, 0, 2000, "lab/momentum");
    mgr.create1D("ee_p_lab", "e^{+}e^{-} momentum (LAB);p [MeV/c];Counts",
                 100, 0, 3000, "lab/momentum");

    // Angular distributions
    mgr.create1D("ep_theta_lab", "e^{+} #theta (LAB);#theta [deg];Counts",
                 90, 0, 90, "lab/angular");
    mgr.create1D("em_theta_lab", "e^{-} #theta (LAB);#theta [deg];Counts",
                 90, 0, 90, "lab/angular");
    mgr.create1D("ee_theta_lab", "e^{+}e^{-} #theta (LAB);#theta [deg];Counts",
                 180, 0, 180, "lab/angular");

    // ========================================================================
    // CMS Kinematics
    // ========================================================================

    // Angular distributions in CMS
    mgr.create1D("cos_theta_ee_cms", "e^{+}e^{-} cos#theta (CMS);cos#theta;Counts",
                 40, -1, 1, "cms/angular");
    mgr.create1D("cos_theta_ep_cms", "e^{+} cos#theta (CMS);cos#theta;Counts",
                 40, -1, 1, "cms/angular");
    mgr.create1D("cos_theta_em_cms", "e^{-} cos#theta (CMS);cos#theta;Counts",
                 40, -1, 1, "cms/angular");

    // CMS momenta
    mgr.create1D("ep_p_cms", "e^{+} momentum (CMS);p [MeV/c];Counts",
                 100, 0, 1500, "cms/momentum");
    mgr.create1D("em_p_cms", "e^{-} momentum (CMS);p [MeV/c];Counts",
                 100, 0, 1500, "cms/momentum");
    mgr.create1D("ee_p_cms", "e^{+}e^{-} momentum (CMS);p [MeV/c];Counts",
                 100, 0, 2000, "cms/momentum");

    // ========================================================================
    // Opening Angle
    // ========================================================================

    mgr.create1D("oa_epem", "Opening angle e^{+}-e^{-};#alpha [deg];Counts",
                 180, 0, 180, "opening_angles");

    // ========================================================================
    // Dilepton Kinematics
    // ========================================================================

    // Rapidity
    mgr.create1D("ee_rapidity", "e^{+}e^{-} rapidity;y;Counts",
                 100, -2, 2, "dilepton");

    // Transverse momentum
    mgr.create1D("ee_pt", "e^{+}e^{-} transverse momentum;p_{T} [MeV/c];Counts",
                 100, 0, 1500, "dilepton");

    // ========================================================================
    // ECAL Photon Histograms
    // ========================================================================

    // Photon multiplicity
    mgr.create1D("gamma_mult", "ECAL photon multiplicity;N_{#gamma};Counts",
                 10, 0, 10, "ecal");

    // Number of good photons (after quality cuts)
    mgr.create1D("gamma_mult_good", "Good photon multiplicity;N_{#gamma}^{good};Counts",
                 10, 0, 10, "ecal");

    // Photon energy (all candidates)
    mgr.create1D("gamma_energy", "Photon energy;E_{#gamma} [MeV];Counts",
                 200, 0, 2000, "ecal");

    // Photon energy (good candidates)
    mgr.create1D("gamma_energy_good", "Good photon energy;E_{#gamma} [MeV];Counts",
                 200, 0, 2000, "ecal");

    // Photon theta
    mgr.create1D("gamma_theta", "Photon #theta;#theta [deg];Counts",
                 90, 0, 90, "ecal");

    // Photon phi
    mgr.create1D("gamma_phi", "Photon #phi;#phi [deg];Counts",
                 360, -180, 180, "ecal");

    // Individual photon energies
    mgr.create1D("gamma1_energy", "#gamma_{1} energy;E [MeV];Counts",
                 200, 0, 2000, "ecal");
    mgr.create1D("gamma2_energy", "#gamma_{2} energy;E [MeV];Counts",
                 200, 0, 2000, "ecal");
    mgr.create1D("gamma3_energy", "#gamma_{3} energy;E [MeV];Counts",
                 200, 0, 2000, "ecal");

    // ========================================================================
    // e+e-γ Invariant Mass Histograms
    // ========================================================================

    // e+e- + gamma1
    mgr.create1D("mass_eeg1", "M(e^{+}e^{-}#gamma_{1});M [GeV/c^{2}];Counts",
                 500, 0.0, 1.5, "eegamma");

    // e+e- + gamma2
    mgr.create1D("mass_eeg2", "M(e^{+}e^{-}#gamma_{2});M [GeV/c^{2}];Counts",
                 500, 0.0, 1.5, "eegamma");

    // e+e- + gamma3
    mgr.create1D("mass_eeg3", "M(e^{+}e^{-}#gamma_{3});M [GeV/c^{2}];Counts",
                 500, 0.0, 1.5, "eegamma");

    // e+e- + any good gamma (combined)
    mgr.create1D("mass_eegamma", "M(e^{+}e^{-}#gamma);M [GeV/c^{2}];Counts",
                 500, 0.0, 1.5, "eegamma");

    // e+e- + gamma (zoomed for pi0 -> e+e-gamma Dalitz)
    mgr.create1D("mass_eegamma_pi0", "M(e^{+}e^{-}#gamma) #pi^{0} region;M [GeV/c^{2}];Counts",
                 200, 0.0, 0.3, "eegamma");

    // e+e- + gamma (zoomed for eta -> e+e-gamma Dalitz)
    mgr.create1D("mass_eegamma_eta", "M(e^{+}e^{-}#gamma) #eta region;M [GeV/c^{2}];Counts",
                 200, 0.3, 0.8, "eegamma");

    // ========================================================================
    // Pi0 Reconstruction (gamma-gamma)
    // ========================================================================

    // Two-gamma invariant mass (pi0 candidate)
    mgr.create1D("mass_gg", "M(#gamma#gamma);M [GeV/c^{2}];Counts",
                 200, 0.0, 0.5, "pi0");

    // Two-gamma invariant mass (zoomed on pi0 peak)
    mgr.create1D("mass_gg_pi0", "M(#gamma#gamma) #pi^{0} region;M [GeV/c^{2}];Counts",
                 100, 0.05, 0.25, "pi0");

    // Gamma-gamma opening angle
    mgr.create1D("oa_gg", "Opening angle #gamma-#gamma;#alpha [deg];Counts",
                 180, 0, 180, "pi0");

    // Pi0 energy
    mgr.create1D("pi0_energy", "#pi^{0} energy;E [MeV];Counts",
                 200, 0, 2000, "pi0");

    // Pi0 theta
    mgr.create1D("pi0_theta", "#pi^{0} #theta;#theta [deg];Counts",
                 90, 0, 90, "pi0");

    // ========================================================================
    // Gamma-Lepton Correlations
    // ========================================================================

    // Opening angle gamma - positron
    mgr.create1D("oa_g_ep", "Opening angle #gamma-e^{+};#alpha [deg];Counts",
                 180, 0, 180, "gamma_lepton");

    // Opening angle gamma - electron
    mgr.create1D("oa_g_em", "Opening angle #gamma-e^{-};#alpha [deg];Counts",
                 180, 0, 180, "gamma_lepton");

    // Opening angle gamma - dilepton
    mgr.create1D("oa_g_ee", "Opening angle #gamma-e^{+}e^{-};#alpha [deg];Counts",
                 180, 0, 180, "gamma_lepton");

    // Minimum opening angle gamma to any lepton
    mgr.create1D("oa_g_lepton_min", "Min opening angle #gamma-lepton;#alpha [deg];Counts",
                 180, 0, 180, "gamma_lepton");

    // ========================================================================
    // ECAL Quality Histograms
    // ========================================================================

    // Gamma beta (velocity)
    mgr.create1D("gamma_beta", "Photon #beta;#beta;Counts",
                 100, 0.5, 1.5, "ecal_quality");

    // Gamma chi2
    mgr.create1D("gamma_chi2", "Photon #chi^{2};#chi^{2};Counts",
                 100, 0, 10, "ecal_quality");

    // Gamma distance to EMC
    mgr.create1D("gamma_dist", "Photon distance to EMC;d [mm];Counts",
                 100, 0, 500, "ecal_quality");

    // Gamma cluster position
    mgr.create2D("gamma_rz", "Photon cluster position;r [mm];z [mm]",
                 100, 0, 600, 200, 200, 800, "ecal_quality");

    // ========================================================================
    // 2D Correlations
    // ========================================================================

    mgr.create2D("mass_vs_costh_ee", "M_{e^{+}e^{-}} vs cos#theta (CMS);M [GeV/c^{2}];cos#theta",
                 100, 0.0, 1.0, 40, -1, 1, "correlations");

    mgr.create2D("theta_ep_vs_em_lab", "#theta_{e^{-}} vs #theta_{e^{+}} (LAB);#theta_{e^{-}} [deg];#theta_{e^{+}} [deg]",
                 90, 0, 90, 90, 0, 90, "correlations");

    mgr.create2D("p_ep_vs_em_lab", "p_{e^{-}} vs p_{e^{+}} (LAB);p_{e^{-}} [MeV/c];p_{e^{+}} [MeV/c]",
                 100, 0, 2000, 100, 0, 2000, "correlations");

    // e+e-gamma correlations
    mgr.create2D("mass_ee_vs_eegamma", "M(e^{+}e^{-}) vs M(e^{+}e^{-}#gamma);M_{ee} [GeV/c^{2}];M_{ee#gamma} [GeV/c^{2}]",
                 100, 0.0, 0.5, 100, 0.0, 1.0, "correlations");

    mgr.create2D("mass_eegamma_vs_gamma_energy", "M(e^{+}e^{-}#gamma) vs E_{#gamma};M [GeV/c^{2}];E_{#gamma} [MeV]",
                 100, 0.0, 1.0, 100, 0, 1000, "correlations");

    // ========================================================================
    // PWA Variables (Partial Wave Analysis)
    // ========================================================================

    // Helicity angles (in dilepton rest frame)
    mgr.create1D("pwa_ep_helicity_ee", "e^{+} helicity in e^{+}e^{-} frame;cos#theta_{H};Counts",
                 40, -1, 1, "pwa/helicity");
    mgr.create1D("pwa_em_helicity_ee", "e^{-} helicity in e^{+}e^{-} frame;cos#theta_{H};Counts",
                 40, -1, 1, "pwa/helicity");

    // Gottfried-Jackson angles
    mgr.create1D("pwa_ep_gj_ee", "e^{+} GJ angle in e^{+}e^{-} frame;cos#theta_{GJ};Counts",
                 40, -1, 1, "pwa/gottfried_jackson");

    std::cout << "✓ Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

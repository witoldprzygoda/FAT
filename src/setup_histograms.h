/**
 * @file setup_histograms.h
 * @brief Histogram definitions for analysis
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
 * @brief Define all histograms for the analysis
 * @param mgr Manager object for histogram creation
 *
 * EDIT THIS FUNCTION to customize histograms for your analysis.
 */
inline void setupHistograms(Manager& mgr) {
    std::cout << "Setting up histogram system...\n";
    
    // ========================================================================
    // Quality Control Histograms
    // ========================================================================
    
    // Missing mass (neutron)
    mgr.create1D("mass_n", "Missing mass (neutron);M [GeV/c^{2}];Counts", 
                 200, 0.5, 1.5, "quality");
    mgr.create1D("mass_n_cut", "Missing mass (after cuts);M [GeV/c^{2}];Counts", 
                 200, 0.5, 1.5, "quality");
    
    // Proton mass check
    mgr.create1D("mass_p", "Proton mass;M [GeV/c^{2}];Counts", 
                 200, 0.5, 1.5, "quality");
    
    // Pion mass check
    mgr.create1D("mass_pip", "Pion mass;M [GeV/c^{2}];Counts", 
                 200, 0.0, 0.5, "quality");
    
    // Event vertex
    mgr.create1D("eVertX", "Event Vertex X;X [mm];Counts", 200, -100, 100, "quality");
    mgr.create1D("eVertY", "Event Vertex Y;Y [mm];Counts", 200, -100, 100, "quality");
    mgr.create1D("eVertZ", "Event Vertex Z;Z [mm];Counts", 500, -400, 100, "quality");
    
    // ========================================================================
    // Composite Particle Masses
    // ========================================================================
    
    mgr.create1D("mass_deltaPP", "#Delta^{++} mass;M [GeV/c^{2}];Counts", 
                 100, 0.8, 1.8, "composites");
    mgr.create1D("mass_deltaP", "#Delta^{+} mass;M [GeV/c^{2}];Counts", 
                 100, 0.8, 1.8, "composites");
    mgr.create1D("mass_ppip", "p#pi^{+} invariant mass;M [GeV/c^{2}];Counts", 
                 100, 1.0, 2.0, "composites");
    mgr.create1D("mass_npip", "n#pi^{+} invariant mass;M [GeV/c^{2}];Counts", 
                 100, 1.0, 2.0, "composites");
    mgr.create1D("mass_pn", "pn invariant mass;M [GeV/c^{2}];Counts", 
                 100, 1.7, 2.5, "composites");
    
    // ========================================================================
    // LAB Frame Kinematics
    // ========================================================================
    
    // Momenta
    mgr.create1D("p_p_lab", "Proton momentum (LAB);p [MeV/c];Counts", 
                 100, 0, 3000, "lab/momentum");
    mgr.create1D("pip_p_lab", "#pi^{+} momentum (LAB);p [MeV/c];Counts", 
                 100, 0, 2000, "lab/momentum");
    mgr.create1D("n_p_lab", "Neutron momentum (LAB);p [MeV/c];Counts", 
                 100, 0, 2000, "lab/momentum");
    
    // Angular distributions
    mgr.create1D("p_theta_lab", "Proton #theta (LAB);#theta [deg];Counts", 
                 90, 0, 90, "lab/angular");
    mgr.create1D("pip_theta_lab", "#pi^{+} #theta (LAB);#theta [deg];Counts", 
                 90, 0, 90, "lab/angular");
    mgr.create1D("n_theta_lab", "Neutron #theta (LAB);#theta [deg];Counts", 
                 180, 0, 180, "lab/angular");
    
    // ========================================================================
    // CMS Kinematics
    // ========================================================================
    
    mgr.create1D("cos_theta_deltaPP_cms", "#Delta^{++} cos#theta (CMS);cos#theta;Counts", 
                 40, -1, 1, "cms/angular");
    mgr.create1D("cos_theta_deltaP_cms", "#Delta^{+} cos#theta (CMS);cos#theta;Counts", 
                 40, -1, 1, "cms/angular");
    mgr.create1D("cos_theta_p_cms", "Proton cos#theta (CMS);cos#theta;Counts", 
                 40, -1, 1, "cms/angular");
    mgr.create1D("cos_theta_pip_cms", "#pi^{+} cos#theta (CMS);cos#theta;Counts", 
                 40, -1, 1, "cms/angular");
    mgr.create1D("cos_theta_n_cms", "Neutron cos#theta (CMS);cos#theta;Counts", 
                 40, -1, 1, "cms/angular");
    
    // CMS momenta
    mgr.create1D("p_p_cms", "Proton momentum (CMS);p [MeV/c];Counts", 
                 100, 0, 1500, "cms/momentum");
    mgr.create1D("pip_p_cms", "#pi^{+} momentum (CMS);p [MeV/c];Counts", 
                 100, 0, 1000, "cms/momentum");
    mgr.create1D("n_p_cms", "Neutron momentum (CMS);p [MeV/c];Counts", 
                 100, 0, 1500, "cms/momentum");
    
    // ========================================================================
    // Opening Angles
    // ========================================================================
    
    mgr.create1D("oa_ppip", "Opening angle p-#pi^{+};#alpha [deg];Counts", 
                 180, 0, 180, "opening_angles");
    mgr.create1D("oa_npip", "Opening angle n-#pi^{+};#alpha [deg];Counts", 
                 180, 0, 180, "opening_angles");
    mgr.create1D("oa_pn", "Opening angle p-n;#alpha [deg];Counts", 
                 180, 0, 180, "opening_angles");
    
    // ========================================================================
    // 2D Correlations
    // ========================================================================
    
    mgr.create2D("dalitz_ppip_npip", "Dalitz plot;M^{2}(p#pi^{+}) [GeV^{2}/c^{4}];M^{2}(n#pi^{+}) [GeV^{2}/c^{4}]",
                 100, 1.0, 4.0, 100, 1.0, 4.0, "correlations");
    
    mgr.create2D("mass_vs_costh_deltaPP", "#Delta^{++}: M vs cos#theta;M [GeV/c^{2}];cos#theta",
                 50, 1.0, 1.6, 40, -1, 1, "correlations");
    
    mgr.create2D("theta_p_vs_pip_lab", "#theta_{p} vs #theta_{#pi^{+}} (LAB);#theta_{#pi^{+}} [deg];#theta_{p} [deg]",
                 90, 0, 90, 90, 0, 90, "correlations");
    
    // ========================================================================
    // PWA Variables (Partial Wave Analysis)
    // ========================================================================
    
    // Helicity angles
    mgr.create1D("pwa_pip_helicity_ppip", "#pi^{+} helicity in p#pi^{+} frame;cos#theta_{H};Counts",
                 40, -1, 1, "pwa/helicity");
    mgr.create1D("pwa_n_helicity_ppip", "n helicity in p#pi^{+} frame;cos#theta_{H};Counts",
                 40, -1, 1, "pwa/helicity");
    
    // Gottfried-Jackson angles  
    mgr.create1D("pwa_pip_gj_ppip", "#pi^{+} GJ angle in p#pi^{+} frame;cos#theta_{GJ};Counts",
                 40, -1, 1, "pwa/gottfried_jackson");
    
    std::cout << "✓ Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

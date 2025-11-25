/**
 * @file setup_cuts.h
 * @brief Cut definitions for analysis
 *
 * This file contains all cut definitions.
 * Edit this file to add/modify cuts for your analysis.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_CUTS_H
#define SETUP_CUTS_H

#include "cut_manager.h"
#include <iostream>

/**
 * @brief Define all cuts for the analysis
 * @param cuts CutManager object for cut definition
 *
 * EDIT THIS FUNCTION to customize cuts for your analysis.
 */
inline void setupCuts(CutManager& cuts) {
    std::cout << "Setting up cuts...\n";
    
    // ========================================================================
    // Range Cuts (min <= value <= max)
    // ========================================================================
    
    // Missing mass (neutron) cut
    cuts.defineRangeCut("neutron_mass", 0.899, 0.986, "Neutron mass window [GeV]");
    
    // Delta++ mass cut
    cuts.defineRangeCut("deltaPP_mass", 0.8, 1.8, "Delta++ mass window [GeV]");
    
    // Example: Momentum cuts (uncomment if needed)
    // cuts.defineRangeCut("proton_momentum", 100.0, 3000.0, "Proton p range [MeV/c]");
    // cuts.defineRangeCut("pion_momentum", 50.0, 2000.0, "Pion p range [MeV/c]");
    
    // Example: Angular cuts (uncomment if needed)
    // cuts.defineRangeCut("proton_theta", 5.0, 85.0, "Proton theta range [deg]");
    
    // ========================================================================
    // Trigger Cuts (uncomment if needed)
    // ========================================================================
    
    // cuts.defineTriggerCut("physics", 4, false, "PT3 trigger");
    // cuts.defineTriggerCut("calibration", 0x0C, true, "Calibration triggers");
    
    // ========================================================================
    // Graphical Cuts (uncomment and edit path if needed)
    // ========================================================================
    
    // try {
    //     cuts.loadGraphicalCut("proton_pid", "cuts/proton_pid.root", "cutg_proton");
    // } catch (const std::exception& e) {
    //     std::cerr << "Warning: Could not load proton PID cut: " << e.what() << "\n";
    // }
    
    cuts.printDefinedCuts();
}

#endif // SETUP_CUTS_H

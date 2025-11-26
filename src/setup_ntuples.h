/**
 * @file setup_ntuples.h
 * @brief Ntuple definitions for the analysis
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
 * @brief Setup all output ntuples
 *
 * Define your analysis ntuples here. Each ntuple:
 * - Uses intermediate TTree (supports dynamic variable addition)
 * - Converts to TNtuple at closeFile() for easy plotting
 * - Variables are alphabetically ordered in final ntuple
 *
 * @param manager Reference to the Manager
 * @param config Reference to AnalysisConfig (for missing_value and keep_intermediate options)
 *
 * Example adding a new ntuple:
 * @code
 *   manager.createDynamicNtuple(
 *       "nt_systematics",           // name
 *       "Systematic checks",        // title
 *       config.getMissingValue(),   // sentinel for missing values
 *       config.getKeepIntermediateTree()  // keep TTree file?
 *   );
 * @endcode
 */
inline void setupNtuples(Manager& manager, const AnalysisConfig& config) {
    
    // ========================================================================
    // Ntuple 1: Basic particle observables (p, π+, n)
    // ========================================================================
    // Contains: momenta, angles, masses of individual particles
    
    manager.createDynamicNtuple(
        "nt_particles",
        "Basic particle observables",
        config.getMissingValue(),
        config.getKeepIntermediateTree()
    );
    
    // ========================================================================
    // Ntuple 2: Compound observables (Δ++, Δ+, pπ+, etc.)
    // ========================================================================
    // Contains: composite masses, CMS angles, opening angles, PWA variables
    
    manager.createDynamicNtuple(
        "nt_compound",
        "Compound particle observables",
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
}

#endif // SETUP_NTUPLES_H

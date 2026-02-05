// ========================================================================
// FAT Framework - Tutorial
// ========================================================================
// Step-by-step analysis development.
//
// Step 0: Empty framework - just reads events, does nothing
// Step 1: Create e+ and e- particles from ntuple variables
//
// Usage:
//   ./ana [config.json]
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2025
// ========================================================================

#include "src/manager.h"
#include "src/pparticle.h"
#include "src/ntuple_reader.h"
#include "src/cut_manager.h"
#include "src/analysis_config.h"
#include "src/setup_histograms.h"
#include "src/setup_ntuples.h"
#include "src/setup_cuts.h"
#include "src/progressbar.h"
#include "src/console_box.h"
#include <iostream>

// Use Physics namespace for mass constants
using namespace Physics;

// ============================================================================
// PROCESS SINGLE EVENT
// ============================================================================
// This function is called for each event in the input data.
// Currently empty - we will build it step by step.
// ============================================================================

void processEvent(NTupleReader& reader, Manager& mgr, CutManager& cuts,
                 const AnalysisConfig& config) {

    // ========================================================================
    // STEP 1: Create e+ and e- particles
    // ========================================================================
    // PParticle is a class that represents a particle with 4-momentum.
    // We create "shells" with the particle mass, then fill kinematics.
    //
    // The EpEm_ID ntuple contains:
    //   ep_p, ep_theta, ep_phi  - positron momentum (MeV/c) and angles (deg)
    //   em_p, em_theta, em_phi  - electron momentum (MeV/c) and angles (deg)

    // Create particle shells (just mass, no kinematics yet)
    PParticle positron(MASS_ELECTRON, "e+");   // e+ has same mass as e-
    PParticle electron(MASS_ELECTRON, "e-");

    // Fill kinematics from ntuple variables
    // setFromSpherical(momentum, theta, phi, kinematic_type)
    positron.setFromSpherical(reader["ep_p"], reader["ep_theta"], reader["ep_phi"],
                              KinematicType::RECONSTRUCTED);

    electron.setFromSpherical(reader["em_p"], reader["em_theta"], reader["em_phi"],
                              KinematicType::RECONSTRUCTED);

    // Now positron and electron have their 4-momenta set.
    // We can access their properties:
    //   positron.momentum()  - |p| in MeV/c
    //   positron.theta()     - polar angle in degrees
    //   positron.phi()       - azimuthal angle in degrees
    //   positron.energy()    - E in MeV
    //   positron.mass()      - invariant mass in MeV/c²

}

// ============================================================================
// MAIN PROGRAM
// ============================================================================

int main(int argc, char* argv[]) {
    // Install signal handler for graceful Ctrl+C termination
    SignalHandler::install();

    ConsoleBox::newLine();
    ConsoleBox::printHeader("FAT Framework - Tutorial",
                           "Step-by-step Analysis Development");
    ConsoleBox::newLine();

    // ========================================================================
    // 1. LOAD CONFIGURATION
    // ========================================================================

    std::string config_file = "config.json";
    if (argc > 1) {
        config_file = argv[1];
    }

    AnalysisConfig config;
    try {
        config.load(config_file);
    } catch (const std::exception& e) {
        std::cerr << "Error loading configuration: " << e.what() << "\n";
        return 1;
    }

    config.print();

    // ========================================================================
    // 2. OPEN INPUT DATA
    // ========================================================================

    NTupleReader reader;

    try {
        std::string tree_name = config.getInputTreeName();

        if (config.isInputFileList()) {
            reader.openFromList(config.getInputSource(), tree_name);
        } else if (config.isInputMultipleRootFiles()) {
            std::vector<std::string> files = config.getInputFiles();
            std::cout << "\nInput: " << files.size() << " ROOT files (chain)\n";
            reader.openChain(files, tree_name);
        } else if (config.isInputRootFile()) {
            reader.open(config.getInputSource(), tree_name);
        } else {
            throw std::runtime_error("Unknown input format");
        }
    } catch (const std::exception& e) {
        std::cerr << "Error opening input: " << e.what() << "\n";
        return 1;
    }

    // ========================================================================
    // 3. OPEN OUTPUT FILE & SETUP
    // ========================================================================

    Manager manager;
    manager.openFile(config.getOutputFilename(), config.getOutputOption());

    // Setup histograms (defined in src/setup_histograms.h)
    setupHistograms(manager);

    // Setup ntuples (defined in src/setup_ntuples.h)
    setupNtuples(manager, config);

    // Setup cuts (defined in src/setup_cuts.h)
    CutManager cuts;
    setupCuts(cuts);

    // ========================================================================
    // 4. EVENT LOOP
    // ========================================================================

    Long64_t total_entries = reader.entries();
    Long64_t start_event = config.getStartEvent();
    Long64_t max_events = config.getMaxEvents();

    Long64_t end_event = total_entries;
    if (max_events > 0) {
        end_event = std::min(start_event + max_events, total_entries);
    }

    Long64_t events_to_process = end_event - start_event;

    ConsoleBox::newLine();
    ConsoleBox::printInfoBox("Press Ctrl+C at any time to stop and save partial results");
    ConsoleBox::newLine();
    std::cout << "Processing events " << start_event << " to " << end_event
              << " (" << events_to_process << " events)...\n\n";

    Long64_t processed = 0;
    bool was_interrupted = false;

    ProgressBar progress(events_to_process);

    for (Long64_t i = start_event; i < end_event; ++i) {
        if (SignalHandler::wasInterrupted()) {
            was_interrupted = true;
            break;
        }

        reader.getEntry(i);
        ++processed;
        progress.update(processed);

        // Process event
        try {
            processEvent(reader, manager, cuts, config);
        } catch (const std::exception& e) {
            continue;
        }
    }

    progress.finish(was_interrupted);

    // ========================================================================
    // 5. FINALIZATION
    // ========================================================================

    std::cout << "\n";
    if (was_interrupted) {
        std::cout << "Processing interrupted by user (Ctrl+C).\n";
    } else {
        std::cout << "Processing complete!\n";
    }
    std::cout << "  Events processed: " << processed << "\n";

    cuts.printCutFlow();

    std::cout << "\nSaving results to " << config.getOutputFilename() << "...\n";
    manager.printSummary();
    manager.closeFile();

    ConsoleBox::newLine();
    ConsoleBox::printStatus("Analysis Complete!");
    ConsoleBox::newLine();

    return 0;
}

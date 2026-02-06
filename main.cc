// ========================================================================
// FAT Framework - Tutorial
// ========================================================================
// Step-by-step analysis development.
//
// Step 0: Empty framework - just reads events, does nothing
// Step 1: Create e+ and e- particles from ntuple variables
// Step 1b: Add CORRECTED kinematics (energy-loss corrected momentum)
// Step 2: Add first histograms (lepton momentum)
// Step 2b: Add 2D histograms (momentum correction vs momentum)
// Step 3: Compound object - dilepton (e+ + e-)
// Step 4a: Event-level cuts (isBest, vertex_z)
// Step 4b: CutSet for particle-level cuts (implemented in cut_manager.h)
// Step 5a: Opening angle cut
//
// Usage:
//   ./ana [config.json]
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2025
// ========================================================================

#include "src/manager.h"
#include "src/pparticle.h"
#include "src/physics_utils.h"
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
    // STEP 4a: Event-level cuts
    // ========================================================================
    // These cuts are applied FIRST, BEFORE creating any particles.
    // If an event fails these cuts, we skip it entirely (return early).
    // This is efficient because we avoid creating unnecessary objects.

    if (!cuts.passValueCut("isBest", reader["isBest"])) return;
    if (!cuts.passMinCut("vertex_z", reader["eVertReco_z"])) return;

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

    // ========================================================================
    // STEP 1b: Add CORRECTED kinematics
    // ========================================================================
    // Energy-loss corrected momentum - same angles, different |p|
    // The correction accounts for energy lost in detector material.
    //
    // PParticle can store MULTIPLE kinematic representations simultaneously!
    // This allows comparing reconstructed vs corrected in the same analysis.

    positron.setFromSpherical(reader["ep_p_corr_ep"], reader["ep_theta"], reader["ep_phi"],
                              KinematicType::CORRECTED);

    electron.setFromSpherical(reader["em_p_corr_em"], reader["em_theta"], reader["em_phi"],
                              KinematicType::CORRECTED);

    // Now each particle has TWO momentum representations:
    //   positron.momentum(KinematicType::RECONSTRUCTED)  - raw measurement
    //   positron.momentum(KinematicType::CORRECTED)      - energy-loss corrected
    //
    // Default (no argument) uses RECONSTRUCTED:
    //   positron.momentum() == positron.momentum(KinematicType::RECONSTRUCTED)

    // ========================================================================
    // STEP 2: Fill histograms
    // ========================================================================
    // mgr.fill(histogram_name, value) - fill histogram with value
    //
    // We use positron.momentum() which returns |p| in MeV/c

    mgr.fill("ep_p", positron.momentum());   // positron momentum
    mgr.fill("em_p", electron.momentum());   // electron momentum

    // ========================================================================
    // STEP 2b: Fill 2D histograms - momentum correction
    // ========================================================================
    // mgr.fill(histogram_name, x_value, y_value) - fill 2D histogram
    //
    // Calculate: delta_p = p_corrected - p_reconstructed

    double ep_p_rec = positron.momentum(KinematicType::RECONSTRUCTED);
    double ep_p_cor = positron.momentum(KinematicType::CORRECTED);
    double ep_dp = ep_p_cor - ep_p_rec;

    double em_p_rec = electron.momentum(KinematicType::RECONSTRUCTED);
    double em_p_cor = electron.momentum(KinematicType::CORRECTED);
    double em_dp = em_p_cor - em_p_rec;

    mgr.fill("ep_dp_vs_p", ep_p_rec, ep_dp);   // positron correction vs p
    mgr.fill("em_dp_vs_p", em_p_rec, em_dp);   // electron correction vs p

    // ========================================================================
    // STEP 5a: Opening angle cut
    // ========================================================================
    // The opening angle is the angle between e+ and e- momentum vectors.
    // Must be calculated BEFORE combining into dilepton (operator+ loses info).
    // Physics::openingAngle() is a symmetric function from physics_utils.h.

    double oa = Physics::openingAngle(positron, electron);
    mgr.fill("opening_angle", oa);

    // ========================================================================
    // STEP 5b: Before/after histograms
    // ========================================================================
    // Calculate mass BEFORE the cut to fill "before" histogram.
    // Then apply cut, and fill "after" histogram only for passing events.
    // Comparing these shows what the opening angle cut removes.

    PParticle dilepton = positron + electron;
    double m_ee = dilepton.massGeV();  // uses RECONSTRUCTED by default

    mgr.fill("mass_ee_before_oa", m_ee);  // BEFORE opening angle cut

    // Apply cut: reject close pairs (e.g., from conversions)
    if (!cuts.passMinCut("opening_angle", oa)) return;

    mgr.fill("mass_ee_after_oa", m_ee);   // AFTER opening angle cut

    // ========================================================================
    // STEP 3: Compound object - dilepton
    // ========================================================================
    // PParticle supports operator+ to create composite particles.
    // The result is a new PParticle with combined 4-momentum.
    //
    // dilepton = e+ + e-  (virtual photon)
    //
    // The invariant mass M = sqrt(E² - p²) gives the dilepton mass.
    //
    // (dilepton already created above for Step 5b)

    // The dilepton inherits BOTH kinematic types from its parents:
    //   dilepton.mass(KinematicType::RECONSTRUCTED)
    //   dilepton.mass(KinematicType::CORRECTED)

    // Get invariant mass in GeV/c² (massGeV divides by 1000)
    // (m_ee already calculated above for Step 5b)

    mgr.fill("mass_ee", m_ee);

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

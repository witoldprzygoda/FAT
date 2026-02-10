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
// Step 5b: Before/after histograms
// Step 6: Output ntuple with dilepton variables
// Step 7: CMS boost - transform dilepton to center of mass frame
// Step 8: ECAL objects - create particle objects from ECAL detector
// Step 8b: ECAL quality cuts (pid==1, 0.8<beta<1.2, cluster_energy>100)
// Step 9: Forward Tracker objects - create particle objects from FT detector
// Step 10: Combinatorial background - variable aliases for like-sign pairs
// Step 11: Compound objects e+e-gamma for pi0 Dalitz decay identification
//
// Usage:
//   ./ana [config.json]
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2025
// ========================================================================

#include "src/manager.h"
#include "src/pparticle.h"
#include "src/pparticle_ecal.h"
#include "src/pparticle_fwd.h"
#include "src/physics_utils.h"
#include "src/boost_frame.h"
#include "src/ntuple_reader.h"
#include "src/cut_manager.h"
#include "src/analysis_config.h"
#include "src/setup_histograms.h"
#include "src/setup_ntuples.h"
#include "src/setup_cuts.h"
#include "src/progressbar.h"
#include "src/console_box.h"
#include <iostream>
#include <vector>

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

    // ========================================================================
    // STEP 7: CMS boost - transform dilepton to center of mass frame
    // ========================================================================
    // The beam-target CMS frame is where total momentum is zero.
    // This is the natural frame for studying the reaction dynamics.
    //
    // EventFrames uses the beam kinetic energy from config to create
    // the proper boost vector. Then we boost a COPY of the dilepton
    // (original stays in LAB frame for other uses).
    //
    // Key CMS observables:
    //   y_cms    - rapidity (centered around 0 in symmetric collisions)
    //   pt       - transverse momentum (same in LAB and CMS)
    //   theta_cms - polar angle in CMS (different from LAB theta)

    EventFrames frames;
    frames.setBeamFrameFromKineticEnergy(config.getBeamKineticEnergy());

    // Boost dilepton to CMS (creates a COPY, original unchanged)
    PParticle dilepton_cms = frames.getFrame("beam").boost(dilepton);

    // Extract CMS quantities
    double y_cms = dilepton_cms.rapidity();
    double pt = dilepton_cms.vec().Pt();       // transverse momentum
    double theta_cms = dilepton_cms.theta();   // polar angle in CMS

    // ========================================================================
    // STEP 6: Fill output ntuple (before applying OA cut)
    // ========================================================================
    // DynamicHNtuple uses operator[] to set variables.
    // Variables are automatically discovered on first use.
    // Fill ntuple BEFORE the OA cut to store all events with cut decisions.

    // Check OA cut but store decision as flag (don't apply yet)
    bool oa_pass = cuts.passMinCut("opening_angle", oa);

    auto& nt = mgr.getDynamicNtuple("dilepton_nt");

    // Positron variables
    nt["ep_p_rec"] = ep_p_rec;
    nt["ep_p_cor"] = ep_p_cor;
    nt["ep_theta"] = positron.theta();
    nt["ep_phi"] = positron.phi();
    nt["ep_theta_rich"] = reader["ep_theta_rich"];
    nt["ep_phi_rich"] = reader["ep_phi_rich"];

    // Electron variables
    nt["em_p_rec"] = em_p_rec;
    nt["em_p_cor"] = em_p_cor;
    nt["em_theta"] = electron.theta();
    nt["em_phi"] = electron.phi();
    nt["em_theta_rich"] = reader["em_theta_rich"];
    nt["em_phi_rich"] = reader["em_phi_rich"];

    // Dilepton variables (LAB frame)
    nt["oa"] = oa;
    nt["m_ee"] = m_ee;

    // CMS variables (Step 7)
    nt["y_cms"] = y_cms;
    nt["pt"] = pt;
    nt["theta_cms"] = theta_cms;

    // Cut decisions (0 = fail, 1 = pass)
    nt["oa_pass"] = oa_pass ? 1.0f : 0.0f;

    nt.fill();

    // Apply cut: reject close pairs (e.g., from conversions)
    if (!oa_pass) return;

    mgr.fill("mass_ee_after_oa", m_ee);   // AFTER opening angle cut

    // ========================================================================
    // STEP 7: Fill CMS histograms (after OA cut)
    // ========================================================================
    // These histograms show the physics observables for good dilepton pairs.

    mgr.fill("rapidity_cms", y_cms);
    mgr.fill("pt_cms", pt);
    mgr.fill("theta_cms", theta_cms);
    mgr.fill("rapidity_vs_mass", m_ee, y_cms);  // 2D: mass vs rapidity

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

    // ========================================================================
    // STEP 8: ECAL objects
    // ========================================================================
    // Create particle objects from ECAL (electromagnetic calorimeter).
    // ECAL can detect various particles: photons, neutrons, pions, electrons, etc.
    // neutr_mult tells us how many ECAL hits are present (0-5).
    //
    // PParticleEcal extends PParticle with ECAL-specific data:
    //   - Works with operator+ for composite particle creation
    //   - Stores all ECAL detector variables (energy, chi2, tof, etc.)
    //   - Uses cluster reconstruction (neutr_cluster_*) for kinematics
    //
    // Controlled by config.json: "ecal": true/false

    if (config.isEcalEnabled()) {
        int neutr_mult = static_cast<int>(reader["neutr_mult"]);

        std::vector<PParticleEcal> ecal_objects;
        ecal_objects.reserve(neutr_mult);

        for (int i = 1; i <= neutr_mult && i <= 5; ++i) {
            // Create ECAL particle with photon hypothesis (mass = 0)
            PParticleEcal ecal_obj(0.0, "ecal" + std::to_string(i));

            // Fill from ntuple variables (neutr_*_N where N = i)
            // Uses cluster_theta/phi for PParticle kinematics
            if (ecal_obj.setFromReader(reader, i)) {
                ecal_objects.push_back(ecal_obj);
            }
        }

        // ecal_objects now contains 0-5 valid ECAL particles

        // STEP 8b: Check ECAL quality cuts for each object
        // Store pass/fail flags for use by both ecal_nt and Step 11 compounds
        std::vector<bool> ecal_pass(ecal_objects.size(), false);
        for (size_t j = 0; j < ecal_objects.size(); ++j) {
            const auto& obj = ecal_objects[j];
            ecal_pass[j] = cuts.passCutSet("ecal_quality", {
                static_cast<double>(obj.ecal_pid),
                obj.ecal_beta,
                obj.cluster_energy
            });
        }

        // Fill ECAL ntuple with detector variables
        auto& ecal_nt = mgr.getDynamicNtuple("ecal_nt");

        ecal_nt["ecal_mult"] = neutr_mult;

        // Helper lambda to fill ntuple for each ECAL hit
        auto fillEcalHit = [&](size_t idx, const std::string& suffix) {
            if (ecal_objects.size() > idx) {
                const auto& obj = ecal_objects[idx];

                ecal_nt["ecal_pass" + suffix] = ecal_pass[idx] ? 1.0f : 0.0f;

                // Cluster reconstruction (used for kinematics)
                ecal_nt["cluster_energy" + suffix] = obj.cluster_energy;
                ecal_nt["cluster_theta" + suffix] = obj.cluster_theta;
                ecal_nt["cluster_phi" + suffix] = obj.cluster_phi;

                // Primary reconstruction
                ecal_nt["ecal_beta" + suffix] = obj.ecal_beta;
                ecal_nt["ecal_pid" + suffix] = obj.ecal_pid;
                ecal_nt["ecal_energy" + suffix] = obj.ecal_energy;
                ecal_nt["ecal_theta" + suffix] = obj.ecal_theta;
                ecal_nt["ecal_phi" + suffix] = obj.ecal_phi;
                ecal_nt["ecal_chi2" + suffix] = obj.ecal_chi2;
                ecal_nt["ecal_tof" + suffix] = obj.ecal_tof;
                ecal_nt["ecal_r" + suffix] = obj.ecal_r;
                ecal_nt["ecal_z" + suffix] = obj.ecal_z;
            }
        };

        // ECAL hits 1-5
        fillEcalHit(0, "_1");
        fillEcalHit(1, "_2");
        fillEcalHit(2, "_3");
        fillEcalHit(3, "_4");
        fillEcalHit(4, "_5");

        ecal_nt.fill();

        // ====================================================================
        // STEP 11: Compound objects e+e-gamma (pi0 Dalitz candidates)
        // ====================================================================
        // For each ECAL photon passing quality cuts, build compound:
        //   epemg = dilepton + gamma
        // One ntuple entry per passing gamma combination.
        // Goal: identify pi0 from Dalitz decay (pi0 -> e+e-gamma)

        auto& epemg_nt = mgr.getDynamicNtuple("epemg_nt");

        // Count passing gammas and build energy-ranked index
        int gamma_pass_mult = std::count(ecal_pass.begin(), ecal_pass.end(), true);

        // Build sorted indices of passing gammas by descending cluster energy
        std::vector<size_t> pass_indices;
        for (size_t j = 0; j < ecal_objects.size(); ++j) {
            if (ecal_pass[j]) pass_indices.push_back(j);
        }
        std::sort(pass_indices.begin(), pass_indices.end(),
                  [&](size_t a, size_t b) {
                      return ecal_objects[a].cluster_energy > ecal_objects[b].cluster_energy;
                  });

        // Assign energy rank: rank_map[index] = rank (1 = highest energy)
        std::map<size_t, int> rank_map;
        for (size_t r = 0; r < pass_indices.size(); ++r) {
            rank_map[pass_indices[r]] = static_cast<int>(r + 1);
        }

        for (size_t j = 0; j < ecal_objects.size(); ++j) {
            if (!ecal_pass[j]) continue;

            // Build compound: e+e-gamma
            PParticle epemg = dilepton + ecal_objects[j];

            // Boost to CMS
            PParticle epemg_cms = frames.getFrame("beam").boost(epemg);

            // Compound variables
            epemg_nt["epemg_mass"] = epemg.massGeV();
            epemg_nt["epemg_p"] = epemg.momentum();
            epemg_nt["epemg_theta"] = epemg.theta();
            epemg_nt["epemg_phi"] = epemg.phi();

            // CMS variables
            epemg_nt["epemg_rapidity_cms"] = epemg_cms.rapidity();
            epemg_nt["epemg_pt_cms"] = epemg_cms.vec().Pt();
            epemg_nt["epemg_theta_cms"] = epemg_cms.theta();

            // Dilepton (e+e-) sub-variables
            epemg_nt["ee_oa"] = oa;
            epemg_nt["ee_mass"] = m_ee;

            // Gamma variables
            epemg_nt["gamma_energy"] = ecal_objects[j].cluster_energy;
            epemg_nt["gamma_theta"] = ecal_objects[j].cluster_theta;
            epemg_nt["gamma_phi"] = ecal_objects[j].cluster_phi;
            epemg_nt["gamma_index"] = static_cast<Float_t>(j + 1);

            // Multiplicity and ranking
            epemg_nt["ecal_mult"] = static_cast<Float_t>(neutr_mult);
            epemg_nt["gamma_pass_mult"] = static_cast<Float_t>(gamma_pass_mult);
            epemg_nt["gamma_rank_energy"] = static_cast<Float_t>(rank_map[j]);

            epemg_nt.fill();
        }
    }

    // ========================================================================
    // STEP 9: Forward Tracker objects
    // ========================================================================
    // Create particle objects from Forward Tracker detector.
    // FT can detect hadrons (protons, pions, etc.) in the forward direction.
    // fwdet_mult tells us how many FT hits are present (0-3).
    //
    // PParticleFwd extends PParticle with FT-specific data:
    //   - Works with operator+ for composite particle creation
    //   - Stores all FT detector variables (beta, chi2, mass, etc.)
    //   - Uses fwdet_p, fwdet_theta, fwdet_phi for kinematics
    //
    // Controlled by config.json: "fwdet": true/false

    if (config.isFwdEnabled()) {
        int fwdet_mult = static_cast<int>(reader["fwdet_mult"]);

        std::vector<PParticleFwd> fwd_objects;
        fwd_objects.reserve(fwdet_mult);

        for (int i = 1; i <= fwdet_mult && i <= 3; ++i) {
            // Create FT particle with proton hypothesis (most common)
            PParticleFwd fwd_obj(MASS_PROTON, "fwd" + std::to_string(i));

            // Fill from ntuple variables (fwdet_*_N where N = i)
            if (fwd_obj.setFromReader(reader, i)) {
                fwd_objects.push_back(fwd_obj);
            }
        }

        // fwd_objects now contains 0-3 valid FT particles
        // These can be combined with other PParticles:
        //   PParticle composite = dilepton + fwd_objects[0];

        // Fill FT ntuple with detector variables
        auto& fwd_nt = mgr.getDynamicNtuple("fwdet_nt");

        fwd_nt["fwd_mult"] = fwdet_mult;

        // Helper lambda to fill ntuple for each FT hit
        auto fillFwdHit = [&](size_t idx, const std::string& suffix) {
            if (fwd_objects.size() > idx) {
                const auto& obj = fwd_objects[idx];

                // Kinematics
                fwd_nt["fwd_p" + suffix] = obj.momentum();
                fwd_nt["fwd_theta" + suffix] = obj.fwd_theta;
                fwd_nt["fwd_phi" + suffix] = obj.fwd_phi;

                // Measured quantities
                fwd_nt["fwd_beta" + suffix] = obj.fwd_beta;
                fwd_nt["fwd_mass" + suffix] = obj.fwd_mass;
                fwd_nt["fwd_mass2" + suffix] = obj.fwd_mass2;
                fwd_nt["fwd_charge" + suffix] = obj.fwd_charge;

                // Quality
                fwd_nt["fwd_chi2" + suffix] = obj.fwd_chi2;
                fwd_nt["fwd_ndf" + suffix] = obj.fwd_ndf;
                fwd_nt["fwd_chi2ndf" + suffix] = obj.fwd_chi2ndf;

                // Position
                fwd_nt["fwd_r" + suffix] = obj.fwd_r;
                fwd_nt["fwd_z" + suffix] = obj.fwd_z;
                fwd_nt["fwd_distToRpc" + suffix] = obj.fwd_distToRpc;
            }
        };

        // FT hits 1-3
        fillFwdHit(0, "_1");
        fillFwdHit(1, "_2");
        fillFwdHit(2, "_3");

        fwd_nt.fill();
    }

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
    // STEP 10: Apply lepton prefix mapping (for combinatorial background)
    // ========================================================================
    // If config has lepton_prefixes, apply them to the reader.
    // This allows the same processEvent() code to work with EpEp_ID and EmEm_ID
    // trees where variable names differ (ep1_p/ep2_p instead of ep_p/em_p).
    // When no prefixes are defined (EpEm_ID), this does nothing.

    auto [prefix1, prefix2] = config.getLeptonPrefixes();
    if (!prefix1.empty()) {
        reader.setLeptonPrefixes(prefix1, prefix2);
        std::string channel = config.getChannelName();
        std::cout << "NTupleReader: Prefix mapping ep_ -> " << prefix1
                  << "_, em_ -> " << prefix2 << "_";
        if (!channel.empty()) std::cout << " (channel: " << channel << ")";
        std::cout << "\n";
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

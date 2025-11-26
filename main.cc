// ========================================================================
// FAT Framework - Modern Analysis Example
// ========================================================================
// Main analysis file - contains only physics-related code:
// - processEvent: Event-by-event physics analysis
// - main: Program entry point
//
// Setup code is in separate files:
// - src/setup_histograms.h: Histogram definitions
// - src/setup_ntuples.h: Ntuple definitions
// - src/setup_cuts.h: Cut definitions
//
// Usage:
//   ./ana [config.json]
//   ./ana                    # Uses default config.json
//   ./ana my_analysis.json   # Uses custom config file
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2025
// ========================================================================

#include "src/manager.h"
#include "src/pparticle.h"
#include "src/boost_frame.h"
#include "src/ntuple_reader.h"
#include "src/cut_manager.h"
#include "src/analysis_config.h"
#include "src/setup_histograms.h"
#include "src/setup_ntuples.h"
#include "src/setup_cuts.h"
#include "src/progressbar.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

// Use Physics namespace for mass constants
using namespace Physics;

// ============================================================================
// PROCESS SINGLE EVENT - Physics Analysis
// ============================================================================
// This is the main physics function. Edit here to customize your analysis.
// 
// Structure:
// 1. Read variables from ntuple
// 2. Create particles
// 3. Fill quality histograms
// 4. Apply cuts
// 5. Boost to CMS
// 6. Fill physics histograms
// ============================================================================

void processEvent(NTupleReader& reader, Manager& mgr, CutManager& cuts,
                 const PParticle& beam, const PParticle& projectile,
                 EventFrames& frames, bool use_corrected) {
    
    // ========================================================================
    // 1. READ KINEMATIC VARIABLES FROM NTUPLE
    // ========================================================================
    // These are the ONLY places where reader["..."] is needed.
    // Variable names must exactly match your tree branch names.
    
    // Proton kinematics
    double p_p = use_corrected ? reader["p_p_corr_p"] : reader["p_p"];
    double p_theta = reader["p_theta"];
    double p_phi = reader["p_phi"];
    
    // Pion kinematics
    double pip_p = use_corrected ? reader["pip_p_corr_pip"] : reader["pip_p"];
    double pip_theta = reader["pip_theta"];
    double pip_phi = reader["pip_phi"];
    
    // Event weight (optional)
    double weight = reader.hasVariable("weight") ? reader["weight"] : 1.0;
    
    // Vertex (optional)
    if (reader.hasVariable("eVertX")) {
        mgr.fill("eVertX", reader["eVertX"]);
        mgr.fill("eVertY", reader["eVertY"]);
        mgr.fill("eVertZ", reader["eVertZ"]);
    }
    
    // ========================================================================
    // 2. CREATE PARTICLES
    // ========================================================================
    // From here on, use regular C++ variables - no more reader["..."]
    
    PParticle proton = ParticleFactory::createProton(p_p, p_theta, p_phi);
    PParticle pion = ParticleFactory::createPiPlus(pip_p, pip_theta, pip_phi);
    
    // Missing mass technique for neutron
    PParticle neutron = beam - proton - pion;
    
    // Composite particles
    PParticle deltaPP = proton + pion;      // Delta++
    PParticle deltaP = beam - proton;        // Delta+ (missing)
    PParticle p_pip = proton + pion;
    PParticle n_pip = neutron + pion;
    PParticle pn = proton + neutron;
    
    // ========================================================================
    // 3. QUALITY HISTOGRAMS (before cuts)
    // ========================================================================
    
    double m_n = neutron.massGeV();
    double m_p = proton.massGeV();
    double m_pip = pion.massGeV();
    
    mgr.fill("mass_n", m_n);
    mgr.fill("mass_p", m_p);
    mgr.fill("mass_pip", m_pip);
    
    // ========================================================================
    // 4. APPLY CUTS
    // ========================================================================
    
    // Neutron mass cut
    if (cuts.hasRangeCut("neutron_mass")) {
        if (!cuts.passRangeCut("neutron_mass", m_n)) return;
    }
    
    // Fill after neutron cut
    mgr.fill("mass_n_cut", m_n);
    
    // Delta++ mass cut
    double m_deltaPP = deltaPP.massGeV();
    if (cuts.hasRangeCut("deltaPP_mass")) {
        if (!cuts.passRangeCut("deltaPP_mass", m_deltaPP)) return;
    }
    
    // ========================================================================
    // 5. BOOST TO CMS
    // ========================================================================
    
    const BoostFrame& beam_frame = frames.getFrame("beam");
    
    PParticle p_cms = beam_frame.boost(proton);
    PParticle pip_cms = beam_frame.boost(pion);
    PParticle n_cms = beam_frame.boost(neutron);
    PParticle deltaPP_cms = beam_frame.boost(deltaPP);
    PParticle deltaP_cms = beam_frame.boost(deltaP);
    
    // ========================================================================
    // 6. FILL HISTOGRAMS
    // ========================================================================
    
    // Composite masses
    mgr.fill("mass_deltaPP", m_deltaPP);
    mgr.fill("mass_deltaP", deltaP.massGeV());
    mgr.fill("mass_ppip", p_pip.massGeV());
    mgr.fill("mass_npip", n_pip.massGeV());
    mgr.fill("mass_pn", pn.massGeV());
    
    // LAB frame kinematics
    mgr.fill("p_p_lab", proton.momentum());
    mgr.fill("pip_p_lab", pion.momentum());
    mgr.fill("n_p_lab", neutron.momentum());
    
    mgr.fill("p_theta_lab", proton.theta());
    mgr.fill("pip_theta_lab", pion.theta());
    mgr.fill("n_theta_lab", neutron.theta());
    
    // CMS kinematics
    mgr.fill("cos_theta_deltaPP_cms", deltaPP_cms.cosTheta());
    mgr.fill("cos_theta_deltaP_cms", deltaP_cms.cosTheta());
    mgr.fill("cos_theta_p_cms", p_cms.cosTheta());
    mgr.fill("cos_theta_pip_cms", pip_cms.cosTheta());
    mgr.fill("cos_theta_n_cms", n_cms.cosTheta());
    
    mgr.fill("p_p_cms", p_cms.momentum());
    mgr.fill("pip_p_cms", pip_cms.momentum());
    mgr.fill("n_p_cms", n_cms.momentum());
    
    // Opening angles
    mgr.fill("oa_ppip", proton.openingAngle(pion));
    mgr.fill("oa_npip", neutron.openingAngle(pion));
    mgr.fill("oa_pn", proton.openingAngle(neutron));
    
    // 2D correlations
    double m2_ppip = p_pip.mass() * p_pip.mass() / 1e6;  // GeV^2
    double m2_npip = n_pip.mass() * n_pip.mass() / 1e6;  // GeV^2
    mgr.fill("dalitz_ppip_npip", m2_ppip, m2_npip);
    
    mgr.fill("mass_vs_costh_deltaPP", m_deltaPP, deltaPP_cms.cosTheta());
    mgr.fill("theta_p_vs_pip_lab", pion.theta(), proton.theta());
    
    // ========================================================================
    // 7. PWA VARIABLES (in composite rest frames)
    // ========================================================================
    
    // Create p+pi+ rest frame
    BoostFrame ppip_frame(p_pip);
    
    PParticle pip_in_ppip = ppip_frame.boost(pion);
    PParticle n_in_ppip = ppip_frame.boost(neutron);
    PParticle proj_in_ppip = ppip_frame.boost(projectile);
    
    // Helicity angle: pion angle relative to beam direction in ppip frame
    mgr.fill("pwa_pip_helicity_ppip", pip_in_ppip.cosTheta());
    mgr.fill("pwa_n_helicity_ppip", n_in_ppip.cosTheta());
    
    // Gottfried-Jackson: angle relative to beam in composite frame
    double gj_angle = pip_in_ppip.vec().Angle(proj_in_ppip.vec().Vect());
    mgr.fill("pwa_pip_gj_ppip", cos(gj_angle));
    
    // ========================================================================
    // 8. FILL OUTPUT NTUPLES
    // ========================================================================
    // Two ntuples demonstrate multiple TTree→TNtuple handling
    
    // --- Ntuple 1: Basic particle observables ---
    DynamicHNtuple& nt_particles = mgr.getDynamicNtuple("nt_particles");
    
    // Proton observables (LAB)
    nt_particles["p_p"] = proton.momentum();
    nt_particles["p_theta"] = proton.theta();
    nt_particles["p_phi"] = proton.phi();
    nt_particles["p_mass"] = m_p;
    
    // Pion observables (LAB)
    nt_particles["pip_p"] = pion.momentum();
    nt_particles["pip_theta"] = pion.theta();
    nt_particles["pip_phi"] = pion.phi();
    nt_particles["pip_mass"] = m_pip;
    
    // Neutron observables (missing mass)
    nt_particles["n_p"] = neutron.momentum();
    nt_particles["n_theta"] = neutron.theta();
    nt_particles["n_phi"] = neutron.phi();
    nt_particles["n_mass"] = m_n;
    
    // Event weight
    nt_particles["weight"] = weight;
    
    nt_particles.fill();
    
    // --- Ntuple 2: Compound observables ---
    DynamicHNtuple& nt_compound = mgr.getDynamicNtuple("nt_compound");
    
    // Composite masses
    nt_compound["m_deltaPP"] = m_deltaPP;
    nt_compound["m_deltaP"] = deltaP.massGeV();
    nt_compound["m_ppip"] = p_pip.massGeV();
    nt_compound["m_npip"] = n_pip.massGeV();
    nt_compound["m_pn"] = pn.massGeV();
    
    // CMS angles (composite particles)
    nt_compound["cos_th_deltaPP_cms"] = deltaPP_cms.cosTheta();
    nt_compound["cos_th_deltaP_cms"] = deltaP_cms.cosTheta();
    nt_compound["cos_th_p_cms"] = p_cms.cosTheta();
    nt_compound["cos_th_pip_cms"] = pip_cms.cosTheta();
    nt_compound["cos_th_n_cms"] = n_cms.cosTheta();
    
    // Opening angles
    nt_compound["oa_ppip"] = proton.openingAngle(pion);
    nt_compound["oa_npip"] = neutron.openingAngle(pion);
    nt_compound["oa_pn"] = proton.openingAngle(neutron);
    
    // PWA variables (helicity and Gottfried-Jackson angles)
    nt_compound["pip_helicity"] = pip_in_ppip.cosTheta();
    nt_compound["pip_gj"] = cos(gj_angle);
    nt_compound["n_helicity"] = n_in_ppip.cosTheta();
    
    // Dalitz plot variables (squared masses)
    nt_compound["m2_ppip"] = m2_ppip;
    nt_compound["m2_npip"] = m2_npip;
    
    // Event weight
    nt_compound["weight"] = weight;
    
    nt_compound.fill();
}

// ============================================================================
// MAIN PROGRAM
// ============================================================================

int main(int argc, char* argv[]) {
    // Install signal handler for graceful Ctrl+C termination
    SignalHandler::install();
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                                                  ║\n";
    std::cout << "║     FAT Framework - Final Analysis Tool                          ║\n";
    std::cout << "║     pp → npπ+ (n missing) Analysis                               ║\n";
    std::cout << "║                                                                  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
    
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
        std::cerr << "Usage: " << argv[0] << " [config.json]\n";
        return 1;
    }
    
    config.print();
    
    // ========================================================================
    // 2. SETUP BEAM
    // ========================================================================
    
    double beam_kinetic_energy = config.getBeamKineticEnergy();
    
    PParticle projectile = ParticleFactory::createBeamProton(beam_kinetic_energy);
    PParticle target = ParticleFactory::createTargetProton();
    PParticle beam = projectile + target;
    
    std::cout << "\nBeam Setup:\n";
    std::cout << "  Kinetic energy: " << beam_kinetic_energy << " MeV\n";
    std::cout << "  Beam beta: " << beam.beta() << "\n";
    std::cout << "  sqrt(s): " << beam.mass() / 1000.0 << " GeV\n";
    
    // Setup event frames
    EventFrames frames;
    frames.setBeamFrame(projectile, target);
    
    // ========================================================================
    // 3. OPEN INPUT DATA
    // ========================================================================
    
    NTupleReader reader;
    
    try {
        std::string source = config.getInputSource();
        std::string tree_name = config.getInputTreeName();
        
        if (config.isInputFileList()) {
            reader.openFromList(source, tree_name);
        } else if (config.isInputRootFile()) {
            reader.open(source, tree_name);
        } else {
            throw std::runtime_error("Unknown input format. Use .root or .list file");
        }
    } catch (const std::exception& e) {
        std::cerr << "Error opening input: " << e.what() << "\n";
        return 1;
    }
    
    // ========================================================================
    // 4. OPEN OUTPUT FILE & SETUP HISTOGRAMS
    // ========================================================================
    
    Manager manager;
    manager.openFile(config.getOutputFilename(), config.getOutputOption());
    
    // Setup histograms (defined in src/setup_histograms.h)
    setupHistograms(manager);
    
    // Setup ntuples (defined in src/setup_ntuples.h)
    setupNtuples(manager, config);
    
    // ========================================================================
    // 5. SETUP CUTS
    // ========================================================================
    
    CutManager cuts;
    // Setup cuts (defined in src/setup_cuts.h)
    setupCuts(cuts);
    
    // ========================================================================
    // 6. EVENT LOOP
    // ========================================================================
    
    Long64_t total_entries = reader.entries();
    Long64_t start_event = config.getStartEvent();
    Long64_t max_events = config.getMaxEvents();
    
    // Calculate end event
    Long64_t end_event = total_entries;
    if (max_events > 0) {
        end_event = std::min(start_event + max_events, total_entries);
    }
    
    Long64_t events_to_process = end_event - start_event;
    
    // User decides which momentum to use in processEvent()
    // Set this flag based on your analysis needs:
    bool use_corrected = true;  // Change to false for raw momentum
    
    std::cout << "\n";
    std::cout << "┌───────────────────────────────────────────────────────────────┐\n";
    std::cout << "│  Press Ctrl+C at any time to stop and save partial results    │\n";
    std::cout << "└───────────────────────────────────────────────────────────────┘\n";
    std::cout << "\n";
    std::cout << "Processing events " << start_event << " to " << end_event 
              << " (" << events_to_process << " events)...\n";
    std::cout << "\n";
    
    Long64_t processed = 0;
    bool was_interrupted = false;
    
    // Progress bar with time estimation
    ProgressBar progress(events_to_process);
    
    for (Long64_t i = start_event; i < end_event; ++i) {
        // Check for Ctrl+C - graceful termination
        if (SignalHandler::wasInterrupted()) {
            was_interrupted = true;
            break;
        }
        
        reader.getEntry(i);
        
        ++processed;
        
        // Update progress bar (updates only on percent change)
        progress.update(processed);
        
        // Process event
        try {
            processEvent(reader, manager, cuts, beam, projectile, frames, use_corrected);
        } catch (const std::exception& e) {
            // Skip events with missing variables
            continue;
        }
    }
    
    // Finish progress bar (shows total elapsed time or interrupted status)
    progress.finish(was_interrupted);
    
    std::cout << "\n";
    if (was_interrupted) {
        std::cout << "Processing interrupted by user (Ctrl+C).\n";
    } else {
        std::cout << "Processing complete!\n";
    }
    std::cout << "  Events processed: " << processed << "\n";
    
    // ========================================================================
    // 7. PRINT CUT FLOW
    // ========================================================================
    
    cuts.printCutFlow();
    
    // ========================================================================
    // 8. SAVE AND CLOSE
    // ========================================================================
    
    std::cout << "\nSaving results to " << config.getOutputFilename() << "...\n";
    manager.printSummary();
    manager.closeFile();
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                     Analysis Complete!                           ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    return 0;
}

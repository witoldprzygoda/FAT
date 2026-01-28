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
#include "src/pparticle_fwd.h"
#include "src/pparticle_ecal.h"
#include "src/boost_frame.h"
#include "src/ntuple_reader.h"
#include "src/cut_manager.h"
#include "src/analysis_config.h"
#include "src/setup_histograms.h"
#include "src/setup_ntuples.h"
#include "src/setup_cuts.h"
#include "src/progressbar.h"
#include "src/console_box.h"
#include "src/reactionvertexfind.h"
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
                 EventFrames& frames, const AnalysisConfig& config) {
    
    // Get which kinematic type to use for analysis (from config)
    KinematicType analysis_type = config.getAnalysisKinematicType();
    
    // ========================================================================
    // 1. CREATE PARTICLE SHELLS (just mass, no kinematics yet)
    // ========================================================================
    
    PParticle proton(Physics::MASS_PROTON, "p");
    PParticle pion(Physics::MASS_PION_PLUS, "pi+");
    
    // Forward Tracker particles (assume proton mass hypothesis)
    PParticleFwd fw_p1(Physics::MASS_PROTON, "fw_p1");
    PParticleFwd fw_p2(Physics::MASS_PROTON, "fw_p2");
    PParticleFwd fw_p3(Physics::MASS_PROTON, "fw_p3");
    
    // ECAL neutral particles (default: photon hypothesis, mass = 0)
    PParticleEcal ecal_n1(0.0, "ecal_n1");
    PParticleEcal ecal_n2(0.0, "ecal_n2");
    PParticleEcal ecal_n3(0.0, "ecal_n3");
    
    // ========================================================================
    // 2. FILL KINEMATIC TYPES BASED ON CONFIG
    // ========================================================================
    // Each PParticle can hold multiple representations (RECONSTRUCTED, CORRECTED, SIMULATED)
    // The config determines which types are populated from the ntuple
    
    // --- RECONSTRUCTED kinematics (spherical: p, theta, phi) ---
    if (config.hasReconstructed()) {
        proton.setFromSpherical(reader["p_p"], reader["p_theta"], reader["p_phi"],
                                KinematicType::RECONSTRUCTED);
        pion.setFromSpherical(reader["pip_p"], reader["pip_theta"], reader["pip_phi"],
                              KinematicType::RECONSTRUCTED);
    }
    
    // --- CORRECTED kinematics (spherical: only momentum corrected, same angles) ---
    if (config.hasCorrected()) {
        proton.setFromSpherical(reader["p_p_corr_p"], reader["p_theta"], reader["p_phi"],
                                KinematicType::CORRECTED);
        pion.setFromSpherical(reader["pip_p_corr_pip"], reader["pip_theta"], reader["pip_phi"],
                              KinematicType::CORRECTED);
    }
    
    // --- SIMULATED kinematics (Cartesian: px, py, pz) ---
    if (config.hasSimulated()) {
        proton.setFromCartesian(reader["p_sim_px"], reader["p_sim_py"], reader["p_sim_pz"],
                                KinematicType::SIMULATED);
        pion.setFromCartesian(reader["pip_sim_px"], reader["pip_sim_py"], reader["pip_sim_pz"],
                              KinematicType::SIMULATED);
    }
    
    // --- Forward Tracker particles (if enabled in config) ---
    // Currently only RECONSTRUCTED; CORRECTED/SIMULATED prepared for future
    if (config.hasFwdet1()) {
        fw_p1.setFromReader(reader, 1, KinematicType::RECONSTRUCTED);
    }
    if (config.hasFwdet2()) {
        fw_p2.setFromReader(reader, 2, KinematicType::RECONSTRUCTED);
    }
    if (config.hasFwdet3()) {
        fw_p3.setFromReader(reader, 3, KinematicType::RECONSTRUCTED);
    }
    
    // --- ECAL particles (if enabled in config) ---
    // Default: photon hypothesis (mass=0, uses energy as momentum)
    // Currently only RECONSTRUCTED; CORRECTED/SIMULATED prepared for future
    if (config.hasEcal1()) {
        ecal_n1.setFromReader(reader, 1, KinematicType::RECONSTRUCTED);
    }
    if (config.hasEcal2()) {
        ecal_n2.setFromReader(reader, 2, KinematicType::RECONSTRUCTED);
    }
    if (config.hasEcal3()) {
        ecal_n3.setFromReader(reader, 3, KinematicType::RECONSTRUCTED);
    }
    
    // ========================================================================
    // WEIGHTING: Combine event weight with luminosity for cross-section scaling
    // ========================================================================
    double event_weight = reader.hasVariable("weight") ? reader["weight"] : 1.0;
    double luminosity = config.getLuminosity();  // From config.json
    double weight = event_weight * luminosity;   // Combined weight for all histograms
    
    // Vertex (optional)
    if (reader.hasVariable("eVertX")) {
        mgr.fill("eVertX", reader["eVertX"]);
        mgr.fill("eVertY", reader["eVertY"]);
        mgr.fill("eVertZ", reader["eVertZ"]);
    }
    
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
    // All getters use analysis_type to select which kinematic representation to use
    
    double m_n = neutron.massGeV(analysis_type);
    double m_p = proton.massGeV(analysis_type);
    double m_pip = pion.massGeV(analysis_type);
    
    mgr.fillw("mass_n", m_n, weight);
    mgr.fillw("mass_p", m_p, weight);
    mgr.fillw("mass_pip", m_pip, weight);
    
    // ========================================================================
    // 4. APPLY CUTS
    // ========================================================================
    
    // Neutron mass cut
    if (cuts.hasRangeCut("neutron_mass")) {
        if (!cuts.passRangeCut("neutron_mass", m_n)) return;
    }
    
    // Fill after neutron cut
    mgr.fillw("mass_n_cut", m_n, weight);
    
    // Delta++ mass cut
    double m_deltaPP = deltaPP.massGeV(analysis_type);
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
    // All accessors use analysis_type to select which kinematic representation
    
    // Composite masses
    mgr.fillw("mass_deltaPP", m_deltaPP, weight);
    mgr.fillw("mass_deltaP", deltaP.massGeV(analysis_type), weight);
    mgr.fillw("mass_ppip", p_pip.massGeV(analysis_type), weight);
    mgr.fillw("mass_npip", n_pip.massGeV(analysis_type), weight);
    mgr.fillw("mass_pn", pn.massGeV(analysis_type), weight);
    
    // LAB frame kinematics
    mgr.fillw("p_p_lab", proton.momentum(analysis_type), weight);
    mgr.fillw("pip_p_lab", pion.momentum(analysis_type), weight);
    mgr.fillw("n_p_lab", neutron.momentum(analysis_type), weight);
    
    mgr.fillw("p_theta_lab", proton.theta(analysis_type), weight);
    mgr.fillw("pip_theta_lab", pion.theta(analysis_type), weight);
    mgr.fillw("n_theta_lab", neutron.theta(analysis_type), weight);
    
    // CMS kinematics
    mgr.fillw("cos_theta_deltaPP_cms", deltaPP_cms.cosTheta(analysis_type), weight);
    mgr.fillw("cos_theta_deltaP_cms", deltaP_cms.cosTheta(analysis_type), weight);
    mgr.fillw("cos_theta_p_cms", p_cms.cosTheta(analysis_type), weight);
    mgr.fillw("cos_theta_pip_cms", pip_cms.cosTheta(analysis_type), weight);
    mgr.fillw("cos_theta_n_cms", n_cms.cosTheta(analysis_type), weight);
    
    mgr.fillw("p_p_cms", p_cms.momentum(analysis_type), weight);
    mgr.fillw("pip_p_cms", pip_cms.momentum(analysis_type), weight);
    mgr.fillw("n_p_cms", n_cms.momentum(analysis_type), weight);
    
    // Opening angles
    mgr.fillw("oa_ppip", proton.openingAngle(pion, analysis_type), weight);
    mgr.fillw("oa_npip", neutron.openingAngle(pion, analysis_type), weight);
    mgr.fillw("oa_pn", proton.openingAngle(neutron, analysis_type), weight);
    
    // 2D correlations
    double m2_ppip = p_pip.mass() * p_pip.mass() / 1e6;  // GeV^2
    double m2_npip = n_pip.mass() * n_pip.mass() / 1e6;  // GeV^2
    mgr.fillw("dalitz_ppip_npip", m2_ppip, m2_npip, weight);
    
    mgr.fillw("mass_vs_costh_deltaPP", m_deltaPP, deltaPP_cms.cosTheta(), weight);
    mgr.fillw("theta_p_vs_pip_lab", pion.theta(), proton.theta(), weight);
    
    // ========================================================================
    // 7. PWA VARIABLES (in composite rest frames)
    // ========================================================================
    
    // Create p+pi+ rest frame
    BoostFrame ppip_frame(p_pip);
    
    PParticle pip_in_ppip = ppip_frame.boost(pion);
    PParticle n_in_ppip = ppip_frame.boost(neutron);
    PParticle proj_in_ppip = ppip_frame.boost(projectile);
    
    // Helicity angle: pion angle relative to beam direction in ppip frame
    mgr.fillw("pwa_pip_helicity_ppip", pip_in_ppip.cosTheta(), weight);
    mgr.fillw("pwa_n_helicity_ppip", n_in_ppip.cosTheta(), weight);
    
    // Gottfried-Jackson: angle relative to beam in composite frame
    double gj_angle = pip_in_ppip.vec().Angle(proj_in_ppip.vec().Vect());
    mgr.fillw("pwa_pip_gj_ppip", cos(gj_angle), weight);
    
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
    
    // ========================================================================
    // 9. REACTION VERTEX FINDING (OPTIONAL - uncomment to use)
    // ========================================================================
    // Calculate vertex from track closest approach to beam axis (x=y=0).
    // Each track has (r, z) = closest distance and z-position to beam axis.
    // ReactionVertexFind combines N tracks (N >= 2) to find the common vertex.
    //
    // UNCOMMENT THE BLOCK BELOW AND ADAPT VARIABLE NAMES TO YOUR NTUPLE:
    //
    // /*
    // // --- Vertex finder setup (can be moved outside event loop for efficiency) ---
    // static ReactionVertexFind vtxFinder;
    // 
    // // Configure cuts (do once, or per-event if needed)
    // // setCuts(minZ, maxZ, maxR) - z and r acceptance
    // vtxFinder.setCuts(-120., 20., 15.);  // Default: z in [-120, 20] mm, r < 15 mm
    // 
    // // setQualityCuts(maxRKChi2, maxSegChi2, minBeta, maxBeta)
    // vtxFinder.setQualityCuts(100., 6., 0., 1.2);  // Default quality cuts
    // 
    // // To disable all cuts (for testing):
    // // vtxFinder.setCuts(-1000., 1000., 1000.);       // z: any, r: any
    // // vtxFinder.setQualityCuts(1e9, 1e9, -10., 10.); // chi2: any, beta: any
    // 
    // // --- Set tracks from ntuple (2-track example: proton + pion) ---
    // // Variable names must match your ntuple branches!
    // vtxFinder.setTracks2(
    //     reader["p_r"],   reader["p_z"],   reader["p_theta"],   reader["p_phi"],
    //     reader["p_rkchi2"],   reader["p_mdcinnerchi2"],   reader["p_beta"],
    //     reader["pip_r"], reader["pip_z"], reader["pip_theta"], reader["pip_phi"],
    //     reader["pip_rkchi2"], reader["pip_mdcinnerchi2"], reader["pip_beta"]
    // );
    // 
    // // --- Alternative: 3 tracks (p, pip, pim) ---
    // // vtxFinder.setTracks3(
    // //     reader["p_r"],   reader["p_z"],   reader["p_theta"],   reader["p_phi"],
    // //     reader["p_rkchi2"],   reader["p_mdcinnerchi2"],   reader["p_beta"],
    // //     reader["pip_r"], reader["pip_z"], reader["pip_theta"], reader["pip_phi"],
    // //     reader["pip_rkchi2"], reader["pip_mdcinnerchi2"], reader["pip_beta"],
    // //     reader["pim_r"], reader["pim_z"], reader["pim_theta"], reader["pim_phi"],
    // //     reader["pim_rkchi2"], reader["pim_mdcinnerchi2"], reader["pim_beta"]
    // // );
    // 
    // // --- Alternative: 4 tracks (pim, pip, em, ep) ---
    // // vtxFinder.setTracks4(
    // //     reader["pim_r"], reader["pim_z"], reader["pim_theta"], reader["pim_phi"],
    // //     reader["pim_rkchi2"], reader["pim_mdcinnerchi2"], reader["pim_beta"],
    // //     reader["pip_r"], reader["pip_z"], reader["pip_theta"], reader["pip_phi"],
    // //     reader["pip_rkchi2"], reader["pip_mdcinnerchi2"], reader["pip_beta"],
    // //     reader["em_r"],  reader["em_z"],  reader["em_theta"],  reader["em_phi"],
    // //     reader["em_rkchi2"],  reader["em_mdcinnerchi2"],  reader["em_beta"],
    // //     reader["ep_r"],  reader["ep_z"],  reader["ep_theta"],  reader["ep_phi"],
    // //     reader["ep_rkchi2"],  reader["ep_mdcinnerchi2"],  reader["ep_beta"]
    // // );
    // 
    // // --- Find vertex and get results ---
    // if (vtxFinder.findVertex()) {
    //     // Vertex position
    //     Float_t vx = vtxFinder.getVx();
    //     Float_t vy = vtxFinder.getVy();
    //     Float_t vz = vtxFinder.getVz();
    //     
    //     // Quality parameters
    //     Float_t chi2     = vtxFinder.getChi2();
    //     Float_t rVertex  = vtxFinder.getRVertex();  // sqrt(vx^2 + vy^2)
    //     Float_t sumW     = vtxFinder.getSumOfWeights();
    //     Float_t zSpread  = vtxFinder.getZSpread();
    //     Int_t   nIter    = vtxFinder.getIterations();
    //     
    //     // Per-track info: track index 0 = proton, 1 = pion (order of setTracks2)
    //     Float_t p_dca   = vtxFinder.getTrackDCA(0);    // Proton distance to vertex
    //     Float_t pip_dca = vtxFinder.getTrackDCA(1);    // Pion distance to vertex
    //     
    //     // Quality check (recommended: chi2 < 2, rVertex < 3 mm)
    //     Bool_t isGood = vtxFinder.isGoodVertex(2.0, 3.0);
    //     
    //     // Fill histograms (create these in setup_histograms.h first!)
    //     // mgr.fillw("vtx_x", vx, weight);
    //     // mgr.fillw("vtx_y", vy, weight);
    //     // mgr.fillw("vtx_z", vz, weight);
    //     // mgr.fillw("vtx_r", rVertex, weight);
    //     // mgr.fillw("vtx_chi2", chi2, weight);
    //     // mgr.fillw("vtx_xy", vx, vy, weight);  // 2D
    //     
    //     // Add to ntuple
    //     // nt_particles["vtx_x"] = vx;
    //     // nt_particles["vtx_y"] = vy;
    //     // nt_particles["vtx_z"] = vz;
    //     // nt_particles["vtx_chi2"] = chi2;
    // }
    // */
}

// ============================================================================
// MAIN PROGRAM
// ============================================================================

int main(int argc, char* argv[]) {
    // Install signal handler for graceful Ctrl+C termination
    SignalHandler::install();
    
    ConsoleBox::newLine();
    ConsoleBox::printHeader("FAT Framework - Final Analysis Tool",
                           "pp → npπ+ (n missing) Analysis");
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
    
    ConsoleBox::newLine();
    ConsoleBox::printInfoBox("Press Ctrl+C at any time to stop and save partial results");
    ConsoleBox::newLine();
    std::cout << "Processing events " << start_event << " to " << end_event 
              << " (" << events_to_process << " events)...\n\n";
    
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
            processEvent(reader, manager, cuts, beam, projectile, frames, config);
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
    
    ConsoleBox::newLine();
    ConsoleBox::printStatus("Analysis Complete!");
    ConsoleBox::newLine();
    
    return 0;
}

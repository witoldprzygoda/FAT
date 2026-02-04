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
#include <vector>

// Use Physics namespace for mass constants
using namespace Physics;

// ============================================================================
// Helper function: Check if photon passes quality cuts
// ============================================================================
bool isGoodPhoton(const PParticleEcal& gamma, CutManager& cuts) {
    if (!gamma.isValid()) return false;

    double energy = gamma.getEnergy();
    double theta = gamma.ecal_theta;

    // Apply energy cut
    if (cuts.hasRangeCut("photon_energy")) {
        if (!cuts.passRangeCut("photon_energy", energy)) return false;
    }

    // Apply theta cut
    if (cuts.hasRangeCut("photon_theta")) {
        if (!cuts.passRangeCut("photon_theta", theta)) return false;
    }

    return true;
}

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

    PParticle positron(Physics::MASS_ELECTRON, "e+");   // e+ (positron)
    PParticle electron(Physics::MASS_ELECTRON, "e-");   // e- (electron)

    // Forward Tracker particles (assume proton mass hypothesis)
    PParticleFwd fw_p1(Physics::MASS_PROTON, "fw_p1");
    PParticleFwd fw_p2(Physics::MASS_PROTON, "fw_p2");
    PParticleFwd fw_p3(Physics::MASS_PROTON, "fw_p3");

    // ECAL photon candidates (mass = 0 for photon hypothesis)
    PParticleEcal gamma1(0.0, "gamma1");
    PParticleEcal gamma2(0.0, "gamma2");
    PParticleEcal gamma3(0.0, "gamma3");

    // ========================================================================
    // 2. FILL KINEMATIC TYPES BASED ON CONFIG
    // ========================================================================
    // Each PParticle can hold multiple representations (RECONSTRUCTED, CORRECTED, SIMULATED)
    // The config determines which types are populated from the ntuple

    // --- RECONSTRUCTED kinematics (spherical: p, theta, phi) ---
    if (config.hasReconstructed()) {
        positron.setFromSpherical(reader["ep_p"], reader["ep_theta"], reader["ep_phi"],
                                KinematicType::RECONSTRUCTED);
        electron.setFromSpherical(reader["em_p"], reader["em_theta"], reader["em_phi"],
                              KinematicType::RECONSTRUCTED);
    }

    // --- CORRECTED kinematics (spherical: only momentum corrected, same angles) ---
    if (config.hasCorrected()) {
        positron.setFromSpherical(reader["ep_p_corr_ep"], reader["ep_theta"], reader["ep_phi"],
                                KinematicType::CORRECTED);
        electron.setFromSpherical(reader["em_p_corr_em"], reader["em_theta"], reader["em_phi"],
                              KinematicType::CORRECTED);
    }

    // --- SIMULATED kinematics (Cartesian: px, py, pz) ---
    if (config.hasSimulated()) {
        positron.setFromCartesian(reader["ep_sim_px"], reader["ep_sim_py"], reader["ep_sim_pz"],
                                KinematicType::SIMULATED);
        electron.setFromCartesian(reader["em_sim_px"], reader["em_sim_py"], reader["em_sim_pz"],
                              KinematicType::SIMULATED);
    }

    // --- Forward Tracker particles (if enabled in config) ---
    if (config.hasFwdet1()) {
        fw_p1.setFromReader(reader, 1, KinematicType::RECONSTRUCTED);
    }
    if (config.hasFwdet2()) {
        fw_p2.setFromReader(reader, 2, KinematicType::RECONSTRUCTED);
    }
    if (config.hasFwdet3()) {
        fw_p3.setFromReader(reader, 3, KinematicType::RECONSTRUCTED);
    }

    // --- ECAL photon candidates (if enabled in config) ---
    if (config.hasEcal1()) {
        gamma1.setFromReader(reader, 1, KinematicType::RECONSTRUCTED);
    }
    if (config.hasEcal2()) {
        gamma2.setFromReader(reader, 2, KinematicType::RECONSTRUCTED);
    }
    if (config.hasEcal3()) {
        gamma3.setFromReader(reader, 3, KinematicType::RECONSTRUCTED);
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

    // ========================================================================
    // DILEPTON (e+e-) - Main physics observable
    // ========================================================================

    // Dilepton = e+ + e- (virtual photon)
    PParticle dilepton = positron + electron;

    // ========================================================================
    // 3. QUALITY HISTOGRAMS (before cuts)
    // ========================================================================

    double m_ep = positron.massGeV(analysis_type);
    double m_em = electron.massGeV(analysis_type);
    double m_ee = dilepton.massGeV(analysis_type);  // Dilepton invariant mass

    mgr.fillw("mass_ep", m_ep, weight);
    mgr.fillw("mass_em", m_em, weight);
    mgr.fillw("mass_ee", m_ee, weight);

    // ========================================================================
    // 4. APPLY CUTS
    // ========================================================================

    // Dilepton mass cut (e.g., to select specific mass regions)
    if (cuts.hasRangeCut("dilepton_mass")) {
        if (!cuts.passRangeCut("dilepton_mass", m_ee)) return;
    }

    // Fill after dilepton mass cut
    mgr.fillw("mass_ee_cut", m_ee, weight);

    // Opening angle cut (optional)
    double oa_epem = positron.openingAngle(electron, analysis_type);
    if (cuts.hasRangeCut("opening_angle")) {
        if (!cuts.passRangeCut("opening_angle", oa_epem)) return;
    }

    // ========================================================================
    // 5. BOOST TO CMS
    // ========================================================================

    const BoostFrame& beam_frame = frames.getFrame("beam");

    PParticle ep_cms = beam_frame.boost(positron);
    PParticle em_cms = beam_frame.boost(electron);
    PParticle ee_cms = beam_frame.boost(dilepton);

    // ========================================================================
    // 6. FILL HISTOGRAMS
    // ========================================================================

    // Dilepton mass
    mgr.fillw("mass_dilepton", m_ee, weight);

    // LAB frame kinematics
    mgr.fillw("ep_p_lab", positron.momentum(analysis_type), weight);
    mgr.fillw("em_p_lab", electron.momentum(analysis_type), weight);
    mgr.fillw("ee_p_lab", dilepton.momentum(analysis_type), weight);

    mgr.fillw("ep_theta_lab", positron.theta(analysis_type), weight);
    mgr.fillw("em_theta_lab", electron.theta(analysis_type), weight);
    mgr.fillw("ee_theta_lab", dilepton.theta(analysis_type), weight);

    // CMS kinematics
    mgr.fillw("cos_theta_ee_cms", ee_cms.cosTheta(analysis_type), weight);
    mgr.fillw("cos_theta_ep_cms", ep_cms.cosTheta(analysis_type), weight);
    mgr.fillw("cos_theta_em_cms", em_cms.cosTheta(analysis_type), weight);

    mgr.fillw("ep_p_cms", ep_cms.momentum(analysis_type), weight);
    mgr.fillw("em_p_cms", em_cms.momentum(analysis_type), weight);
    mgr.fillw("ee_p_cms", ee_cms.momentum(analysis_type), weight);

    // Opening angle (e+ - e-)
    mgr.fillw("oa_epem", oa_epem, weight);

    // 2D correlations
    mgr.fillw("mass_vs_costh_ee", m_ee, ee_cms.cosTheta(analysis_type), weight);
    mgr.fillw("theta_ep_vs_em_lab", electron.theta(analysis_type), positron.theta(analysis_type), weight);
    mgr.fillw("p_ep_vs_em_lab", electron.momentum(analysis_type), positron.momentum(analysis_type), weight);

    // Rapidity distribution
    mgr.fillw("ee_rapidity", dilepton.rapidity(analysis_type), weight);

    // Transverse momentum
    double pt_ee = sqrt(dilepton.px(analysis_type)*dilepton.px(analysis_type) +
                        dilepton.py(analysis_type)*dilepton.py(analysis_type));
    mgr.fillw("ee_pt", pt_ee, weight);

    // ========================================================================
    // 7. ECAL PHOTON PROCESSING
    // ========================================================================

    // Count photon multiplicity
    int n_gamma = 0;
    int n_gamma_good = 0;

    // Store good photons for e+e-γ combinations
    std::vector<PParticleEcal*> good_gammas;

    // Process gamma1
    if (gamma1.isValid()) {
        n_gamma++;
        mgr.fillw("gamma_energy", gamma1.getEnergy(), weight);
        mgr.fillw("gamma_theta", gamma1.ecal_theta, weight);
        mgr.fillw("gamma_phi", gamma1.ecal_phi, weight);
        mgr.fillw("gamma1_energy", gamma1.getEnergy(), weight);

        if (isGoodPhoton(gamma1, cuts)) {
            n_gamma_good++;
            mgr.fillw("gamma_energy_good", gamma1.getEnergy(), weight);
            good_gammas.push_back(&gamma1);
        }
    }

    // Process gamma2
    if (gamma2.isValid()) {
        n_gamma++;
        mgr.fillw("gamma_energy", gamma2.getEnergy(), weight);
        mgr.fillw("gamma_theta", gamma2.ecal_theta, weight);
        mgr.fillw("gamma_phi", gamma2.ecal_phi, weight);
        mgr.fillw("gamma2_energy", gamma2.getEnergy(), weight);

        if (isGoodPhoton(gamma2, cuts)) {
            n_gamma_good++;
            mgr.fillw("gamma_energy_good", gamma2.getEnergy(), weight);
            good_gammas.push_back(&gamma2);
        }
    }

    // Process gamma3
    if (gamma3.isValid()) {
        n_gamma++;
        mgr.fillw("gamma_energy", gamma3.getEnergy(), weight);
        mgr.fillw("gamma_theta", gamma3.ecal_theta, weight);
        mgr.fillw("gamma_phi", gamma3.ecal_phi, weight);
        mgr.fillw("gamma3_energy", gamma3.getEnergy(), weight);

        if (isGoodPhoton(gamma3, cuts)) {
            n_gamma_good++;
            mgr.fillw("gamma_energy_good", gamma3.getEnergy(), weight);
            good_gammas.push_back(&gamma3);
        }
    }

    // Fill multiplicity histograms
    mgr.fillw("gamma_mult", n_gamma, weight);
    mgr.fillw("gamma_mult_good", n_gamma_good, weight);

    // ========================================================================
    // 8. GAMMA QUALITY HISTOGRAMS
    // ========================================================================

    // Fill gamma quality histograms for all valid gammas
    auto fillGammaQuality = [&](const PParticleEcal& g, double w) {
        if (!g.isValid()) return;
        if (g.ecal_beta > 0) mgr.fillw("gamma_beta", g.ecal_beta, w);
        if (g.ecal_chi2 > 0) mgr.fillw("gamma_chi2", g.ecal_chi2, w);
        if (g.ecal_dist > 0) mgr.fillw("gamma_dist", g.ecal_dist, w);
        if (g.ecal_r > 0 && g.ecal_z > 0) mgr.fillw("gamma_rz", g.ecal_r, g.ecal_z, w);
    };

    fillGammaQuality(gamma1, weight);
    fillGammaQuality(gamma2, weight);
    fillGammaQuality(gamma3, weight);

    // ========================================================================
    // 9. GAMMA-LEPTON OPENING ANGLES
    // ========================================================================

    // Calculate opening angles between gammas and leptons
    for (const auto* g : good_gammas) {
        double oa_g_ep = positron.openingAngle(*g, analysis_type);
        double oa_g_em = electron.openingAngle(*g, analysis_type);
        double oa_g_ee = dilepton.openingAngle(*g, analysis_type);
        double oa_g_min = std::min(oa_g_ep, oa_g_em);

        mgr.fillw("oa_g_ep", oa_g_ep, weight);
        mgr.fillw("oa_g_em", oa_g_em, weight);
        mgr.fillw("oa_g_ee", oa_g_ee, weight);
        mgr.fillw("oa_g_lepton_min", oa_g_min, weight);
    }

    // ========================================================================
    // 10. PI0 RECONSTRUCTION (gamma-gamma)
    // ========================================================================

    double m_gg_best = -1.0;
    double oa_gg_best = -1.0;
    double pi0_energy_best = -1.0;
    double pi0_theta_best = -1.0;

    // Try all gamma-gamma combinations
    if (good_gammas.size() >= 2) {
        for (size_t i = 0; i < good_gammas.size(); ++i) {
            for (size_t j = i + 1; j < good_gammas.size(); ++j) {
                const PParticleEcal& g1 = *good_gammas[i];
                const PParticleEcal& g2 = *good_gammas[j];

                // Create pi0 candidate
                PParticle pi0_cand = g1 + g2;
                double m_gg = pi0_cand.massGeV(analysis_type);
                double oa_gg = g1.openingAngle(g2, analysis_type);
                double pi0_E = pi0_cand.energy(analysis_type);
                double pi0_th = pi0_cand.theta(analysis_type);

                // Fill all combinations
                mgr.fillw("mass_gg", m_gg, weight);
                mgr.fillw("mass_gg_pi0", m_gg, weight);
                mgr.fillw("oa_gg", oa_gg, weight);
                mgr.fillw("pi0_energy", pi0_E, weight);
                mgr.fillw("pi0_theta", pi0_th, weight);

                // Keep best pi0 candidate (closest to pi0 mass)
                if (m_gg_best < 0 || std::abs(m_gg - 0.135) < std::abs(m_gg_best - 0.135)) {
                    m_gg_best = m_gg;
                    oa_gg_best = oa_gg;
                    pi0_energy_best = pi0_E;
                    pi0_theta_best = pi0_th;
                }
            }
        }
    }

    // ========================================================================
    // 11. e+e-γ INVARIANT MASS COMBINATIONS
    // ========================================================================

    // e+e- + gamma1 (if valid and passes cuts)
    double m_eeg1 = -1.0;
    if (gamma1.isValid() && isGoodPhoton(gamma1, cuts)) {
        PParticle eeg1 = dilepton + gamma1;
        m_eeg1 = eeg1.massGeV(analysis_type);
        mgr.fillw("mass_eeg1", m_eeg1, weight);
        mgr.fillw("mass_eegamma", m_eeg1, weight);
        mgr.fillw("mass_eegamma_pi0", m_eeg1, weight);
        mgr.fillw("mass_eegamma_eta", m_eeg1, weight);
        mgr.fillw("mass_ee_vs_eegamma", m_ee, m_eeg1, weight);
        mgr.fillw("mass_eegamma_vs_gamma_energy", m_eeg1, gamma1.getEnergy(), weight);
    }

    // e+e- + gamma2 (if valid and passes cuts)
    double m_eeg2 = -1.0;
    if (gamma2.isValid() && isGoodPhoton(gamma2, cuts)) {
        PParticle eeg2 = dilepton + gamma2;
        m_eeg2 = eeg2.massGeV(analysis_type);
        mgr.fillw("mass_eeg2", m_eeg2, weight);
        mgr.fillw("mass_eegamma", m_eeg2, weight);
        mgr.fillw("mass_eegamma_pi0", m_eeg2, weight);
        mgr.fillw("mass_eegamma_eta", m_eeg2, weight);
        mgr.fillw("mass_ee_vs_eegamma", m_ee, m_eeg2, weight);
        mgr.fillw("mass_eegamma_vs_gamma_energy", m_eeg2, gamma2.getEnergy(), weight);
    }

    // e+e- + gamma3 (if valid and passes cuts)
    double m_eeg3 = -1.0;
    if (gamma3.isValid() && isGoodPhoton(gamma3, cuts)) {
        PParticle eeg3 = dilepton + gamma3;
        m_eeg3 = eeg3.massGeV(analysis_type);
        mgr.fillw("mass_eeg3", m_eeg3, weight);
        mgr.fillw("mass_eegamma", m_eeg3, weight);
        mgr.fillw("mass_eegamma_pi0", m_eeg3, weight);
        mgr.fillw("mass_eegamma_eta", m_eeg3, weight);
        mgr.fillw("mass_ee_vs_eegamma", m_ee, m_eeg3, weight);
        mgr.fillw("mass_eegamma_vs_gamma_energy", m_eeg3, gamma3.getEnergy(), weight);
    }

    // ========================================================================
    // 12. PWA VARIABLES (in dilepton rest frame)
    // ========================================================================

    // Create e+e- rest frame (dilepton rest frame)
    BoostFrame ee_frame(dilepton);

    PParticle ep_in_ee = ee_frame.boost(positron);
    PParticle em_in_ee = ee_frame.boost(electron);
    PParticle proj_in_ee = ee_frame.boost(projectile);

    // Helicity angle: positron angle relative to beam direction in dilepton frame
    mgr.fillw("pwa_ep_helicity_ee", ep_in_ee.cosTheta(analysis_type), weight);
    mgr.fillw("pwa_em_helicity_ee", em_in_ee.cosTheta(analysis_type), weight);

    // Gottfried-Jackson: angle relative to beam in dilepton frame
    double gj_angle = ep_in_ee.vec(analysis_type).Angle(proj_in_ee.vec(analysis_type).Vect());
    mgr.fillw("pwa_ep_gj_ee", cos(gj_angle), weight);

    // ========================================================================
    // 13. FILL OUTPUT NTUPLES
    // ========================================================================

    // --- Ntuple 1: Basic particle observables ---
    DynamicHNtuple& nt_particles = mgr.getDynamicNtuple("nt_particles");

    // Positron observables (LAB)
    nt_particles["ep_p"] = positron.momentum(analysis_type);
    nt_particles["ep_theta"] = positron.theta(analysis_type);
    nt_particles["ep_phi"] = positron.phi(analysis_type);
    nt_particles["ep_mass"] = m_ep;

    // Electron observables (LAB)
    nt_particles["em_p"] = electron.momentum(analysis_type);
    nt_particles["em_theta"] = electron.theta(analysis_type);
    nt_particles["em_phi"] = electron.phi(analysis_type);
    nt_particles["em_mass"] = m_em;

    // Dilepton observables
    nt_particles["ee_p"] = dilepton.momentum(analysis_type);
    nt_particles["ee_theta"] = dilepton.theta(analysis_type);
    nt_particles["ee_phi"] = dilepton.phi(analysis_type);
    nt_particles["ee_mass"] = m_ee;

    // Opening angle
    nt_particles["oa_epem"] = oa_epem;

    // Photon multiplicity
    nt_particles["n_gamma"] = n_gamma;
    nt_particles["n_gamma_good"] = n_gamma_good;

    // Photon observables (gamma1)
    if (gamma1.isValid()) {
        nt_particles["g1_energy"] = gamma1.getEnergy();
        nt_particles["g1_theta"] = gamma1.ecal_theta;
        nt_particles["g1_phi"] = gamma1.ecal_phi;
        nt_particles["g1_good"] = isGoodPhoton(gamma1, cuts) ? 1.0 : 0.0;
    }

    // Photon observables (gamma2)
    if (gamma2.isValid()) {
        nt_particles["g2_energy"] = gamma2.getEnergy();
        nt_particles["g2_theta"] = gamma2.ecal_theta;
        nt_particles["g2_phi"] = gamma2.ecal_phi;
        nt_particles["g2_good"] = isGoodPhoton(gamma2, cuts) ? 1.0 : 0.0;
    }

    // Photon observables (gamma3)
    if (gamma3.isValid()) {
        nt_particles["g3_energy"] = gamma3.getEnergy();
        nt_particles["g3_theta"] = gamma3.ecal_theta;
        nt_particles["g3_phi"] = gamma3.ecal_phi;
        nt_particles["g3_good"] = isGoodPhoton(gamma3, cuts) ? 1.0 : 0.0;
    }

    // Event weight
    nt_particles["weight"] = weight;

    nt_particles.fill();

    // --- Ntuple 2: Compound observables ---
    DynamicHNtuple& nt_compound = mgr.getDynamicNtuple("nt_compound");

    // Dilepton mass and kinematics
    nt_compound["m_ee"] = m_ee;
    nt_compound["ee_pt"] = pt_ee;
    nt_compound["ee_rapidity"] = dilepton.rapidity(analysis_type);

    // CMS angles
    nt_compound["cos_th_ee_cms"] = ee_cms.cosTheta(analysis_type);
    nt_compound["cos_th_ep_cms"] = ep_cms.cosTheta(analysis_type);
    nt_compound["cos_th_em_cms"] = em_cms.cosTheta(analysis_type);

    // Opening angle
    nt_compound["oa_epem"] = oa_epem;

    // e+e-γ invariant masses
    nt_compound["m_eeg1"] = m_eeg1;
    nt_compound["m_eeg2"] = m_eeg2;
    nt_compound["m_eeg3"] = m_eeg3;

    // PWA variables (helicity and Gottfried-Jackson angles)
    nt_compound["ep_helicity"] = ep_in_ee.cosTheta(analysis_type);
    nt_compound["ep_gj"] = cos(gj_angle);
    nt_compound["em_helicity"] = em_in_ee.cosTheta(analysis_type);

    // Pi0 reconstruction (gamma-gamma)
    nt_compound["m_gg"] = m_gg_best;
    nt_compound["oa_gg"] = oa_gg_best;
    nt_compound["pi0_energy"] = pi0_energy_best;
    nt_compound["pi0_theta"] = pi0_theta_best;

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

    ConsoleBox::newLine();
    ConsoleBox::printHeader("FAT Framework - Final Analysis Tool",
                           "pp → pp e+e- (dilepton) + ECAL γ Analysis");
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
        std::string tree_name = config.getInputTreeName();

        if (config.isInputFileList()) {
            // File list (.list)
            reader.openFromList(config.getInputSource(), tree_name);
        } else if (config.isInputMultipleRootFiles()) {
            // Multiple ROOT files (comma-separated or JSON array)
            std::vector<std::string> files = config.getInputFiles();
            std::cout << "\nInput: " << files.size() << " ROOT files (chain)\n";
            for (const auto& f : files) {
                std::cout << "  - " << f << "\n";
            }
            reader.openChain(files, tree_name);
        } else if (config.isInputRootFile()) {
            // Single ROOT file
            reader.open(config.getInputSource(), tree_name);
        } else {
            throw std::runtime_error("Unknown input format. Use .root, .list, or comma-separated .root files");
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

// ========================================================================
// FAT Framework - Proton (HADES) + Proton (FT) Analysis
// ========================================================================
// p + p(FT) analysis for HADES experiment (pp @ 1.58 GeV).
// Processes one HADES proton and Forward Tracker proton candidates
// with pi0 as missing particle.
//
// Usage:
//   ./ana [config.json]
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2025
// ========================================================================

#include "src/manager.h"
#include "src/pparticle.h"
#include "src/pparticle_fwd.h"
#include "src/boost_frame.h"
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
#include <vector>

// Use Physics namespace for mass constants
using namespace Physics;

// ============================================================================
// PROCESS SINGLE EVENT
// ============================================================================

void processEvent(NTupleReader& reader, Manager& mgr, CutManager& cuts,
                 const AnalysisConfig& config) {

    // Event-level cuts (applied before particle creation)
    if (!cuts.passValueCut("isBest", reader["isBest"])) return;
    if (!cuts.passMinCut("vertex_z", reader["eVertReco_z"])) return;

    // Create proton1 (HADES) with reconstructed kinematics (p_ prefix)
    PParticle proton1(MASS_PROTON, "p1");

    proton1.setFromSpherical(reader["p_p"], reader["p_theta"], reader["p_phi"],
                             KinematicType::RECONSTRUCTED);

    // Energy-loss corrected kinematics (same angles, corrected momentum)
    proton1.setFromSpherical(reader["p_p_corr_p"], reader["p_theta"], reader["p_phi"],
                             KinematicType::CORRECTED);

    // Beam + target for missing mass calculations
    PParticle beam = ParticleFactory::createBeamProton(config.getBeamKineticEnergy());
    PParticle target = ParticleFactory::createTargetProton();
    PParticle initial = beam + target;

    // Forward Tracker objects (forward hadrons, up to 3 hits)
    if (!config.isFwdEnabled()) return;

    int fwdet_mult = static_cast<int>(reader["fwdet_mult"]);
    if (fwdet_mult < 0) fwdet_mult = 0;

    std::vector<PParticleFwd> fwd_objects;
    fwd_objects.reserve(fwdet_mult);

    for (int i = 1; i <= fwdet_mult && i <= 3; ++i) {
        PParticleFwd fwd_obj(MASS_PROTON, "fwd" + std::to_string(i));
        if (fwd_obj.setFromReader(reader, i) && fwd_obj.fwd_theta > 0) {
            fwd_objects.push_back(fwd_obj);
        }
    }

    // Fill fwd ntuple
    auto& fwd_nt = mgr.getDynamicNtuple("fwdet_nt");

    fwd_nt["fwd_mult"] = fwdet_mult;

    auto fillFwdHit = [&](size_t idx, const std::string& suffix) {
        if (fwd_objects.size() > idx) {
            const auto& obj = fwd_objects[idx];

            fwd_nt["fwd_p" + suffix] = obj.momentum();
            fwd_nt["fwd_theta" + suffix] = obj.fwd_theta;
            fwd_nt["fwd_phi" + suffix] = obj.fwd_phi;

            fwd_nt["fwd_beta" + suffix] = obj.fwd_beta;
            fwd_nt["fwd_mass" + suffix] = obj.fwd_mass;
            fwd_nt["fwd_mass2" + suffix] = obj.fwd_mass2;
            fwd_nt["fwd_charge" + suffix] = obj.fwd_charge;

            fwd_nt["fwd_chi2" + suffix] = obj.fwd_chi2;
            fwd_nt["fwd_ndf" + suffix] = obj.fwd_ndf;
            fwd_nt["fwd_chi2ndf" + suffix] = obj.fwd_chi2ndf;

            fwd_nt["fwd_r" + suffix] = obj.fwd_r;
            fwd_nt["fwd_z" + suffix] = obj.fwd_z;
            fwd_nt["fwd_distToRpc" + suffix] = obj.fwd_distToRpc;

            fwd_nt["fwd_tof" + suffix] = obj.fwd_tof;
        }
    };

    fillFwdHit(0, "_1");
    fillFwdHit(1, "_2");
    fillFwdHit(2, "_3");

    fwd_nt.fill();

    // Setup CMS frame
    EventFrames frames;
    frames.setBeamFrameFromKineticEnergy(config.getBeamKineticEnergy());

    // Loop over FWD proton candidates (proton2)
    auto& pp_nt = mgr.getDynamicNtuple("pp_nt");

    for (size_t j = 0; j < fwd_objects.size(); ++j) {
        // FWD proton beta quality cut
        /////if (!cuts.passMinCut("proton_fwd_cut", fwd_objects[j].fwd_beta)) continue;

        PParticle proton2 = fwd_objects[j];  // FT proton candidate

        // Anti-elastic cut (reject elastic events)
        double tantan = TMath::Tan(reader["p_theta"] * TMath::DegToRad()) *
                        TMath::Tan(proton2.theta() * TMath::DegToRad());
        double dphi = TMath::Abs(reader["p_phi"] - proton2.phi());
        if (cuts.passGraphicalCut("anti_elastic", dphi, tantan)) continue;  // reject if inside

        // Missing mass cuts (reject events with MM(p) < 1.05 GeV)
        double mm_p1 = (initial - proton1).massGeV();  // MM of proton1
        double mm_p2 = (initial - proton2).massGeV();  // MM of proton2
        if (!cuts.passMinCut("mm_p1", mm_p1)) continue;
        if (!cuts.passMinCut("mm_p2", mm_p2)) continue;

        // Build compound: proton1 + proton2
        PParticle pp = proton1 + proton2;
        double pp_mass = pp.massGeV();

        // Boost to CMS frame
        PParticle pp_cms = frames.getFrame("beam").boost(pp);

        // Missing mass: MM(pp) = beam + target - proton1 - proton2 (pi0)
        PParticle mm_pp = initial - proton1 - proton2;
        double mm_pp_mass = mm_pp.massGeV();

        // Boost missing mass to CMS frame
        PParticle mm_pp_cms = frames.getFrame("beam").boost(mm_pp);

        // Delta systems (proton + pi0)
        PParticle deltaP1 = initial - proton2;  // p1 + pi0
        PParticle deltaP2 = initial - proton1;  // p2 + pi0

        // Fill pp ntuple - LAB frame quantities
        pp_nt["pp_mass"] = pp_mass;
        pp_nt["pp_p"] = pp.momentum();
        pp_nt["pp_theta"] = pp.theta();
        pp_nt["pp_phi"] = pp.phi();
        pp_nt["pp_rapidity"] = pp.rapidity();
        pp_nt["pp_pt"] = pp.vec().Pt();

        // CMS frame quantities (with _cms suffix)
        pp_nt["pp_mass_cms"] = pp_cms.massGeV();
        pp_nt["pp_p_cms"] = pp_cms.momentum();
        pp_nt["pp_theta_cms"] = pp_cms.theta();
        pp_nt["pp_phi_cms"] = pp_cms.phi();
        pp_nt["pp_rapidity_cms"] = pp_cms.rapidity();
        pp_nt["pp_pt_cms"] = pp_cms.vec().Pt();
        pp_nt["pp_costheta_cms"] = pp_cms.cosTheta();

        // Missing mass - LAB frame (pi0)
        pp_nt["mm_pp_mass"] = mm_pp_mass;
        pp_nt["mm_pp_p"] = mm_pp.momentum();
        pp_nt["mm_pp_theta"] = mm_pp.theta();
        pp_nt["mm_pp_phi"] = mm_pp.phi();
        pp_nt["mm_pp_rapidity"] = mm_pp.rapidity();
        pp_nt["mm_pp_pt"] = mm_pp.vec().Pt();

        // Missing mass - CMS frame
        pp_nt["mm_pp_mass_cms"] = mm_pp_cms.massGeV();
        pp_nt["mm_pp_p_cms"] = mm_pp_cms.momentum();
        pp_nt["mm_pp_theta_cms"] = mm_pp_cms.theta();
        pp_nt["mm_pp_phi_cms"] = mm_pp_cms.phi();
        pp_nt["mm_pp_rapidity_cms"] = mm_pp_cms.rapidity();
        pp_nt["mm_pp_pt_cms"] = mm_pp_cms.vec().Pt();
        pp_nt["mm_pp_costheta_cms"] = mm_pp_cms.cosTheta();

        // Proton1 info (HADES)
        pp_nt["p1_p_rec"] = proton1.momentum(KinematicType::RECONSTRUCTED);
        pp_nt["p1_p_cor"] = proton1.momentum(KinematicType::CORRECTED);
        pp_nt["p1_theta"] = proton1.theta();
        pp_nt["p1_phi"] = proton1.phi();

        // Proton2 info (FT)
        pp_nt["p2_p"] = proton2.momentum();
        pp_nt["p2_theta"] = proton2.theta();
        pp_nt["p2_phi"] = proton2.phi();
        pp_nt["fwd_beta"] = fwd_objects[j].fwd_beta;
        pp_nt["fwd_mass2"] = fwd_objects[j].fwd_mass2;
        pp_nt["fwd_tof"] = fwd_objects[j].fwd_tof;
        pp_nt["fwd_index"] = static_cast<Float_t>(j + 1);
        pp_nt["fwd_mult"] = static_cast<Float_t>(fwdet_mult);

        // Check pi0 mass cut and store result
        bool passes_pion0_cut = cuts.passRangeCut("pion0_mass_cut", mm_pp_mass);
        pp_nt["pion0_mass_cut"] = passes_pion0_cut ? 1 : 0;

        pp_nt.fill();

        // Fill histograms
        // ppi0 invariant mass (delta system) - fill twice with 0.5 weight for symmetry
        mgr.fillw("ppi0_inv_mass", deltaP1.massGeV(), 0.5);
        mgr.fillw("ppi0_inv_mass", deltaP2.massGeV(), 0.5);

        // pp invariant mass
        mgr.fill("pp_inv_mass", pp_mass);

        // pi0 missing mass
        mgr.fill("pi0_miss_mass", mm_pp_mass);

        // pi0 missing mass squared (from 4-vector, can be negative)
        mgr.fill("pi0_miss_mass2", mm_pp.vec().M2() / 1e6);

        // 2D debug: dphi vs tantan
        mgr.fill("dphi_vs_tantan", dphi, tantan);

        // 2D debug: missing mass p1 vs missing mass p2
        // deltaP2 = initial - proton1 (MM of p1), deltaP1 = initial - proton2 (MM of p2)
        mgr.fill("mm_p1_vs_mm_p2", deltaP2.massGeV(), deltaP1.massGeV());

        // Pi0 mass cut - mandatory for PWA histogram filling
        if (!passes_pion0_cut) continue;

        // ================================================================
        // PWA (Partial Wave Analysis)
        // ================================================================

        // Save LAB frame copies
        PParticle p1_LAB = proton1;
        PParticle p2_LAB = proton2;
        PParticle pi0_LAB = mm_pp;  // pi0 is the missing mass

        // Boost to CMS frame (beam direction)
        PParticle p1_CMS = frames.getFrame("beam").boost(proton1);
        PParticle p2_CMS = frames.getFrame("beam").boost(proton2);
        PParticle pi0_CMS = frames.getFrame("beam").boost(mm_pp);
        PParticle deltaP1_CMS = frames.getFrame("beam").boost(deltaP1);
        PParticle deltaP2_CMS = frames.getFrame("beam").boost(deltaP2);

        // 2D: M(ppi0) vs cos_theta_cms - fill twice with 0.5 weight for symmetry
        mgr.fillw("mppi0_vs_costh_cms", deltaP1.massGeV(), deltaP1_CMS.cosTheta(), 0.5);
        mgr.fillw("mppi0_vs_costh_cms", deltaP2.massGeV(), deltaP2_CMS.cosTheta(), 0.5);

        // Group A: cos(theta) in CMS - fill twice with 0.5 weight for symmetry
        mgr.fillw("pwa_pi0_costh", pi0_CMS.cosTheta(), 0.5);
        mgr.fillw("pwa_pi0_costh", pi0_CMS.cosTheta(), 0.5);
        mgr.fillw("pwa_p_costh", p1_CMS.cosTheta(), 0.5);
        mgr.fillw("pwa_p_costh", p2_CMS.cosTheta(), 0.5);

        // Group B: Momenta in LAB frame - fill twice with 0.5 weight
        mgr.fillw("pwa_pi0_p", pi0_LAB.momentum() / 1000.0, 0.5);  // Convert MeV to GeV
        mgr.fillw("pwa_pi0_p", pi0_LAB.momentum() / 1000.0, 0.5);
        mgr.fillw("pwa_p_p", p1_LAB.momentum() / 1000.0, 0.5);
        mgr.fillw("pwa_p_p", p2_LAB.momentum() / 1000.0, 0.5);

        // Group C: Invariant masses
        mgr.fillw("pwa_ppi0_m", deltaP1.massGeV(), 0.5);  // p1 + pi0
        mgr.fillw("pwa_ppi0_m", deltaP2.massGeV(), 0.5);  // p2 + pi0
        mgr.fillw("pwa_pp_m", pp.massGeV(), 0.5);
        mgr.fillw("pwa_pp_m", pp.massGeV(), 0.5);

        // Helicity frames (boost particles to rest frame of parent)
        // For p + pi0 systems (deltaP1 and deltaP2)
        BoostFrame deltaP1_frame(deltaP1);  // p1 + pi0 rest frame
        BoostFrame deltaP2_frame(deltaP2);  // p2 + pi0 rest frame
        BoostFrame pp_frame(pp);            // p1 + p2 rest frame

        // Boost p1 to deltaP2 rest frame (p2 + pi0 system)
        PParticle p1_P2PI0 = deltaP2_frame.boost(proton1);

        // Boost p2 to deltaP1 rest frame (p1 + pi0 system)
        PParticle p2_P1PI0 = deltaP1_frame.boost(proton2);

        // Boost pi0 to both delta frames
        PParticle pi0_P1PI0 = deltaP1_frame.boost(mm_pp);
        PParticle pi0_P2PI0 = deltaP2_frame.boost(mm_pp);

        // Boost particles to pp rest frame
        PParticle p1_pp = pp_frame.boost(proton1);
        PParticle p2_pp = pp_frame.boost(proton2);
        PParticle pi0_pp = pp_frame.boost(mm_pp);

        // Group D: Helicity distributions (opening angles in rest frames)
        // pi0 helicity in p+pi0 system - fill twice with 0.25 weight (0.5/2)
        double helicity_pi0_1 = Physics::openingAngle(pi0_P1PI0, p2_P1PI0);
        double helicity_pi0_2 = Physics::openingAngle(pi0_P2PI0, p1_P2PI0);
        mgr.fillw("pwa_pi0_helicity", cos(helicity_pi0_1 * M_PI / 180.0), 0.25);
        mgr.fillw("pwa_pi0_helicity", cos(helicity_pi0_2 * M_PI / 180.0), 0.25);

        // p helicity in pp+pi0 system
        double helicity_p1 = Physics::openingAngle(p1_pp, pi0_pp);
        double helicity_p2 = Physics::openingAngle(p2_pp, pi0_pp);
        mgr.fillw("pwa_p_helicity", cos(helicity_p1 * M_PI / 180.0), 0.5);
        mgr.fillw("pwa_p_helicity", cos(helicity_p2 * M_PI / 180.0), 0.5);

        // Gottfried-Jackson frames (with projectile)
        // Boost projectile (beam) to same rest frames
        PParticle proj_P1PI0 = deltaP1_frame.boost(beam);
        PParticle proj_P2PI0 = deltaP2_frame.boost(beam);
        PParticle proj_pp = pp_frame.boost(beam);

        // Group E: GJ distributions (angles with respect to beam in rest frames)
        // pi0 GJ angle - fill twice with 0.25 weight
        double gj_pi0_1 = Physics::openingAngle(pi0_P1PI0, proj_P1PI0);
        double gj_pi0_2 = Physics::openingAngle(pi0_P2PI0, proj_P2PI0);
        mgr.fillw("pwa_pi0_gj", cos(gj_pi0_1 * M_PI / 180.0), 0.25);
        mgr.fillw("pwa_pi0_gj", cos(gj_pi0_2 * M_PI / 180.0), 0.25);

        // p GJ angle
        double gj_p1 = Physics::openingAngle(p1_pp, proj_pp);
        double gj_p2 = Physics::openingAngle(p2_pp, proj_pp);
        mgr.fillw("pwa_p_gj", cos(gj_p1 * M_PI / 180.0), 0.5);
        mgr.fillw("pwa_p_gj", cos(gj_p2 * M_PI / 180.0), 0.5);

        // Additional histograms (D+)
        // Delta (p+pi0) mass and cos(theta)
        mgr.fillw("mass_deltaP", deltaP1.massGeV(), 0.5);
        mgr.fillw("mass_deltaP", deltaP2.massGeV(), 0.5);

        mgr.fillw("cos_theta_deltaP", deltaP1_CMS.cosTheta(), 0.5);
        mgr.fillw("cos_theta_deltaP", deltaP2_CMS.cosTheta(), 0.5);

        // pi0 mass and cos(theta)
        mgr.fillw("cos_theta_pi0", pi0_CMS.cosTheta(), 0.5);
        mgr.fillw("cos_theta_pi0", pi0_CMS.cosTheta(), 0.5);

        mgr.fillw("mass_pi0", mm_pp.massGeV(), 0.5);
        mgr.fillw("mass_pi0", mm_pp.massGeV(), 0.5);

    }  // End FWD proton loop

}

// ============================================================================
// MAIN PROGRAM
// ============================================================================

int main(int argc, char* argv[]) {
    // Install signal handler for graceful Ctrl+C termination
    SignalHandler::install();

    ConsoleBox::newLine();
    ConsoleBox::printHeader("FAT Framework",
                           "Proton (HADES) + Proton (FT) Analysis");
    ConsoleBox::newLine();

    // Load configuration
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

    // Open input data
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

    // Open output file and setup analysis components
    Manager manager;
    manager.openFile(config.getOutputFilename(), config.getOutputOption());

    setupHistograms(manager);
    setupNtuples(manager, config);

    CutManager cuts;
    setupCuts(cuts);

    // Event loop
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

        try {
            processEvent(reader, manager, cuts, config);
        } catch (const std::exception& e) {
            continue;
        }
    }

    progress.finish(was_interrupted);

    // Finalization
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

// ========================================================================
// FAT Framework - Pion Analysis
// ========================================================================
// π+ pion analysis for HADES experiment (pp @ 4.5 GeV).
// Processes positive pions and Forward Tracker hadrons.
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

    // Create π+ (positive pion) with reconstructed kinematics
    PParticle pion(MASS_PION_PLUS, "pi+");

    pion.setFromSpherical(reader["pip_p"], reader["pip_theta"], reader["pip_phi"],
                          KinematicType::RECONSTRUCTED);

    // Energy-loss corrected kinematics (same angles, corrected momentum)
    pion.setFromSpherical(reader["pip_p_corr_pip"], reader["pip_theta"], reader["pip_phi"],
                          KinematicType::CORRECTED);

    // Setup CMS frame (beam + target)
    EventFrames frames;
    frames.setBeamFrameFromKineticEnergy(config.getBeamKineticEnergy());

    // Beam + target for missing mass calculations
    PParticle beam = ParticleFactory::createBeamProton(config.getBeamKineticEnergy());
    PParticle target = ParticleFactory::createTargetProton();
    PParticle initial = beam + target;

    // Forward Tracker objects (forward hadrons, up to 3 hits)
    if (config.isFwdEnabled()) {
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

        // Compound ppip objects (one entry per FWD proton candidate)
        auto& ppip_nt = mgr.getDynamicNtuple("ppip_nt");

        for (size_t j = 0; j < fwd_objects.size(); ++j) {
            // Apply neutron cut (FWD beta quality)
            if (!cuts.passMinCut("neutron_cut", fwd_objects[j].fwd_beta)) continue;

            // Build compound: FWD proton + pi+
            PParticle ppip = fwd_objects[j] + pion;
            double ppip_mass = ppip.massGeV();

            // Boost to CMS frame
            PParticle ppip_cms = frames.getFrame("beam").boost(ppip);

            // Missing mass: MM(ppip) = beam + target - proton - pi+
            PParticle mm_ppip = initial - fwd_objects[j] - pion;
            double mm_ppip_mass = mm_ppip.massGeV();

            // Boost missing mass to CMS frame
            PParticle mm_ppip_cms = frames.getFrame("beam").boost(mm_ppip);

            // Fill ppip ntuple - LAB frame quantities
            ppip_nt["ppip_mass"] = ppip_mass;
            ppip_nt["ppip_p"] = ppip.momentum();
            ppip_nt["ppip_theta"] = ppip.theta();
            ppip_nt["ppip_phi"] = ppip.phi();
            ppip_nt["ppip_rapidity"] = ppip.rapidity();
            ppip_nt["ppip_pt"] = ppip.vec().Pt();

            // CMS frame quantities (with _cms suffix)
            ppip_nt["ppip_mass_cms"] = ppip_cms.massGeV();
            ppip_nt["ppip_p_cms"] = ppip_cms.momentum();
            ppip_nt["ppip_theta_cms"] = ppip_cms.theta();
            ppip_nt["ppip_phi_cms"] = ppip_cms.phi();
            ppip_nt["ppip_rapidity_cms"] = ppip_cms.rapidity();
            ppip_nt["ppip_pt_cms"] = ppip_cms.vec().Pt();
            ppip_nt["ppip_costheta_cms"] = ppip_cms.cosTheta();

            // Missing mass - LAB frame
            ppip_nt["mm_ppip_mass"] = mm_ppip_mass;
            ppip_nt["mm_ppip_p"] = mm_ppip.momentum();
            ppip_nt["mm_ppip_theta"] = mm_ppip.theta();
            ppip_nt["mm_ppip_phi"] = mm_ppip.phi();
            ppip_nt["mm_ppip_rapidity"] = mm_ppip.rapidity();
            ppip_nt["mm_ppip_pt"] = mm_ppip.vec().Pt();

            // Missing mass - CMS frame
            ppip_nt["mm_ppip_mass_cms"] = mm_ppip_cms.massGeV();
            ppip_nt["mm_ppip_p_cms"] = mm_ppip_cms.momentum();
            ppip_nt["mm_ppip_theta_cms"] = mm_ppip_cms.theta();
            ppip_nt["mm_ppip_phi_cms"] = mm_ppip_cms.phi();
            ppip_nt["mm_ppip_rapidity_cms"] = mm_ppip_cms.rapidity();
            ppip_nt["mm_ppip_pt_cms"] = mm_ppip_cms.vec().Pt();
            ppip_nt["mm_ppip_costheta_cms"] = mm_ppip_cms.cosTheta();

            // FWD proton info
            ppip_nt["fwd_p"] = fwd_objects[j].momentum();
            ppip_nt["fwd_theta"] = fwd_objects[j].fwd_theta;
            ppip_nt["fwd_phi"] = fwd_objects[j].fwd_phi;
            ppip_nt["fwd_beta"] = fwd_objects[j].fwd_beta;
            ppip_nt["fwd_mass2"] = fwd_objects[j].fwd_mass2;
            ppip_nt["fwd_tof"] = fwd_objects[j].fwd_tof;
            ppip_nt["fwd_index"] = static_cast<Float_t>(j + 1);
            ppip_nt["fwd_mult"] = static_cast<Float_t>(fwdet_mult);

            // Pion info
            ppip_nt["pip_p_rec"] = pion.momentum(KinematicType::RECONSTRUCTED);
            ppip_nt["pip_p_cor"] = pion.momentum(KinematicType::CORRECTED);
            ppip_nt["pip_theta"] = pion.theta();
            ppip_nt["pip_phi"] = pion.phi();

            ppip_nt.fill();

            // Fill histograms (with fwd_time cut)
            if (/*true ||*/cuts.passMaxCut("fwd_time", fwd_objects[j].fwd_tof)) {
                mgr.fill("ppip_inv_mass", ppip_mass);
                mgr.fill("ppip_miss_mass", mm_ppip_mass);

                // ================================================================
                // PWA (Partial Wave Analysis)
                // ================================================================
                // Apply neutron mass cut for PWA analysis
                if (/*true ||*/ cuts.passRangeCut("neutron_mass_cut", mm_ppip_mass)) {

                    // Save LAB frame copies
                    PParticle p_LAB = fwd_objects[j];
                PParticle pip_LAB = pion;
                PParticle n_LAB = mm_ppip;  // neutron is the missing mass

                // Additional compound systems
                PParticle npip = mm_ppip + pion;  // n + pip
                PParticle pn = fwd_objects[j] + mm_ppip;  // p + n

                // Boost to CMS frame for Group A (cos theta distributions)
                PParticle p_CMS = frames.getFrame("beam").boost(fwd_objects[j]);
                PParticle pip_CMS = frames.getFrame("beam").boost(pion);
                PParticle n_CMS = frames.getFrame("beam").boost(mm_ppip);

                // Group A: cos(theta) in CMS
                mgr.fill("pwa_pip_costh", pip_CMS.cosTheta());
                mgr.fill("pwa_p_costh", p_CMS.cosTheta());
                mgr.fill("pwa_n_costh", n_CMS.cosTheta());

                // Group B: Momenta in LAB frame
                mgr.fill("pwa_pip_p", pip_LAB.momentum() / 1000.0);  // Convert MeV to GeV
                mgr.fill("pwa_p_p", p_LAB.momentum() / 1000.0);
                mgr.fill("pwa_n_p", n_LAB.momentum() / 1000.0);

                // Group C: Invariant masses
                mgr.fill("pwa_ppip_m", ppip.massGeV());
                mgr.fill("pwa_npip_m", npip.massGeV());
                mgr.fill("pwa_pn_m", pn.massGeV());

                // Helicity frames (boost particles to rest frame of parent)
                // Create boost frames for each compound system
                BoostFrame ppip_frame(ppip);  // p + pip rest frame
                BoostFrame npip_frame(npip);  // n + pip rest frame
                BoostFrame pn_frame(pn);      // p + n rest frame

                // Boost particles to ppip rest frame
                PParticle pip_PPIP = ppip_frame.boost(pion);
                PParticle n_PPIP = ppip_frame.boost(mm_ppip);

                // Boost particles to npip rest frame
                PParticle pip_NPIP = npip_frame.boost(pion);
                PParticle p_NPIP = npip_frame.boost(fwd_objects[j]);

                // Boost particles to pn rest frame
                PParticle n_PN = pn_frame.boost(mm_ppip);
                PParticle pip_PN = pn_frame.boost(pion);

                // Group D: Helicity distributions (opening angles in rest frames)
                double helicity_pip = Physics::openingAngle(pip_PPIP, n_PPIP);
                double helicity_pipn = Physics::openingAngle(pip_NPIP, p_NPIP);
                double helicity_n = Physics::openingAngle(n_PN, pip_PN);

                mgr.fill("pwa_pip_helicity", cos(helicity_pip * M_PI / 180.0));
                mgr.fill("pwa_pipn_helicity", cos(helicity_pipn * M_PI / 180.0));
                mgr.fill("pwa_n_helicity", cos(helicity_n * M_PI / 180.0));

                // Gottfried-Jackson frames (with projectile)
                // Boost projectile to same rest frames
                PParticle proj_PPIP = ppip_frame.boost(beam);
                PParticle proj_NPIP = npip_frame.boost(beam);
                PParticle proj_PN = pn_frame.boost(beam);

                // Group E: GJ distributions (angles with respect to beam in rest frames)
                double gj_pip = Physics::openingAngle(pip_PPIP, proj_PPIP);
                double gj_pipn = Physics::openingAngle(pip_NPIP, proj_NPIP);
                double gj_n = Physics::openingAngle(n_PN, proj_PN);

                mgr.fill("pwa_pip_gj", cos(gj_pip * M_PI / 180.0));
                mgr.fill("pwa_pipn_gj", cos(gj_pipn * M_PI / 180.0));
                mgr.fill("pwa_n_gj", cos(gj_n * M_PI / 180.0));
                }  // End neutron_mass_cut
            }
        }
    }

}

// ============================================================================
// MAIN PROGRAM
// ============================================================================

int main(int argc, char* argv[]) {
    // Install signal handler for graceful Ctrl+C termination
    SignalHandler::install();

    ConsoleBox::newLine();
    ConsoleBox::printHeader("FAT Framework",
                           "Pion Analysis");
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

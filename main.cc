// ========================================================================
// FAT Framework - Dilepton Analysis
// ========================================================================
// e+e- dilepton analysis for HADES experiment (pp @ 4.5 GeV).
// Processes lepton pairs, ECAL photons, and Forward Tracker hadrons.
// Supports EpEm, EpEp, EmEm channels via config-driven prefix mapping.
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

void processEvent(NTupleReader& reader, Manager& mgr, CutManager& cuts,
                 const AnalysisConfig& config) {

    // Event-level cuts (applied before particle creation)
    if (!cuts.passValueCut("isBest", reader["isBest"])) return;
    if (!cuts.passMinCut("vertex_z", reader["eVertReco_z"])) return;

    // Create e+ and e- with reconstructed kinematics
    PParticle positron(MASS_ELECTRON, "e+");
    PParticle electron(MASS_ELECTRON, "e-");

    positron.setFromSpherical(reader["ep_p"], reader["ep_theta"], reader["ep_phi"],
                              KinematicType::RECONSTRUCTED);

    electron.setFromSpherical(reader["em_p"], reader["em_theta"], reader["em_phi"],
                              KinematicType::RECONSTRUCTED);

    // Energy-loss corrected kinematics (same angles, corrected momentum)
    positron.setFromSpherical(reader["ep_p_corr_ep"], reader["ep_theta"], reader["ep_phi"],
                              KinematicType::CORRECTED);

    electron.setFromSpherical(reader["em_p_corr_em"], reader["em_theta"], reader["em_phi"],
                              KinematicType::CORRECTED);

    // Lepton momentum histograms
    mgr.fill("ep_p", positron.momentum());
    mgr.fill("em_p", electron.momentum());

    // Momentum correction: delta_p = p_corrected - p_reconstructed
    double ep_p_rec = positron.momentum(KinematicType::RECONSTRUCTED);
    double ep_p_cor = positron.momentum(KinematicType::CORRECTED);
    double ep_dp = ep_p_cor - ep_p_rec;

    double em_p_rec = electron.momentum(KinematicType::RECONSTRUCTED);
    double em_p_cor = electron.momentum(KinematicType::CORRECTED);
    double em_dp = em_p_cor - em_p_rec;

    mgr.fill("ep_dp_vs_p", ep_p_rec, ep_dp);
    mgr.fill("em_dp_vs_p", em_p_rec, em_dp);

    // Opening angle between e+ and e- (must be computed before operator+)
    double oa = Physics::openingAngle(positron, electron);
    mgr.fill("opening_angle", oa);

    // Dilepton (e+ + e-) invariant mass
    PParticle dilepton = positron + electron;
    double m_ee = dilepton.massGeV();

    mgr.fill("mass_ee_before_oa", m_ee);

    // Boost dilepton to beam-target CMS frame
    EventFrames frames;
    frames.setBeamFrameFromKineticEnergy(config.getBeamKineticEnergy());

    PParticle dilepton_cms = frames.getFrame("beam").boost(dilepton);

    double y_cms = dilepton_cms.rapidity();
    double pt = dilepton_cms.vec().Pt();
    double theta_cms = dilepton_cms.theta();

    // Fill dilepton ntuple (before OA cut, store cut decision as flag)
    bool oa_pass = cuts.passMinCut("opening_angle", oa);

    auto& nt = mgr.getDynamicNtuple("dilepton_nt");

    nt["ep_p_rec"] = ep_p_rec;
    nt["ep_p_cor"] = ep_p_cor;
    nt["ep_theta"] = positron.theta();
    nt["ep_phi"] = positron.phi();
    nt["ep_theta_rich"] = reader["ep_theta_rich"];
    nt["ep_phi_rich"] = reader["ep_phi_rich"];

    nt["em_p_rec"] = em_p_rec;
    nt["em_p_cor"] = em_p_cor;
    nt["em_theta"] = electron.theta();
    nt["em_phi"] = electron.phi();
    nt["em_theta_rich"] = reader["em_theta_rich"];
    nt["em_phi_rich"] = reader["em_phi_rich"];

    nt["oa"] = oa;
    nt["m_ee"] = m_ee;

    nt["y_cms"] = y_cms;
    nt["pt"] = pt;
    nt["theta_cms"] = theta_cms;

    nt["oa_pass"] = oa_pass ? 1.0f : 0.0f;  // cut decision flag

    nt.fill();

    // Apply opening angle cut (reject close pairs)
    if (!oa_pass) return;

    mgr.fill("mass_ee_after_oa", m_ee);

    // CMS histograms (after OA cut)
    mgr.fill("rapidity_cms", y_cms);
    mgr.fill("pt_cms", pt);
    mgr.fill("theta_cms", theta_cms);
    mgr.fill("rapidity_vs_mass", m_ee, y_cms);

    // Dilepton invariant mass (after all cuts)
    mgr.fill("mass_ee", m_ee);

    // Beam + target for missing mass calculations
    PParticle beam = ParticleFactory::createBeamProton(config.getBeamKineticEnergy());
    PParticle target = ParticleFactory::createTargetProton();
    PParticle initial = beam + target;

    // Missing mass of e+e- system: MM(e+e-) = beam + target - e+ - e-
    PParticle miss_epem = initial - positron - electron;
    double mm_epem_mass = miss_epem.massGeV();
    double mm_epem_mass2 = miss_epem.vec().M2() / 1e6;  // GeV²/c⁴

    // ECAL objects (electromagnetic calorimeter, up to 5 hits)
    if (config.isEcalEnabled()) {
        int neutr_mult = static_cast<int>(reader["neutr_mult"]);

        std::vector<PParticleEcal> ecal_objects;
        ecal_objects.reserve(neutr_mult);

        for (int i = 1; i <= neutr_mult && i <= 5; ++i) {
            PParticleEcal ecal_obj(0.0, "ecal" + std::to_string(i));
            if (ecal_obj.setFromReader(reader, i)) {
                ecal_objects.push_back(ecal_obj);
            }
        }

        // ECAL quality cuts (pass/fail flags reused for compound building)
        std::vector<bool> ecal_pass(ecal_objects.size(), false);
        for (size_t j = 0; j < ecal_objects.size(); ++j) {
            const auto& obj = ecal_objects[j];
            ecal_pass[j] = cuts.passCutSet("ecal_quality", {
                static_cast<double>(obj.ecal_pid),
                obj.ecal_beta,
                obj.cluster_energy
            });
        }

        auto& ecal_nt = mgr.getDynamicNtuple("ecal_nt");

        ecal_nt["ecal_mult"] = neutr_mult;

        auto fillEcalHit = [&](size_t idx, const std::string& suffix) {
            if (ecal_objects.size() > idx) {
                const auto& obj = ecal_objects[idx];

                ecal_nt["ecal_pass" + suffix] = ecal_pass[idx] ? 1.0f : 0.0f;

                ecal_nt["cluster_energy" + suffix] = obj.cluster_energy;
                ecal_nt["cluster_theta" + suffix] = obj.cluster_theta;
                ecal_nt["cluster_phi" + suffix] = obj.cluster_phi;

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

        fillEcalHit(0, "_1");
        fillEcalHit(1, "_2");
        fillEcalHit(2, "_3");
        fillEcalHit(3, "_4");
        fillEcalHit(4, "_5");

        ecal_nt.fill();

        // Compound e+e-gamma objects (one entry per passing gamma)
        auto& epemg_nt = mgr.getDynamicNtuple("epemg_nt");

        int gamma_pass_mult = std::count(ecal_pass.begin(), ecal_pass.end(), true);

        // Energy-ranked indices of passing gammas (rank 1 = highest energy)
        std::vector<size_t> pass_indices;
        for (size_t j = 0; j < ecal_objects.size(); ++j) {
            if (ecal_pass[j]) pass_indices.push_back(j);
        }
        std::sort(pass_indices.begin(), pass_indices.end(),
                  [&](size_t a, size_t b) {
                      return ecal_objects[a].cluster_energy > ecal_objects[b].cluster_energy;
                  });

        std::map<size_t, int> rank_map;
        for (size_t r = 0; r < pass_indices.size(); ++r) {
            rank_map[pass_indices[r]] = static_cast<int>(r + 1);
        }

        for (size_t j = 0; j < ecal_objects.size(); ++j) {
            if (!ecal_pass[j]) continue;

            PParticle epemg = dilepton + ecal_objects[j];
            PParticle epemg_cms = frames.getFrame("beam").boost(epemg);

            epemg_nt["epemg_mass"] = epemg.massGeV();
            epemg_nt["epemg_p"] = epemg.momentum();
            epemg_nt["epemg_theta"] = epemg.theta();
            epemg_nt["epemg_phi"] = epemg.phi();

            epemg_nt["epemg_rapidity_cms"] = epemg_cms.rapidity();
            epemg_nt["epemg_pt_cms"] = epemg_cms.vec().Pt();
            epemg_nt["epemg_theta_cms"] = epemg_cms.theta();

            epemg_nt["ee_oa"] = oa;
            epemg_nt["ee_mass"] = m_ee;

            epemg_nt["gamma_energy"] = ecal_objects[j].cluster_energy;
            epemg_nt["gamma_theta"] = ecal_objects[j].cluster_theta;
            epemg_nt["gamma_phi"] = ecal_objects[j].cluster_phi;
            epemg_nt["gamma_index"] = static_cast<Float_t>(j + 1);

            epemg_nt["ecal_mult"] = static_cast<Float_t>(neutr_mult);
            epemg_nt["gamma_pass_mult"] = static_cast<Float_t>(gamma_pass_mult);
            epemg_nt["gamma_rank_energy"] = static_cast<Float_t>(rank_map[j]);

            // Missing masses
            epemg_nt["mm_epem_mass"] = mm_epem_mass;
            epemg_nt["mm_epem_mass2"] = mm_epem_mass2;

            PParticle miss_epemg = initial - positron - electron - ecal_objects[j];
            epemg_nt["mm_epemg_mass"] = miss_epemg.massGeV();
            epemg_nt["mm_epemg_mass2"] = miss_epemg.vec().M2() / 1e6;  // GeV²/c⁴

            epemg_nt.fill();
        }
    }

    // Forward Tracker objects (forward hadrons, up to 3 hits)
    if (config.isFwdEnabled()) {
        int fwdet_mult = static_cast<int>(reader["fwdet_mult"]);

        std::vector<PParticleFwd> fwd_objects;
        fwd_objects.reserve(fwdet_mult);

        for (int i = 1; i <= fwdet_mult && i <= 3; ++i) {
            PParticleFwd fwd_obj(MASS_PROTON, "fwd" + std::to_string(i));
            if (fwd_obj.setFromReader(reader, i)) {
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
            }
        };

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
    ConsoleBox::printHeader("FAT Framework",
                           "Dilepton Analysis");
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

    // Lepton prefix mapping (for like-sign channels: EpEp, EmEm)
    auto [prefix1, prefix2] = config.getLeptonPrefixes();
    if (!prefix1.empty()) {
        reader.setLeptonPrefixes(prefix1, prefix2);
        std::string channel = config.getChannelName();
        std::cout << "NTupleReader: Prefix mapping ep_ -> " << prefix1
                  << "_, em_ -> " << prefix2 << "_";
        if (!channel.empty()) std::cout << " (channel: " << channel << ")";
        std::cout << "\n";
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

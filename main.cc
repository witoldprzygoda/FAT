// ========================================================================
// FAT Framework - Pi+Pi- Hadronic Analysis
// ========================================================================
// pi+ pi- analysis for HADES experiment (pp @ 4.5 GeV).
// Tree:    PipPim_ID
// Optional ECAL: builds gamma+gamma compound (when ecal_mult == 2 and both
//                pass ecal_quality + pi0 mass window) and pippim+gg = "pippimgg"
//                as candidate for eta -> pi+ pi- pi0(-> gamma gamma).
//
// All compound observables are filled in two flavors:
//   - RECONSTRUCTED (default)
//   - CORRECTED  (KinematicType::CORRECTED)
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

using namespace Physics;

// ============================================================================
// PROCESS SINGLE EVENT
// ============================================================================

void processEvent(NTupleReader& reader, Manager& mgr, CutManager& cuts,
                 const AnalysisConfig& config) {

    // Event-level cuts (applied first)
    if (!cuts.passValueCut("isBest", reader["isBest"])) return;
    if (!cuts.passMinCut("vertex_z", reader["eVertReco_z"])) return;
    if (!cuts.passValueCut("trigger_PT3", reader["trigbit"])) return;
    if (!cuts.passCutSet("start_detector", { reader["start_iteration"] })) return;

    // ========================================================================
    // pi+ and pi- with reconstructed and corrected kinematics
    // ========================================================================
    PParticle piplus(MASS_PION_PLUS, "pi+");
    PParticle piminus(MASS_PION_MINUS, "pi-");

    piplus.setFromSpherical(reader["pip_p"], reader["pip_theta"], reader["pip_phi"],
                            KinematicType::RECONSTRUCTED);
    piminus.setFromSpherical(reader["pim_p"], reader["pim_theta"], reader["pim_phi"],
                             KinematicType::RECONSTRUCTED);

    piplus.setFromSpherical(reader["pip_p_corr_pip"], reader["pip_theta"], reader["pip_phi"],
                            KinematicType::CORRECTED);
    piminus.setFromSpherical(reader["pim_p_corr_pim"], reader["pim_theta"], reader["pim_phi"],
                             KinematicType::CORRECTED);

    // Pion momentum histograms (REC + COR)
    mgr.fill("pip_p", piplus.momentum());
    mgr.fill("pim_p", piminus.momentum());
    mgr.fill("pip_p_cor", piplus.momentum(KinematicType::CORRECTED));
    mgr.fill("pim_p_cor", piminus.momentum(KinematicType::CORRECTED));

    // Per-particle correction scatter
    double pip_p_rec = piplus.momentum(KinematicType::RECONSTRUCTED);
    double pip_p_cor = piplus.momentum(KinematicType::CORRECTED);
    double pim_p_rec = piminus.momentum(KinematicType::RECONSTRUCTED);
    double pim_p_cor = piminus.momentum(KinematicType::CORRECTED);
    mgr.fill("pip_dp_vs_p", pip_p_rec, pip_p_cor - pip_p_rec);
    mgr.fill("pim_dp_vs_p", pim_p_rec, pim_p_cor - pim_p_rec);

    // ========================================================================
    // Compound: pippim
    // ========================================================================
    PParticle pippim = piplus + piminus;
    double m_pippim     = pippim.massGeV();
    double m_pippim_cor = pippim.massGeV(KinematicType::CORRECTED);
    mgr.fill("mass_pippim", m_pippim);
    mgr.fill("mass_pippim_cor", m_pippim_cor);

    // Opening angle pi+ pi- (LAB) — direction-only, identical for REC and COR
    double oa_pip_pim = Physics::openingAngle(piplus, piminus);
    mgr.fill("oa_pip_pim", oa_pip_pim);

    // Beam + target for missing-mass calculations (synthetic — REC == COR)
    PParticle beam = ParticleFactory::createBeamProton(config.getBeamKineticEnergy());
    PParticle target = ParticleFactory::createTargetProton();
    PParticle initial = beam + target;

    PParticle miss_pippim = initial - pippim;
    double mm_pippim     = miss_pippim.massGeV();
    double mm_pippim_cor = miss_pippim.massGeV(KinematicType::CORRECTED);
    mgr.fill("mm_pippim", mm_pippim);
    mgr.fill("mm_pippim_cor", mm_pippim_cor);

    // Boost pippim to beam-target CMS for kinematic observables
    EventFrames frames;
    frames.setBeamFrameFromKineticEnergy(config.getBeamKineticEnergy());
    PParticle pippim_cms = frames.getFrame("beam").boost(pippim);

    double y_cms     = pippim_cms.rapidity();
    double pt        = pippim_cms.vec().Pt();
    double theta_cms = pippim_cms.theta();

    double y_cms_cor     = pippim_cms.rapidity(KinematicType::CORRECTED);
    double pt_cor        = pippim_cms.vec(KinematicType::CORRECTED).Pt();
    double theta_cms_cor = pippim_cms.theta(KinematicType::CORRECTED);

    mgr.fill("rapidity_cms", y_cms);
    mgr.fill("pt_cms", pt);
    mgr.fill("theta_cms", theta_cms);
    mgr.fill("rapidity_cms_cor", y_cms_cor);
    mgr.fill("pt_cms_cor", pt_cor);
    mgr.fill("theta_cms_cor", theta_cms_cor);

    // ========================================================================
    // ECAL: gamma+gamma -> pi0 hypothesis (only when ecal_mult == 2)
    // Both photons must pass ecal_quality, M(gg) must be in pi0 window.
    // Then build pippimgg = pippim + gg as eta candidate.
    // ========================================================================
    double m_gg          = -1.0,  m_gg_cor          = -1.0;
    double m_pippimgg    = -1.0,  m_pippimgg_cor    = -1.0;
    double mm_pippimgg   = -1.0,  mm_pippimgg_cor   = -1.0;
    double oa_gg         = -1.0;
    bool   gg_pass        = false;   // wide  pi0_mass_window        — gates ntuple values
    bool   gg_pass_narrow = false;   // narrow pi0_mass_window_narrow — gates histograms

    if (config.isEcalEnabled() && static_cast<int>(reader["neutr_mult"]) == 2) {
        PParticleEcal g1(0.0, "g1");
        PParticleEcal g2(0.0, "g2");
        bool ok1 = g1.setFromReader(reader, 1);
        bool ok2 = g2.setFromReader(reader, 2);

        if (ok1 && ok2) {
            bool q1 = cuts.passCutSet("ecal_quality", {
                static_cast<double>(g1.ecal_pid), g1.ecal_beta, g1.cluster_energy
            });
            bool q2 = cuts.passCutSet("ecal_quality", {
                static_cast<double>(g2.ecal_pid), g2.ecal_beta, g2.cluster_energy
            });
            if (q1 && q2) {
                // Mirror RECONSTRUCTED to CORRECTED for the photons (no measured correction)
                g1.setFromSpherical(g1.cluster_energy, g1.cluster_theta, g1.cluster_phi,
                                    KinematicType::CORRECTED);
                g2.setFromSpherical(g2.cluster_energy, g2.cluster_theta, g2.cluster_phi,
                                    KinematicType::CORRECTED);

                oa_gg = Physics::openingAngle(g1, g2);
                mgr.fill("oa_gg", oa_gg);

                PParticle gg = g1 + g2;
                m_gg     = gg.massGeV();
                m_gg_cor = gg.massGeV(KinematicType::CORRECTED);
                mgr.fill("mass_gg", m_gg);
                mgr.fill("mass_gg_cor", m_gg_cor);

                // pi0 mass window: WIDE gates pippimgg construction & ntuple values;
                // NARROW (subset of wide) gates histogram filling — tighter signal region
                // for plotted spectra now that statistics is large.
                if (cuts.passRangeCut("pi0_mass_window", m_gg)) {
                    PParticle pippimgg = pippim + gg;
                    m_pippimgg     = pippimgg.massGeV();
                    m_pippimgg_cor = pippimgg.massGeV(KinematicType::CORRECTED);

                    PParticle miss_pippimgg = initial - pippimgg;
                    mm_pippimgg     = miss_pippimgg.massGeV();
                    mm_pippimgg_cor = miss_pippimgg.massGeV(KinematicType::CORRECTED);

                    gg_pass = true;

                    // Narrow window: fill histograms only for the tight signal region
                    if (cuts.passRangeCut("pi0_mass_window_narrow", m_gg)) {
                        mgr.fill("mass_pippimgg", m_pippimgg);
                        mgr.fill("mass_pippimgg_cor", m_pippimgg_cor);
                        mgr.fill("mm_pippimgg", mm_pippimgg);
                        mgr.fill("mm_pippimgg_cor", mm_pippimgg_cor);
                        gg_pass_narrow = true;
                    }
                }
            }
        }
    }

    // ========================================================================
    // Output ntuples (RECONSTRUCTED-derived and CORRECTED-derived)
    // ========================================================================
    auto& nt = mgr.getDynamicNtuple("pippim_nt");

    nt["pip_p_rec"] = pip_p_rec;
    nt["pip_p_cor"] = pip_p_cor;
    nt["pip_theta"] = piplus.theta();
    nt["pip_phi"]   = piplus.phi();

    nt["pim_p_rec"] = pim_p_rec;
    nt["pim_p_cor"] = pim_p_cor;
    nt["pim_theta"] = piminus.theta();
    nt["pim_phi"]   = piminus.phi();

    nt["oa_pip_pim"] = oa_pip_pim;
    nt["oa_gg"]      = oa_gg;

    nt["m_pippim"]    = m_pippim;
    nt["m_gg"]        = m_gg;
    nt["m_pippimgg"]  = m_pippimgg;

    nt["mm_pippim"]   = mm_pippim;
    nt["mm_pippimgg"] = mm_pippimgg;

    nt["y_cms"]     = y_cms;
    nt["pt"]        = pt;
    nt["theta_cms"] = theta_cms;

    nt["gg_pass"]        = gg_pass        ? 1.0f : 0.0f;  // wide   M(gg) in [0.10, 0.18]
    nt["gg_pass_narrow"] = gg_pass_narrow ? 1.0f : 0.0f;  // narrow M(gg) in [0.125, 0.145]
    nt.fill();

    // Mirror with CORRECTED-derived compound observables
    auto& nt_cor = mgr.getDynamicNtuple("pippim_nt_cor");

    nt_cor["pip_p_rec"] = pip_p_rec;
    nt_cor["pip_p_cor"] = pip_p_cor;
    nt_cor["pip_theta"] = piplus.theta();
    nt_cor["pip_phi"]   = piplus.phi();

    nt_cor["pim_p_rec"] = pim_p_rec;
    nt_cor["pim_p_cor"] = pim_p_cor;
    nt_cor["pim_theta"] = piminus.theta();
    nt_cor["pim_phi"]   = piminus.phi();

    nt_cor["oa_pip_pim"] = oa_pip_pim;       // direction-only — same as REC
    nt_cor["oa_gg"]      = oa_gg;            // direction-only — same as REC

    nt_cor["m_pippim"]    = m_pippim_cor;
    nt_cor["m_gg"]        = m_gg_cor;
    nt_cor["m_pippimgg"]  = m_pippimgg_cor;

    nt_cor["mm_pippim"]   = mm_pippim_cor;
    nt_cor["mm_pippimgg"] = mm_pippimgg_cor;

    nt_cor["y_cms"]     = y_cms_cor;
    nt_cor["pt"]        = pt_cor;
    nt_cor["theta_cms"] = theta_cms_cor;

    nt_cor["gg_pass"]        = gg_pass        ? 1.0f : 0.0f;  // wide
    nt_cor["gg_pass_narrow"] = gg_pass_narrow ? 1.0f : 0.0f;  // narrow
    nt_cor.fill();

    // ========================================================================
    // ECAL detector ntuple (full per-hit info, regardless of selection)
    // ========================================================================
    if (config.isEcalEnabled()) {
        int neutr_mult = static_cast<int>(reader["neutr_mult"]);
        if (neutr_mult < 0 || neutr_mult > 5) neutr_mult = (neutr_mult < 0) ? 0 : 5;

        std::vector<PParticleEcal> ecal_objects;
        ecal_objects.reserve(neutr_mult);
        for (int i = 1; i <= neutr_mult && i <= 5; ++i) {
            PParticleEcal ecal_obj(0.0, "ecal" + std::to_string(i));
            if (ecal_obj.setFromReader(reader, i)) {
                ecal_objects.push_back(ecal_obj);
            }
        }

        // ECAL quality flags per hit (for ntuple inspection)
        std::vector<bool> ecal_pass(ecal_objects.size(), false);
        for (size_t j = 0; j < ecal_objects.size(); ++j) {
            const auto& obj = ecal_objects[j];
            ecal_pass[j] = cuts.passCutSet("ecal_quality", {
                static_cast<double>(obj.ecal_pid), obj.ecal_beta, obj.cluster_energy
            });
        }

        auto& ecal_nt = mgr.getDynamicNtuple("ecal_nt");
        ecal_nt["ecal_mult"] = neutr_mult;

        auto fillEcalHit = [&](size_t idx, const std::string& suffix) {
            if (ecal_objects.size() > idx) {
                const auto& obj = ecal_objects[idx];
                ecal_nt["ecal_pass" + suffix] = ecal_pass[idx] ? 1.0f : 0.0f;
                ecal_nt["cluster_energy" + suffix] = obj.cluster_energy;
                ecal_nt["cluster_theta" + suffix]  = obj.cluster_theta;
                ecal_nt["cluster_phi" + suffix]    = obj.cluster_phi;
                ecal_nt["ecal_beta" + suffix]      = obj.ecal_beta;
                ecal_nt["ecal_pid" + suffix]       = obj.ecal_pid;
                ecal_nt["ecal_energy" + suffix]    = obj.ecal_energy;
                ecal_nt["ecal_theta" + suffix]     = obj.ecal_theta;
                ecal_nt["ecal_phi" + suffix]       = obj.ecal_phi;
                ecal_nt["ecal_chi2" + suffix]      = obj.ecal_chi2;
                ecal_nt["ecal_tof" + suffix]       = obj.ecal_tof;
                ecal_nt["ecal_r" + suffix]         = obj.ecal_r;
                ecal_nt["ecal_z" + suffix]         = obj.ecal_z;
            }
        };
        fillEcalHit(0, "_1");
        fillEcalHit(1, "_2");
        fillEcalHit(2, "_3");
        fillEcalHit(3, "_4");
        fillEcalHit(4, "_5");
        ecal_nt.fill();
    }

    // ========================================================================
    // Forward Tracker ntuple (kept for hadron-side completeness)
    // ========================================================================
    if (config.isFwdEnabled()) {
        int fwdet_mult = static_cast<int>(reader["fwdet_mult"]);
        if (fwdet_mult < 0 || fwdet_mult > 3) fwdet_mult = (fwdet_mult < 0) ? 0 : 3;
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
                fwd_nt["fwd_p" + suffix]        = obj.momentum();
                fwd_nt["fwd_theta" + suffix]    = obj.fwd_theta;
                fwd_nt["fwd_phi" + suffix]      = obj.fwd_phi;
                fwd_nt["fwd_beta" + suffix]     = obj.fwd_beta;
                fwd_nt["fwd_mass" + suffix]     = obj.fwd_mass;
                fwd_nt["fwd_mass2" + suffix]    = obj.fwd_mass2;
                fwd_nt["fwd_charge" + suffix]   = obj.fwd_charge;
                fwd_nt["fwd_chi2" + suffix]     = obj.fwd_chi2;
                fwd_nt["fwd_ndf" + suffix]      = obj.fwd_ndf;
                fwd_nt["fwd_chi2ndf" + suffix]  = obj.fwd_chi2ndf;
                fwd_nt["fwd_r" + suffix]        = obj.fwd_r;
                fwd_nt["fwd_z" + suffix]        = obj.fwd_z;
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
    SignalHandler::install();

    ConsoleBox::newLine();
    ConsoleBox::printHeader("FAT Framework",
                           "Pi+Pi- Hadronic Analysis");
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

    // Optional: lepton-prefix mapping (used by like-sign CB channels in other analyses).
    // For PipPim_ID there are no leptons, so this is typically empty in config.
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

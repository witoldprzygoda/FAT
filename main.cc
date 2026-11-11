// ========================================================================
// FAT Framework - Dilepton Analysis
// ========================================================================
// e+e- dilepton analysis for HADES experiment (pp @ 4.53 GeV, T_kin = 4530 MeV).
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
#include <TError.h>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// Use Physics namespace for mass constants
using namespace Physics;

// ============================================================================
// FAILURE HANDLING
// ============================================================================
// Exit codes of ./ana (run_parallel.sh does not retry 2, 3 and 130):
//     0  output written and verified
//     1  configuration or input cannot be opened
//     2  input variables missing in the tree, or bound only inside the event loop
//     3  exception in processEvent
//     5  read error, ROOT error, or output file not written correctly
//   130  interrupted (Ctrl+C); the partial output is kept
// In the cases 2, 3 and 5 no output file is left behind.

// ROOT reports I/O trouble (input file that cannot be opened, read error, full
// disk) only as a printed message and carries on. This handler counts the
// messages of level Error and above; main() treats a non-zero count as failure.
static long g_root_errors = 0;

void countingErrorHandler(int level, Bool_t abort_bool, const char* location, const char* msg) {
    if (level >= kError) ++g_root_errors;
    DefaultErrorHandler(level, abort_bool, location, msg);
}

// End the program after a failure without leaving a valid-looking output file:
// close the output if asked (this also converts the ntuples and removes their
// intermediate files), remove it, remove intermediate ntuple files that are
// still there, and exit with `code`. std::_Exit skips the destructors: after a
// failed closeFile() the Manager destructor would finalize the ntuples again.
[[noreturn]] void abandonOutput(Manager& manager, const std::string& output, int code,
                                bool close_first) {
    if (close_first) {
        try {
            manager.closeFile();
        } catch (const std::exception& e) {
            std::cerr << "ERROR while closing the output file: " << e.what() << "\n";
        }
    }
    if (std::remove(output.c_str()) == 0) {
        std::cerr << "ERROR: incomplete output file removed: " << output << "\n";
    } else if (errno != ENOENT) {
        std::cerr << "ERROR: could not remove the incomplete output file: " << output << "\n";
    }
    // Same naming rule as in DynamicHNtuple
    const size_t dot = output.rfind('.');
    const std::string base = (dot != std::string::npos) ? output.substr(0, dot) : output;
    for (const auto& name : outputNtupleNames()) {
        std::remove((base + "_" + name + "_tree.root").c_str());
    }
    std::cerr << "ERROR: job failed (exit code " << code << ").\n";
    std::cout.flush();
    std::cerr.flush();
    std::_Exit(code);
}

// Reopen the closed output file and compare every ntuple with the number of
// fills. Returns an empty string if the file is good, otherwise the reason.
// (A write error does not stop ROOT: the file can end up truncated.)
std::string verifyOutputFile(const std::string& output, Manager& manager) {
    std::unique_ptr<TFile> file(TFile::Open(output.c_str(), "READ"));
    if (!file || file->IsZombie()) return "the file cannot be reopened";
    if (file->TestBit(TFile::kRecovered)) return "the file was not closed properly (ROOT had to recover it)";
    for (const auto& name : outputNtupleNames()) {
        if (!manager.hasDynamicNtuple(name)) continue;
        const DynamicHNtuple& nt = manager.getDynamicNtuple(name);
        if (nt.getVariableCount() == 0) continue;   // an ntuple without columns is not written
        TTree* tree = dynamic_cast<TTree*>(file->Get(name.c_str()));
        if (!tree) return "ntuple '" + name + "' is missing in the file";
        if (tree->GetEntries() != nt.getFillCount()) {
            return "ntuple '" + name + "' has " + std::to_string(tree->GetEntries()) +
                   " rows in the file, but " + std::to_string(nt.getFillCount()) + " were filled";
        }
    }
    return "";
}

// ============================================================================
// READER VARIABLES USED BY processEvent
// ============================================================================
// Every PAT branch that processEvent reads (directly, or through
// PParticleEcal / PParticleFwd::setFromReader) is listed here. main() binds
// the whole list BEFORE the output file is opened: a missing branch (wrong
// tree name or lepton prefix in the config, different PAT production) then
// stops the program with a message, instead of throwing on every row.
// Keep this list in step with processEvent.

std::vector<std::string> requiredReaderVariables(const AnalysisConfig& config) {
    // Event level
    std::vector<std::string> vars = {
        "isBest", "trigbit", "runnumber", "event",
        "eVertReco_x", "eVertReco_y", "eVertReco_z", "eVertReco_chi2",
        "neutr_mult", "fwdet_mult", "evtStart", "start_iteration",
        "ep_hyp_chi2_combined"
    };

    // All four particles ("ep_"/"em_" go through the lepton prefix mapping)
    for (const std::string prefix : {"ep", "em", "pip", "pim"}) {
        for (const char* name : {"p", "theta", "phi", "beta_new", "system", "sector",
                                 "r", "z", "rkchi2", "isUsedVertex", "oa_lept",
                                 "oa_hadr", "tof_new", "tof_mom"}) {
            vars.push_back(prefix + "_" + name);
        }
        vars.push_back(prefix + "_p_corr_" + prefix);   // energy-loss corrected momentum
    }
    // Leptons only
    for (const std::string prefix : {"ep", "em"}) {
        for (const char* name : {"theta_rich", "phi_rich", "rich_ind", "richmatchquality",
                                 "rich_padnum", "ringchi2", "kIsLepton"}) {
            vars.push_back(prefix + "_" + name);
        }
    }
    // Pions only
    for (const std::string prefix : {"pip", "pim"}) {
        for (const char* name : {"isring", "dedx_mdc"}) {
            vars.push_back(prefix + "_" + name);
        }
    }

    // ECAL neutral candidates, slots 1-5 (PParticleEcal::setFromReader)
    if (config.isEcalEnabled()) {
        vars.push_back("neutr_counter");
        for (int i = 1; i <= 5; ++i) {
            for (const char* name : {"cluster_theta", "cluster_phi", "cluster_energy", "p",
                                     "theta", "phi", "pid", "clusterid", "q", "energy",
                                     "p_pid", "mass", "mass2", "beta", "tof", "tofrec",
                                     "dist", "r", "z", "chi2", "phi2", "theta2", "r2", "z2"}) {
                vars.push_back(std::string("neutr_") + name + "_" + std::to_string(i));
            }
        }
    }

    // Forward detector candidates, slots 1-3 (PParticleFwd::setFromReader)
    if (config.isFwdEnabled()) {
        vars.push_back("fwdet_counter");
        for (int i = 1; i <= 3; ++i) {
            for (const char* name : {"p", "theta", "phi", "beta", "mass", "mass2", "q",
                                     "tofrec", "tof", "r", "z", "phi2", "theta2", "r2",
                                     "z2", "chi2", "ndf", "dist_tofhit"}) {
                vars.push_back(std::string("fwdet_") + name + "_" + std::to_string(i));
            }
        }
    }

    return vars;
}

// ============================================================================
// FILE NUMBER OF THE CURRENT INPUT FILE
// ============================================================================
// The last group of digits in the base name of the input file being read
// (lep_034.root -> 34); -1 if the name contains no digit. The name is parsed
// once per file (cached by the tree number of the chain).

float currentFileId(NTupleReader& reader) {
    static int   cached_tree = -1;
    static float cached_id   = -1.0f;

    const int tree_number = reader.getCurrentTreeNumber();
    if (tree_number != cached_tree) {
        cached_tree = tree_number;
        cached_id   = -1.0f;

        TFile* file = reader.getTree() ? reader.getTree()->GetCurrentFile() : nullptr;
        if (file) {
            std::string name = file->GetName();
            size_t slash = name.find_last_of('/');
            if (slash != std::string::npos) name = name.substr(slash + 1);
            size_t dot = name.find_last_of('.');
            if (dot != std::string::npos) name = name.substr(0, dot);

            size_t last = name.find_last_of("0123456789");
            if (last != std::string::npos) {
                size_t first = last;
                while (first > 0 && std::isdigit(static_cast<unsigned char>(name[first - 1]))) --first;
                cached_id = static_cast<float>(std::stol(name.substr(first, last - first + 1)));
            }
        }
    }
    return cached_id;
}

// Entry index of the current row inside its own PAT file. For a chain of
// several files this differs from the chain index reader.currentEntry(), so
// that (file_id, pat_entry) addresses a PAT row however the list was split.
Long64_t currentFileEntry(NTupleReader& reader) {
    TTree* tree = reader.getTree();
    if (reader.isChain() && tree && tree->GetTree()) {
        return tree->GetTree()->GetReadEntry();
    }
    return reader.currentEntry();
}

// ============================================================================
// PROCESS SINGLE EVENT
// ============================================================================

void processEvent(NTupleReader& reader, Manager& mgr, CutManager& cuts,
                 const AnalysisConfig& config) {

    // NOTE: with isBestAll (isBest >= 0) one call = one stored COMBINATION, and
    // one event can give several rows (same e+e- pair with another pion, ...).
    // Every histogram entry and every ntuple row below is therefore per
    // combination, not per event.

    // Event-level cuts (applied before particle creation)
    // All stored combinations are used; the "isBest" (== 1) cut is not applied.
    if (!cuts.passMinCut("isBestAll", reader["isBest"])) return;
    if (!cuts.passMinCut("vertex_z", reader["eVertReco_z"])) return;

    // Trigger selection (driven by config.trigger.selection):
    //   "PT3"  -> trigbit == 8192   (default; 2-particle leptonic with bias)
    //   "PT2"  -> trigbit == 4096   (unbiased 2-particle hadronic, downscale 64)
    //   "none" -> no trigger filter (debug / no-bias ≥3-particle topologies)
    {
        const std::string trig_cut = config.getTriggerCutName();
        if (!trig_cut.empty() &&
            !cuts.passValueCut(trig_cut, reader["trigbit"])) return;
    }

    // No start-detector requirement: the betanew productions reconstruct the
    // relative time (beta_new) without START. start_iteration is only stored
    // in the ntuples.

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

    // Create pi+ and pi- with reconstructed kinematics
    PParticle piplus(MASS_PION_PLUS, "pi+");
    PParticle piminus(MASS_PION_MINUS, "pi-");

    piplus.setFromSpherical(reader["pip_p"], reader["pip_theta"], reader["pip_phi"],
                            KinematicType::RECONSTRUCTED);

    piminus.setFromSpherical(reader["pim_p"], reader["pim_theta"], reader["pim_phi"],
                             KinematicType::RECONSTRUCTED);

    // Energy-loss corrected kinematics for pions (same angles, corrected momentum)
    piplus.setFromSpherical(reader["pip_p_corr_pip"], reader["pip_theta"], reader["pip_phi"],
                            KinematicType::CORRECTED);

    piminus.setFromSpherical(reader["pim_p_corr_pim"], reader["pim_theta"], reader["pim_phi"],
                             KinematicType::CORRECTED);

    // Lepton momentum histograms (RECONSTRUCTED + CORRECTED)
    mgr.fill("ep_p", positron.momentum());
    mgr.fill("em_p", electron.momentum());
    mgr.fill("ep_p_cor", positron.momentum(KinematicType::CORRECTED));
    mgr.fill("em_p_cor", electron.momentum(KinematicType::CORRECTED));

    // Pion momentum histograms (RECONSTRUCTED + CORRECTED)
    mgr.fill("pip_p", piplus.momentum());
    mgr.fill("pim_p", piminus.momentum());
    mgr.fill("pip_p_cor", piplus.momentum(KinematicType::CORRECTED));
    mgr.fill("pim_p_cor", piminus.momentum(KinematicType::CORRECTED));

    // Momentum correction: delta_p = p_corrected - p_reconstructed
    double ep_p_rec = positron.momentum(KinematicType::RECONSTRUCTED);
    double ep_p_cor = positron.momentum(KinematicType::CORRECTED);
    double ep_dp = ep_p_cor - ep_p_rec;

    double em_p_rec = electron.momentum(KinematicType::RECONSTRUCTED);
    double em_p_cor = electron.momentum(KinematicType::CORRECTED);
    double em_dp = em_p_cor - em_p_rec;

    double pip_p_rec = piplus.momentum(KinematicType::RECONSTRUCTED);
    double pip_p_cor = piplus.momentum(KinematicType::CORRECTED);
    double pip_dp = pip_p_cor - pip_p_rec;

    double pim_p_rec = piminus.momentum(KinematicType::RECONSTRUCTED);
    double pim_p_cor = piminus.momentum(KinematicType::CORRECTED);
    double pim_dp = pim_p_cor - pim_p_rec;

    mgr.fill("ep_dp_vs_p", ep_p_rec, ep_dp);
    mgr.fill("em_dp_vs_p", em_p_rec, em_dp);
    mgr.fill("pip_dp_vs_p", pip_p_rec, pip_dp);
    mgr.fill("pim_dp_vs_p", pim_p_rec, pim_dp);

    // Opening angle between e+ and e- (must be computed before operator+)
    double oa = Physics::openingAngle(positron, electron);
    mgr.fill("opening_angle", oa);

    // Compound objects
    PParticle epem = positron + electron;                   // e+ e-
    PParticle pippim = piplus + piminus;                    // pi+ pi-
    PParticle pippimepem = pippim + epem;                   // pi+ pi- e+ e-

    double m_ee = epem.massGeV();
    double m_pippim = pippim.massGeV();
    double m_pippimepem = pippimepem.massGeV();

    // CORRECTED counterparts (compound PParticles already carry both kinematics)
    double m_ee_cor = epem.massGeV(KinematicType::CORRECTED);
    double m_pippim_cor = pippim.massGeV(KinematicType::CORRECTED);
    double m_pippimepem_cor = pippimepem.massGeV(KinematicType::CORRECTED);

    mgr.fill("mass_ee_before_oa", m_ee);
    mgr.fill("mass_pippim", m_pippim);
    mgr.fill("mass_pippimepem", m_pippimepem);

    mgr.fill("mass_ee_before_oa_cor", m_ee_cor);
    mgr.fill("mass_pippim_cor", m_pippim_cor);
    mgr.fill("mass_pippimepem_cor", m_pippimepem_cor);

    // Beam + target for missing mass calculations (synthetic — REC == COR)
    PParticle beam = ParticleFactory::createBeamProton(config.getBeamKineticEnergy());
    PParticle target = ParticleFactory::createTargetProton();
    PParticle initial = beam + target;

    // Missing masses: MM(X) = beam + target - X (operator- propagates both kinematics)
    PParticle miss_epem = initial - epem;
    PParticle miss_pippim = initial - pippim;
    PParticle miss_pippimepem = initial - pippimepem;

    double mm_epem = miss_epem.massGeV();
    double mm_pippim = miss_pippim.massGeV();
    double mm_pippimepem = miss_pippimepem.massGeV();

    double mm_epem_cor = miss_epem.massGeV(KinematicType::CORRECTED);
    double mm_pippim_cor = miss_pippim.massGeV(KinematicType::CORRECTED);
    double mm_pippimepem_cor = miss_pippimepem.massGeV(KinematicType::CORRECTED);

    mgr.fill("mm_epem", mm_epem);
    mgr.fill("mm_pippim", mm_pippim);
    mgr.fill("mm_pippimepem", mm_pippimepem);
    mgr.fill("mm_vs_m_pippimepem", mm_pippimepem, m_pippimepem);

    mgr.fill("mm_epem_cor", mm_epem_cor);
    mgr.fill("mm_pippim_cor", mm_pippim_cor);
    mgr.fill("mm_pippimepem_cor", mm_pippimepem_cor);
    mgr.fill("mm_vs_m_pippimepem_cor", mm_pippimepem_cor, m_pippimepem_cor);

    // Boost compound systems to beam-target CMS frame.
    EventFrames frames;
    frames.setBeamFrameFromKineticEnergy(config.getBeamKineticEnergy());

    PParticle epem_cms = frames.getFrame("beam").boost(epem);

    double y_cms = epem_cms.rapidity();
    double pt = epem_cms.vec().Pt();
    double theta_cms = epem_cms.theta();

    // CORRECTED CMS-frame observables: same beam frame (synthetic, no correction),
    // read CORRECTED data of the boosted epem.
    PParticle epem_cms_cor = frames.getFrame("beam").boost(epem);  // boosts both kinematics
    double y_cms_cor = epem_cms_cor.rapidity(KinematicType::CORRECTED);
    double pt_cor = epem_cms_cor.vec(KinematicType::CORRECTED).Pt();
    double theta_cms_cor = epem_cms_cor.theta(KinematicType::CORRECTED);

    // Lab-frame opening angle between the (π+π-) and (e+e-) momenta — REC + COR.
    double oa_pippim_epem_lab     = Physics::openingAngle(pippim, epem);
    double oa_pippim_epem_lab_cor = Physics::openingAngle(pippim, epem,
                                                          KinematicType::CORRECTED);

    // OA between (pippim) and (epem) in the rest frame of a HYPOTHETICAL η whose
    // 3-momentum matches (pippim+epem) but whose mass is forced to m_η. When the
    // measured M(pippim+epem) ≈ m_η this approaches 180° (back-to-back); when it
    // differs the OA distribution spreads out — useful as a near-η selection.
    double oa_pippim_epem_eta_rest     = Physics::openingAngleInMassConstrainedRestFrame(
                                            pippim, epem, MASS_ETA);
    double oa_pippim_epem_eta_rest_cor = Physics::openingAngleInMassConstrainedRestFrame(
                                            pippim, epem, MASS_ETA,
                                            KinematicType::CORRECTED);

    // Control histograms — filled BEFORE the 4-body selection so the cut
    // boundaries (oa_pippim_epem_lab < 50, oa_pippim_epem_eta_rest > 140) can
    // be verified directly on these spectra.
    mgr.fill("oa_pippim_epem_lab",          oa_pippim_epem_lab);
    mgr.fill("oa_pippim_epem_eta_rest",     oa_pippim_epem_eta_rest);
    mgr.fill("oa_pippim_epem_lab_cor",      oa_pippim_epem_lab_cor);
    mgr.fill("oa_pippim_epem_eta_rest_cor", oa_pippim_epem_eta_rest_cor);

    // 4-body selection as standalone cuts: OA_LAB < 50, M(pi+pi-) < 0.420,
    // eta-rest OA > 140. Each cut is evaluated on EVERY row (no short-circuit)
    // and its decision is kept as a flag. RECONSTRUCTED and CORRECTED values go
    // through the same cut objects, so the cut-flow counters of these three
    // cuts mix REC and COR (each row is tested twice).
    bool pass_oa_lab      = cuts.passMaxCut("oa_pippim_epem_lab", oa_pippim_epem_lab);
    bool pass_m_pippim    = cuts.passMaxCut("m_pippim", m_pippim);
    bool pass_oa_eta_rest = cuts.passMinCut("oa_pippim_epem_eta_rest", oa_pippim_epem_eta_rest);
    // Same cuts with CORRECTED inputs
    bool pass_oa_lab_cor      = cuts.passMaxCut("oa_pippim_epem_lab", oa_pippim_epem_lab_cor);
    bool pass_m_pippim_cor    = cuts.passMaxCut("m_pippim", m_pippim_cor);
    bool pass_oa_eta_rest_cor = cuts.passMinCut("oa_pippim_epem_eta_rest", oa_pippim_epem_eta_rest_cor);

    // sel_pass keeps its meaning: the AND of the three cuts (REC / COR)
    bool sel_pass     = pass_oa_lab && pass_m_pippim && pass_oa_eta_rest;
    bool sel_pass_cor = pass_oa_lab_cor && pass_m_pippim_cor && pass_oa_eta_rest_cor;

    // Plain e+e- opening angle from below and from above (direction only: REC == COR)
    bool pass_oa_min = cuts.passMinCut("oa_min", oa);
    bool pass_oa_max = cuts.passMaxCut("oa_max", oa);

    // Four-body system for scan_nt: transverse momentum in the lab [MeV/c] and
    // rapidity in the beam-target centre-of-mass frame (REC + COR).
    double pt4     = pippimepem.vec().Pt();
    double pt4_cor = pippimepem.vec(KinematicType::CORRECTED).Pt();
    PParticle pippimepem_cms = frames.getFrame("beam").boost(pippimepem);  // boosts both kinematics
    double y4_cms     = pippimepem_cms.rapidity();
    double y4_cms_cor = pippimepem_cms.rapidity(KinematicType::CORRECTED);

    if (sel_pass) {
        mgr.fill("mass_pippimepem_selected", m_pippimepem);
        mgr.fill("mm_vs_m_pippimepem_selected", mm_pippimepem, m_pippimepem);

        // MM(pi+pi-e+e-) slice scan (selection + one slice window)
        if (cuts.passRangeCut("mm_slice_20_22", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_20_22", m_pippimepem);
        if (cuts.passRangeCut("mm_slice_22_24", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_22_24", m_pippimepem);
        if (cuts.passRangeCut("mm_slice_24_26", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_24_26", m_pippimepem);
        if (cuts.passRangeCut("mm_slice_26_28", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_26_28", m_pippimepem);
        if (cuts.passRangeCut("mm_slice_28_30", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_28_30", m_pippimepem);
    }

    // === CORRECTED selection block (mirror of the above, using *_cor values) ===
    if (sel_pass_cor) {
        mgr.fill("mass_pippimepem_selected_cor", m_pippimepem_cor);
        mgr.fill("mm_vs_m_pippimepem_selected_cor", mm_pippimepem_cor, m_pippimepem_cor);

        if (cuts.passRangeCut("mm_slice_20_22", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_20_22_cor", m_pippimepem_cor);
        if (cuts.passRangeCut("mm_slice_22_24", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_22_24_cor", m_pippimepem_cor);
        if (cuts.passRangeCut("mm_slice_24_26", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_24_26_cor", m_pippimepem_cor);
        if (cuts.passRangeCut("mm_slice_26_28", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_26_28_cor", m_pippimepem_cor);
        if (cuts.passRangeCut("mm_slice_28_30", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_28_30_cor", m_pippimepem_cor);
    }

    // === ETA / f1 candidates from ECAL (REC + COR) ===============================
    // ECAL N_gamma==1 case:
    //   eta -> pi+pi- pi0 -> pi+pi- e+e-gamma (Dalitz pi0)  — pippimepemg_pass flag
    //   f1  -> pi+pi- eta(e+e-gamma)  (Dalitz eta)           — eta_dalitz_pass flag
    //   The compound M(pi+pi-e+e-gamma) is the same in both cases — only the cut
    //   on M(e+e-gamma) differs (pi0 narrow vs eta window).
    // ECAL N_gamma==2 case:
    //   f1 -> pi+pi- eta(gamma gamma) — eta_gg_pass flag.
    //   New compound M(pi+pi-gamma gamma) is the f1 candidate from the gg branch.
    //
    // Photon has no measured correction (gamma is mirrored REC≡COR), so M(gg)
    // is identical between REC and COR. M(epemg) and the compounds with pippim,
    // however, DO differ between REC and COR via the lepton/pion momenta.
    double m_epemg_rec        = -1.0,   m_epemg_cor        = -1.0;
    double m_pippimepemg_rec  = -1.0,   m_pippimepemg_cor  = -1.0;
    double mm_pippimepemg_rec = -1.0,   mm_pippimepemg_cor = -1.0;  // beam+target − pippimepemg
    bool   pippimepemg_pass_rec = false, pippimepemg_pass_cor = false;
    bool   pippimepemg_pass_narrow_rec = false, pippimepemg_pass_narrow_cor = false;
    bool   eta_dalitz_pass_rec  = false, eta_dalitz_pass_cor  = false;

    double m_gg             = -1.0;     // gamma is REC≡COR — single value
    double m_pippim_gg_rec  = -1.0,   m_pippim_gg_cor   = -1.0;
    bool   eta_gg_pass      = false;    // M(gg) in eta window — REC and COR agree

    if (config.isEcalEnabled()) {
        int neutr_mult_for_eta = static_cast<int>(reader["neutr_mult"]);

        // -- mult==1: epemg compound + pi0/eta cuts on M(epemg) --
        if (neutr_mult_for_eta == 1) {
            PParticleEcal gamma(0.0, "gamma");
            if (gamma.setFromReader(reader, 1)) {
                bool gamma_quality = cuts.passCutSet("ecal_quality", {
                    static_cast<double>(gamma.ecal_pid),
                    gamma.ecal_beta,
                    gamma.cluster_energy
                });
                if (gamma_quality) {
                    // Mirror RECONSTRUCTED to CORRECTED (photon has no measured correction)
                    gamma.setFromSpherical(gamma.cluster_energy,
                                           gamma.cluster_theta,
                                           gamma.cluster_phi,
                                           KinematicType::CORRECTED);
                    PParticle epemg      = epem   + gamma;
                    PParticle pippimepemg = pippim + epemg;
                    m_epemg_rec       = epemg.massGeV(KinematicType::RECONSTRUCTED);
                    m_epemg_cor       = epemg.massGeV(KinematicType::CORRECTED);
                    m_pippimepemg_rec = pippimepemg.massGeV(KinematicType::RECONSTRUCTED);
                    m_pippimepemg_cor = pippimepemg.massGeV(KinematicType::CORRECTED);

                    // Missing mass to the full 5-body: MM(pi+pi-e+e-gamma)
                    PParticle miss_pippimepemg = initial - pippimepemg;
                    mm_pippimepemg_rec = miss_pippimepemg.massGeV(KinematicType::RECONSTRUCTED);
                    mm_pippimepemg_cor = miss_pippimepemg.massGeV(KinematicType::CORRECTED);

                    pippimepemg_pass_rec        = cuts.passRangeCut("pi0_mass_window",        m_epemg_rec);
                    pippimepemg_pass_cor        = cuts.passRangeCut("pi0_mass_window",        m_epemg_cor);
                    pippimepemg_pass_narrow_rec = cuts.passRangeCut("pi0_mass_window_narrow", m_epemg_rec);
                    pippimepemg_pass_narrow_cor = cuts.passRangeCut("pi0_mass_window_narrow", m_epemg_cor);
                    eta_dalitz_pass_rec         = cuts.passRangeCut("eta_mass_window",        m_epemg_rec);
                    eta_dalitz_pass_cor         = cuts.passRangeCut("eta_mass_window",        m_epemg_cor);
                }
            }
        }

        // -- mult==2: gg compound + eta cut on M(gg) (f1 -> pi+pi- eta(gg)) --
        if (neutr_mult_for_eta == 2) {
            PParticleEcal g1(0.0, "g1");
            PParticleEcal g2(0.0, "g2");
            bool ok1 = g1.setFromReader(reader, 1);
            bool ok2 = g2.setFromReader(reader, 2);
            if (ok1 && ok2) {
                bool q1 = cuts.passCutSet("ecal_quality", {
                    static_cast<double>(g1.ecal_pid), g1.ecal_beta, g1.cluster_energy});
                bool q2 = cuts.passCutSet("ecal_quality", {
                    static_cast<double>(g2.ecal_pid), g2.ecal_beta, g2.cluster_energy});
                if (q1 && q2) {
                    g1.setFromSpherical(g1.cluster_energy, g1.cluster_theta, g1.cluster_phi,
                                        KinematicType::CORRECTED);
                    g2.setFromSpherical(g2.cluster_energy, g2.cluster_theta, g2.cluster_phi,
                                        KinematicType::CORRECTED);
                    PParticle gg         = g1 + g2;
                    PParticle pippim_gg  = pippim + gg;
                    m_gg            = gg.massGeV(KinematicType::CORRECTED);   // REC≡COR
                    m_pippim_gg_rec = pippim_gg.massGeV(KinematicType::RECONSTRUCTED);
                    m_pippim_gg_cor = pippim_gg.massGeV(KinematicType::CORRECTED);

                    mgr.fill("mass_gg_cor", m_gg);              // control plot
                    eta_gg_pass = cuts.passRangeCut("eta_mass_window", m_gg);
                }
            }
        }
    }

    // Fill pi+pi-e+e- ntuple (no OA cut applied; the decision is stored as a flag)
    // oa_pass: 2 deg <= OA(e+e-) <= 6 deg (opening_angle_2_6).
    bool oa_pass = cuts.passRangeCut("opening_angle_2_6", oa);

    auto& nt = mgr.getDynamicNtuple("pippimepem_nt");

    nt["ep_p_rec"] = ep_p_rec;
    nt["ep_p_cor"] = ep_p_cor;
    nt["ep_theta"] = positron.theta();
    nt["ep_phi"] = positron.phi();
    nt["ep_theta_rich"] = reader["ep_theta_rich"];
    nt["ep_phi_rich"] = reader["ep_phi_rich"];
    nt["ep_beta_new"] = reader["ep_beta_new"];

    nt["em_p_rec"] = em_p_rec;
    nt["em_p_cor"] = em_p_cor;
    nt["em_theta"] = electron.theta();
    nt["em_phi"] = electron.phi();
    nt["em_theta_rich"] = reader["em_theta_rich"];
    nt["em_phi_rich"] = reader["em_phi_rich"];
    nt["em_beta_new"] = reader["em_beta_new"];

    nt["pip_p_rec"] = pip_p_rec;
    nt["pip_p_cor"] = pip_p_cor;
    nt["pip_theta"] = piplus.theta();
    nt["pip_phi"] = piplus.phi();
    nt["pip_beta_new"] = reader["pip_beta_new"];

    nt["pim_p_rec"] = pim_p_rec;
    nt["pim_p_cor"] = pim_p_cor;
    nt["pim_theta"] = piminus.theta();
    nt["pim_phi"] = piminus.phi();
    nt["pim_beta_new"] = reader["pim_beta_new"];

    nt["start_iteration"] = reader["start_iteration"];  // stored only, not cut on

    // Combination keys: isBest (0/1; exactly one row with 1 per event), run, event.
    // (runnumber, event) is NOT unique; rows of one event are adjacent.
    nt["isBest"] = reader["isBest"];
    nt["runnumber"] = reader["runnumber"];
    nt["event"] = reader["event"];

    nt["oa"] = oa;

    nt["m_ee"] = m_ee;                       // M(e+e-)
    nt["m_pippim"] = m_pippim;               // M(pi+pi-)
    nt["m_pippimepem"] = m_pippimepem;       // M(pi+pi-e+e-)

    nt["mm_epem"] = mm_epem;                 // MM(e+e-)
    nt["mm_pippim"] = mm_pippim;             // MM(pi+pi-)
    nt["mm_pippimepem"] = mm_pippimepem;     // MM(pi+pi-e+e-)

    nt["oa_pippim_epem_lab"]      = oa_pippim_epem_lab;
    nt["oa_pippim_epem_eta_rest"] = oa_pippim_epem_eta_rest;

    nt["y_cms"] = y_cms;
    nt["pt"] = pt;
    nt["theta_cms"] = theta_cms;

    nt["oa_pass"] = oa_pass ? 1.0f : 0.0f;       // 2 <= OA <= 6 deg decision
    nt["sel_pass"] = sel_pass ? 1.0f : 0.0f;     // 4-body selection decision (AND of three cuts)

    // ECAL-derived eta / f1 fields — RECONSTRUCTED.
    //   mult==1: shared compound M(pi+pi-e+e-gamma); pi0 vs eta window flags
    //   mult==2: M(pi+pi-gamma gamma); M(gg) is REC≡COR so the gg flag is shared
    nt["m_epemg"]                = m_epemg_rec;
    nt["m_pippimepemg"]          = m_pippimepemg_rec;
    nt["mm_pippimepemg"]         = mm_pippimepemg_rec;
    nt["pippimepemg_pass"]        = pippimepemg_pass_rec        ? 1.0f : 0.0f;  // wide  [0.10, 0.18]
    nt["pippimepemg_pass_narrow"] = pippimepemg_pass_narrow_rec ? 1.0f : 0.0f;  // ACTIVE [0.125, 0.145]
    nt["eta_dalitz_pass"]        = eta_dalitz_pass_rec  ? 1.0f : 0.0f;

    nt["m_gg"]               = m_gg;
    nt["m_pippim_gg"]        = m_pippim_gg_rec;
    nt["eta_gg_pass"]        = eta_gg_pass         ? 1.0f : 0.0f;

    nt.fill();

    // === Mirror ntuple: same field names, compound observables from CORRECTED ===
    auto& nt_cor = mgr.getDynamicNtuple("pippimepem_nt_cor");

    // Per-particle fields are identical (raw _rec/_cor and angles don't change)
    nt_cor["ep_p_rec"] = ep_p_rec;
    nt_cor["ep_p_cor"] = ep_p_cor;
    nt_cor["ep_theta"] = positron.theta();
    nt_cor["ep_phi"] = positron.phi();
    nt_cor["ep_theta_rich"] = reader["ep_theta_rich"];
    nt_cor["ep_phi_rich"] = reader["ep_phi_rich"];
    nt_cor["ep_beta_new"] = reader["ep_beta_new"];

    nt_cor["em_p_rec"] = em_p_rec;
    nt_cor["em_p_cor"] = em_p_cor;
    nt_cor["em_theta"] = electron.theta();
    nt_cor["em_phi"] = electron.phi();
    nt_cor["em_theta_rich"] = reader["em_theta_rich"];
    nt_cor["em_phi_rich"] = reader["em_phi_rich"];
    nt_cor["em_beta_new"] = reader["em_beta_new"];

    nt_cor["pip_p_rec"] = pip_p_rec;
    nt_cor["pip_p_cor"] = pip_p_cor;
    nt_cor["pip_theta"] = piplus.theta();
    nt_cor["pip_phi"] = piplus.phi();
    nt_cor["pip_beta_new"] = reader["pip_beta_new"];

    nt_cor["pim_p_rec"] = pim_p_rec;
    nt_cor["pim_p_cor"] = pim_p_cor;
    nt_cor["pim_theta"] = piminus.theta();
    nt_cor["pim_phi"] = piminus.phi();
    nt_cor["pim_beta_new"] = reader["pim_beta_new"];

    nt_cor["start_iteration"] = reader["start_iteration"];  // stored only, not cut on

    // Combination keys (same as in pippimepem_nt)
    nt_cor["isBest"] = reader["isBest"];
    nt_cor["runnumber"] = reader["runnumber"];
    nt_cor["event"] = reader["event"];

    nt_cor["oa"] = oa;  // single-pair OA is direction-only — same as REC

    // Compound observables — computed from CORRECTED kinematics
    nt_cor["m_ee"] = m_ee_cor;
    nt_cor["m_pippim"] = m_pippim_cor;
    nt_cor["m_pippimepem"] = m_pippimepem_cor;

    nt_cor["mm_epem"] = mm_epem_cor;
    nt_cor["mm_pippim"] = mm_pippim_cor;
    nt_cor["mm_pippimepem"] = mm_pippimepem_cor;

    nt_cor["oa_pippim_epem_lab"]      = oa_pippim_epem_lab_cor;
    nt_cor["oa_pippim_epem_eta_rest"] = oa_pippim_epem_eta_rest_cor;

    nt_cor["y_cms"] = y_cms_cor;
    nt_cor["pt"] = pt_cor;
    nt_cor["theta_cms"] = theta_cms_cor;

    nt_cor["oa_pass"] = oa_pass ? 1.0f : 0.0f;             // OA single-pair: same decision
    nt_cor["sel_pass"] = sel_pass_cor ? 1.0f : 0.0f;       // CORRECTED selection chain

    // ECAL-derived eta / f1 fields — CORRECTED.
    //   mult==1: shared compound M(pi+pi-e+e-gamma); two flags select interpretation
    //            (pi0 Dalitz inside epemg vs eta Dalitz)
    //   mult==2: separate compound M(pi+pi-gamma gamma) for f1 -> pi+pi-eta(gg);
    //            M(gg) itself is REC≡COR (photon mirrored), so eta_gg_pass shared.
    nt_cor["m_epemg"]                = m_epemg_cor;
    nt_cor["m_pippimepemg"]          = m_pippimepemg_cor;
    nt_cor["mm_pippimepemg"]         = mm_pippimepemg_cor;
    nt_cor["pippimepemg_pass"]        = pippimepemg_pass_cor        ? 1.0f : 0.0f;  // wide  [0.10, 0.18]
    nt_cor["pippimepemg_pass_narrow"] = pippimepemg_pass_narrow_cor ? 1.0f : 0.0f;  // ACTIVE [0.125, 0.145]
    nt_cor["eta_dalitz_pass"]        = eta_dalitz_pass_cor  ? 1.0f : 0.0f;

    nt_cor["m_gg"]               = m_gg;
    nt_cor["m_pippim_gg"]        = m_pippim_gg_cor;
    nt_cor["eta_gg_pass"]        = eta_gg_pass         ? 1.0f : 0.0f;

    nt_cor.fill();

    // No OA cut on the event: only the 2-6 deg window spectrum is filled separately
    if (oa_pass) {
        mgr.fill("mass_ee_after_oa", m_ee);
        mgr.fill("mass_ee_after_oa_cor", m_ee_cor);
    }
    // CMS histograms (after OA cut)
    mgr.fill("rapidity_cms", y_cms);
    mgr.fill("pt_cms", pt);
    mgr.fill("theta_cms", theta_cms);
    mgr.fill("rapidity_vs_mass", m_ee, y_cms);

    mgr.fill("rapidity_cms_cor", y_cms_cor);
    mgr.fill("pt_cms_cor", pt_cor);
    mgr.fill("theta_cms_cor", theta_cms_cor);
    mgr.fill("rapidity_vs_mass_cor", m_ee_cor, y_cms_cor);

    // Dilepton invariant mass (after all cuts)
    mgr.fill("mass_ee", m_ee);
    mgr.fill("mass_ee_cor", m_ee_cor);

    // Photon and forward-candidate values for scan_nt (written at the end of
    // processEvent). -1 = absent slot; n_gam stays -1 when ECAL is off in the config.
    double scan_n_gam = -1.0;
    double scan_g_e[3]      = {-1.0, -1.0, -1.0};
    double scan_g_theta[3]  = {-1.0, -1.0, -1.0};
    double scan_g_phi[3]    = {-1.0, -1.0, -1.0};
    double scan_fw_p[2]     = {-1.0, -1.0};
    double scan_fw_theta[2] = {-1.0, -1.0};
    double scan_fw_phi[2]   = {-1.0, -1.0};

    // ECAL objects (electromagnetic calorimeter, up to 5 hits)
    if (config.isEcalEnabled()) {
        int neutr_mult = static_cast<int>(reader["neutr_mult"]);
        // PAT stores -100 when there is no ECAL cluster; clamp before reserve()
        if (neutr_mult < 0) neutr_mult = 0;

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
            } else {
                // Absent hit: still write every column (as -1), so ntuples from
                // separate run_parallel chunks have identical columns for hadd.
                for (const char* name : {"ecal_pass", "cluster_energy", "cluster_theta",
                                         "cluster_phi", "ecal_beta", "ecal_pid",
                                         "ecal_energy", "ecal_theta", "ecal_phi",
                                         "ecal_chi2", "ecal_tof", "ecal_r", "ecal_z"}) {
                    ecal_nt[name + suffix] = -1.0f;
                }
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

        // scan_nt: number of clusters passing ecal_quality (among the five
        // slots) and the three most energetic of them
        scan_n_gam = static_cast<double>(pass_indices.size());
        for (size_t r = 0; r < pass_indices.size() && r < 3; ++r) {
            const auto& obj = ecal_objects[pass_indices[r]];
            scan_g_e[r]     = obj.cluster_energy;
            scan_g_theta[r] = obj.cluster_theta;
            scan_g_phi[r]   = obj.cluster_phi;
        }

        for (size_t j = 0; j < ecal_objects.size(); ++j) {
            if (!ecal_pass[j]) continue;

            PParticle epemg = epem + ecal_objects[j];
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
            epemg_nt["mm_epem_mass"] = mm_epem;
            epemg_nt["mm_epem_mass2"] = miss_epem.vec().M2() / 1e6;  // GeV²/c⁴

            PParticle miss_epemg = initial - epem - ecal_objects[j];
            epemg_nt["mm_epemg_mass"] = miss_epemg.massGeV();
            epemg_nt["mm_epemg_mass2"] = miss_epemg.vec().M2() / 1e6;  // GeV²/c⁴

            // fill() resets every column and Manager::setEventWeight() sets
            // trigger_corr only once per PAT row: set it again for each photon
            // (it was -1 for the second and later photons of a combination).
            epemg_nt["trigger_corr"] = static_cast<Float_t>(mgr.getEventWeight());

            epemg_nt.fill();
        }
    }

    // Forward Tracker objects (forward hadrons, up to 3 hits)
    if (config.isFwdEnabled()) {
        int fwdet_mult = static_cast<int>(reader["fwdet_mult"]);
        // PAT stores -100 when there is no FwDet track; clamp before reserve()
        if (fwdet_mult < 0) fwdet_mult = 0;

        std::vector<PParticleFwd> fwd_objects;
        fwd_objects.reserve(fwdet_mult);

        for (int i = 1; i <= fwdet_mult && i <= 3; ++i) {
            PParticleFwd fwd_obj(MASS_PROTON, "fwd" + std::to_string(i));
            if (fwd_obj.setFromReader(reader, i)) {
                fwd_objects.push_back(fwd_obj);
            }
        }

        // scan_nt: the first two forward candidates
        for (size_t k = 0; k < fwd_objects.size() && k < 2; ++k) {
            scan_fw_p[k]     = fwd_objects[k].momentum();
            scan_fw_theta[k] = fwd_objects[k].fwd_theta;
            scan_fw_phi[k]   = fwd_objects[k].fwd_phi;
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
            } else {
                // Absent hit: write every column as -1 (same reason as in ecal_nt)
                for (const char* name : {"fwd_p", "fwd_theta", "fwd_phi", "fwd_beta",
                                         "fwd_mass", "fwd_mass2", "fwd_charge",
                                         "fwd_chi2", "fwd_ndf", "fwd_chi2ndf",
                                         "fwd_r", "fwd_z", "fwd_distToRpc"}) {
                    fwd_nt[name + suffix] = -1.0f;
                }
            }
        };

        fillFwdHit(0, "_1");
        fillFwdHit(1, "_2");
        fillFwdHit(2, "_3");

        fwd_nt.fill();
    }

    // ========================================================================
    // scan_nt: cut-scan ntuple, one row per combination
    // ========================================================================
    // Filled for every row that passed isBestAll + vertex + trigger; there is
    // no other pre-selection (no mass window). EVERY column is assigned on
    // EVERY fill, all of them in this one block: a DynamicHNtuple column first
    // written after N fills would be shifted by N rows. Do not write a scan_nt
    // column anywhere else. (trigger_corr is added by Manager::setEventWeight.)
    {
        auto& sn = mgr.getDynamicNtuple("scan_nt");

        // --- Keys and event ---
        // file_id: number in the name of the input file (lep_034.root -> 34).
        // pat_entry: entry index of the row in the tree of that file (equal to
        // the index in the chain of the job when the job reads one file).
        // Float columns hold integers exactly below 2^24.
        sn["pat_entry"] = static_cast<Float_t>(currentFileEntry(reader));
        sn["file_id"]   = currentFileId(reader);
        sn["runnumber"] = reader["runnumber"];
        sn["event"]     = reader["event"];
        sn["isBest"]    = reader["isBest"];
        sn["trigbit"]   = reader["trigbit"];
        // Ranking chi2 of the hypothesis (lowest in the event <=> isBest == 1);
        // PAT stores the same value on all four particles.
        sn["hyp_chi2"]  = reader["ep_hyp_chi2_combined"];

        sn["eVertReco_x"]    = reader["eVertReco_x"];
        sn["eVertReco_y"]    = reader["eVertReco_y"];
        sn["eVertReco_z"]    = reader["eVertReco_z"];
        sn["eVertReco_chi2"] = reader["eVertReco_chi2"];

        // PAT stores -100 when there is no ECAL cluster / forward candidate: 0 here
        sn["neutr_mult"] = std::max(0.0f, reader["neutr_mult"]);
        sn["fwdet_mult"] = std::max(0.0f, reader["fwdet_mult"]);

        sn["evtStart"]        = reader["evtStart"];         // stored only, not cut on
        sn["start_iteration"] = reader["start_iteration"];  // stored only, not cut on

        // --- Tracks ---
        // p_rec, p_cor, theta, phi as in pippimepem_nt (phi in -180..180 deg);
        // all other columns are copied from PAT (phi_rich in 0..360 deg).
        auto fillTrack = [&](const std::string& prefix, const PParticle& particle,
                             double p_rec, double p_cor) {
            sn[prefix + "p_rec"] = p_rec;
            sn[prefix + "p_cor"] = p_cor;
            sn[prefix + "theta"] = particle.theta();
            sn[prefix + "phi"]   = particle.phi();
            for (const char* name : {"beta_new", "system", "sector", "r", "z", "rkchi2",
                                     "isUsedVertex", "oa_lept", "oa_hadr", "tof_new",
                                     "tof_mom"}) {
                sn[prefix + name] = reader[prefix + name];
            }
        };

        fillTrack("ep_",  positron, ep_p_rec,  ep_p_cor);
        fillTrack("em_",  electron, em_p_rec,  em_p_cor);
        fillTrack("pip_", piplus,   pip_p_rec, pip_p_cor);
        fillTrack("pim_", piminus,  pim_p_rec, pim_p_cor);

        // cor_ok: 1 if all four energy-loss corrected momenta are finite. PAT
        // gives NaN for some pi- (pim_p_corr_pim); every _cor quantity with
        // that pion is then NaN and the _cor flags are 0. The values are stored
        // as they come: no replacement by the reconstructed momentum.
        sn["cor_ok"] = (std::isfinite(ep_p_cor) && std::isfinite(em_p_cor) &&
                        std::isfinite(pip_p_cor) && std::isfinite(pim_p_cor)) ? 1.0f : 0.0f;

        // Leptons only: RICH ring and lepton flag
        for (const std::string prefix : {"ep_", "em_"}) {
            for (const char* name : {"theta_rich", "phi_rich", "rich_ind", "richmatchquality",
                                     "rich_padnum", "ringchi2", "kIsLepton"}) {
                sn[prefix + name] = reader[prefix + name];
            }
        }

        // Pions only: ring flag, MDC dE/dx and mass squared.
        // mass2 = p_rec^2 (1/beta_new^2 - 1) in (MeV/c^2)^2 is computed here: all
        // PAT mass2 branches (mass2, mass2_pip, mass2_corr_pip, ...) are built
        // from the START-based beta, none from beta_new.
        auto mass2FromBetaNew = [](double p, double beta_new) {
            return beta_new > 0.0 ? p * p * (1.0 / (beta_new * beta_new) - 1.0) : -1.0e6;
        };
        for (const std::string prefix : {"pip_", "pim_"}) {
            sn[prefix + "isring"]   = reader[prefix + "isring"];
            sn[prefix + "dedx_mdc"] = reader[prefix + "dedx_mdc"];
        }
        sn["pip_mass2"] = mass2FromBetaNew(pip_p_rec, reader["pip_beta_new"]);
        sn["pim_mass2"] = mass2FromBetaNew(pim_p_rec, reader["pim_beta_new"]);

        // --- Cut variables as values: RECONSTRUCTED (plain name) and CORRECTED (_cor) ---
        sn["oa"] = oa;                                   // OA(e+e-), direction only

        sn["m_pippimepem"]     = m_pippimepem;           // M(pi+pi-e+e-) [GeV/c^2]
        sn["m_pippimepem_cor"] = m_pippimepem_cor;
        sn["m_ee"]             = m_ee;                   // M(e+e-)
        sn["m_ee_cor"]         = m_ee_cor;
        sn["m_pippim"]         = m_pippim;               // M(pi+pi-)
        sn["m_pippim_cor"]     = m_pippim_cor;

        sn["oa_pippim_epem_lab"]          = oa_pippim_epem_lab;
        sn["oa_pippim_epem_lab_cor"]      = oa_pippim_epem_lab_cor;
        sn["oa_pippim_epem_eta_rest"]     = oa_pippim_epem_eta_rest;
        sn["oa_pippim_epem_eta_rest_cor"] = oa_pippim_epem_eta_rest_cor;

        sn["mm_pippimepem"]     = mm_pippimepem;         // MM(pi+pi-e+e-) [GeV/c^2]
        sn["mm_pippimepem_cor"] = mm_pippimepem_cor;
        sn["mm_pippim"]         = mm_pippim;             // MM(pi+pi-)
        sn["mm_pippim_cor"]     = mm_pippim_cor;
        sn["mm_epem"]           = mm_epem;               // MM(e+e-)
        sn["mm_epem_cor"]       = mm_epem_cor;

        sn["pt4"]        = pt4;                          // pT of pi+pi-e+e-, lab [MeV/c]
        sn["pt4_cor"]    = pt4_cor;
        sn["y4_cms"]     = y4_cms;                       // rapidity of pi+pi-e+e-, beam-target CMS
        sn["y4_cms_cor"] = y4_cms_cor;

        // --- Cut decisions (1/0), each evaluated on every row ---
        sn["pass_oa_min"] = pass_oa_min ? 1.0f : 0.0f;   // OA(e+e-) > 2 deg
        sn["pass_oa_max"] = pass_oa_max ? 1.0f : 0.0f;   // OA(e+e-) < 6 deg
        sn["oa_pass"]     = oa_pass     ? 1.0f : 0.0f;   // 2 <= OA(e+e-) <= 6 deg (legacy flag)

        sn["pass_oa_lab"]          = pass_oa_lab          ? 1.0f : 0.0f;
        sn["pass_oa_lab_cor"]      = pass_oa_lab_cor      ? 1.0f : 0.0f;
        sn["pass_m_pippim"]        = pass_m_pippim        ? 1.0f : 0.0f;
        sn["pass_m_pippim_cor"]    = pass_m_pippim_cor    ? 1.0f : 0.0f;
        sn["pass_oa_eta_rest"]     = pass_oa_eta_rest     ? 1.0f : 0.0f;
        sn["pass_oa_eta_rest_cor"] = pass_oa_eta_rest_cor ? 1.0f : 0.0f;
        sn["sel_pass"]             = sel_pass             ? 1.0f : 0.0f;
        sn["sel_pass_cor"]         = sel_pass_cor         ? 1.0f : 0.0f;

        // --- ECAL: clusters passing ecal_quality, the three most energetic ---
        sn["n_gam"] = scan_n_gam;
        for (int k = 0; k < 3; ++k) {
            const std::string prefix = "g" + std::to_string(k + 1) + "_";
            sn[prefix + "e"]     = scan_g_e[k];          // cluster energy [MeV]
            sn[prefix + "theta"] = scan_g_theta[k];
            sn[prefix + "phi"]   = scan_g_phi[k];
        }

        // --- Forward detector: the first two candidates ---
        for (int k = 0; k < 2; ++k) {
            const std::string prefix = "fw" + std::to_string(k + 1) + "_";
            sn[prefix + "p"]     = scan_fw_p[k];
            sn[prefix + "theta"] = scan_fw_theta[k];
            sn[prefix + "phi"]   = scan_fw_phi[k];
        }

        sn.fill();
    }

}

// ============================================================================
// MAIN PROGRAM
// ============================================================================

int main(int argc, char* argv[]) {
    // Install signal handler for graceful Ctrl+C termination
    SignalHandler::install();

    // Count ROOT error messages (see FAILURE HANDLING at the top of this file)
    SetErrorHandler(countingErrorHandler);

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

    // A file of the list that is missing, unreadable or has no such tree gives
    // only a ROOT error message; the chain would skip the file silently.
    if (g_root_errors > 0) {
        std::cerr << "\nERROR: ROOT reported " << g_root_errors
                  << " error(s) while the input was opened (messages above):"
                  << " a file of the input list is missing, unreadable or has no tree '"
                  << config.getInputTreeName() << "'. No output file was written.\n";
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

    // Resolve every PAT branch used in processEvent now, before the output
    // file exists. A missing branch stops the program here (exit code 2).
    {
        const std::vector<std::string> required = requiredReaderVariables(config);
        std::vector<std::string> missing;
        for (const auto& var : required) {
            try {
                reader.bindVariables({var});
            } catch (const std::exception& e) {
                missing.push_back(e.what());
            }
        }
        if (!missing.empty()) {
            std::cerr << "\nERROR: " << missing.size() << " of " << required.size()
                      << " input variables used by the analysis are missing in tree '"
                      << config.getInputTreeName() << "':\n";
            for (const auto& msg : missing) std::cerr << "  " << msg << "\n";
            std::cerr << "Check input.tree_name and input.lepton_prefixes in "
                      << config_file << ". No output file was written.\n";
            return 2;
        }
        if (g_root_errors > 0) {
            std::cerr << "\nERROR: ROOT reported " << g_root_errors
                      << " error(s) while the input variables were bound (messages above)."
                      << " No output file was written.\n";
            return 2;
        }
        std::cout << "NTupleReader: all " << required.size()
                  << " input variables used by the analysis are bound.\n";
    }
    const size_t n_bound_before_loop = reader.boundVariableCount();

    // Open output file and setup analysis components
    const std::string output_filename = config.getOutputFilename();
    Manager manager;
    CutManager cuts;
    try {
        manager.openFile(output_filename, config.getOutputOption());
    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << "\n";
        return 5;   // nothing was created
    }
    try {
        setupHistograms(manager);
        setupNtuples(manager, config);
        setupCuts(cuts);
    } catch (const std::exception& e) {
        std::cerr << "\nERROR while setting up histograms, ntuples or cuts: " << e.what() << "\n";
        abandonOutput(manager, output_filename, 5, true);
    }

    // Event loop
    Long64_t total_entries = reader.entries();
    Long64_t start_event = config.getStartEvent();
    Long64_t max_events = config.getMaxEvents();

    Long64_t end_event = total_entries;
    if (max_events > 0) {
        end_event = std::min(start_event + max_events, total_entries);
    }

    Long64_t events_to_process = end_event - start_event;

    // ------------------------------------------------------------------
    // Pass 1 (only if PT3 trigger-bias correction is enabled).
    //   Counts trigbit==8192 (PT3) and trigbit==4096 (PT2) per input
    //   ROOT file (= per HADES run) — beam conditions and trigger rates
    //   drift between runs, so a single global w would carry a few-%
    //   systematic. Per-file:
    //
    //       w[file] = (63 · N_PT2[file]) / N_PT3[file]
    //
    //   Pass 2 then sets the per-event weight to w[current_file]; every
    //   histogram fill is multiplied by it and the trigger_corr branch
    //   on every dynamic ntuple records the same number per event.
    //   When the correction is disabled, this pass is skipped and every
    //   event keeps weight 1.0.
    // ------------------------------------------------------------------
    const bool do_trigger_corr = config.isTriggerBiasCorrectionEnabled();
    const int  n_files = reader.getNTrees();
    std::vector<double> w_per_file(n_files, 1.0);

    // Entries that could not be read (getEntry returns <= 0): the variables
    // would keep the values of the previous row. Counted here, fatal after the loop.
    Long64_t n_read_errors = 0;
    const Long64_t max_errors_printed = 5;
    auto readEntry = [&](Long64_t entry) {
        if (reader.getEntry(entry) > 0) return true;
        ++n_read_errors;
        if (n_read_errors <= max_errors_printed) {
            std::cerr << "\nERROR: entry " << entry << " could not be read from the input.\n";
        }
        return false;
    };

    if (do_trigger_corr) {
        ConsoleBox::newLine();
        std::cout << "Pass 1/2: per-file trigger-bias counts (PT3=8192, PT2=4096)...\n\n";

        std::vector<Long64_t> n_PT3(n_files, 0), n_PT2(n_files, 0);
        ProgressBar pre_progress(events_to_process);
        Long64_t pre_processed = 0;
        for (Long64_t i = start_event; i < end_event; ++i) {
            if (SignalHandler::wasInterrupted()) break;
            const bool entry_ok = readEntry(i);
            ++pre_processed;
            pre_progress.update(pre_processed);
            if (!entry_ok) continue;
            const int fi      = reader.getCurrentTreeNumber();
            const int trigbit = static_cast<int>(reader["trigbit"]);
            if      (trigbit == 8192) ++n_PT3[fi];
            else if (trigbit == 4096) ++n_PT2[fi];
        }
        pre_progress.finish(SignalHandler::wasInterrupted());

        std::cout << "\nPer-file PT3 trigger-bias correction"
                  << " w = (63·N_PT2)/N_PT3:\n";
        for (int f = 0; f < n_files; ++f) {
            if (n_PT3[f] > 0 && n_PT2[f] > 0) {
                w_per_file[f] = (63.0 * static_cast<double>(n_PT2[f]))
                                       / static_cast<double>(n_PT3[f]);
            }
            std::cout << "  file " << std::setw(3) << f
                      << "  N_PT3=" << std::setw(10) << n_PT3[f]
                      << "  N_PT2=" << std::setw(8)  << n_PT2[f]
                      << "  w=" << std::fixed << std::setprecision(5)
                      << w_per_file[f] << "\n";
        }
    } else {
        std::cout << "PT3 trigger-bias correction DISABLED"
                  << " (config.trigger.bias_correction == false)."
                  << " All event weights = 1.\n";
    }

    ConsoleBox::newLine();
    ConsoleBox::printInfoBox("Press Ctrl+C at any time to stop and save partial results");
    ConsoleBox::newLine();
    std::cout << (do_trigger_corr ? "Pass 2/2: " : "")
              << "Processing events " << start_event << " to " << end_event
              << " (" << events_to_process << " events)...\n\n";

    Long64_t processed = 0;
    Long64_t n_exceptions = 0;          // entries on which processEvent threw
    const Long64_t max_exceptions_printed = max_errors_printed;
    bool was_interrupted = false;

    ProgressBar progress(events_to_process);

    for (Long64_t i = start_event; i < end_event; ++i) {
        if (SignalHandler::wasInterrupted()) {
            was_interrupted = true;
            break;
        }

        const bool entry_ok = readEntry(i);
        ++processed;
        progress.update(processed);
        if (!entry_ok) continue;

        // Set per-event weight (1.0 when correction disabled).
        manager.setEventWeight(w_per_file[reader.getCurrentTreeNumber()]);

        // An exception leaves the row half-written (some ntuples filled, others
        // not), so it is never ignored: counted here, fatal after the loop.
        try {
            processEvent(reader, manager, cuts, config);
        } catch (const std::exception& e) {
            ++n_exceptions;
            if (n_exceptions <= max_exceptions_printed) {
                std::cerr << "\nERROR: exception in processEvent at entry " << i
                          << ": " << e.what() << "\n";
            }
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

    // Failures found in the loop: report, remove the output file (closeFile()
    // also removes the intermediate ntuple files) and exit with a non-zero code.
    int fail_code = 0;

    // Entries that could not be read: the output is incomplete (code 5).
    if (n_read_errors > 0) {
        std::cerr << "\nERROR: " << n_read_errors << " entries could not be read from the input"
                  << " (the first " << std::min(n_read_errors, max_errors_printed)
                  << " are printed above).\n";
        fail_code = 5;
    }

    // A variable bound only inside the loop is missing in requiredReaderVariables():
    // the check of the input branches before the loop was not complete (code 2).
    if (reader.boundVariableCount() != n_bound_before_loop) {
        std::cerr << "\nERROR: " << (reader.boundVariableCount() - n_bound_before_loop)
                  << " input variable(s) were bound inside the event loop;"
                  << " add them to requiredReaderVariables() in main.cc.\n";
        fail_code = 2;
    }

    // Exceptions in the event loop: rows are half-written (code 3).
    if (n_exceptions > 0) {
        std::cerr << "\nERROR: processEvent threw an exception on " << n_exceptions
                  << " of " << processed << " entries (the first "
                  << std::min(n_exceptions, max_exceptions_printed)
                  << " are printed above).\n";
        fail_code = 3;
    }

    if (fail_code != 0) {
        abandonOutput(manager, output_filename, fail_code, true);
    }

    std::cout << "\nSaving results to " << output_filename << "...\n";
    manager.printSummary();
    try {
        manager.closeFile();
    } catch (const std::exception& e) {
        std::cerr << "\nERROR while writing the output file: " << e.what() << "\n";
        abandonOutput(manager, output_filename, 5, false);
    }

    // ROOT errors during the run (read error, full disk, ...) do not stop ROOT.
    if (g_root_errors > 0) {
        std::cerr << "\nERROR: ROOT reported " << g_root_errors
                  << " error(s) during the run (messages above); the output cannot be trusted.\n";
        abandonOutput(manager, output_filename, 5, false);
    }

    // Final check of the written file
    const std::string output_problem = verifyOutputFile(output_filename, manager);
    if (!output_problem.empty() || g_root_errors > 0) {
        std::cerr << "\nERROR: check of the output file failed: "
                  << (output_problem.empty() ? "ROOT error while reading it back" : output_problem)
                  << "\n";
        abandonOutput(manager, output_filename, 5, false);
    }
    std::cout << "Output file checked: every ntuple has as many rows as were filled.\n";

    // Interrupted: the partial output is kept, but the job is not a success.
    // ("Analysis Complete!" is not printed: monitor_parallel.sh reads it as done.)
    if (was_interrupted) {
        std::cerr << "\nINTERRUPTED: the output file holds only the first " << processed
                  << " of " << events_to_process << " entries (exit code 130).\n";
        return 130;
    }

    ConsoleBox::newLine();
    ConsoleBox::printStatus("Analysis Complete!");
    ConsoleBox::newLine();

    return 0;
}

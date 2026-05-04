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
    if (!cuts.passValueCut("trigger_PT3", reader["trigbit"])) return;
    if (!cuts.passCutSet("start_detector", { reader["start_iteration"] })) return;

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

    // Lepton momentum histograms (RECONSTRUCTED + CORRECTED)
    mgr.fill("ep_p", positron.momentum());
    mgr.fill("em_p", electron.momentum());
    mgr.fill("ep_p_cor", positron.momentum(KinematicType::CORRECTED));
    mgr.fill("em_p_cor", electron.momentum(KinematicType::CORRECTED));

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

    // Dilepton (e+ + e-) — operator+ propagates both REC and COR kinematics
    PParticle dilepton = positron + electron;
    double m_ee = dilepton.massGeV();
    double m_ee_cor = dilepton.massGeV(KinematicType::CORRECTED);

    mgr.fill("mass_ee_before_oa", m_ee);
    mgr.fill("mass_ee_before_oa_cor", m_ee_cor);

    // Boost dilepton to beam-target CMS frame (boost propagates both kinematics)
    EventFrames frames;
    frames.setBeamFrameFromKineticEnergy(config.getBeamKineticEnergy());

    PParticle dilepton_cms = frames.getFrame("beam").boost(dilepton);

    double y_cms = dilepton_cms.rapidity();
    double pt = dilepton_cms.vec().Pt();
    double theta_cms = dilepton_cms.theta();

    double y_cms_cor = dilepton_cms.rapidity(KinematicType::CORRECTED);
    double pt_cor = dilepton_cms.vec(KinematicType::CORRECTED).Pt();
    double theta_cms_cor = dilepton_cms.theta(KinematicType::CORRECTED);

    // Fill dilepton ntuple (before OA cut, store cut decision as flag).
    // OA > 4 deg is the active analysis selection.
    bool oa_pass = cuts.passMinCut("opening_angle_4", oa);

    // Helper to fill either dilepton_nt (REC) or dilepton_nt_cor (CORRECTED) —
    // identical field layout, only compound observables (m_ee, CMS) differ.
    auto fillDileptonNt = [&](const char* nt_name,
                              double m_ee_v, double y_v, double pt_v, double theta_v) {
        auto& nt = mgr.getDynamicNtuple(nt_name);
        nt["ep_p_rec"]      = ep_p_rec;
        nt["ep_p_cor"]      = ep_p_cor;
        nt["ep_theta"]      = positron.theta();
        nt["ep_phi"]        = positron.phi();
        nt["ep_theta_rich"] = reader["ep_theta_rich"];
        nt["ep_phi_rich"]   = reader["ep_phi_rich"];

        nt["em_p_rec"]      = em_p_rec;
        nt["em_p_cor"]      = em_p_cor;
        nt["em_theta"]      = electron.theta();
        nt["em_phi"]        = electron.phi();
        nt["em_theta_rich"] = reader["em_theta_rich"];
        nt["em_phi_rich"]   = reader["em_phi_rich"];

        nt["oa"]            = oa;
        nt["m_ee"]          = m_ee_v;

        nt["y_cms"]         = y_v;
        nt["pt"]            = pt_v;
        nt["theta_cms"]     = theta_v;

        nt["oa_pass"]       = oa_pass ? 1.0f : 0.0f;
        nt.fill();
    };

    fillDileptonNt("dilepton_nt",     m_ee,     y_cms,     pt,     theta_cms);
    fillDileptonNt("dilepton_nt_cor", m_ee_cor, y_cms_cor, pt_cor, theta_cms_cor);

    // Apply opening angle cut (reject close pairs)
    // if (!oa_pass) return;
    if (oa_pass) {
        // Still fill histograms for failed OA cut for comparison
        mgr.fill("mass_ee_after_oa", m_ee);
        mgr.fill("mass_ee_after_oa_cor", m_ee_cor);
    }
    // CMS histograms (after OA cut)
    mgr.fill("rapidity_cms", y_cms);
    mgr.fill("pt_cms", pt);
    mgr.fill("theta_cms", theta_cms);
    mgr.fill("rapidity_vs_mass", m_ee, y_cms);

    // Dilepton invariant mass (after all cuts)
    mgr.fill("mass_ee", m_ee);
    mgr.fill("mass_ee_cor", m_ee_cor);

    // Beam + target for missing mass calculations
    PParticle beam = ParticleFactory::createBeamProton(config.getBeamKineticEnergy());
    PParticle target = ParticleFactory::createTargetProton();
    PParticle initial = beam + target;

    // Missing mass of e+e- system: MM(e+e-) = beam + target - e+ - e-
    // operator- propagates both REC and COR (synthetic beam/target carry both).
    PParticle miss_epem = initial - positron - electron;
    double mm_epem_mass      = miss_epem.massGeV();
    double mm_epem_mass2     = miss_epem.vec().M2() / 1e6;  // GeV²/c⁴
    double mm_epem_mass_cor  = miss_epem.massGeV(KinematicType::CORRECTED);
    double mm_epem_mass2_cor = miss_epem.vec(KinematicType::CORRECTED).M2() / 1e6;

    // ECAL objects (electromagnetic calorimeter, up to 5 hits)
    if (config.isEcalEnabled()) {
        int neutr_mult = static_cast<int>(reader["neutr_mult"]);

        std::vector<PParticleEcal> ecal_objects;
        ecal_objects.reserve(neutr_mult);

        for (int i = 1; i <= neutr_mult && i <= 5; ++i) {
            PParticleEcal ecal_obj(0.0, "ecal" + std::to_string(i));
            if (ecal_obj.setFromReader(reader, i)) {
                // Photon has no measured momentum correction — mirror REC into
                // CORRECTED so dilepton+gamma compounds propagate both kinematics
                // (the difference between REC/COR comes solely from the leptons).
                ecal_obj.setFromSpherical(ecal_obj.vec().P(),
                                          ecal_obj.cluster_theta,
                                          ecal_obj.cluster_phi,
                                          KinematicType::CORRECTED);
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

        // Compound e+e-gamma — ecal_mult == 1, single neutral hit available.
        // Two fills:
        //   meson_dalitz_nt      — UNCONDITIONAL on ecal_pass; quality stored
        //                          as a flag, so RICH/ECAL ID effects can be
        //                          studied also for events that fail quality.
        //   epemg_nt(_cor) +     — only when gamma passes ecal_quality
        //   mass_epemg histos      (existing pi0/eta Dalitz pipeline).
        if (neutr_mult == 1 && !ecal_objects.empty()) {
            const auto& gamma = ecal_objects[0];

            PParticle epemg = dilepton + gamma;
            PParticle epemg_cms = frames.getFrame("beam").boost(epemg);
            PParticle miss_epemg = initial - positron - electron - gamma;

            double m_epemg     = epemg.massGeV();
            double m_epemg_cor = epemg.massGeV(KinematicType::CORRECTED);
            double oa_epem_g   = Physics::openingAngle(dilepton, gamma);

            // === meson_dalitz_nt — Dalitz pi0/eta study with RICH+ECAL ========
            // Mult==1 only; no ecal_quality gating — quality is recorded as
            // a flag. Lepton kinematics here are RECONSTRUCTED.
            {
                auto& nt = mgr.getDynamicNtuple("meson_dalitz_nt");

                // Compound observables (REC kinematics)
                nt["m_ee"]               = m_ee;
                nt["m_epemg"]            = m_epemg;
                nt["epemg_theta"]        = epemg.theta();
                nt["epemg_phi"]          = epemg.phi();

                // Lepton kinematics — momentum REC + COR, angles (REC ≡ COR
                // for theta/phi: energy-loss correction adjusts |p| only).
                nt["ep_p_rec"]           = ep_p_rec;
                nt["ep_p_cor"]           = ep_p_cor;
                nt["ep_theta"]           = positron.theta();
                nt["ep_phi"]             = positron.phi();

                nt["em_p_rec"]           = em_p_rec;
                nt["em_p_cor"]           = em_p_cor;
                nt["em_theta"]           = electron.theta();
                nt["em_phi"]             = electron.phi();

                // Mass-window flags applied to M(epemg)
                nt["pi0_pass"]           = cuts.passRangeCut("pi0_mass_window",        m_epemg) ? 1.0f : 0.0f;
                nt["pi0_pass_narrow"]    = cuts.passRangeCut("pi0_mass_window_narrow", m_epemg) ? 1.0f : 0.0f;
                nt["eta_pass"]           = cuts.passRangeCut("eta_mass_window",        m_epemg) ? 1.0f : 0.0f;

                // Opening angles
                nt["oa_epem"]            = oa;
                nt["oa_epem_g"]          = oa_epem_g;

                // ECAL quality cut decision (1 = passed quality)
                nt["ecal_quality_pass"]  = ecal_pass[0] ? 1.0f : 0.0f;

                // Neutral-cluster diagnostics from input ntuple (mult==1: _1 suffix
                // stripped on output, since this ntuple is mult==1-only).
                nt["neutr_mult"]            = reader["neutr_mult"];
                nt["neutr_counter"]         = reader["neutr_counter"];
                nt["neutr_pid"]             = reader["neutr_pid_1"];
                nt["neutr_tof"]             = reader["neutr_tof_1"];
                nt["neutr_tofrec"]          = reader["neutr_tofrec_1"];
                nt["neutr_dist"]            = reader["neutr_dist_1"];
                nt["neutr_clusterid"]       = reader["neutr_clusterid_1"];
                nt["neutr_cluster_theta"]   = reader["neutr_cluster_theta_1"];
                nt["neutr_cluster_phi"]     = reader["neutr_cluster_phi_1"];
                nt["neutr_cluster_energy"]  = reader["neutr_cluster_energy_1"];
                nt["neutr_cluster_ncells"]  = reader["neutr_cluster_ncells_1"];
                nt["neutr_beta"]            = reader["neutr_beta_1"];
                nt["neutr_p"]               = reader["neutr_p_1"];
                nt["neutr_p_pid"]           = reader["neutr_p_pid_1"];
                nt["neutr_mass"]            = reader["neutr_mass_1"];
                nt["neutr_mass2"]           = reader["neutr_mass2_1"];
                nt["neutr_q"]               = reader["neutr_q_1"];
                nt["neutr_phi"]             = reader["neutr_phi_1"];
                nt["neutr_theta"]           = reader["neutr_theta_1"];
                nt["neutr_r"]               = reader["neutr_r_1"];
                nt["neutr_z"]               = reader["neutr_z_1"];
                nt["neutr_chi2"]            = reader["neutr_chi2_1"];
                nt["neutr_phi2"]            = reader["neutr_phi2_1"];
                nt["neutr_theta2"]          = reader["neutr_theta2_1"];
                nt["neutr_r2"]              = reader["neutr_r2_1"];
                nt["neutr_z2"]              = reader["neutr_z2_1"];
                nt["neutr_energy"]          = reader["neutr_energy_1"];

                // RICH ring parameters — separately for ep and em
                nt["ep_rich_amp"]              = reader["ep_rich_amp"];
                nt["ep_rich_avg_ringcharge"]   = reader["ep_rich_avg_ringcharge"];
                nt["ep_rich_padnum"]           = reader["ep_rich_padnum"];
                nt["ep_rich_centr"]            = reader["ep_rich_centr"];
                nt["ep_rich_patmat"]           = reader["ep_rich_patmat"];
                nt["ep_rich_houtra"]           = reader["ep_rich_houtra"];
                nt["ep_rich_radius"]           = reader["ep_rich_radius"];

                nt["em_rich_amp"]              = reader["em_rich_amp"];
                nt["em_rich_avg_ringcharge"]   = reader["em_rich_avg_ringcharge"];
                nt["em_rich_padnum"]           = reader["em_rich_padnum"];
                nt["em_rich_centr"]            = reader["em_rich_centr"];
                nt["em_rich_patmat"]           = reader["em_rich_patmat"];
                nt["em_rich_houtra"]           = reader["em_rich_houtra"];
                nt["em_rich_radius"]           = reader["em_rich_radius"];

                nt.fill();
            }

            // === Existing pi0/eta Dalitz pipeline — quality-gated ===========
            if (ecal_pass[0]) {
                mgr.fill("mass_epemg", m_epemg);
                mgr.fill("mass_epemg_cor", m_epemg_cor);

                // Helper — fills epemg_nt (REC) or epemg_nt_cor (CORRECTED).
                auto fillEpemgNt = [&](const char* nt_name, KinematicType k) {
                    auto& nt = mgr.getDynamicNtuple(nt_name);
                    nt["epemg_mass"]         = epemg.massGeV(k);
                    nt["epemg_p"]            = epemg.momentum(k);
                    nt["epemg_theta"]        = epemg.theta(k);
                    nt["epemg_phi"]          = epemg.phi(k);

                    nt["epemg_rapidity_cms"] = epemg_cms.rapidity(k);
                    nt["epemg_pt_cms"]       = epemg_cms.vec(k).Pt();
                    nt["epemg_theta_cms"]    = epemg_cms.theta(k);

                    nt["ee_oa"]              = oa;
                    nt["ee_mass"]            = (k == KinematicType::CORRECTED) ? m_ee_cor : m_ee;

                    nt["gamma_energy"]       = gamma.cluster_energy;
                    nt["gamma_theta"]        = gamma.cluster_theta;
                    nt["gamma_phi"]          = gamma.cluster_phi;

                    nt["ecal_mult"]          = static_cast<Float_t>(neutr_mult);

                    if (k == KinematicType::CORRECTED) {
                        nt["mm_epem_mass"]   = mm_epem_mass_cor;
                        nt["mm_epem_mass2"]  = mm_epem_mass2_cor;
                    } else {
                        nt["mm_epem_mass"]   = mm_epem_mass;
                        nt["mm_epem_mass2"]  = mm_epem_mass2;
                    }

                    nt["mm_epemg_mass"]      = miss_epemg.massGeV(k);
                    nt["mm_epemg_mass2"]     = miss_epemg.vec(k).M2() / 1e6;
                    nt.fill();
                };

                fillEpemgNt("epemg_nt",     KinematicType::RECONSTRUCTED);
                fillEpemgNt("epemg_nt_cor", KinematicType::CORRECTED);
            }
        }

        // ===========================================================
        // MULT == 2 — pi0 candidate from gg, omega candidate from epemgg
        // M(gg) is REC≡COR (photon mirrored). M(epemgg) differs via leptons.
        // ===========================================================
        if (neutr_mult == 2 && ecal_objects.size() >= 2 &&
            ecal_pass[0] && ecal_pass[1]) {

            PParticle gg     = ecal_objects[0] + ecal_objects[1];
            PParticle epemgg = dilepton + gg;

            double m_gg         = gg.massGeV();              // REC≡COR
            double m_epemgg     = epemgg.massGeV();
            double m_epemgg_cor = epemgg.massGeV(KinematicType::CORRECTED);

            mgr.fill("mass_gg", m_gg);

            bool pi0_pass_narrow = cuts.passRangeCut("pi0_mass_window_narrow", m_gg);

            // REC ntuple
            auto& nt = mgr.getDynamicNtuple("epemgg_nt");
            nt["m_ee"]            = m_ee;
            nt["m_gg"]            = m_gg;
            nt["m_epemgg"]        = m_epemgg;
            nt["pi0_pass_narrow"] = pi0_pass_narrow ? 1.0f : 0.0f;
            nt.fill();

            // CORRECTED mirror — same field names, m_ee/m_epemgg from COR
            auto& nt_cor = mgr.getDynamicNtuple("epemgg_nt_cor");
            nt_cor["m_ee"]            = m_ee_cor;
            nt_cor["m_gg"]            = m_gg;
            nt_cor["m_epemgg"]        = m_epemgg_cor;
            nt_cor["pi0_pass_narrow"] = pi0_pass_narrow ? 1.0f : 0.0f;
            nt_cor.fill();
        }

        // ===========================================================
        // MULT == 3 — 3 rotational combinations
        //   epemg_i = dilepton + g_i,  gg_jk = g_j + g_k,  full = epemg_i + gg_jk
        // One ntuple row per rotation (3 rows per qualifying event).
        // M(gg) is REC≡COR; M(epemg) and M(epemggg) differ between REC and COR.
        // Cut decisions evaluated separately for REC and COR (m_epemg differs).
        // ===========================================================
        if (neutr_mult == 3 && ecal_objects.size() >= 3 &&
            ecal_pass[0] && ecal_pass[1] && ecal_pass[2]) {

            const int rot[3][3] = {{0,1,2}, {1,2,0}, {2,0,1}};
            auto& nt     = mgr.getDynamicNtuple("epemggg_nt");
            auto& nt_cor = mgr.getDynamicNtuple("epemggg_nt_cor");

            for (int r = 0; r < 3; ++r) {
                int i = rot[r][0], j = rot[r][1], k = rot[r][2];

                PParticle epemg_i  = dilepton + ecal_objects[i];
                PParticle gg_jk    = ecal_objects[j] + ecal_objects[k];
                PParticle epemggg  = epemg_i + gg_jk;

                double m_epemg       = epemg_i.massGeV();
                double m_epemg_cor   = epemg_i.massGeV(KinematicType::CORRECTED);
                double m_gg_jk       = gg_jk.massGeV();                       // REC≡COR
                double m_epemggg     = epemggg.massGeV();
                double m_epemggg_cor = epemggg.massGeV(KinematicType::CORRECTED);

                bool pi0_pass_narrow     = cuts.passRangeCut("pi0_mass_window_narrow", m_gg_jk);
                bool eta_pass_rec        = cuts.passRangeCut("eta_mass_window",        m_epemg);
                bool eta_pass_cor        = cuts.passRangeCut("eta_mass_window",        m_epemg_cor);

                // REC row
                nt["m_ee"]            = m_ee;
                nt["m_epemg"]         = m_epemg;
                nt["m_gg"]            = m_gg_jk;
                nt["m_epemggg"]       = m_epemggg;
                nt["rot_idx"]         = static_cast<Float_t>(r);
                nt["pi0_pass_narrow"] = pi0_pass_narrow ? 1.0f : 0.0f;
                nt["eta_pass"]        = eta_pass_rec    ? 1.0f : 0.0f;
                nt.fill();

                // CORRECTED row
                nt_cor["m_ee"]            = m_ee_cor;
                nt_cor["m_epemg"]         = m_epemg_cor;
                nt_cor["m_gg"]            = m_gg_jk;
                nt_cor["m_epemggg"]       = m_epemggg_cor;
                nt_cor["rot_idx"]         = static_cast<Float_t>(r);
                nt_cor["pi0_pass_narrow"] = pi0_pass_narrow ? 1.0f : 0.0f;
                nt_cor["eta_pass"]        = eta_pass_cor    ? 1.0f : 0.0f;
                nt_cor.fill();
            }
        }

        // ===========================================================
        // MULT == 4 — 4-gamma combinatorics under (pi0, pi0) constraint
        //   3 pairings: (01)(23), (02)(13), (03)(12). Both gg pairs must
        //   sit in the narrow pi0 window. epem is irrelevant here — this
        //   is purely an ECAL combinatorial signal, so plotted directly
        //   from output_epem_exp.root without CB extraction.
        // ===========================================================
        if (neutr_mult == 4 && ecal_objects.size() >= 4 &&
            ecal_pass[0] && ecal_pass[1] && ecal_pass[2] && ecal_pass[3]) {

            const int pairings[3][4] = {{0,1, 2,3}, {0,2, 1,3}, {0,3, 1,2}};

            for (int p = 0; p < 3; ++p) {
                PParticle gg_a = ecal_objects[pairings[p][0]] + ecal_objects[pairings[p][1]];
                PParticle gg_b = ecal_objects[pairings[p][2]] + ecal_objects[pairings[p][3]];

                bool a_in_pi0 = cuts.passRangeCut("pi0_mass_window_narrow", gg_a.massGeV());
                bool b_in_pi0 = cuts.passRangeCut("pi0_mass_window_narrow", gg_b.massGeV());

                if (a_in_pi0 && b_in_pi0) {
                    PParticle gggg = gg_a + gg_b;
                    mgr.fill("mass_gggg_pi0pi0", gggg.massGeV());
                }
            }
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

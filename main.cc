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

    // 2D graphical cut on (mm_pippimepem, m_pippimepem) — separate decisions for REC and COR
    bool cut2d_pass     = cuts.passGraphicalCut("cut_2d", mm_pippimepem,     m_pippimepem);
    bool cut2d_pass_cor = cuts.passGraphicalCut("cut_2d", mm_pippimepem_cor, m_pippimepem_cor);
    if (cut2d_pass) {
        mgr.fill("mm_epem_cut2d", mm_epem);
        mgr.fill("mm_pippim_cut2d", mm_pippim);
        mgr.fill("mm_pippimepem_cut2d", mm_pippimepem);
        mgr.fill("mm_vs_m_pippimepem_cut2d", mm_pippimepem, m_pippimepem);
    }
    if (cut2d_pass_cor) {
        mgr.fill("mm_epem_cut2d_cor", mm_epem_cor);
        mgr.fill("mm_pippim_cut2d_cor", mm_pippim_cor);
        mgr.fill("mm_pippimepem_cut2d_cor", mm_pippimepem_cor);
        mgr.fill("mm_vs_m_pippimepem_cut2d_cor", mm_pippimepem_cor, m_pippimepem_cor);
    }

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

    // Control histograms — filled BEFORE pippimepem_selection so the cut
    // boundaries (oa_pippim_epem_lab < 50, oa_pippim_epem_eta_rest > 140) can
    // be verified directly on these spectra.
    mgr.fill("oa_pippim_epem_lab",          oa_pippim_epem_lab);
    mgr.fill("oa_pippim_epem_eta_rest",     oa_pippim_epem_eta_rest);
    mgr.fill("oa_pippim_epem_lab_cor",      oa_pippim_epem_lab_cor);
    mgr.fill("oa_pippim_epem_eta_rest_cor", oa_pippim_epem_eta_rest_cor);

    // 4-body selection (RECONSTRUCTED): OA_LAB < 50, M(pi+pi-) < 0.420, eta-rest OA > 140
    bool sel_pass = cuts.passCutSet("pippimepem_selection", {
        oa_pippim_epem_lab,
        m_pippim,
        oa_pippim_epem_eta_rest
    });
    // Same selection with CORRECTED inputs
    bool sel_pass_cor = cuts.passCutSet("pippimepem_selection", {
        oa_pippim_epem_lab_cor,
        m_pippim_cor,
        oa_pippim_epem_eta_rest_cor
    });

    if (sel_pass) {
        mgr.fill("mass_pippimepem_selected", m_pippimepem);
        mgr.fill("mm_vs_m_pippimepem_selected", mm_pippimepem, m_pippimepem);

        // MM(pi+pi-e+e-) slice scan (selection + one slice window)
        if (cuts.passRangeCut("mm_slice_20_22", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_20_22", m_pippimepem);
        if (cuts.passRangeCut("mm_slice_22_24", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_22_24", m_pippimepem);
        if (cuts.passRangeCut("mm_slice_24_26", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_24_26", m_pippimepem);
        if (cuts.passRangeCut("mm_slice_26_28", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_26_28", m_pippimepem);
        if (cuts.passRangeCut("mm_slice_28_30", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_28_30", m_pippimepem);

        // Same slices, additionally requiring cut_2d
        if (cut2d_pass) {
            if (cuts.passRangeCut("mm_slice_20_22", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_20_22_cut2d", m_pippimepem);
            if (cuts.passRangeCut("mm_slice_22_24", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_22_24_cut2d", m_pippimepem);
            if (cuts.passRangeCut("mm_slice_24_26", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_24_26_cut2d", m_pippimepem);
            if (cuts.passRangeCut("mm_slice_26_28", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_26_28_cut2d", m_pippimepem);
            if (cuts.passRangeCut("mm_slice_28_30", mm_pippimepem)) mgr.fill("mass_pippimepem_slice_28_30_cut2d", m_pippimepem);
        }
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

        if (cut2d_pass_cor) {
            if (cuts.passRangeCut("mm_slice_20_22", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_20_22_cut2d_cor", m_pippimepem_cor);
            if (cuts.passRangeCut("mm_slice_22_24", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_22_24_cut2d_cor", m_pippimepem_cor);
            if (cuts.passRangeCut("mm_slice_24_26", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_24_26_cut2d_cor", m_pippimepem_cor);
            if (cuts.passRangeCut("mm_slice_26_28", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_26_28_cut2d_cor", m_pippimepem_cor);
            if (cuts.passRangeCut("mm_slice_28_30", mm_pippimepem_cor)) mgr.fill("mass_pippimepem_slice_28_30_cut2d_cor", m_pippimepem_cor);
        }
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

    // Fill pi+pi-e+e- ntuple (before OA cut, store cut decision as flag)
    // opening_angle_4 is the ACTIVE cut; opening_angle_9 is kept defined for future use.
    bool oa_pass = cuts.passMinCut("opening_angle_4", oa);

    auto& nt = mgr.getDynamicNtuple("pippimepem_nt");

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

    nt["pip_p_rec"] = pip_p_rec;
    nt["pip_p_cor"] = pip_p_cor;
    nt["pip_theta"] = piplus.theta();
    nt["pip_phi"] = piplus.phi();

    nt["pim_p_rec"] = pim_p_rec;
    nt["pim_p_cor"] = pim_p_cor;
    nt["pim_theta"] = piminus.theta();
    nt["pim_phi"] = piminus.phi();

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

    nt["oa_pass"] = oa_pass ? 1.0f : 0.0f;       // opening_angle_4 decision
    nt["sel_pass"] = sel_pass ? 1.0f : 0.0f;     // pippimepem_selection chain decision
    nt["cut2d_pass"] = cut2d_pass ? 1.0f : 0.0f; // TCutG cut_2d on (MM, M) of pi+pi-e+e-

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

    nt_cor["em_p_rec"] = em_p_rec;
    nt_cor["em_p_cor"] = em_p_cor;
    nt_cor["em_theta"] = electron.theta();
    nt_cor["em_phi"] = electron.phi();
    nt_cor["em_theta_rich"] = reader["em_theta_rich"];
    nt_cor["em_phi_rich"] = reader["em_phi_rich"];

    nt_cor["pip_p_rec"] = pip_p_rec;
    nt_cor["pip_p_cor"] = pip_p_cor;
    nt_cor["pip_theta"] = piplus.theta();
    nt_cor["pip_phi"] = piplus.phi();

    nt_cor["pim_p_rec"] = pim_p_rec;
    nt_cor["pim_p_cor"] = pim_p_cor;
    nt_cor["pim_theta"] = piminus.theta();
    nt_cor["pim_phi"] = piminus.phi();

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
    nt_cor["cut2d_pass"] = cut2d_pass_cor ? 1.0f : 0.0f;   // CORRECTED cut_2d

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

    mgr.fill("rapidity_cms_cor", y_cms_cor);
    mgr.fill("pt_cms_cor", pt_cor);
    mgr.fill("theta_cms_cor", theta_cms_cor);
    mgr.fill("rapidity_vs_mass_cor", m_ee_cor, y_cms_cor);

    // Dilepton invariant mass (after all cuts)
    mgr.fill("mass_ee", m_ee);
    mgr.fill("mass_ee_cor", m_ee_cor);

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

/**
 * @file setup_histograms.h
 * @brief Histogram definitions
 *
 * This file defines all histograms for the analysis.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef SETUP_HISTOGRAMS_H
#define SETUP_HISTOGRAMS_H

#include "manager.h"
#include <iostream>

/**
 * @brief Define histograms for the analysis
 * @param mgr Manager object for histogram creation
 */
inline void setupHistograms(Manager& mgr) {
    std::cout << "Setting up histograms...\n";

    // Lepton momentum
    mgr.create1D("ep_p", "e^{+} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "leptons");
    mgr.create1D("em_p", "e^{-} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "leptons");

    // Pion momentum
    mgr.create1D("pip_p", "#pi^{+} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "pions");
    mgr.create1D("pim_p", "#pi^{-} momentum;p [MeV/c];Counts",
                 100, 0, 2000, "pions");

    // Momentum correction: delta_p vs p_reconstructed
    mgr.create2D("ep_dp_vs_p", "e^{+}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");
    mgr.create2D("em_dp_vs_p", "e^{-}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");
    mgr.create2D("pip_dp_vs_p", "#pi^{+}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");
    mgr.create2D("pim_dp_vs_p", "#pi^{-}: #Deltap vs p;p_{rec} [MeV/c];p_{corr} - p_{rec} [MeV/c]",
                 100, 0, 2000, 100, 0, 10, "corrections");

    // Compound invariant masses
    mgr.create1D("mass_ee", "e^{+}e^{-} invariant mass;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "compound");
    mgr.create1D("mass_pippim", "#pi^{+}#pi^{-} invariant mass;M_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 2.0, "compound");
    mgr.create1D("mass_pippimepem", "#pi^{+}#pi^{-}e^{+}e^{-} invariant mass;M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 2.0, "compound");
    // With pippimepem_selection cut chain applied
    mgr.create1D("mass_pippimepem_selected",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) after selection;M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 2.0, "compound");

    // M(pi+pi-e+e-) in MM(pi+pi-e+e-) slice windows (selection chain + slice cut)
    mgr.create1D("mass_pippimepem_slice_20_22",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel, MM in [2.0,2.2] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");
    mgr.create1D("mass_pippimepem_slice_22_24",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel, MM in [2.2,2.4] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");
    mgr.create1D("mass_pippimepem_slice_24_26",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel, MM in [2.4,2.6] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");
    mgr.create1D("mass_pippimepem_slice_26_28",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel, MM in [2.6,2.8] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");
    mgr.create1D("mass_pippimepem_slice_28_30",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel, MM in [2.8,3.0] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");

    // Same slice histograms with the 2D graphical cut (cut_2d) additionally applied
    mgr.create1D("mass_pippimepem_slice_20_22_cut2d",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d, MM in [2.0,2.2] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");
    mgr.create1D("mass_pippimepem_slice_22_24_cut2d",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d, MM in [2.2,2.4] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");
    mgr.create1D("mass_pippimepem_slice_24_26_cut2d",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d, MM in [2.4,2.6] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");
    mgr.create1D("mass_pippimepem_slice_26_28_cut2d",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d, MM in [2.6,2.8] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");
    mgr.create1D("mass_pippimepem_slice_28_30_cut2d",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d, MM in [2.8,3.0] GeV/c^{2};M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices");

    // Missing masses: MM(X) = beam + target - X
    mgr.create1D("mm_epem", "MM(e^{+}e^{-});MM(e^{+}e^{-}) [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses");
    mgr.create1D("mm_pippim", "MM(#pi^{+}#pi^{-});MM(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses");
    mgr.create1D("mm_pippimepem", "MM(#pi^{+}#pi^{-}e^{+}e^{-});MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses");

    // 2D: MM (X) vs M (Y) of pi+pi-e+e-
    mgr.create2D("mm_vs_m_pippimepem",
                 "MM(#pi^{+}#pi^{-}e^{+}e^{-}) vs M(#pi^{+}#pi^{-}e^{+}e^{-});MM [GeV/c^{2}];M [GeV/c^{2}]",
                 1000, 0.0, 4.0, 1000, 0.0, 2.0, "missing_masses");
    // Same 2D with pippimepem_selection cut chain applied
    mgr.create2D("mm_vs_m_pippimepem_selected",
                 "MM vs M of #pi^{+}#pi^{-}e^{+}e^{-} after selection;MM [GeV/c^{2}];M [GeV/c^{2}]",
                 1000, 0.0, 4.0, 1000, 0.0, 2.0, "missing_masses");

    // Versions with the 2D graphical cut (cut_2d) applied
    mgr.create1D("mm_epem_cut2d", "MM(e^{+}e^{-}) after cut_2d;MM(e^{+}e^{-}) [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses");
    mgr.create1D("mm_pippim_cut2d", "MM(#pi^{+}#pi^{-}) after cut_2d;MM(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses");
    mgr.create1D("mm_pippimepem_cut2d", "MM(#pi^{+}#pi^{-}e^{+}e^{-}) after cut_2d;MM(#pi^{+}#pi^{-}e^{+}e^{-}) [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses");
    mgr.create2D("mm_vs_m_pippimepem_cut2d",
                 "MM vs M of #pi^{+}#pi^{-}e^{+}e^{-} after cut_2d;MM [GeV/c^{2}];M [GeV/c^{2}]",
                 1000, 0.0, 4.0, 1000, 0.0, 2.0, "missing_masses");

    // Opening angle and before/after OA cut mass spectra
    mgr.create1D("opening_angle", "e^{+}e^{-} opening angle;#theta_{open} [deg];Counts",
                 180, 0, 180, "dilepton");
    mgr.create1D("mass_ee_before_oa", "M_{ee} before OA cut;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");
    mgr.create1D("mass_ee_after_oa", "M_{ee} after OA cut;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton");

    // CMS frame observables (after boost to beam-target CMS)
    mgr.create1D("rapidity_cms", "Dilepton rapidity in CMS;y_{CMS};Counts",
                 100, -2.0, 2.0, "cms");
    mgr.create1D("pt_cms", "Dilepton transverse momentum;p_{T} [MeV/c];Counts",
                 100, 0, 1500, "cms");
    mgr.create1D("theta_cms", "Dilepton polar angle in CMS;#theta_{CMS} [deg];Counts",
                 90, 0, 180, "cms");
    mgr.create2D("rapidity_vs_mass", "y_{CMS} vs M_{ee};M_{e^{+}e^{-}} [GeV/c^{2}];y_{CMS}",
                 100, 0.0, 1.0, 100, -2.0, 2.0, "cms");

    // ========================================================================
    // CORRECTED kinematics: parallel set of histograms (suffix "_cor")
    // Per-particle theta/phi are direction-only, so opening_angle (single-pair)
    // and the dp_vs_p plots are NOT duplicated. Everything else is mirrored.
    // ========================================================================

    // Lepton momentum (CORRECTED)
    mgr.create1D("ep_p_cor", "e^{+} momentum (corrected);p [MeV/c];Counts",
                 100, 0, 2000, "leptons_cor");
    mgr.create1D("em_p_cor", "e^{-} momentum (corrected);p [MeV/c];Counts",
                 100, 0, 2000, "leptons_cor");

    // Pion momentum (CORRECTED)
    mgr.create1D("pip_p_cor", "#pi^{+} momentum (corrected);p [MeV/c];Counts",
                 100, 0, 2000, "pions_cor");
    mgr.create1D("pim_p_cor", "#pi^{-} momentum (corrected);p [MeV/c];Counts",
                 100, 0, 2000, "pions_cor");

    // Compound invariant masses (CORRECTED)
    mgr.create1D("mass_ee_cor", "e^{+}e^{-} invariant mass (cor);M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "compound_cor");
    mgr.create1D("mass_pippim_cor", "#pi^{+}#pi^{-} invariant mass (cor);M_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 2.0, "compound_cor");
    mgr.create1D("mass_pippimepem_cor", "#pi^{+}#pi^{-}e^{+}e^{-} invariant mass (cor);M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 2.0, "compound_cor");
    mgr.create1D("mass_pippimepem_selected_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) after selection (cor);M_{#pi^{+}#pi^{-}e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 2.0, "compound_cor");

    // Slice histograms (CORRECTED)
    mgr.create1D("mass_pippimepem_slice_20_22_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel cor, MM in [2.0,2.2];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    mgr.create1D("mass_pippimepem_slice_22_24_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel cor, MM in [2.2,2.4];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    mgr.create1D("mass_pippimepem_slice_24_26_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel cor, MM in [2.4,2.6];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    mgr.create1D("mass_pippimepem_slice_26_28_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel cor, MM in [2.6,2.8];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    mgr.create1D("mass_pippimepem_slice_28_30_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel cor, MM in [2.8,3.0];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    // Same slices, additionally with cut_2d (CORRECTED)
    mgr.create1D("mass_pippimepem_slice_20_22_cut2d_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d cor, MM in [2.0,2.2];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    mgr.create1D("mass_pippimepem_slice_22_24_cut2d_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d cor, MM in [2.2,2.4];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    mgr.create1D("mass_pippimepem_slice_24_26_cut2d_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d cor, MM in [2.4,2.6];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    mgr.create1D("mass_pippimepem_slice_26_28_cut2d_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d cor, MM in [2.6,2.8];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");
    mgr.create1D("mass_pippimepem_slice_28_30_cut2d_cor",
                 "M(#pi^{+}#pi^{-}e^{+}e^{-}) sel+cut_2d cor, MM in [2.8,3.0];M [GeV/c^{2}];Counts",
                 100, 0.0, 2.0, "slices_cor");

    // Missing masses (CORRECTED)
    mgr.create1D("mm_epem_cor", "MM(e^{+}e^{-}) (cor);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses_cor");
    mgr.create1D("mm_pippim_cor", "MM(#pi^{+}#pi^{-}) (cor);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses_cor");
    mgr.create1D("mm_pippimepem_cor", "MM(#pi^{+}#pi^{-}e^{+}e^{-}) (cor);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses_cor");
    mgr.create2D("mm_vs_m_pippimepem_cor",
                 "MM vs M of #pi^{+}#pi^{-}e^{+}e^{-} (cor);MM [GeV/c^{2}];M [GeV/c^{2}]",
                 1000, 0.0, 4.0, 1000, 0.0, 2.0, "missing_masses_cor");
    mgr.create2D("mm_vs_m_pippimepem_selected_cor",
                 "MM vs M of #pi^{+}#pi^{-}e^{+}e^{-} after selection (cor);MM [GeV/c^{2}];M [GeV/c^{2}]",
                 1000, 0.0, 4.0, 1000, 0.0, 2.0, "missing_masses_cor");
    mgr.create1D("mm_epem_cut2d_cor", "MM(e^{+}e^{-}) after cut_2d (cor);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses_cor");
    mgr.create1D("mm_pippim_cut2d_cor", "MM(#pi^{+}#pi^{-}) after cut_2d (cor);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses_cor");
    mgr.create1D("mm_pippimepem_cut2d_cor", "MM(#pi^{+}#pi^{-}e^{+}e^{-}) after cut_2d (cor);MM [GeV/c^{2}];Counts",
                 200, 0.0, 4.0, "missing_masses_cor");
    mgr.create2D("mm_vs_m_pippimepem_cut2d_cor",
                 "MM vs M of #pi^{+}#pi^{-}e^{+}e^{-} after cut_2d (cor);MM [GeV/c^{2}];M [GeV/c^{2}]",
                 1000, 0.0, 4.0, 1000, 0.0, 2.0, "missing_masses_cor");

    // Dilepton OA-cut spectra (CORRECTED). opening_angle itself is direction-only,
    // identical to RECONSTRUCTED, so it is NOT duplicated.
    mgr.create1D("mass_ee_before_oa_cor", "M_{ee} before OA cut (cor);M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton_cor");
    mgr.create1D("mass_ee_after_oa_cor", "M_{ee} after OA cut (cor);M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 200, 0.0, 1.0, "dilepton_cor");

    // CMS frame observables (CORRECTED)
    mgr.create1D("rapidity_cms_cor", "Dilepton rapidity in CMS (cor);y_{CMS};Counts",
                 100, -2.0, 2.0, "cms_cor");
    mgr.create1D("pt_cms_cor", "Dilepton transverse momentum (cor);p_{T} [MeV/c];Counts",
                 100, 0, 1500, "cms_cor");
    mgr.create1D("theta_cms_cor", "Dilepton polar angle in CMS (cor);#theta_{CMS} [deg];Counts",
                 90, 0, 180, "cms_cor");
    mgr.create2D("rapidity_vs_mass_cor", "y_{CMS} vs M_{ee} (cor);M_{e^{+}e^{-}} [GeV/c^{2}];y_{CMS}",
                 100, 0.0, 1.0, 100, -2.0, 2.0, "cms_cor");

    std::cout << "  Created " << mgr.histogramCount() << " histograms\n";
}

#endif // SETUP_HISTOGRAMS_H

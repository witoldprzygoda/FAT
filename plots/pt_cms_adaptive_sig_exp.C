// pt_cms_adaptive_sig_exp.C — adaptive-binning CB-subtracted dilepton p_T (CMS).
// Thin wrapper over AdaptiveSigPlot::run; see AdaptiveSigPlot.h for details.
//
// Branch: dilepton_nt::pt   (transverse momentum in MeV/c, in beam-target CMS)
//
// Usage:
//   root -l -b -q plots/pt_cms_adaptive_sig_exp.C
//   root -l -b -q 'plots/pt_cms_adaptive_sig_exp.C(500, 200)'

#include "AdaptiveSigPlot.h"

void pt_cms_adaptive_sig_exp(double      N_min_pt2         = 200.0,
                             double      max_bin_width     = 150.0,   // MeV/c (10% of range)
                             double      N_min_pt2_ratio   = 2000.0,
                             const char* nt_name           = "dilepton_nt") {
    AdaptiveSigPlot::Config cfg;
    cfg.var_expr          = "pt";
    cfg.x_axis_label      = "p_{T} [MeV/c]";
    cfg.y_axis_unit       = "MeV/c";
    cfg.xmin              = 0.0;
    cfg.xmax              = 1500.0;
    cfg.fine_nbins        = 1500;             // 1 MeV/c seed
    cfg.max_bin_width     = max_bin_width;
    cfg.N_min_pt2_spectra = N_min_pt2;
    cfg.N_min_pt2_ratio   = N_min_pt2_ratio;
    cfg.base_filename     = "pt_cms";
    cfg.canvas_title      = "p_{T} (CMS)";
    cfg.nt_name           = nt_name;
    cfg.log_y_spectra     = true;
    AdaptiveSigPlot::run(cfg);
}

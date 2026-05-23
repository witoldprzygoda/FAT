// theta_cms_adaptive_sig_exp.C — adaptive-binning CB-subtracted dilepton θ (CMS).
// Thin wrapper over AdaptiveSigPlot::run; see AdaptiveSigPlot.h for details.
//
// Branch: dilepton_nt::theta_cms   (polar angle in degrees, beam-target CMS)
//
// Usage:
//   root -l -b -q plots/theta_cms_adaptive_sig_exp.C
//   root -l -b -q 'plots/theta_cms_adaptive_sig_exp.C(500, 30)'

#include "AdaptiveSigPlot.h"

void theta_cms_adaptive_sig_exp(double      N_min_pt2         = 200.0,
                                double      max_bin_width     = 20.0,   // deg (~11% of range)
                                double      N_min_pt2_ratio   = 2000.0,
                                const char* nt_name           = "dilepton_nt") {
    AdaptiveSigPlot::Config cfg;
    cfg.var_expr          = "theta_cms";
    cfg.x_axis_label      = "#theta_{CMS} [deg]";
    cfg.y_axis_unit       = "deg";
    cfg.xmin              =   0.0;
    cfg.xmax              = 180.0;
    cfg.fine_nbins        = 1800;              // 0.1 deg/bin seed
    cfg.max_bin_width     = max_bin_width;
    cfg.N_min_pt2_spectra = N_min_pt2;
    cfg.N_min_pt2_ratio   = N_min_pt2_ratio;
    cfg.base_filename     = "theta_cms";
    cfg.canvas_title      = "#theta_{CMS}";
    cfg.nt_name           = nt_name;
    cfg.log_y_spectra     = false;
    AdaptiveSigPlot::run(cfg);
}

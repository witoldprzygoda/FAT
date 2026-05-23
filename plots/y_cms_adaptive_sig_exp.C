// y_cms_adaptive_sig_exp.C — adaptive-binning CB-subtracted dilepton rapidity (CMS).
// Thin wrapper over AdaptiveSigPlot::run; see AdaptiveSigPlot.h for details.
//
// Branch: dilepton_nt::y_cms   (dimensionless rapidity, beam-target CMS)
//
// Usage:
//   root -l -b -q plots/y_cms_adaptive_sig_exp.C
//   root -l -b -q 'plots/y_cms_adaptive_sig_exp.C(500, 0.5)'

#include "AdaptiveSigPlot.h"

void y_cms_adaptive_sig_exp(double      N_min_pt2         = 200.0,
                            double      max_bin_width     = 0.4,   // 10% of range
                            double      N_min_pt2_ratio   = 2000.0,
                            const char* nt_name           = "dilepton_nt") {
    AdaptiveSigPlot::Config cfg;
    cfg.var_expr          = "y_cms";
    cfg.x_axis_label      = "y_{CMS}";
    cfg.y_axis_unit       = "";                // dimensionless — label "Counts / bin"
    cfg.xmin              = -2.0;
    cfg.xmax              =  2.0;
    cfg.fine_nbins        = 400;               // 0.01 / bin seed
    cfg.max_bin_width     = max_bin_width;
    cfg.N_min_pt2_spectra = N_min_pt2;
    cfg.N_min_pt2_ratio   = N_min_pt2_ratio;
    cfg.base_filename     = "y_cms";
    cfg.canvas_title      = "y_{CMS}";
    cfg.nt_name           = nt_name;
    cfg.log_y_spectra     = false;             // rapidity is bell-shaped → linear Y
    AdaptiveSigPlot::run(cfg);
}

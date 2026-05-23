// oa_adaptive_sig_exp.C — adaptive-binning CB-subtracted e+e- opening angle.
// Thin wrapper over AdaptiveSigPlot::run; see AdaptiveSigPlot.h for details.
//
// Branch: dilepton_nt::oa   (lab-frame opening angle in degrees)
//
// Usage:
//   root -l -b -q plots/oa_adaptive_sig_exp.C
//   root -l -b -q 'plots/oa_adaptive_sig_exp.C(500, 30)'

#include "AdaptiveSigPlot.h"

void oa_adaptive_sig_exp(double      N_min_pt2         = 200.0,
                         double      max_bin_width     = 20.0,    // deg
                         double      N_min_pt2_ratio   = 2000.0,
                         const char* nt_name           = "dilepton_nt") {
    AdaptiveSigPlot::Config cfg;
    cfg.var_expr          = "oa";
    cfg.x_axis_label      = "#theta_{open}  [deg]";
    cfg.y_axis_unit       = "deg";
    cfg.xmin              =   0.0;
    cfg.xmax              = 180.0;
    cfg.fine_nbins        = 1800;              // 0.1 deg/bin seed
    cfg.max_bin_width     = max_bin_width;
    cfg.N_min_pt2_spectra = N_min_pt2;
    cfg.N_min_pt2_ratio   = N_min_pt2_ratio;
    cfg.base_filename     = "oa";
    cfg.canvas_title      = "Opening angle";
    cfg.nt_name           = nt_name;
    cfg.log_y_spectra     = true;              // OA distribution falls steeply
    AdaptiveSigPlot::run(cfg);
}

// mass_ee_adaptive_sig_exp.C — adaptive-binning CB-subtracted M_ee.
// Thin wrapper over AdaptiveSigPlot::run — see AdaptiveSigPlot.h for
// algorithm, layout, and output format details.
//
// Outputs (in plots/output/):
//   mass_ee_adaptive_sig.{pdf,png}                 spectra + fine ratio (sig)
//   mass_ee_adaptive_sig_ratios_coarse.{pdf,png}   3-panel coarse sig ratios
//
// Usage:
//   root -l -b -q plots/mass_ee_adaptive_sig_exp.C
//   root -l -b -q 'plots/mass_ee_adaptive_sig_exp.C(500)'              // wider spectra bins
//   root -l -b -q 'plots/mass_ee_adaptive_sig_exp.C(200, 0.2, 5000)'   // coarser ratio bins

#include "AdaptiveSigPlot.h"

void mass_ee_adaptive_sig_exp(double      N_min_pt2         = 200.0,
                              double      max_bin_width     = 0.2,
                              double      N_min_pt2_ratio   = 2000.0,
                              const char* nt_name           = "dilepton_nt") {
    AdaptiveSigPlot::Config cfg;
    cfg.var_expr          = "m_ee";
    cfg.x_axis_label      = "M_{e^{+}e^{-}} [GeV/c^{2}]";
    cfg.y_axis_unit       = "GeV/c^{2}";
    cfg.xmin              = 0.0;
    cfg.xmax              = 1.4;
    cfg.fine_nbins        = 1400;             // 1 MeV/bin seed
    cfg.max_bin_width     = max_bin_width;
    cfg.N_min_pt2_spectra = N_min_pt2;
    cfg.N_min_pt2_ratio   = N_min_pt2_ratio;
    cfg.base_filename     = "mass_ee";
    cfg.canvas_title      = "M_{e^{+}e^{-}}";
    cfg.nt_name           = nt_name;
    cfg.log_y_spectra     = true;
    AdaptiveSigPlot::run(cfg);
}

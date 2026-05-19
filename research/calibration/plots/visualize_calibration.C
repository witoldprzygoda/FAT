// research/calibration/plots/visualize_calibration.C
//
// Phase (c) of the trigger-bias calibration pipeline: visualisation only.
//
// Reads BOTH artefacts produced by phases (a) and (b):
//   ../../output_<channel>_cal.root       (trigger_cal_nt + trigger_cal_files)
//   ../../pt3_calibration_<channel>.root  (pt3_calibration TTree, in repo root)
//
// and produces per-channel and overlay diagnostic plots in
// research/calibration/plots/output/.
//
// The lookup path for the calibration ROOT is intentionally THE SAME path
// the main analysis (./ana on pp45_epem) will use — if the visualisation
// agrees with what you expect, main.cc will see exactly the same weights.
//
// Per-channel canvas, two pads:
//
//   Pad 1  — cumulative N_PT3 (red) and N_PT2 (blue) vs chain event index,
//            log-y. Shows how trigger counts accumulate across the run;
//            kinks visible here often correspond to segment boundaries.
//
//   Pad 2  — weight w plot:
//              * cumulative w     (red line, growing from the chain start)
//              * segment w        (blue thick horizontal lines spanning
//                                  each segment's chain range; each line
//                                  has a thin BLACK vertical error bar at
//                                  its midpoint showing ±σ_seg, the
//                                  Poisson-propagated uncertainty of the
//                                  calibration value)
//              * file boundaries  (light gray dotted verticals)
//
// Cross-channel overlay canvas:
//   * Segment-only piecewise-constant w for epem, epep, emem on one axis,
//     to spot whether the three channels share segment boundaries / drift
//     coherently.
//
// Usage (from research/calibration/):
//   root -l -b -q plots/visualize_calibration.C
//
// Per-channel pad-2 y-ranges and pad-1 fixed [1e2, 1e8] log range are
// hardcoded in visualize_calibration() at the bottom of this file —
// edit there to tune.
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TColor.h>
#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TAxis.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>

namespace {

constexpr double kK = 63.0;

// Pastel band colours, registered once. Done as explicit RGB rather than
// named ROOT colour indices because (a) named "kAzure-9 / kOrange-9" can
// render unexpectedly saturated in some PDF viewers, and (b) PDF/ROOT
// alpha rendering is unreliable across viewer combinations. Using a flat
// pale fill avoids both pitfalls — the bands are always faint by
// construction, no transparency required.
int paleColor(float r, float g, float b) {
    static std::map<int, int> cache;
    const int key = (int(r*255)) | (int(g*255)<<8) | (int(b*255)<<16);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    const int idx = TColor::GetFreeColorIndex();
    new TColor(idx, r, g, b);
    cache[key] = idx;
    return idx;
}

// One file's metadata + its cumulative offset in the chain.
struct FileInfo {
    std::string path;
    Long64_t    n_events;        // total events in the file's tree
    Long64_t    chain_offset;    // sum of n_events of all preceding files
};

// One per-segment per-file row (matches calibration TTree schema).
struct CalRow {
    int          seg_idx;
    std::string  file_path;
    Long64_t     event_lo;
    Long64_t     event_hi;
    double       w;
    double       sigma;
};

// Aggregated by seg_idx — for overlay drawing of segment lines.
struct SegmentSpan {
    int       seg_idx;
    Long64_t  global_lo;     // chain event coord of segment start
    Long64_t  global_hi;     // chain event coord of segment end (inclusive)
    double    w;
    double    sigma;
};

// ----- channel data loaded once and passed around ----------------------
struct ChannelData {
    std::string                label;
    std::string                scan_path;
    std::string                cal_path;
    std::vector<FileInfo>      files;          // ordered, with cum offsets
    std::map<std::string, int> path2idx;
    // Streaming accumulators built from trigger_events:
    std::vector<Long64_t>      cum_x;          // chain event idx checkpoints
    std::vector<long long>     cum_pt3;        // cumulative N_PT3 at checkpoints
    std::vector<long long>     cum_pt2;        // cumulative N_PT2 at checkpoints
    std::vector<CalRow>        cal_rows;       // calibration table rows
    std::vector<SegmentSpan>   seg_spans;      // segments mapped to chain coords
    Long64_t                   total_chain_events = 0;
};

bool loadFiles(TTree* t_files, ChannelData& cd) {
    std::string  fpath;
    std::string* pp = &fpath;
    Long64_t     fn = 0;
    t_files->SetBranchAddress("file_path",      &pp);
    t_files->SetBranchAddress("n_events_total", &fn);
    const Long64_t N = t_files->GetEntries();
    cd.files.reserve(N);
    Long64_t cum = 0;
    for (Long64_t i = 0; i < N; ++i) {
        t_files->GetEntry(i);
        cd.files.push_back({fpath, fn, cum});
        cd.path2idx[fpath] = (int) cd.files.size() - 1;
        cum += fn;
    }
    cd.total_chain_events = cum;
    return !cd.files.empty();
}

// Stream trigger_cal_nt with the same quality cuts trigger_calibration.C
// applies, accumulating cumulative N_PT3 and N_PT2 at periodic checkpoints
// for the cumulative-counts pad and the cumulative-w line. Cuts are mirrored
// here so the visualisation reflects what the calibration actually saw —
// not the raw 113.6M-entry stream.
//
// File-boundary detection mirrors trigger_calibration.C: local_event_idx
// monotonically non-decreases within a file and resets to 0 at each file
// change (including post-hadd chunk boundaries).
bool loadEvents(TTree* t_evts, ChannelData& cd) {
    Float_t f_file_idx = 0, f_local_idx = 0, f_trigbit = 0;
    Float_t f_oa = 0, f_isBest = 0, f_vz = 0, f_si = 0;
    t_evts->SetBranchAddress("file_idx",        &f_file_idx);
    t_evts->SetBranchAddress("local_event_idx", &f_local_idx);
    t_evts->SetBranchAddress("trigbit",         &f_trigbit);
    t_evts->SetBranchAddress("oa",              &f_oa);
    t_evts->SetBranchAddress("isBest",          &f_isBest);
    t_evts->SetBranchAddress("eVertReco_z",     &f_vz);
    t_evts->SetBranchAddress("start_iteration", &f_si);

    const Long64_t N = t_evts->GetEntries();
    if (N == 0) return false;

    long long cum_p3 = 0, cum_p2 = 0;
    int      current_file_idx = 0;
    Long64_t prev_local_idx   = -1;

    // Cumulative checkpoint cadence: keep the curve at ~few thousand
    // points regardless of input size.
    const Long64_t CHECK_EVERY = std::max<Long64_t>(1, N / 5000);

    for (Long64_t i = 0; i < N; ++i) {
        t_evts->GetEntry(i);

        const Long64_t li = static_cast<Long64_t>(f_local_idx);
        if (i > 0 && li < prev_local_idx) ++current_file_idx;
        prev_local_idx = li;

        // Same cuts as trigger_calibration.C — NO OA cut, segmentation
        // edges are OA-independent (see trigger_calibration.C header).
        if (f_isBest != 1.0f)  continue;
        if (f_vz   <= -500.0f) continue;
        if (f_si   != 3.0f)    continue;
        const int tb = static_cast<int>(f_trigbit);
        if (tb != 8192 && tb != 4096) continue;

        if (tb == 8192) ++cum_p3;
        else            ++cum_p2;

        if (i % CHECK_EVERY == 0 || i == N - 1) {
            if (current_file_idx >= 0 &&
                current_file_idx < (int) cd.files.size()) {
                const Long64_t global_idx =
                    cd.files[current_file_idx].chain_offset + li;
                cd.cum_x.push_back(global_idx);
                cd.cum_pt3.push_back(cum_p3);
                cd.cum_pt2.push_back(cum_p2);
            }
        }
    }
    return true;
}

bool loadCalibration(const std::string& path, ChannelData& cd) {
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "  WARNING: cannot open " << path
                  << " — segments will be empty\n";
        if (f) delete f;
        return false;
    }
    TTree* t = dynamic_cast<TTree*>(f->Get("pt3_calibration"));
    if (!t) {
        std::cerr << "  WARNING: pt3_calibration tree missing in " << path << "\n";
        f->Close(); delete f;
        return false;
    }
    Int_t        seg = 0;
    std::string  pth;
    std::string* ppth = &pth;
    Long64_t lo = 0, hi = 0;
    Double_t w = 0, sigma = 0;
    t->SetBranchAddress("seg_idx",   &seg);
    t->SetBranchAddress("file_path", &ppth);
    t->SetBranchAddress("event_lo",  &lo);
    t->SetBranchAddress("event_hi",  &hi);
    t->SetBranchAddress("w",         &w);
    t->SetBranchAddress("sigma",     &sigma);
    const Long64_t N = t->GetEntries();
    cd.cal_rows.reserve(N);
    for (Long64_t i = 0; i < N; ++i) {
        t->GetEntry(i);
        cd.cal_rows.push_back({seg, pth, lo, hi, w, sigma});
    }
    f->Close(); delete f;

    // Aggregate into segments: lo of segment = lowest global_lo across its
    // rows; hi = highest global_hi; w/sigma constant across rows of a seg.
    std::map<int, SegmentSpan> by_seg;
    for (const CalRow& r : cd.cal_rows) {
        auto fi = cd.path2idx.find(r.file_path);
        if (fi == cd.path2idx.end()) continue;
        const Long64_t off = cd.files[fi->second].chain_offset;
        const Long64_t glo = off + r.event_lo;
        const Long64_t ghi = off + r.event_hi;
        auto it = by_seg.find(r.seg_idx);
        if (it == by_seg.end()) {
            by_seg[r.seg_idx] = {r.seg_idx, glo, ghi, r.w, r.sigma};
        } else {
            it->second.global_lo = std::min(it->second.global_lo, glo);
            it->second.global_hi = std::max(it->second.global_hi, ghi);
        }
    }
    cd.seg_spans.reserve(by_seg.size());
    for (auto& p : by_seg) cd.seg_spans.push_back(p.second);
    return true;
}

bool loadChannel(const std::string& chan, ChannelData& cd) {
    cd.label     = chan;
    // Paths are relative to the cwd from which root is invoked. Standard
    // usage is `cd research/calibration && root -l -b -q plots/...` so the
    // scan ROOTs live next to the cwd (./) and the calibration ROOTs in
    // the repo root (../../) per the project convention.
    cd.scan_path = "../../output_" + chan + "_cal.root";
    cd.cal_path  = "../../pt3_calibration_" + chan + ".root";

    TFile* fscan = TFile::Open(cd.scan_path.c_str(), "READ");
    if (!fscan || fscan->IsZombie()) {
        std::cerr << "  WARNING: cannot open " << cd.scan_path << " — skipping\n";
        if (fscan) delete fscan;
        return false;
    }
    TTree* t_files = dynamic_cast<TTree*>(fscan->Get("trigger_cal_files"));
    TTree* t_evts  = dynamic_cast<TTree*>(fscan->Get("trigger_cal_nt"));
    if (!t_files || !t_evts) {
        std::cerr << "  WARNING: missing trees in " << cd.scan_path << "\n";
        fscan->Close(); delete fscan;
        return false;
    }
    if (!loadFiles(t_files, cd)) {
        fscan->Close(); delete fscan;
        return false;
    }
    if (!loadEvents(t_evts, cd)) {
        fscan->Close(); delete fscan;
        return false;
    }
    fscan->Close(); delete fscan;

    loadCalibration(cd.cal_path, cd);   // OK if missing — segments just absent

    std::cout << "  loaded " << chan
              << ": files=" << cd.files.size()
              << "  segments=" << cd.seg_spans.size() << "\n";
    return true;
}

// ----- drawing ---------------------------------------------------------

// Identify "outlier" segments: short-duration AND deviating significantly
// from the local baseline (mean of left + right neighbour weights).
// Returns one bool per segment in cd.seg_spans order.
//   min_jump_pct       — required relative deviation from neighbour
//                        baseline (default 0.05 = 5%)
//   max_duration_frac  — segment span must be SHORTER than this fraction
//                        of the channel's total chain to qualify
//                        (default 0.01 = 1% of chain)
// Both conditions ANDed: a segment is flagged only if it's both brief AND
// far from its neighbours. A long-duration segment that genuinely shifts
// the baseline is NOT an outlier — that's a real trend change.
std::vector<bool> tagOutliers(const ChannelData& cd,
                              double min_jump_pct,
                              double max_duration_frac) {
    const size_t N = cd.seg_spans.size();
    std::vector<bool> out(N, false);
    if (N < 2 || cd.total_chain_events <= 0) return out;
    for (size_t s = 0; s < N; ++s) {
        const SegmentSpan& seg = cd.seg_spans[s];
        const Long64_t span = seg.global_hi - seg.global_lo + 1;
        const double dur_frac = (double) span / (double) cd.total_chain_events;
        if (dur_frac >= max_duration_frac) continue;       // too long
        double baseline = 0; int nb = 0;
        if (s > 0)     { baseline += cd.seg_spans[s-1].w; ++nb; }
        if (s + 1 < N) { baseline += cd.seg_spans[s+1].w; ++nb; }
        if (nb == 0)   continue;
        baseline /= nb;
        if (baseline <= 0) continue;
        const double rel_jump = std::abs(seg.w - baseline) / baseline;
        if (rel_jump >= min_jump_pct) out[s] = true;
    }
    return out;
}

void drawChannelCanvas(const ChannelData& cd, double y_lo, double y_hi,
                       double min_jump_pct      = 0.05,
                       double max_duration_frac = 0.01) {
    if (cd.cum_x.empty()) return;

    // Outlier accounting — count segments whose w falls outside the fixed
    // y range. They will be drawn off-canvas (gaps in the visible step
    // function); user needs to know they exist.
    int  n_outliers = 0;
    double w_min =  std::numeric_limits<double>::infinity();
    double w_max = -std::numeric_limits<double>::infinity();
    for (const SegmentSpan& s : cd.seg_spans) {
        w_min = std::min(w_min, s.w);
        w_max = std::max(w_max, s.w);
        if (s.w < y_lo || s.w > y_hi) ++n_outliers;
    }
    std::cout << "  [" << cd.label << "]  w range observed: ["
              << std::fixed << std::setprecision(3) << w_min << ", "
              << w_max << "]   off-pad: " << n_outliers
              << " / " << cd.seg_spans.size()
              << " segments outside [" << y_lo << ", " << y_hi << "]\n";

    // -- Pad 1: cumulative N_PT3, N_PT2 vs chain event idx (log-y) --------
    std::vector<double> x(cd.cum_x.begin(), cd.cum_x.end());
    std::vector<double> y3(cd.cum_pt3.begin(), cd.cum_pt3.end());
    std::vector<double> y2(cd.cum_pt2.begin(), cd.cum_pt2.end());

    auto* g_p3 = new TGraph((int) x.size(), x.data(), y3.data());
    auto* g_p2 = new TGraph((int) x.size(), x.data(), y2.data());
    g_p3->SetLineColor(kRed+1);  g_p3->SetLineWidth(2);
    g_p2->SetLineColor(kBlue+1); g_p2->SetLineWidth(2);

    // -- Pad 2: cumulative w (red line), segments (blue with black error bars)
    // Cumulative w(x) computed from cum_pt3, cum_pt2 at checkpoints.
    std::vector<double> cwx, cwy;
    cwx.reserve(cd.cum_x.size()); cwy.reserve(cd.cum_x.size());
    for (size_t i = 0; i < cd.cum_x.size(); ++i) {
        if (cd.cum_pt3[i] <= 0 || cd.cum_pt2[i] <= 0) continue;
        cwx.push_back((double) cd.cum_x[i]);
        cwy.push_back(kK * (double) cd.cum_pt2[i] / (double) cd.cum_pt3[i]);
    }
    auto* g_cw = new TGraph((int) cwx.size(), cwx.data(), cwy.data());
    g_cw->SetLineColor(kRed+1); g_cw->SetLineWidth(2);

    const TString cname = TString::Format("c_calib_%s", cd.label.c_str());
    auto* c = new TCanvas(cname, cname, 1500, 950);
    c->Divide(1, 2, 0.001, 0.001);

    // Common chain x-extent reused across pad 1 / pad 2 / band drawing.
    const double xmax_chain = cd.files.empty()
        ? 1.0
        : (double) (cd.files.back().chain_offset + cd.files.back().n_events);

    // Shared band styling (alternating per file across the full pad height).
    // Very pale pastel RGB — same hue family as kAzure/kOrange but flat
    // light so they don't depend on alpha rendering at all.
    const Color_t kBandA = paleColor(0.92f, 0.95f, 0.99f);  // pale azure
    const Color_t kBandB = paleColor(0.99f, 0.95f, 0.92f);  // pale linen / orange
    const double  kLabelYNDC  = 0.83;     // inside data area, near top, below legend strip
    const double  kLabelSize  = 0.022;

    // -- Pad 1 ----
    auto* p1 = (TPad*) c->cd(1);
    p1->SetLogy(); p1->SetGridx(); p1->SetGridy();
    p1->SetLeftMargin(0.10); p1->SetRightMargin(0.04);
    p1->SetTopMargin(0.13);  p1->SetBottomMargin(0.13);

    // Pad 1 y-range: pinned to a fixed log decade [1e2, 1e8] across all
    // three channels so the cumulative-counts curves are directly visually
    // comparable. The dynamic-fit version made each channel's pad 1
    // self-similar but mutually incomparable.
    const double y_min_p1 = 1e2;
    const double y_max_p1 = 1e8;

    auto* frame1 = p1->DrawFrame(0.0, y_min_p1, xmax_chain, y_max_p1);
    frame1->SetTitle(TString::Format(
        "cumulative trigger counts  (%s);chain event index;count",
        cd.label.c_str()));

    {
        const double lm = p1->GetLeftMargin(), rm = p1->GetRightMargin();
        for (size_t i = 0; i < cd.files.size(); ++i) {
            const double xb_lo = (double) cd.files[i].chain_offset;
            const double xb_hi = xb_lo + (double) cd.files[i].n_events;
            auto* box = new TBox(xb_lo, y_min_p1, xb_hi, y_max_p1);
            box->SetFillColor((i % 2 == 0) ? kBandA : kBandB);
            box->SetLineWidth(0);
            box->Draw();
            const double ndc_x = lm + (0.5 * (xb_lo + xb_hi) / xmax_chain)
                                       * (1.0 - lm - rm);
            auto* lbl = new TLatex();
            lbl->SetNDC();
            lbl->SetTextSize(kLabelSize);
            lbl->SetTextAlign(22);
            lbl->SetTextColor(kGray + 3);
            lbl->DrawLatex(ndc_x, kLabelYNDC, TString::Format("%zu", i));
        }
    }
    // Force tick marks back on top of the bands so they remain visible
    // even when alpha rendering is unreliable in the user's PDF viewer.
    p1->RedrawAxis();

    g_p3->Draw("L");
    g_p2->Draw("L");

    auto* leg1 = new TLegend(0.55, 0.86, 0.96, 0.93);
    leg1->SetTextSize(0.026);
    leg1->SetBorderSize(0);
    leg1->SetFillColorAlpha(kWhite, 0.7);
    leg1->SetNColumns(2);
    leg1->AddEntry(g_p3, "cum N_{PT3}", "l");
    leg1->AddEntry(g_p2, "cum N_{PT2}", "l");
    leg1->Draw();

    // -- Pad 2 ----
    auto* p2 = (TPad*) c->cd(2);
    p2->SetGridx(); p2->SetGridy();
    p2->SetLeftMargin(0.10); p2->SetRightMargin(0.04);
    p2->SetTopMargin(0.13);  p2->SetBottomMargin(0.16);

    // Frame manually so axis stack is: axes < bands < cumulative-line/segments.
    auto* frame2 = p2->DrawFrame(0.0, y_lo, xmax_chain, y_hi);
    frame2->SetTitle(TString::Format(
        "weight w = 63 #upoint N_{PT2}/N_{PT3}  (%s);chain event index;w",
        cd.label.c_str()));

    // Alternating semi-transparent file bands spanning the full pad height,
    // with the file index labelled at the top of each band (NDC). Replaces
    // the old gray-dashed file-boundary verticals with something that
    // simultaneously identifies WHICH file you're looking at.
    {
        const double lm = p2->GetLeftMargin(), rm = p2->GetRightMargin();
        for (size_t i = 0; i < cd.files.size(); ++i) {
            const double xb_lo = (double) cd.files[i].chain_offset;
            const double xb_hi = xb_lo + (double) cd.files[i].n_events;
            auto* box = new TBox(xb_lo, y_lo, xb_hi, y_hi);
            box->SetFillColor((i % 2 == 0) ? kBandA : kBandB);
            box->SetLineWidth(0);
            box->Draw();
            const double ndc_x = lm + (0.5 * (xb_lo + xb_hi) / xmax_chain)
                                       * (1.0 - lm - rm);
            auto* lbl = new TLatex();
            lbl->SetNDC();
            lbl->SetTextSize(kLabelSize);
            lbl->SetTextAlign(22);
            lbl->SetTextColor(kGray + 3);
            lbl->DrawLatex(ndc_x, kLabelYNDC, TString::Format("%zu", i));
        }
    }
    // Force tick marks back on top of the bands.
    p2->RedrawAxis();

    // Cumulative w drawn ON TOP of bands.
    g_cw->Draw("L");

    // Each segment drawn as: a thick horizontal line (the value w_seg
    // spanning the segment's chain range) PLUS a thin BLACK vertical error
    // bar at the segment midpoint (±σ_seg). Colour is BLUE for "trend"
    // segments and RED for short-duration outliers (per tagOutliers
    // criteria — short segment far from local baseline). The red band
    // marks segments that the user may want to drop entirely from the
    // physics analysis via a downstream weight-acceptance cut in main.cc.
    const std::vector<bool> outlier = tagOutliers(cd, min_jump_pct,
                                                  max_duration_frac);
    int n_brief_outliers = std::count(outlier.begin(), outlier.end(), true);
    std::cout << "  [" << cd.label << "]  brief-outliers: "
              << n_brief_outliers << " / " << cd.seg_spans.size()
              << " (jump > " << std::fixed << std::setprecision(1)
              << 100*min_jump_pct << "% AND duration < "
              << 100*max_duration_frac << "% of chain)\n";
    for (size_t s = 0; s < cd.seg_spans.size(); ++s) {
        if (!outlier[s]) continue;
        const SegmentSpan& seg = cd.seg_spans[s];
        const Long64_t span = seg.global_hi - seg.global_lo + 1;
        std::cout << "    seg " << std::setw(3) << seg.seg_idx
                  << "  w=" << std::fixed << std::setprecision(3) << seg.w
                  << " +/-" << seg.sigma
                  << "  span=" << std::setw(8) << span
                  << " (" << std::setprecision(2)
                  << 100.0*span/cd.total_chain_events << "%)\n";
    }

    const Color_t kSegColor    = kBlue + 1;
    const Color_t kSegOutColor = kRed  + 1;
    TLine* seg_legend_normal = nullptr;
    TLine* seg_legend_outlier = nullptr;
    for (size_t s = 0; s < cd.seg_spans.size(); ++s) {
        const SegmentSpan& sp = cd.seg_spans[s];
        const double xlo  = (double) sp.global_lo;
        const double xhi  = (double) sp.global_hi;
        const double xmid = 0.5 * (xlo + xhi);
        const Color_t col = outlier[s] ? kSegOutColor : kSegColor;

        auto* lh = new TLine(xlo, sp.w, xhi, sp.w);
        lh->SetLineColor(col); lh->SetLineWidth(3);
        lh->Draw();

        if (sp.sigma > 0) {
            auto* lv = new TLine(xmid, sp.w - sp.sigma, xmid, sp.w + sp.sigma);
            lv->SetLineColor(kBlack); lv->SetLineWidth(1);
            lv->Draw();
        }

        if (outlier[s] && !seg_legend_outlier) seg_legend_outlier = lh;
        if (!outlier[s] && !seg_legend_normal)  seg_legend_normal  = lh;
    }

    // Compact legend in the freed top-margin strip (above the data area).
    auto* leg2 = new TLegend(0.45, 0.86, 0.96, 0.93);
    leg2->SetTextSize(0.026);
    leg2->SetBorderSize(0);
    leg2->SetFillColorAlpha(kWhite, 0.7);
    leg2->SetNColumns(3);
    leg2->AddEntry(g_cw, "cumulative w_{cum}", "l");
    if (seg_legend_normal) {
        leg2->AddEntry(seg_legend_normal,
                       TString::Format("trend (N=%zu)",
                                       cd.seg_spans.size() - n_brief_outliers),
                       "l");
    }
    if (seg_legend_outlier) {
        leg2->AddEntry(seg_legend_outlier,
                       TString::Format("outliers (N=%d)", n_brief_outliers),
                       "l");
    }
    leg2->Draw();

    gSystem->mkdir("plots/output", kTRUE);
    const TString stem = TString::Format(
        "plots/output/visualize_calibration_%s", cd.label.c_str());
    c->SaveAs(stem + ".pdf");
    c->SaveAs(stem + ".png");
    std::cout << "  saved " << stem << ".{pdf,png}\n";
}

void drawOverlay(const std::vector<ChannelData*>& chans,
                 double y_lo, double y_hi,
                 double min_jump_pct,
                 double max_duration_frac) {
    auto* c = new TCanvas("c_calib_overlay", "c_calib_overlay", 1500, 700);
    c->SetGridx(); c->SetGridy();
    // Wider right margin so the axis-overflow indicator (×10⁶ or similar)
    // for the chain event index doesn't get clipped at the right edge.
    c->SetLeftMargin(0.11); c->SetRightMargin(0.07);
    c->SetTopMargin(0.13);  c->SetBottomMargin(0.13);

    // X-axis = REFERENCE channel's chain event index. Reference is the
    // channel with the most events overall — typically epem at pp45.
    // Bands have widths = reference channel's per-file event counts (so
    // wider bands = files with more data in the reference channel = the
    // "important" days). Segments from non-reference channels are
    // re-projected per file: a segment covering events [a, b] in file F
    // of channel C (which has Nc events in F) maps onto the reference
    // axis as
    //   x_lo_ref = ref.chain_offset[F] + (a / Nc) * Nref
    //   x_hi_ref = ref.chain_offset[F] + (b / Nc) * Nref
    // i.e. the fractional position inside file F is preserved when
    // crossing channels (the file is the same physical .root, just
    // selected differently). Files unique to a non-reference channel
    // (none expected in our setup but defensive) are skipped.
    ChannelData* ref = nullptr;
    Long64_t     ref_total = -1;
    for (auto* cd : chans) {
        if (!cd || cd->files.empty()) continue;
        const Long64_t total = cd->files.back().chain_offset
                             + cd->files.back().n_events;
        if (total > ref_total) { ref_total = total; ref = cd; }
    }
    if (!ref) return;

    const double xmax_ref = (double) ref_total;

    auto* frame = c->DrawFrame(0.0, y_lo, xmax_ref, y_hi);
    frame->SetTitle(TString::Format(
        "PT3 trigger-bias segment weights, all channels overlaid;"
        "chain event index (%s reference);"
        "w = 63 #upoint N_{PT2}/N_{PT3}",
        ref->label.c_str()));
    // Smaller X-axis tick labels so the trailing "×10⁶" doesn't overflow
    // the right edge for our 100M-event chains.
    frame->GetXaxis()->SetLabelSize(0.030);

    // Alternating per-file bands across the full pad height — widths
    // proportional to the reference channel's per-file event counts.
    {
        const Color_t kBandA = paleColor(0.92f, 0.95f, 0.99f);
        const Color_t kBandB = paleColor(0.99f, 0.95f, 0.92f);
        const double  lm = c->GetLeftMargin(), rm = c->GetRightMargin();
        for (size_t i = 0; i < ref->files.size(); ++i) {
            const double xb_lo = (double) ref->files[i].chain_offset;
            const double xb_hi = xb_lo + (double) ref->files[i].n_events;
            auto* box = new TBox(xb_lo, y_lo, xb_hi, y_hi);
            box->SetFillColor((i % 2 == 0) ? kBandA : kBandB);
            box->SetLineWidth(0);
            box->Draw();
            const double ndc_x = lm + (0.5 * (xb_lo + xb_hi) / xmax_ref)
                                       * (1.0 - lm - rm);
            auto* lbl = new TLatex();
            lbl->SetNDC();
            lbl->SetTextSize(0.022);
            lbl->SetTextAlign(22);
            lbl->SetTextColor(kGray + 3);
            lbl->DrawLatex(ndc_x, 0.83, TString::Format("%zu", i));
        }
    }
    // Force tick marks back on top of the bands.
    c->RedrawAxis();

    auto* leg = new TLegend(0.45, 0.86, 0.96, 0.93);
    leg->SetTextSize(0.026);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(kWhite, 0.7);
    leg->SetNColumns(4);
    const int colors[3] = { kBlack, kBlue+1, kRed+1 };
    const Color_t kOutColor = kRed + 1;   // outliers are RED across all channels
    int ci = 0;
    TLine* legend_outlier = nullptr;
    int total_outliers = 0;
    for (auto* cd : chans) {
        if (!cd || cd->cal_rows.empty()) { ++ci; continue; }

        // Tag outliers in this channel using the same criteria as pad 2.
        // Build a set of outlier seg_idx values for fast lookup per row.
        const std::vector<bool> outlier =
            tagOutliers(*cd, min_jump_pct, max_duration_frac);
        std::set<int> outlier_seg_ids;
        for (size_t s = 0; s < cd->seg_spans.size(); ++s) {
            if (outlier[s]) outlier_seg_ids.insert(cd->seg_spans[s].seg_idx);
        }
        total_outliers += (int) outlier_seg_ids.size();

        TLine* legend_line = nullptr;
        for (const CalRow& r : cd->cal_rows) {
            // Look up SAME physical file in BOTH this channel and the
            // reference channel, then project event_lo/event_hi from this
            // channel's file-local indexing onto the reference's chain
            // event coordinates by preserving fractional in-file position.
            auto fi_cd  = cd->path2idx.find(r.file_path);
            auto fi_ref = ref->path2idx.find(r.file_path);
            if (fi_cd  == cd->path2idx.end())  continue;
            if (fi_ref == ref->path2idx.end()) continue;
            const Long64_t Nc = cd->files[fi_cd->second].n_events;
            const Long64_t Nr = ref->files[fi_ref->second].n_events;
            if (Nc <= 0 || Nr <= 0) continue;
            const double off_ref = (double) ref->files[fi_ref->second].chain_offset;
            const double x_lo = off_ref + (double) r.event_lo / (double) Nc * (double) Nr;
            const double x_hi = off_ref + (double) r.event_hi / (double) Nc * (double) Nr;
            const bool is_outlier = (outlier_seg_ids.count(r.seg_idx) > 0);
            const Color_t col = is_outlier ? kOutColor : (Color_t) colors[ci];

            auto* lc = new TLine(x_lo, r.w, x_hi, r.w);
            lc->SetLineColor(col); lc->SetLineWidth(2);
            lc->Draw();
            if (is_outlier && !legend_outlier) legend_outlier = lc;
            if (!is_outlier && !legend_line)   legend_line    = lc;
        }
        if (legend_line) leg->AddEntry(legend_line, cd->label.c_str(), "l");
        ++ci;
    }
    if (legend_outlier) {
        leg->AddEntry(legend_outlier,
                      TString::Format("outliers (N=%d)", total_outliers), "l");
    }
    leg->Draw();

    gSystem->mkdir("plots/output", kTRUE);
    c->SaveAs("plots/output/visualize_calibration_overlay.pdf");
    c->SaveAs("plots/output/visualize_calibration_overlay.png");
    std::cout << "  saved plots/output/visualize_calibration_overlay.{pdf,png}\n";
}

}  // anonymous namespace

// Per-channel y-ranges for the weight pad. Hardcoded because the
// observed dynamic range is genuinely different between channels:
//   epem — chunkier sample, real run-to-run drift in [1.4, 2.6]
//   epep, emem — flat near unity (range ~[1.0, 1.5])
// Same scale across {epep, emem} so they're mutually comparable; epem
// gets its own band to keep its detail visible without compressing the
// other two. Overlay uses a span that covers all three.
//
// Outlier-tagging parameters (defaults: 6% jump, 0.5% chain duration —
// validated by the user against the epem PDF as catching the brief
// excursions cleanly without polluting the trend):
//   min_jump_pct       — required relative deviation from neighbour
//                        baseline to flag a segment as outlier
//   max_duration_frac  — segment must be SHORTER than this fraction of
//                        chain to qualify; long deviations are real
//                        trend changes, not outliers
// Outlier segments are drawn RED on pad 2 (and on the cross-channel
// overlay) instead of the channel's normal colour, so the user can
// immediately spot brief excursions across all channels at a glance.
// These red-tagged segments are exactly the candidates for a downstream
// weight-acceptance cut in main.cc on pp45_epem.
void visualize_calibration(double min_jump_pct      = 0.06,
                           double max_duration_frac = 0.005) {
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    std::cout << "=== visualize_calibration ===\n";

    ChannelData epem, epep, emem;
    const bool ok_e = loadChannel("epem", epem);
    const bool ok_p = loadChannel("epep", epep);
    const bool ok_m = loadChannel("emem", emem);

    if (!ok_e && !ok_p && !ok_m) {
        std::cerr << "ERROR: no input files found. Make sure "
                  << "output_<channel>_cal.root and "
                  << "pt3_calibration_<channel>.root exist in repo root.\n";
        return;
    }

    if (ok_e) drawChannelCanvas(epem, 1.0, 4.0, min_jump_pct, max_duration_frac);
    if (ok_p) drawChannelCanvas(epep, 1.0, 2.0, min_jump_pct, max_duration_frac);
    if (ok_m) drawChannelCanvas(emem, 1.0, 2.0, min_jump_pct, max_duration_frac);

    std::vector<ChannelData*> chans;
    if (ok_e) chans.push_back(&epem);
    if (ok_p) chans.push_back(&epep);
    if (ok_m) chans.push_back(&emem);
    if (chans.size() >= 2)
        drawOverlay(chans, 1.0, 4.0, min_jump_pct, max_duration_frac);

    std::cout << "Done.\n";
}

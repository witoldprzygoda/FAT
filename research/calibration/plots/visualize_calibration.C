// research/calibration/plots/visualize_calibration.C
//
// Phase (c) of the trigger-bias calibration pipeline: visualisation only.
//
// Reads BOTH artefacts produced by phases (a) and (b):
//   trigger_scan_<channel>.root           (trigger_events + files TTrees)
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
//              * per-window w   (black markers; each window contains a
//                                fixed number of PT2 events — controlled
//                                by `window_pt2` argument, default 1000;
//                                Poisson-propagated σ_w bars)
//              * cumulative w   (red line, growing from the chain start)
//              * segment w      (blue piecewise-constant step function,
//                                ±1σ thin dashed band per segment)
//              * file boundaries (light gray dashed verticals)
//              * segment boundaries (blue dashed verticals)
//
// Cross-channel overlay canvas:
//   * Segment-only piecewise-constant w for epem, epep, emem on one axis,
//     to spot whether the three channels share segment boundaries / drift
//     coherently.
//
// Usage (from research/calibration/):
//   root -l -b -q plots/visualize_calibration.C
//   root -l -b -q 'plots/visualize_calibration.C(500)'    # tighter (~4.5% σ)
//   root -l -b -q 'plots/visualize_calibration.C(2000)'   # smoother (~2.2% σ)
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
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TAxis.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

namespace {

constexpr double kK = 63.0;

// One file's metadata + its cumulative offset in the chain.
struct FileInfo {
    std::string path;
    Long64_t    n_events;        // total events in the file's tree
    Long64_t    chain_offset;    // sum of n_events of all preceding files
};

// One per-window aggregation for the per-window w plot.
struct Window {
    Long64_t  global_event_mid;   // chain event index of window centre
    long long n_pt3;
    long long n_pt2;
    double w()     const { return (n_pt3 > 0 && n_pt2 > 0)
                                      ? kK * (double) n_pt2 / (double) n_pt3 : 0.0; }
    double sigma() const {
        if (n_pt3 <= 0 || n_pt2 <= 0) return 0.0;
        const double a = (double) n_pt3, b = (double) n_pt2;
        return w() * std::sqrt(1.0/a + 1.0/b);
    }
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
    std::vector<Window>        windows;        // per-window aggregates
    std::vector<CalRow>        cal_rows;       // calibration table rows
    std::vector<SegmentSpan>   seg_spans;      // segments mapped to chain coords
    Long64_t                   total_chain_events = 0;
};

bool loadFiles(TTree* t_files, ChannelData& cd) {
    std::string  fpath;
    std::string* pp = &fpath;
    Long64_t     fn = 0, fp3 = 0, fp2 = 0;
    t_files->SetBranchAddress("file_path",      &pp);
    t_files->SetBranchAddress("n_events_total", &fn);
    t_files->SetBranchAddress("n_pt3",          &fp3);
    t_files->SetBranchAddress("n_pt2",          &fp2);
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

bool loadEventsAndWindows(TTree* t_evts, ChannelData& cd, Long64_t window_pt2) {
    std::string  epath;
    std::string* pp = &epath;
    Long64_t e_local = 0;
    Int_t    e_trig  = 0;
    t_evts->SetBranchAddress("file_path",       &pp);
    t_evts->SetBranchAddress("local_event_idx", &e_local);
    t_evts->SetBranchAddress("trigbit",         &e_trig);

    const Long64_t N = t_evts->GetEntries();
    if (N == 0) return false;

    long long cum_p3 = 0, cum_p2 = 0;
    long long win_p3 = 0, win_p2 = 0;
    Long64_t  win_first_global = -1;

    // Cumulative checkpoint cadence: every CHECK_EVERY trigger events we
    // record (chain_idx, cum_pt3, cum_pt2) — keeps the cumulative graph
    // tractable to draw (~few thousand points) for ~5M trigger events.
    const Long64_t CHECK_EVERY = std::max<Long64_t>(1, N / 5000);

    for (Long64_t i = 0; i < N; ++i) {
        t_evts->GetEntry(i);
        auto it = cd.path2idx.find(epath);
        if (it == cd.path2idx.end()) continue;
        const Long64_t global_idx = cd.files[it->second].chain_offset + e_local;

        if (e_trig == 8192) ++cum_p3;
        else if (e_trig == 4096) ++cum_p2;

        if (i % CHECK_EVERY == 0 || i == N - 1) {
            cd.cum_x.push_back(global_idx);
            cd.cum_pt3.push_back(cum_p3);
            cd.cum_pt2.push_back(cum_p2);
        }

        // Window aggregation: close when we have collected `window_pt2`
        // PT2 events (the rare class — controls the statistics floor).
        // PT3 in the same window is whatever it happens to be (~64×PT2
        // on average given the PT2 downscale).
        if (win_first_global < 0) win_first_global = global_idx;
        if (e_trig == 8192) ++win_p3;
        else if (e_trig == 4096) ++win_p2;
        const bool window_done = (win_p2 >= window_pt2) || (i == N - 1);
        if (window_done && (win_p3 > 0 || win_p2 > 0)) {
            Window w;
            w.global_event_mid = (win_first_global + global_idx) / 2;
            w.n_pt3            = win_p3;
            w.n_pt2            = win_p2;
            cd.windows.push_back(w);
            win_p3 = win_p2 = 0;
            win_first_global = -1;
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

bool loadChannel(const std::string& chan, ChannelData& cd, Long64_t window_pt2) {
    cd.label     = chan;
    // Paths are relative to the cwd from which root is invoked. Standard
    // usage is `cd research/calibration && root -l -b -q plots/...` so the
    // scan ROOTs live next to the cwd (./) and the calibration ROOTs in
    // the repo root (../../) per the project convention.
    cd.scan_path = "trigger_scan_"   + chan + ".root";
    cd.cal_path  = "../../pt3_calibration_" + chan + ".root";

    TFile* fscan = TFile::Open(cd.scan_path.c_str(), "READ");
    if (!fscan || fscan->IsZombie()) {
        std::cerr << "  WARNING: cannot open " << cd.scan_path << " — skipping\n";
        if (fscan) delete fscan;
        return false;
    }
    TTree* t_files = dynamic_cast<TTree*>(fscan->Get("files"));
    TTree* t_evts  = dynamic_cast<TTree*>(fscan->Get("trigger_events"));
    if (!t_files || !t_evts) {
        std::cerr << "  WARNING: missing trees in " << cd.scan_path << "\n";
        fscan->Close(); delete fscan;
        return false;
    }
    if (!loadFiles(t_files, cd)) {
        fscan->Close(); delete fscan;
        return false;
    }
    if (!loadEventsAndWindows(t_evts, cd, window_pt2)) {
        fscan->Close(); delete fscan;
        return false;
    }
    fscan->Close(); delete fscan;

    loadCalibration(cd.cal_path, cd);   // OK if missing — segments just absent

    std::cout << "  loaded " << chan
              << ": files=" << cd.files.size()
              << "  trig_events=" << (cd.windows.empty() ? 0 : cd.windows.back().global_event_mid + 1)
              << "  windows=" << cd.windows.size()
              << "  segments=" << cd.seg_spans.size() << "\n";
    return true;
}

// ----- drawing ---------------------------------------------------------

void drawChannelCanvas(const ChannelData& cd) {
    if (cd.cum_x.empty()) return;

    // -- Pad 1: cumulative N_PT3, N_PT2 vs chain event idx (log-y) --------
    std::vector<double> x(cd.cum_x.begin(), cd.cum_x.end());
    std::vector<double> y3(cd.cum_pt3.begin(), cd.cum_pt3.end());
    std::vector<double> y2(cd.cum_pt2.begin(), cd.cum_pt2.end());

    auto* g_p3 = new TGraph((int) x.size(), x.data(), y3.data());
    auto* g_p2 = new TGraph((int) x.size(), x.data(), y2.data());
    g_p3->SetLineColor(kRed+1);  g_p3->SetLineWidth(2);
    g_p2->SetLineColor(kBlue+1); g_p2->SetLineWidth(2);

    // -- Pad 2: per-window w (markers), cumulative w (line), segments -----
    std::vector<double> wx, wy, wex, wey;
    wx.reserve(cd.windows.size());
    for (const Window& w : cd.windows) {
        if (w.n_pt3 == 0 || w.n_pt2 == 0) continue;
        wx.push_back((double) w.global_event_mid);
        wy.push_back(w.w());
        wex.push_back(0.0);
        wey.push_back(w.sigma());
    }
    auto* g_w = new TGraphErrors((int) wx.size(),
                                 wx.data(), wy.data(),
                                 wex.data(), wey.data());
    g_w->SetMarkerStyle(20); g_w->SetMarkerColor(kBlack);
    g_w->SetLineColor(kBlack); g_w->SetMarkerSize(0.5);

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

    // -- Pad 1 ----
    auto* p1 = (TPad*) c->cd(1);
    p1->SetLogy(); p1->SetGridx(); p1->SetGridy();
    p1->SetLeftMargin(0.10); p1->SetRightMargin(0.04);
    p1->SetTopMargin(0.08);  p1->SetBottomMargin(0.13);

    auto* mg1 = new TMultiGraph();
    mg1->SetTitle(TString::Format(
        "cumulative trigger counts  (%s);chain event index;count",
        cd.label.c_str()));
    mg1->Add(g_p3, "L");
    mg1->Add(g_p2, "L");
    mg1->Draw("A");
    auto* leg1 = new TLegend(0.78, 0.78, 0.95, 0.93);
    leg1->AddEntry(g_p3, "cum N_{PT3}", "l");
    leg1->AddEntry(g_p2, "cum N_{PT2}", "l");
    leg1->Draw();

    // File boundaries on pad 1 (light gray dashed verticals).
    for (size_t i = 1; i < cd.files.size(); ++i) {
        const double xb = (double) cd.files[i].chain_offset;
        auto* lv = new TLine(xb,
                             p1->GetUymin(), xb, p1->GetUymax());
        lv->SetLineColor(kGray + 1);
        lv->SetLineStyle(3);
        lv->Draw();
    }

    // -- Pad 2 ----
    auto* p2 = (TPad*) c->cd(2);
    p2->SetGridx(); p2->SetGridy();
    p2->SetLeftMargin(0.10); p2->SetRightMargin(0.04);
    p2->SetTopMargin(0.08);  p2->SetBottomMargin(0.16);

    auto* mg2 = new TMultiGraph();
    mg2->SetTitle(TString::Format(
        "weight w = 63 #upoint N_{PT2}/N_{PT3}  (%s);chain event index;w",
        cd.label.c_str()));
    mg2->Add(g_w,  "P");
    mg2->Add(g_cw, "L");
    mg2->Draw("A");

    auto* leg2 = new TLegend(0.74, 0.78, 0.95, 0.95);
    leg2->AddEntry(g_w,  "per window  w (N_{PT2}=K)  #pm #sigma", "lp");
    leg2->AddEntry(g_cw, "cumulative w_{cum}",         "l");

    // Segments overlaid as blue horizontal lines per segment (with ±σ
    // band as thin dotted siblings). Blue chosen because the original
    // green lay on top of green vertical boundaries and was unreadable.
    const Color_t kSegColor = kBlue + 1;
    TLine* seg_legend_line = nullptr;
    for (const SegmentSpan& s : cd.seg_spans) {
        const double xlo = (double) s.global_lo;
        const double xhi = (double) s.global_hi;
        if (s.sigma > 0) {
            for (int sgn : {-1, 1}) {
                auto* lb = new TLine(xlo, s.w + sgn * s.sigma,
                                     xhi, s.w + sgn * s.sigma);
                lb->SetLineColor(kSegColor);
                lb->SetLineStyle(3); lb->SetLineWidth(1);
                lb->Draw();
            }
        }
        auto* lc = new TLine(xlo, s.w, xhi, s.w);
        lc->SetLineColor(kSegColor); lc->SetLineWidth(3);
        lc->Draw();
        if (!seg_legend_line) seg_legend_line = lc;
    }
    if (seg_legend_line) {
        leg2->AddEntry(seg_legend_line,
                       TString::Format("segments (N=%zu)", cd.seg_spans.size()),
                       "l");
    }
    leg2->Draw();

    // File boundaries (light gray) and segment boundaries (blue) on pad 2.
    p2->Update();
    const double yymin = p2->GetUymin();
    const double yymax = p2->GetUymax();
    for (size_t i = 1; i < cd.files.size(); ++i) {
        const double xb = (double) cd.files[i].chain_offset;
        auto* lv = new TLine(xb, yymin, xb, yymax);
        lv->SetLineColor(kGray + 1);
        lv->SetLineStyle(3);
        lv->Draw();
    }
    for (size_t i = 1; i < cd.seg_spans.size(); ++i) {
        const double xb = (double) cd.seg_spans[i].global_lo;
        auto* lv = new TLine(xb, yymin, xb, yymax);
        lv->SetLineColor(kSegColor);
        lv->SetLineStyle(2);
        lv->SetLineWidth(2);
        lv->Draw();
    }

    gSystem->mkdir("plots/output", kTRUE);
    const TString stem = TString::Format(
        "plots/output/visualize_calibration_%s", cd.label.c_str());
    c->SaveAs(stem + ".pdf");
    c->SaveAs(stem + ".png");
    std::cout << "  saved " << stem << ".{pdf,png}\n";
}

void drawOverlay(const std::vector<ChannelData*>& chans) {
    auto* c = new TCanvas("c_calib_overlay", "c_calib_overlay", 1500, 700);
    c->SetGridx(); c->SetGridy();
    c->SetLeftMargin(0.08); c->SetRightMargin(0.04);
    c->SetTopMargin(0.08);  c->SetBottomMargin(0.13);

    // Build a dummy frame so segments draw within proper axes.
    double xmin = +1e18, xmax = -1e18, ymin = +1e18, ymax = -1e18;
    for (auto* cd : chans) {
        for (const SegmentSpan& s : cd->seg_spans) {
            xmin = std::min(xmin, (double) s.global_lo);
            xmax = std::max(xmax, (double) s.global_hi);
            ymin = std::min(ymin, s.w - s.sigma);
            ymax = std::max(ymax, s.w + s.sigma);
        }
    }
    if (xmin >= xmax) { xmin = 0; xmax = 1; }
    const double pad = (ymax - ymin) * 0.1 + 1e-3;
    auto* frame = c->DrawFrame(xmin, ymin - pad, xmax, ymax + pad);
    frame->SetTitle("PT3 trigger-bias segment weights, all channels overlaid;"
                    "chain event index;w = 63 #upoint N_{PT2}/N_{PT3}");

    auto* leg = new TLegend(0.83, 0.78, 0.97, 0.95);
    const int colors[3]   = { kBlack, kBlue+1, kRed+1 };
    int ci = 0;
    for (auto* cd : chans) {
        if (!cd || cd->seg_spans.empty()) { ++ci; continue; }
        TLine* legend_line = nullptr;
        for (const SegmentSpan& s : cd->seg_spans) {
            auto* lc = new TLine((double) s.global_lo, s.w,
                                 (double) s.global_hi, s.w);
            lc->SetLineColor(colors[ci]); lc->SetLineWidth(2);
            lc->Draw();
            if (!legend_line) legend_line = lc;
        }
        if (legend_line) leg->AddEntry(legend_line, cd->label.c_str(), "l");
        ++ci;
    }
    leg->Draw();

    gSystem->mkdir("plots/output", kTRUE);
    c->SaveAs("plots/output/visualize_calibration_overlay.pdf");
    c->SaveAs("plots/output/visualize_calibration_overlay.png");
    std::cout << "  saved plots/output/visualize_calibration_overlay.{pdf,png}\n";
}

}  // anonymous namespace

// `window_pt2` is the number of PT2 events per visualisation window (PT2
// is the rare class — ~1.5% of trigger events at pp45 leptons due to the
// 64× downscale, so it controls the statistics floor). Each window
// therefore contains exactly `window_pt2` PT2 events plus ~64·window_pt2
// PT3 events. Per-window σ_w/w ≈ √(1/N_PT2 + 1/N_PT3) ≈ 1/√N_PT2.
//   window_pt2 = 1000  →  ~3% error per marker
//   window_pt2 =  500  →  ~4.5% error per marker
//   window_pt2 =  200  →  ~7% error per marker
// Pure visualisation parameter — does NOT affect the segmenter (which is
// PT2/PT3 stat-driven, not window-driven).
void visualize_calibration(Long64_t window_pt2 = 1000) {
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    std::cout << "=== visualize_calibration  window=" << window_pt2
              << " PT2 events per marker ===\n";

    ChannelData epem, epep, emem;
    const bool ok_e = loadChannel("epem", epem, window_pt2);
    const bool ok_p = loadChannel("epep", epep, window_pt2);
    const bool ok_m = loadChannel("emem", emem, window_pt2);

    if (!ok_e && !ok_p && !ok_m) {
        std::cerr << "ERROR: no scan files found. Run "
                  << "./run_parallel_scan.sh and trigger_calibration.C first.\n";
        return;
    }

    if (ok_e) drawChannelCanvas(epem);
    if (ok_p) drawChannelCanvas(epep);
    if (ok_m) drawChannelCanvas(emem);

    std::vector<ChannelData*> chans;
    if (ok_e) chans.push_back(&epem);
    if (ok_p) chans.push_back(&epep);
    if (ok_m) chans.push_back(&emem);
    if (chans.size() >= 2) drawOverlay(chans);

    std::cout << "Done.\n";
}

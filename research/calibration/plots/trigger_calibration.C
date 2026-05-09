// research/calibration/plots/trigger_calibration.C
//
// Reads per-file PT3 / PT2 trigger counts from trigger_scan_<channel>.root
// (produced by ../trigger_scan or ../run_parallel_scan.sh), runs an online
// z-test segmenter, visualises the run-by-run drift, and writes the final
// calibration tables as JSON in the REPO ROOT — that is where the main FAT
// analysis (pp45_epem branch) reads them from at Pass 2.
//
// Per channel (epem, epep, emem) one canvas with three stacked pads:
//
//   Pad 1  — N_PT3(red) and N_PT2(blue) vs file index (log-y)
//            ±√N error bars; spotting outright drops in beam/trigger rate.
//
//   Pad 2  — per-file weight  w[f]   = 63·N_PT2[f] / N_PT3[f]   (markers)
//            running cumulative w_c[f] = 63·ΣN_PT2 / ΣN_PT3      (red, ±1σ band)
//            segmenter pooled w_seg per segment                   (green
//            horizontal lines with ±1σ thin band).
//
//   Pad 3  — segmenter's z-score (vs the segment open at file f, before f
//            is added):
//                z[f] = (w[f] - w_seg) / √( σ_w[f]² + σ_seg² )
//            ±z_threshold lines (red dashed) and vertical green dashed
//            lines at every segment boundary.
//
// Segmenter algorithm — single forward pass, greedy z-test against the
// currently-open segment. For each file f compute z vs (cur_pt3, cur_pt2);
// if |z| > z_thr close the current segment at f-1 and start a new one at f;
// else add f to the current segment. No explicit min-segment-size — large
// σ_seg on small pools naturally suppresses early splits. Tighten z_thr to
// merge more aggressively, loosen it to react to local fluctuations.
//
// Output:
//   research/calibration/plots/output/trigger_calibration_<channel>.{pdf,png}
//   research/calibration/plots/output/trigger_calibration_overlay.{pdf,png}
//   <repo-root>/pt3_segments_<channel>.json    (LIVES IN REPO ROOT)
//
// Usage (from research/calibration/):
//   root -l -b -q plots/trigger_calibration.C                # default z_thr=3
//   root -l -b -q 'plots/trigger_calibration.C(2.5)'         # softer threshold
//   root -l -b -q 'plots/trigger_calibration.C(5.0)'         # stricter merging
//
// Inputs are searched relative to the research/calibration/ directory:
//   ../trigger_scan_epem.root
//   ../trigger_scan_epep.root
//   ../trigger_scan_emem.root
// Missing files are skipped with a warning (so you can run after just one
// channel pass too). The JSON output is suppressed for any channel whose
// scan ROOT is absent.
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
#include <TLatex.h>
#include <TString.h>
#include <TAxis.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

namespace {

// One channel's per-file payload.
struct ChannelData {
    std::string label;                       // "epem" / "epep" / "emem"
    std::string source_path;                 // path to the trigger_scan_*.root
    std::vector<int>          idx;
    std::vector<long long>    n_pt3;
    std::vector<long long>    n_pt2;
    std::vector<std::string>  file_path;     // populated from TTree branch
};

bool loadChannel(const std::string& path, ChannelData& cd) {
    cd.source_path = path;
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "  WARNING: cannot open " << path << " — skipping\n";
        if (f) delete f;
        return false;
    }
    TTree* t = dynamic_cast<TTree*>(f->Get("pt3_perfile"));
    if (!t) {
        std::cerr << "  WARNING: 'pt3_perfile' tree missing in " << path << "\n";
        f->Close(); delete f;
        return false;
    }
    Int_t        b_idx = 0;
    Long64_t     b_pt3 = 0, b_pt2 = 0;
    std::string  b_path;
    std::string* p_path = &b_path;
    t->SetBranchAddress("file_idx",  &b_idx);
    t->SetBranchAddress("n_pt3",     &b_pt3);
    t->SetBranchAddress("n_pt2",     &b_pt2);
    t->SetBranchAddress("file_path", &p_path);
    const Long64_t n = t->GetEntries();
    cd.idx.reserve(n);
    cd.n_pt3.reserve(n);
    cd.n_pt2.reserve(n);
    cd.file_path.reserve(n);
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        // Use the global entry index, not per-chunk b_idx (which repeats
        // across chunks after hadd). Entry order matches the .list order
        // because split is line-aligned and chunks are merged lexically.
        cd.idx.push_back(static_cast<int>(i));
        cd.n_pt3.push_back(b_pt3);
        cd.n_pt2.push_back(b_pt2);
        cd.file_path.push_back(b_path);
    }
    f->Close(); delete f;
    std::cout << "  loaded " << cd.label << ": " << n << " files from " << path << "\n";
    return true;
}

constexpr double kK = 63.0;

void perFileWeights(const ChannelData& cd,
                    std::vector<double>& w,
                    std::vector<double>& sigma_w) {
    const size_t n = cd.idx.size();
    w.assign(n, std::nan("")); sigma_w.assign(n, std::nan(""));
    for (size_t i = 0; i < n; ++i) {
        const double a = (double) cd.n_pt3[i];
        const double b = (double) cd.n_pt2[i];
        if (a <= 0 || b <= 0) continue;
        const double wi  = kK * b / a;
        w[i]       = wi;
        sigma_w[i] = wi * std::sqrt(1.0/a + 1.0/b);
    }
}

void cumulativeWeights(const ChannelData& cd,
                       std::vector<double>& wc,
                       std::vector<double>& sigma_wc) {
    const size_t n = cd.idx.size();
    wc.assign(n, std::nan("")); sigma_wc.assign(n, std::nan(""));
    long long S3 = 0, S2 = 0;
    for (size_t i = 0; i < n; ++i) {
        S3 += cd.n_pt3[i];
        S2 += cd.n_pt2[i];
        if (S3 <= 0 || S2 <= 0) continue;
        const double a = (double) S3, b = (double) S2;
        const double wi  = kK * b / a;
        wc[i]       = wi;
        sigma_wc[i] = wi * std::sqrt(1.0/a + 1.0/b);
    }
}

struct Segment {
    int       file_lo;     // inclusive (entry index)
    int       file_hi;     // inclusive (entry index)
    long long sum_pt3;
    long long sum_pt2;
    int       n_files;
    double    w()     const {
        return (sum_pt3 > 0 && sum_pt2 > 0)
                   ? kK * (double) sum_pt2 / (double) sum_pt3 : 0.0;
    }
    double    sigma() const {
        if (sum_pt3 <= 0 || sum_pt2 <= 0) return 0.0;
        const double a = (double) sum_pt3, b = (double) sum_pt2;
        return w() * std::sqrt(1.0/a + 1.0/b);
    }
};

struct SegmentationResult {
    std::vector<Segment> segs;
    std::vector<double>  z_seg;     // per-file z vs current open segment
    std::vector<int>     seg_of;    // per-file index of containing segment
};

SegmentationResult segmentChannel(const ChannelData& cd, double z_thr) {
    SegmentationResult res;
    if (cd.idx.empty()) return res;
    const size_t n = cd.idx.size();
    res.z_seg.assign(n, std::nan(""));
    res.seg_of.assign(n, -1);

    int       cur_lo  = cd.idx.front();
    int       cur_lo_i = 0;
    long long cur3    = 0, cur2 = 0;

    for (size_t i = 0; i < n; ++i) {
        const long long a_f = cd.n_pt3[i];
        const long long b_f = cd.n_pt2[i];

        if (cur3 > 0 && cur2 > 0 && a_f > 0 && b_f > 0) {
            const double a = (double) a_f, b = (double) b_f;
            const double w_f  = kK * b / a;
            const double sw_f = w_f * std::sqrt(1.0/a + 1.0/b);
            const double A = (double) cur3, B = (double) cur2;
            const double w_s  = kK * B / A;
            const double sw_s = w_s * std::sqrt(1.0/A + 1.0/B);
            const double s2 = sw_f*sw_f + sw_s*sw_s;
            if (s2 > 0) {
                const double z = (w_f - w_s) / std::sqrt(s2);
                res.z_seg[i] = z;
                if (std::abs(z) > z_thr) {
                    Segment seg;
                    seg.file_lo = cur_lo;
                    seg.file_hi = cd.idx[i-1];
                    seg.sum_pt3 = cur3;
                    seg.sum_pt2 = cur2;
                    seg.n_files = (int)(i - cur_lo_i);
                    res.segs.push_back(seg);
                    cur_lo_i = (int) i;
                    cur_lo   = cd.idx[i];
                    cur3 = 0; cur2 = 0;
                }
            }
        }
        cur3 += a_f;
        cur2 += b_f;
    }

    Segment seg;
    seg.file_lo = cur_lo;
    seg.file_hi = cd.idx.back();
    seg.sum_pt3 = cur3;
    seg.sum_pt2 = cur2;
    seg.n_files = (int)(n - cur_lo_i);
    res.segs.push_back(seg);

    for (size_t s = 0; s < res.segs.size(); ++s) {
        for (size_t i = 0; i < n; ++i) {
            if (cd.idx[i] >= res.segs[s].file_lo &&
                cd.idx[i] <= res.segs[s].file_hi) {
                res.seg_of[i] = (int) s;
            }
        }
    }
    return res;
}

void printSegmentation(const ChannelData& cd,
                       const SegmentationResult& sr,
                       double z_thr) {
    long long tot3 = 0, tot2 = 0;
    for (size_t i = 0; i < cd.idx.size(); ++i) {
        tot3 += cd.n_pt3[i]; tot2 += cd.n_pt2[i];
    }
    const double w_global = (tot3 > 0 && tot2 > 0)
                                ? kK * (double) tot2 / (double) tot3 : 0.0;
    std::cout << "  segmentation [" << cd.label << "]"
              << "  z_thr=" << z_thr
              << "  files=" << cd.idx.size()
              << "  N_PT3=" << tot3
              << "  N_PT2=" << tot2
              << "  w_global=" << std::fixed << std::setprecision(4)
              << w_global
              << "  segments=" << sr.segs.size() << "\n";
    for (size_t s = 0; s < sr.segs.size(); ++s) {
        const Segment& seg = sr.segs[s];
        std::cout << "    seg " << std::setw(3) << s
                  << "  files [" << std::setw(4) << seg.file_lo
                  << " .. "      << std::setw(4) << seg.file_hi << "]"
                  << " ("        << std::setw(3) << seg.n_files << " files)"
                  << "  N_PT3="  << std::setw(10) << seg.sum_pt3
                  << "  N_PT2="  << std::setw(8)  << seg.sum_pt2
                  << "  w="      << std::fixed << std::setprecision(4) << seg.w()
                  << " +/- "     << seg.sigma() << "\n";
    }
}

// Write the calibration JSON to <out_path>. Format:
//
//   { "channel": "...", "z_threshold": 3.0, "scan_source_file": "...",
//     "n_files_total": N, "n_segments": K,
//     "segments": [
//        { "seg_idx": 0, "n_files": 42, "n_pt3_total": ..., "n_pt2_total": ...,
//          "w": ..., "sigma": ...,
//          "files": [ "/path/to/A.root", ... ] },
//        ...
//     ]
//   }
//
// Hand-rolled writer (no JSON dep). Each file in a segment gets the same w
// when read by the main analysis: main.cc builds a map { file_path -> w }
// by iterating segments × files. The "sigma", "n_pt3_total", "n_pt2_total"
// fields are informational; only "files" and "w" are load-bearing.
bool writeSegmentsJson(const ChannelData& cd,
                       const SegmentationResult& sr,
                       double z_thr,
                       const std::string& out_path) {
    std::ofstream ofs(out_path);
    if (!ofs) {
        std::cerr << "  ERROR: cannot open for writing: " << out_path << "\n";
        return false;
    }
    auto fmt = [](double x) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.6f", x);
        return std::string(buf);
    };

    ofs << "{\n";
    ofs << "    \"channel\":          \"" << cd.label << "\",\n";
    ofs << "    \"z_threshold\":      " << fmt(z_thr) << ",\n";
    ofs << "    \"scan_source_file\": \"" << cd.source_path << "\",\n";
    ofs << "    \"n_files_total\":    " << cd.idx.size() << ",\n";
    ofs << "    \"n_segments\":       " << sr.segs.size() << ",\n";
    ofs << "    \"segments\": [\n";
    for (size_t s = 0; s < sr.segs.size(); ++s) {
        const Segment& seg = sr.segs[s];
        ofs << "        {\n";
        ofs << "            \"seg_idx\":     " << s << ",\n";
        ofs << "            \"n_files\":     " << seg.n_files << ",\n";
        ofs << "            \"n_pt3_total\": " << seg.sum_pt3 << ",\n";
        ofs << "            \"n_pt2_total\": " << seg.sum_pt2 << ",\n";
        ofs << "            \"w\":           " << fmt(seg.w())     << ",\n";
        ofs << "            \"sigma\":       " << fmt(seg.sigma()) << ",\n";
        ofs << "            \"files\": [\n";
        // Find which entries belong to this segment by entry-index range.
        // For each, emit the absolute file_path stored in the TTree branch.
        bool first = true;
        for (size_t i = 0; i < cd.idx.size(); ++i) {
            if (cd.idx[i] < seg.file_lo || cd.idx[i] > seg.file_hi) continue;
            if (!first) ofs << ",\n";
            ofs << "                \"" << cd.file_path[i] << "\"";
            first = false;
        }
        ofs << "\n            ]\n";
        ofs << "        }" << (s + 1 < sr.segs.size() ? "," : "") << "\n";
    }
    ofs << "    ]\n";
    ofs << "}\n";

    std::cout << "  wrote " << out_path
              << "  (" << sr.segs.size() << " segments, "
              << cd.idx.size() << " files)\n";
    return true;
}

TGraphErrors* makeGraph(const std::vector<int>& idx,
                        const std::vector<double>& y,
                        const std::vector<double>& yerr) {
    const size_t n = idx.size();
    std::vector<double> xv, yv, exv, eyv;
    xv.reserve(n); yv.reserve(n); exv.reserve(n); eyv.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        if (std::isnan(y[i])) continue;
        xv.push_back((double) idx[i]);
        yv.push_back(y[i]);
        exv.push_back(0.0);
        eyv.push_back(std::isnan(yerr[i]) ? 0.0 : yerr[i]);
    }
    return new TGraphErrors((Int_t) xv.size(),
                            xv.data(), yv.data(),
                            exv.data(), eyv.data());
}

void drawChannelCanvas(const ChannelData& cd,
                       const SegmentationResult& sr,
                       double z_thr) {
    if (cd.idx.empty()) return;

    std::vector<double> w, sw, wc, swc;
    perFileWeights(cd, w, sw);
    cumulativeWeights(cd, wc, swc);

    // ----- counts pad data -----------------------------------------------
    std::vector<double> npt3d(cd.n_pt3.begin(), cd.n_pt3.end());
    std::vector<double> npt2d(cd.n_pt2.begin(), cd.n_pt2.end());
    std::vector<double> e3d(npt3d.size()), e2d(npt2d.size());
    for (size_t i = 0; i < npt3d.size(); ++i) {
        e3d[i] = (npt3d[i] > 0) ? std::sqrt(npt3d[i]) : 0.0;
        e2d[i] = (npt2d[i] > 0) ? std::sqrt(npt2d[i]) : 0.0;
    }
    auto* g_pt3 = makeGraph(cd.idx, npt3d, e3d);
    auto* g_pt2 = makeGraph(cd.idx, npt2d, e2d);
    g_pt3->SetMarkerStyle(20); g_pt3->SetMarkerColor(kRed+1);
    g_pt3->SetLineColor(kRed+1); g_pt3->SetMarkerSize(0.6);
    g_pt2->SetMarkerStyle(21); g_pt2->SetMarkerColor(kBlue+1);
    g_pt2->SetLineColor(kBlue+1); g_pt2->SetMarkerSize(0.6);

    // ----- weight pad data -----------------------------------------------
    auto* g_w  = makeGraph(cd.idx, w, sw);
    g_w->SetMarkerStyle(20); g_w->SetMarkerColor(kBlack);
    g_w->SetLineColor(kBlack); g_w->SetMarkerSize(0.6);

    auto* g_wc = makeGraph(cd.idx, wc, swc);
    g_wc->SetLineColor(kRed+1); g_wc->SetLineWidth(2);
    g_wc->SetFillColorAlpha(kRed-9, 0.35);
    g_wc->SetFillStyle(1001);

    // ----- z-score pad data ---------------------------------------------
    std::vector<double> z_zero_err(sr.z_seg.size(), 0.0);
    auto* g_z  = makeGraph(cd.idx, sr.z_seg, z_zero_err);
    g_z->SetMarkerStyle(20); g_z->SetMarkerColor(kBlack);
    g_z->SetLineColor(kGray+1); g_z->SetMarkerSize(0.5);

    // ----- canvas --------------------------------------------------------
    const TString cname = TString::Format("c_pt3_%s", cd.label.c_str());
    auto* c = new TCanvas(cname, cname, 1400, 1100);
    c->Divide(1, 3, 0.001, 0.001);

    // -- pad 1: counts ----
    auto* p1 = (TPad*) c->cd(1);
    p1->SetLogy(); p1->SetGridx(); p1->SetGridy();
    p1->SetLeftMargin(0.10); p1->SetRightMargin(0.04);
    p1->SetTopMargin(0.08); p1->SetBottomMargin(0.13);

    auto* mg1 = new TMultiGraph();
    mg1->SetTitle(TString::Format(
        "PT3 / PT2 per-file counts  (%s);file index;counts",
        cd.label.c_str()));
    mg1->Add(g_pt3, "P");
    mg1->Add(g_pt2, "P");
    mg1->Draw("A");
    auto* leg1 = new TLegend(0.78, 0.80, 0.95, 0.95);
    leg1->AddEntry(g_pt3, "N_{PT3}", "lp");
    leg1->AddEntry(g_pt2, "N_{PT2}", "lp");
    leg1->Draw();

    // -- pad 2: per-file w + cumulative w ± 1σ band + segments ----
    auto* p2 = (TPad*) c->cd(2);
    p2->SetGridx(); p2->SetGridy();
    p2->SetLeftMargin(0.10); p2->SetRightMargin(0.04);
    p2->SetTopMargin(0.08); p2->SetBottomMargin(0.13);

    auto* mg2 = new TMultiGraph();
    mg2->SetTitle(TString::Format(
        "weight w = 63 #upoint N_{PT2}/N_{PT3}   per-file (black) "
        "vs cumulative (red) vs segmented (green)   (%s);file index;w",
        cd.label.c_str()));
    mg2->Add(g_wc, "L3");
    mg2->Add(g_w,  "P");
    mg2->Draw("A");
    auto* leg2 = new TLegend(0.74, 0.78, 0.95, 0.95);
    leg2->AddEntry(g_w,  "per file w[f] #pm #sigma", "lp");
    leg2->AddEntry(g_wc, "cumulative w_{cum} #pm 1#sigma", "lf");
    TLine* seg_legend_line = nullptr;
    for (const Segment& seg : sr.segs) {
        const double xlo = seg.file_lo - 0.5;
        const double xhi = seg.file_hi + 0.5;
        const double y   = seg.w();
        const double dy  = seg.sigma();
        if (dy > 0) {
            auto* lup = new TLine(xlo, y + dy, xhi, y + dy);
            auto* ldn = new TLine(xlo, y - dy, xhi, y - dy);
            for (TLine* l : {lup, ldn}) {
                l->SetLineColor(kGreen+2); l->SetLineStyle(3);
                l->SetLineWidth(1); l->Draw();
            }
        }
        auto* lc = new TLine(xlo, y, xhi, y);
        lc->SetLineColor(kGreen+2); lc->SetLineWidth(3);
        lc->Draw();
        if (!seg_legend_line) seg_legend_line = lc;
    }
    if (seg_legend_line) {
        leg2->AddEntry(seg_legend_line,
                       TString::Format("segments (z>%.1f)  N=%zu",
                                       z_thr, sr.segs.size()),
                       "l");
    }
    leg2->Draw();

    // -- pad 3: segmenter z-score ----
    auto* p3 = (TPad*) c->cd(3);
    p3->SetGridx(); p3->SetGridy();
    p3->SetLeftMargin(0.10); p3->SetRightMargin(0.04);
    p3->SetTopMargin(0.08); p3->SetBottomMargin(0.16);

    g_z->SetTitle(TString::Format(
        "segmenter z = (w[f] - w_{seg}) / #sqrt{#sigma_{w}^{2} + "
        "#sigma_{seg}^{2}}    (%s, z_{thr}=%.1f, segments=%zu);"
        "file index;z",
        cd.label.c_str(), z_thr, sr.segs.size()));
    g_z->Draw("AP");
    if (g_z->GetN() > 0) {
        const double xlo = g_z->GetXaxis()->GetXmin();
        const double xhi = g_z->GetXaxis()->GetXmax();
        for (double s : {-z_thr, +z_thr}) {
            auto* ln = new TLine(xlo, s, xhi, s);
            ln->SetLineColor(kRed+1); ln->SetLineStyle(2); ln->SetLineWidth(2);
            ln->Draw();
        }
        auto* ln0 = new TLine(xlo, 0, xhi, 0);
        ln0->SetLineColor(kGray+2); ln0->SetLineStyle(1); ln0->Draw();
        const double ylo = g_z->GetYaxis()->GetXmin();
        const double yhi = g_z->GetYaxis()->GetXmax();
        for (size_t s = 1; s < sr.segs.size(); ++s) {
            const double xb = sr.segs[s].file_lo - 0.5;
            auto* lv = new TLine(xb, ylo, xb, yhi);
            lv->SetLineColor(kGreen+2); lv->SetLineStyle(2);
            lv->SetLineWidth(2); lv->Draw();
        }
    }

    // ----- save ----------------------------------------------------------
    gSystem->mkdir("plots/output", kTRUE);
    const TString stem = TString::Format(
        "plots/output/trigger_calibration_%s", cd.label.c_str());
    c->SaveAs(stem + ".pdf");
    c->SaveAs(stem + ".png");
    std::cout << "  saved " << stem << ".{pdf,png}\n";
}

void drawOverlayCanvas(const std::vector<ChannelData*>& chans) {
    auto* c = new TCanvas("c_pt3_overlay", "c_pt3_overlay", 1400, 700);
    c->SetGridx(); c->SetGridy();
    c->SetLeftMargin(0.08); c->SetRightMargin(0.04);
    c->SetTopMargin(0.08);  c->SetBottomMargin(0.13);

    auto* mg = new TMultiGraph();
    mg->SetTitle("per-file weight w = 63 #upoint N_{PT2}/N_{PT3}, "
                 "all channels overlaid;file index;w");
    auto* leg = new TLegend(0.83, 0.78, 0.97, 0.95);

    const int colors[3]   = { kBlack, kBlue+1, kRed+1 };
    const int markers[3]  = { 20, 21, 22 };
    int ci = 0;
    for (auto* cd : chans) {
        if (!cd || cd->idx.empty()) { ++ci; continue; }
        std::vector<double> w, sw;
        perFileWeights(*cd, w, sw);
        auto* g = makeGraph(cd->idx, w, sw);
        g->SetMarkerStyle(markers[ci]);
        g->SetMarkerColor(colors[ci]);
        g->SetLineColor(colors[ci]);
        g->SetMarkerSize(0.6);
        mg->Add(g, "P");
        leg->AddEntry(g, cd->label.c_str(), "lp");
        ++ci;
    }
    mg->Draw("A");
    leg->Draw();

    gSystem->mkdir("plots/output", kTRUE);
    c->SaveAs("plots/output/trigger_calibration_overlay.pdf");
    c->SaveAs("plots/output/trigger_calibration_overlay.png");
    std::cout << "  saved plots/output/trigger_calibration_overlay.{pdf,png}\n";
}

}  // anonymous namespace

void trigger_calibration(double z_threshold = 3.0) {
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    std::cout << "=== trigger_calibration  (z_threshold=" << z_threshold << ") ===\n";

    ChannelData epem{"epem", "", {}, {}, {}, {}};
    ChannelData epep{"epep", "", {}, {}, {}, {}};
    ChannelData emem{"emem", "", {}, {}, {}, {}};
    const bool ok_e = loadChannel("../trigger_scan_epem.root", epem);
    const bool ok_p = loadChannel("../trigger_scan_epep.root", epep);
    const bool ok_m = loadChannel("../trigger_scan_emem.root", emem);

    if (!ok_e && !ok_p && !ok_m) {
        std::cerr << "ERROR: no scan files found. Run "
                  << "./run_parallel_scan.sh config_<channel>.json first.\n";
        return;
    }

    SegmentationResult sr_e, sr_p, sr_m;
    if (ok_e) { sr_e = segmentChannel(epem, z_threshold);
                printSegmentation(epem, sr_e, z_threshold); }
    if (ok_p) { sr_p = segmentChannel(epep, z_threshold);
                printSegmentation(epep, sr_p, z_threshold); }
    if (ok_m) { sr_m = segmentChannel(emem, z_threshold);
                printSegmentation(emem, sr_m, z_threshold); }

    if (ok_e) drawChannelCanvas(epem, sr_e, z_threshold);
    if (ok_p) drawChannelCanvas(epep, sr_p, z_threshold);
    if (ok_m) drawChannelCanvas(emem, sr_m, z_threshold);

    std::vector<ChannelData*> chans;
    if (ok_e) chans.push_back(&epem);
    if (ok_p) chans.push_back(&epep);
    if (ok_m) chans.push_back(&emem);
    if (chans.size() >= 2) drawOverlayCanvas(chans);

    // ----- write calibration JSONs into REPO ROOT (../../) ----------------
    // This is intentional: the JSON is consumed by the main FAT analysis
    // (./ana on pp45_epem branch) which expects to find it next to the
    // main config*.json files, NOT in research/calibration/.
    if (ok_e) writeSegmentsJson(epem, sr_e, z_threshold,
                                "../../pt3_segments_epem.json");
    if (ok_p) writeSegmentsJson(epep, sr_p, z_threshold,
                                "../../pt3_segments_epep.json");
    if (ok_m) writeSegmentsJson(emem, sr_m, z_threshold,
                                "../../pt3_segments_emem.json");

    std::cout << "Done.\n";
}

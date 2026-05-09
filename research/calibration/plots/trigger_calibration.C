// research/calibration/plots/trigger_calibration.C
//
// Phase (b) of the trigger-bias calibration pipeline: compute-only.
//
// Reads trigger_scan_<channel>.root produced by trigger_scan (phase a),
// streams through trigger_events in chain order, runs an online z-test
// change-point segmenter to discover where the PT3/PT2 ratio shifts, and
// writes the calibration table as a ROOT TTree to the REPO ROOT (../../).
//
// NO PLOTS. Visualisation lives in plots/visualize_calibration.C.
//
// Algorithm (online, single forward pass over trigger events, two-phase):
//   Phase 1 — open segment is "ripening": its rel_err = sqrt(1/N_PT3 +
//             1/N_PT2) is still > max_rel_err. Add events to it without
//             considering splits. Without this gate the very first segment
//             would start with a single event and any second event would
//             trigger a meaningless z-test on noise.
//   Phase 2 — open is statistically determined (rel_err ≤ max_rel_err).
//             New events go to candidate. Once candidate ALSO reaches
//             rel_err ≤ max_rel_err, run the z-test:
//               z = (w_cand - w_seg) / sqrt(σ_cand² + σ_seg²)
//             - |z| > z_threshold:  SPLIT — close open at the trigger event
//                                   right before candidate started; the
//                                   candidate becomes the new open
//                                   (already ripe by construction);
//                                   candidate resets.
//             - else:               MERGE — fold candidate into open;
//                                   candidate resets.
//   At end of stream: any leftover candidate that didn't ripen is folded
//   into open; final open is closed at the last trigger event.
//
// max_rel_err controls BOTH the validity of the z-test AND the minimum
// segment precision (every closed segment has σ_w/w ≤ max_rel_err by
// construction, possibly much better if it merged for a long time).
// Tightening max_rel_err (e.g. 0.02) requires more PT2 stats per pool;
// loosening it (e.g. 0.10) accepts noisier candidates and produces more
// segments. At pp45 (PT3:PT2 ≈ 64:1), the rule of thumb is
//   max_rel_err = 0.03  →  N_PT2 ≳ 1100 per pool   (default)
//   max_rel_err = 0.05  →  N_PT2 ≳  400 per pool
//   max_rel_err = 0.10  →  N_PT2 ≳  100 per pool
//
// Output ROOT contains one TTree:
//
//   pt3_calibration:
//     - seg_idx     /I    monotonic, 0-based
//     - file_path   /string  one entry per (segment × file) the segment
//                            covers; multiple files in one segment ⇒
//                            multiple entries with the same seg_idx and w
//     - event_lo    /L    inclusive, in this file's local indexing
//     - event_hi    /L    inclusive, in this file's local indexing
//     - w           /D    weight (63 · N_PT2 / N_PT3)
//     - sigma       /D    Poisson-propagated 1-σ on w
//     - n_pt3       /L    PT3 count of the segment (whole, repeats per file)
//     - n_pt2       /L    PT2 count of the segment (whole, repeats per file)
//
// Convention: segments tile the chain with no gaps. Segment N covers from
// its first event up to and including the event right before segment N+1's
// first event. The first segment starts at event 0 of the first file; the
// last segment runs to end of the last file. Per-file ranges are computed
// from these endpoints using the `files` TTree of the scan output.
//
// Usage (from research/calibration/):
//   root -l -b -q plots/trigger_calibration.C                  # default
//   root -l -b -q 'plots/trigger_calibration.C(3.0, 0.02)'     # tighter pool
//   root -l -b -q 'plots/trigger_calibration.C(4.0, 0.03)'     # stricter z
//   root -l -b -q 'plots/trigger_calibration.C(3.0, 0.05)'     # looser pool
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <string>
#include <cmath>

namespace {

constexpr double kK = 63.0;

// Per-file metadata loaded from the scan output's `files` TTree.
struct FileInfo {
    std::string path;
    Long64_t    n_events;
    Long64_t    n_pt3;
    Long64_t    n_pt2;
};

// Open segment + candidate state during streaming. Both track an interval
// in chain coordinates by recording first and last trigger event positions
// (file_path + local_event_idx).
struct PoolState {
    long long n_pt3 = 0;
    long long n_pt2 = 0;
    std::string first_file;
    long long   first_local = 0;
    std::string last_file;
    long long   last_local = 0;
    bool empty() const { return n_pt3 == 0 && n_pt2 == 0; }
    void reset() { *this = PoolState{}; }
    void addTrigger(const std::string& fp, long long lidx, int trigbit,
                    int trig_pt3, int trig_pt2) {
        if (empty()) {
            first_file  = fp;
            first_local = lidx;
        }
        last_file  = fp;
        last_local = lidx;
        if      (trigbit == trig_pt3) ++n_pt3;
        else if (trigbit == trig_pt2) ++n_pt2;
    }
    double w() const {
        return (n_pt3 > 0 && n_pt2 > 0)
                   ? kK * (double) n_pt2 / (double) n_pt3 : 0.0;
    }
    double sigma() const {
        if (n_pt3 <= 0 || n_pt2 <= 0) return 0.0;
        const double a = (double) n_pt3, b = (double) n_pt2;
        return w() * std::sqrt(1.0/a + 1.0/b);
    }
};

// Final closed segment (logical), still in chain coordinates. Will be
// expanded to per-file ranges before writing to the output TTree.
struct ClosedSegment {
    int        seg_idx = 0;
    long long  n_pt3   = 0;
    long long  n_pt2   = 0;
    std::string first_file;
    long long   first_local = 0;
    // The "end" of a segment is implicit: it runs up to right before the
    // next segment's first event. The very last segment runs to end-of-chain.
    // Actual per-file ranges are computed in expandSegments().
    double w()     const {
        return (n_pt3 > 0 && n_pt2 > 0) ? kK * (double) n_pt2 / (double) n_pt3 : 0.0;
    }
    double sigma() const {
        if (n_pt3 <= 0 || n_pt2 <= 0) return 0.0;
        const double a = (double) n_pt3, b = (double) n_pt2;
        return w() * std::sqrt(1.0/a + 1.0/b);
    }
};

// Per-file portion of a segment; one entry per (seg, file) — what we write.
struct OutRange {
    int          seg_idx;
    std::string  file_path;
    Long64_t     event_lo;
    Long64_t     event_hi;
    double       w;
    double       sigma;
    Long64_t     n_pt3;
    Long64_t     n_pt2;
};

// Find the file index in `files` for a given path. O(N) — N is small (~30).
int findFileIdx(const std::vector<FileInfo>& files, const std::string& path) {
    for (size_t i = 0; i < files.size(); ++i)
        if (files[i].path == path) return (int) i;
    return -1;
}

// Convert closed-segment list (chain-coord boundaries) into per-file ranges
// covering the chain end-to-end.
std::vector<OutRange> expandSegments(const std::vector<ClosedSegment>& segs,
                                     const std::vector<FileInfo>& files) {
    std::vector<OutRange> out;
    if (segs.empty() || files.empty()) return out;

    // For each segment N, its range runs from (segs[N].first_file, first_local)
    // up to (segs[N+1].first_file, first_local - 1). For the last segment,
    // run to (files.back().path, files.back().n_events - 1).
    for (size_t s = 0; s < segs.size(); ++s) {
        const ClosedSegment& seg = segs[s];
        const int fi_start = findFileIdx(files, seg.first_file);
        if (fi_start < 0) {
            std::cerr << "  WARNING: segment first_file not in files TTree: "
                      << seg.first_file << "\n";
            continue;
        }
        // End coordinates
        std::string end_file;
        Long64_t    end_local;
        if (s + 1 < segs.size()) {
            const ClosedSegment& nxt = segs[s + 1];
            const int fi_next = findFileIdx(files, nxt.first_file);
            if (fi_next == fi_start) {
                end_file  = seg.first_file;
                end_local = nxt.first_local - 1;
            } else if (fi_next > fi_start) {
                // End at last event of file just before nxt.first_file
                end_file  = files[fi_next - 1].path;
                end_local = files[fi_next - 1].n_events - 1;
            } else {
                std::cerr << "  WARNING: segment ordering inconsistency between "
                          << seg.first_file << " and " << nxt.first_file << "\n";
                continue;
            }
        } else {
            end_file  = files.back().path;
            end_local = files.back().n_events - 1;
        }
        const int fi_end = findFileIdx(files, end_file);
        if (fi_end < 0) continue;

        // Emit one OutRange per file in [fi_start..fi_end].
        for (int fi = fi_start; fi <= fi_end; ++fi) {
            OutRange r;
            r.seg_idx   = seg.seg_idx;
            r.file_path = files[fi].path;
            r.event_lo  = (fi == fi_start) ? seg.first_local : 0;
            r.event_hi  = (fi == fi_end)
                              ? end_local
                              : (files[fi].n_events - 1);
            r.w     = seg.w();
            r.sigma = seg.sigma();
            r.n_pt3 = seg.n_pt3;
            r.n_pt2 = seg.n_pt2;
            out.push_back(r);
        }
    }
    return out;
}

// Relative Poisson error of w = K·N_PT2/N_PT3 propagated as
//   sigma_w / w  =  sqrt(1/N_PT3 + 1/N_PT2).
// Returns +inf if either count is zero (treat as "not yet stable").
double relErr(long long n_pt3, long long n_pt2) {
    if (n_pt3 <= 0 || n_pt2 <= 0)
        return std::numeric_limits<double>::infinity();
    return std::sqrt(1.0/(double) n_pt3 + 1.0/(double) n_pt2);
}

void processChannel(const std::string& chan,
                    const std::string& in_path,
                    const std::string& out_path,
                    double      z_thr,
                    double      max_rel_err,
                    int         trig_pt3 = 8192,
                    int         trig_pt2 = 4096) {
    std::cout << "\n--- channel " << chan << " ---\n"
              << "  input:  " << in_path  << "\n"
              << "  output: " << out_path << "\n";

    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "  WARNING: cannot open " << in_path << " — skipping channel\n";
        if (fin) delete fin;
        return;
    }

    TTree* t_files = dynamic_cast<TTree*>(fin->Get("files"));
    TTree* t_evts  = dynamic_cast<TTree*>(fin->Get("trigger_events"));
    if (!t_files || !t_evts) {
        std::cerr << "  WARNING: missing trigger_events or files TTree in "
                  << in_path << " — skipping\n";
        fin->Close(); delete fin;
        return;
    }

    // Load files metadata.
    std::vector<FileInfo> files;
    {
        std::string  f_path;
        std::string* p_path = &f_path;
        Long64_t f_n = 0, f_p3 = 0, f_p2 = 0;
        t_files->SetBranchAddress("file_path",      &p_path);
        t_files->SetBranchAddress("n_events_total", &f_n);
        t_files->SetBranchAddress("n_pt3",          &f_p3);
        t_files->SetBranchAddress("n_pt2",          &f_p2);
        const Long64_t N = t_files->GetEntries();
        files.reserve(N);
        for (Long64_t i = 0; i < N; ++i) {
            t_files->GetEntry(i);
            files.push_back({f_path, f_n, f_p3, f_p2});
        }
    }
    std::cout << "  files: " << files.size()
              << "  trigger events: " << t_evts->GetEntries() << "\n";

    if (t_evts->GetEntries() == 0) {
        std::cerr << "  WARNING: no trigger events — channel skipped\n";
        fin->Close(); delete fin;
        return;
    }

    // Stream through trigger_events.
    std::string  e_path;
    std::string* p_e_path = &e_path;
    Long64_t e_local = 0;
    Int_t    e_trigbit = 0;
    t_evts->SetBranchAddress("file_path",       &p_e_path);
    t_evts->SetBranchAddress("local_event_idx", &e_local);
    t_evts->SetBranchAddress("trigbit",         &e_trigbit);

    PoolState open_seg;
    PoolState cand;
    std::vector<ClosedSegment> segs;
    int next_seg_idx = 0;

    auto closeOpen = [&](const std::string& reason) {
        ClosedSegment cs;
        cs.seg_idx     = next_seg_idx++;
        cs.n_pt3       = open_seg.n_pt3;
        cs.n_pt2       = open_seg.n_pt2;
        cs.first_file  = open_seg.first_file;
        cs.first_local = open_seg.first_local;
        segs.push_back(cs);
        std::cout << "    [seg " << std::setw(3) << cs.seg_idx
                  << "  " << reason << "]"
                  << "  N_PT3=" << std::setw(8) << cs.n_pt3
                  << "  N_PT2=" << std::setw(6) << cs.n_pt2
                  << "  w="     << std::fixed << std::setprecision(4) << cs.w()
                  << " +/- "    << cs.sigma()
                  << "  start=" << cs.first_local << " in "
                  << cs.first_file.substr(cs.first_file.find_last_of('/') + 1)
                  << "\n";
    };

    // Two-phase streaming:
    //  Phase 1 — open segment is "ripening": its rel_err > max_rel_err so
    //            we keep adding events to it without testing for change.
    //            Without this gate the very first segment starts with 1
    //            event and any second event triggers a meaningless z-test.
    //  Phase 2 — open is statistically determined (rel_err ≤ max_rel_err);
    //            new events go to candidate. Once candidate ALSO reaches
    //            rel_err ≤ max_rel_err, run the z-test and split or merge.
    // After a split the new open is the old candidate (already ripe by
    // construction), so we re-enter phase 2 immediately.
    const Long64_t N_evts = t_evts->GetEntries();
    for (Long64_t i = 0; i < N_evts; ++i) {
        t_evts->GetEntry(i);

        if (relErr(open_seg.n_pt3, open_seg.n_pt2) > max_rel_err) {
            // Phase 1: keep growing open until it's well-determined.
            open_seg.addTrigger(e_path, e_local, e_trigbit, trig_pt3, trig_pt2);
            continue;
        }

        // Phase 2: open is stable. Build up candidate.
        cand.addTrigger(e_path, e_local, e_trigbit, trig_pt3, trig_pt2);
        if (relErr(cand.n_pt3, cand.n_pt2) > max_rel_err) continue;

        // Both pools are statistically ripe — run the z-test.
        const double w_o = open_seg.w(), s_o = open_seg.sigma();
        const double w_c = cand.w(),     s_c = cand.sigma();
        const double s2  = s_o*s_o + s_c*s_c;
        if (s2 <= 0) {
            // Defensive: shouldn't happen given the relErr guard above.
            cand.reset();
            continue;
        }
        const double z = (w_c - w_o) / std::sqrt(s2);

        if (std::abs(z) > z_thr) {
            // SPLIT — close open segment, candidate becomes new open.
            closeOpen("SPLIT |z|=" + std::to_string(std::abs(z)));
            open_seg = cand;
            cand.reset();
        } else {
            // MERGE — absorb candidate into open.
            open_seg.n_pt3 += cand.n_pt3;
            open_seg.n_pt2 += cand.n_pt2;
            open_seg.last_file  = cand.last_file;
            open_seg.last_local = cand.last_local;
            cand.reset();
        }
    }

    // Absorb any leftover candidate (didn't reach min stats) into open.
    if (!cand.empty()) {
        open_seg.n_pt3 += cand.n_pt3;
        open_seg.n_pt2 += cand.n_pt2;
        open_seg.last_file  = cand.last_file;
        open_seg.last_local = cand.last_local;
        cand.reset();
    }
    closeOpen("END  ");

    std::cout << "  -> " << segs.size() << " segments closed\n";

    // Expand to per-file ranges.
    auto out_ranges = expandSegments(segs, files);
    std::cout << "  -> " << out_ranges.size() << " (segment x file) ranges\n";

    // Write output TTree.
    TFile* fout = TFile::Open(out_path.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "  ERROR: cannot open output " << out_path << "\n";
        if (fout) delete fout;
        fin->Close(); delete fin;
        return;
    }
    fout->cd();
    TTree* tcal = new TTree("pt3_calibration", "PT3 trigger-bias calibration table");
    Int_t       o_seg = 0;
    std::string o_path;
    Long64_t    o_lo = 0, o_hi = 0;
    Double_t    o_w = 0, o_sigma = 0;
    Long64_t    o_pt3 = 0, o_pt2 = 0;
    tcal->Branch("seg_idx",   &o_seg,   "seg_idx/I");
    tcal->Branch("file_path", &o_path);
    tcal->Branch("event_lo",  &o_lo,    "event_lo/L");
    tcal->Branch("event_hi",  &o_hi,    "event_hi/L");
    tcal->Branch("w",         &o_w,     "w/D");
    tcal->Branch("sigma",     &o_sigma, "sigma/D");
    tcal->Branch("n_pt3",     &o_pt3,   "n_pt3/L");
    tcal->Branch("n_pt2",     &o_pt2,   "n_pt2/L");
    for (const OutRange& r : out_ranges) {
        o_seg   = r.seg_idx;
        o_path  = r.file_path;
        o_lo    = r.event_lo;
        o_hi    = r.event_hi;
        o_w     = r.w;
        o_sigma = r.sigma;
        o_pt3   = r.n_pt3;
        o_pt2   = r.n_pt2;
        tcal->Fill();
    }
    tcal->Write();
    fout->Close();
    delete fout;
    std::cout << "  wrote " << out_path << "\n";

    fin->Close(); delete fin;
}

}  // anonymous namespace

void trigger_calibration(double z_threshold = 3.0,
                         double max_rel_err = 0.03) {
    std::cout << "=== trigger_calibration"
              << "  z_threshold=" << z_threshold
              << "  max_rel_err=" << max_rel_err
              << " ===\n";

    // Paths are relative to the cwd from which root is invoked. Standard
    // usage is `cd research/calibration && root -l -b -q plots/...` so the
    // scan ROOTs live next to the cwd (./) and the calibration ROOTs go
    // to the repo root (../../) per the project convention.
    const std::vector<std::string> chans = {"epem", "epep", "emem"};
    for (const std::string& chan : chans) {
        const std::string in_path  = "trigger_scan_"   + chan + ".root";
        const std::string out_path = "../../pt3_calibration_" + chan + ".root";
        processChannel(chan, in_path, out_path, z_threshold, max_rel_err);
    }
    std::cout << "\nDone.\n";
}

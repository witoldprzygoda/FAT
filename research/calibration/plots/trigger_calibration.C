// research/calibration/plots/trigger_calibration.C
//
// Phase (b) of the trigger-bias calibration pipeline: compute-only.
//
// Reads output_<channel>_cal.root from the repo root — that is the output
// of ./ana on pp45_epem run with config_<channel>_cal.json (mode =
// "trigger_calibration"). The cal-mode ./ana writes:
//   trigger_cal_nt    one entry per event (113.6M for epem) with raw
//                     branches: file_idx, local_event_idx, trigbit, oa,
//                     isBest, eVertReco_z, start_iteration. NO cuts at
//                     fill time — both PT3 and PT2 events present.
//   trigger_cal_files one entry per input file: file_idx, n_events_total,
//                     file_path. Used here for the segment-to-file mapping.
//
// This macro applies the SAME quality cuts the main analysis uses
// (isBest == 1, eVertReco_z > -500, start_iteration == 3).
//
// NO opening-angle cut is applied: empirically, the segmentation
// boundaries (the WHEN of trigger-bias changes) are independent of the
// physics selection — only the absolute correction values w(seg) depend
// on it. Since the segmentation defines the calibration's validity
// ranges (i.e. the schedule of bias changes in real time) rather than
// the correction magnitudes, leaving OA inclusive maximises statistics
// per pool and tightens the segment-boundary precision. Correction
// magnitudes ARE consumer-specific, so the main analysis on pp45_epem
// must recompute / re-use a calibration matching its own physics
// selection — but the segment EDGES from this run are reusable across
// any OA threshold.
//
// Surviving events with trigbit ∈ {8192, 4096} feed the online segmenter.
//
// Streams through trigger_cal_nt in chain order, runs an online z-test
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
//   max_rel_err = 0.005 →  N_PT2 ≳ 40000 per pool   (only deep trend changes)
//   max_rel_err = 0.01  →  N_PT2 ≳ 10000 per pool   (default — trend-only;
//                                                    brief excursions
//                                                    naturally absorb into
//                                                    the surrounding segment)
//   max_rel_err = 0.02  →  N_PT2 ≳  2500 per pool   ⚠ resurfaces outliers
//   max_rel_err = 0.03  →  N_PT2 ≳  1100 per pool   ⚠ resurfaces outliers
//
// Anything looser than 0.01 will recover brief outlier-segments because
// the candidate threshold becomes shorter than typical outlier durations.
// Use 0.02 / 0.03 only AFTER outliers have been filtered upstream by a
// dedicated macro — otherwise they pollute the calibration trend.
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
// from these endpoints using the `trigger_cal_files` TTree.
//
// Usage (from research/calibration/):
//   root -l -b -q plots/trigger_calibration.C                  # default
//   root -l -b -q 'plots/trigger_calibration.C(3.0, 0.02)'     # tighter pool
//   root -l -b -q 'plots/trigger_calibration.C(4.0, 0.03)'     # stricter z
//   root -l -b -q 'plots/trigger_calibration.C(3.0, 0.05)'     # looser pool
//
// Cuts applied (mirror main-analysis convention from src/setup_cuts.h):
//   isBest          == 1      best-candidate selection
//   eVertReco_z     > -500    vertex quality [mm]
//   start_iteration == 3      start-detector quality
// (no OA cut — see header note above on why segmentation boundaries are
//  OA-independent while correction magnitudes are not)
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

// Per-file metadata loaded from the cal output's `trigger_cal_files` TTree.
// (n_pt3 / n_pt2 used to live here — the legacy raw scanner pre-counted
// them per file. The cal-mode pipeline no longer pre-counts; segmenter
// computes counts from the streamed events instead.)
struct FileInfo {
    std::string path;
    Long64_t    n_events;
};

// One stored trigger event inside the candidate's per-event buffer.
// Populated only when track == true at addTrigger() — i.e. for the
// candidate pool only. Used by findOptimalSplit to refine where exactly
// the regime change happened within the candidate's span.
struct CandEvent {
    int       file_idx;   // index into the FileInfo vector
    Long64_t  local_idx;  // event index inside that file
    bool      is_pt3;     // true=PT3, false=PT2
};

// Open segment + candidate state during streaming. Both track an interval
// in chain coordinates by recording first and last trigger event positions
// (file_path + local_event_idx). Candidate additionally stores its per-
// event sequence so the split point can be sharpened by likelihood scan
// after the z-test fires.
struct PoolState {
    long long n_pt3 = 0;
    long long n_pt2 = 0;
    std::string first_file;
    long long   first_local = 0;
    std::string last_file;
    long long   last_local = 0;
    std::vector<CandEvent> events;   // populated only when track==true
    bool empty() const { return n_pt3 == 0 && n_pt2 == 0; }
    void reset() { *this = PoolState{}; }
    void addTrigger(const std::string& fp, long long lidx, int trigbit,
                    int trig_pt3, int trig_pt2,
                    int  file_idx_global = -1,
                    bool track           = false) {
        if (empty()) {
            first_file  = fp;
            first_local = lidx;
        }
        last_file  = fp;
        last_local = lidx;
        if      (trigbit == trig_pt3) ++n_pt3;
        else if (trigbit == trig_pt2) ++n_pt2;
        if (track && file_idx_global >= 0) {
            events.push_back({file_idx_global, lidx, trigbit == trig_pt3});
        }
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

// Find the within-candidate position k* that maximises the z² test
// statistic for "merged-left vs right" split:
//     left  = (open + cand[0..k])
//     right = cand[k+1..end]
//     z²(k) = (w_left - w_right)² / (σ²_left + σ²_right)
// The crude algorithm sets the new segment boundary at cand.events[0]
// (i.e. the first candidate event), which is at most a candidate-span
// of slack from the true change point. The refinement narrows this to
// the single event with maximum left/right separation.
//
// Returns the index k* in cand.events (0 ≤ k* < cand.events.size()-1).
// Caller should treat cand.events[0..k*] as belonging to the closing
// segment, and cand.events[k*+1..end] as the new open segment.
size_t findOptimalSplit(const PoolState& open, const PoolState& cand) {
    const size_t N = cand.events.size();
    if (N < 2) return 0;
    long long L_pt3 = open.n_pt3;
    long long L_pt2 = open.n_pt2;
    double   best_z2 = -1.0;
    size_t   best_k  = 0;
    for (size_t k = 0; k + 1 < N; ++k) {
        if (cand.events[k].is_pt3) ++L_pt3; else ++L_pt2;
        const long long R_pt3 = (open.n_pt3 + cand.n_pt3) - L_pt3;
        const long long R_pt2 = (open.n_pt2 + cand.n_pt2) - L_pt2;
        if (L_pt3 <= 0 || L_pt2 <= 0 || R_pt3 <= 0 || R_pt2 <= 0) continue;
        const double w_L = kK * (double) L_pt2 / (double) L_pt3;
        const double w_R = kK * (double) R_pt2 / (double) R_pt3;
        const double s2_L = w_L * w_L * (1.0/L_pt3 + 1.0/L_pt2);
        const double s2_R = w_R * w_R * (1.0/R_pt3 + 1.0/R_pt2);
        const double s2   = s2_L + s2_R;
        if (s2 <= 0) continue;
        const double z2 = (w_L - w_R) * (w_L - w_R) / s2;
        if (z2 > best_z2) { best_z2 = z2; best_k = k; }
    }
    return best_k;
}

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
// covering the chain END-TO-END WITHOUT GAPS. Each event in the chain
// belongs to exactly one segment.
//
// Tiling rules:
//   * Segment 0 always starts at (files[0], event 0). Any non-trigger
//     events before the first PT3/PT2 trigger belong to segment 0.
//   * Boundary between segments N and N+1: just before segment N+1's
//     first trigger event. Concretely:
//        - if segs[N+1].first_local > 0:
//             segment N ends at (segs[N+1].first_file,
//                                segs[N+1].first_local - 1).
//             segment N may extend INTO segs[N+1].first_file (covering
//             non-trigger events at the start of that file).
//        - if segs[N+1].first_local == 0:
//             segment N ends at the LAST event of the file immediately
//             before segs[N+1].first_file.
//   * Last segment always ends at (files.back(), n_events - 1).
//
// Without these rules: segment 0 starts wherever the first trigger sits
// (loses leading events), and segment N stops at the end of the file
// containing its last trigger (loses non-trigger events at the start of
// the next file before segment N+1's first trigger). Either creates
// holes in coverage that the user correctly flagged as bugs.
std::vector<OutRange> expandSegments(const std::vector<ClosedSegment>& segs,
                                     const std::vector<FileInfo>& files) {
    std::vector<OutRange> out;
    if (segs.empty() || files.empty()) return out;

    for (size_t s = 0; s < segs.size(); ++s) {
        const ClosedSegment& seg = segs[s];

        // -- Start coordinates -------------------------------------------
        int      fi_start;
        Long64_t start_local;
        if (s == 0) {
            fi_start    = 0;
            start_local = 0;
        } else {
            fi_start = findFileIdx(files, seg.first_file);
            if (fi_start < 0) {
                std::cerr << "  WARNING: segment first_file not in files: "
                          << seg.first_file << "\n";
                continue;
            }
            start_local = seg.first_local;
        }

        // -- End coordinates --------------------------------------------
        int      fi_end;
        Long64_t end_local;
        if (s + 1 < segs.size()) {
            const ClosedSegment& nxt = segs[s + 1];
            const int fi_nxt = findFileIdx(files, nxt.first_file);
            if (fi_nxt < 0) {
                std::cerr << "  WARNING: next segment first_file not in files: "
                          << nxt.first_file << "\n";
                continue;
            }
            if (nxt.first_local == 0) {
                // Boundary at end of the file before nxt.first_file.
                if (fi_nxt == 0) {
                    std::cerr << "  WARNING: segment "
                              << (s + 1) << " starts at file 0 local 0 — "
                              << "would leave segment " << s << " empty.\n";
                    continue;
                }
                fi_end    = fi_nxt - 1;
                end_local = files[fi_end].n_events - 1;
            } else {
                fi_end    = fi_nxt;
                end_local = nxt.first_local - 1;
            }
        } else {
            fi_end    = (int) files.size() - 1;
            end_local = files[fi_end].n_events - 1;
        }

        // -- Emit one OutRange per file in [fi_start..fi_end] -----------
        for (int fi = fi_start; fi <= fi_end; ++fi) {
            OutRange r;
            r.seg_idx   = seg.seg_idx;
            r.file_path = files[fi].path;
            r.event_lo  = (fi == fi_start) ? start_local : 0;
            r.event_hi  = (fi == fi_end)   ? end_local
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

    TTree* t_files = dynamic_cast<TTree*>(fin->Get("trigger_cal_files"));
    TTree* t_evts  = dynamic_cast<TTree*>(fin->Get("trigger_cal_nt"));
    if (!t_files || !t_evts) {
        std::cerr << "  WARNING: missing trigger_cal_nt or trigger_cal_files in "
                  << in_path << " — skipping\n";
        fin->Close(); delete fin;
        return;
    }

    // Load files metadata. file_idx in trigger_cal_files is local to each
    // run_parallel chunk and after hadd often degenerates to all-zero (1
    // file per chunk), so we use ENTRY ORDER instead — split is line-aligned
    // and hadd merges chunks lexically, which means trigger_cal_files entry
    // order matches the original .list order.
    std::vector<FileInfo> files;
    {
        std::string  f_path;
        std::string* p_path = &f_path;
        Long64_t f_n = 0;
        t_files->SetBranchAddress("file_path",      &p_path);
        t_files->SetBranchAddress("n_events_total", &f_n);
        const Long64_t N = t_files->GetEntries();
        files.reserve(N);
        for (Long64_t i = 0; i < N; ++i) {
            t_files->GetEntry(i);
            files.push_back({f_path, f_n});
        }
    }
    std::cout << "  files: " << files.size()
              << "  events in trigger_cal_nt: " << t_evts->GetEntries() << "\n";

    if (t_evts->GetEntries() == 0) {
        std::cerr << "  WARNING: trigger_cal_nt is empty — channel skipped\n";
        fin->Close(); delete fin;
        return;
    }

    // Stream through trigger_cal_nt. All branches are Float_t (TNtuple).
    Float_t f_file_idx = 0, f_local_idx = 0, f_trigbit = 0;
    Float_t f_oa = 0, f_isBest = 0, f_vz = 0, f_si = 0;
    t_evts->SetBranchAddress("file_idx",        &f_file_idx);
    t_evts->SetBranchAddress("local_event_idx", &f_local_idx);
    t_evts->SetBranchAddress("trigbit",         &f_trigbit);
    t_evts->SetBranchAddress("oa",              &f_oa);
    t_evts->SetBranchAddress("isBest",          &f_isBest);
    t_evts->SetBranchAddress("eVertReco_z",     &f_vz);
    t_evts->SetBranchAddress("start_iteration", &f_si);

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

    // Two-phase streaming with cuts applied INLINE per event. Cuts mirror
    // src/setup_cuts.h on pp45_epem (isBest, vertex_z, start_detector).
    // OA is NOT cut here — see header docstring for the rationale.
    // Events failing any cut — or carrying a trigbit other than PT3 / PT2
    // (e.g. trigbit == 12288 = PT3+PT2 simultaneously) — are skipped
    // entirely: they don't update the segmenter pool and don't trigger the
    // z-test. Tile-end-to-end of segments still holds because the
    // segmenter only records WHERE the first kept trigger of each segment
    // sits; cut events fall inside whatever segment was open at their
    // chain position.
    //
    // File boundary in trigger_cal_nt is detected via local_event_idx
    // resets (it's monotonic non-decreasing within a file, drops to 0 at
    // every file change — including chunk boundaries after hadd, since
    // each chunk's events were written in chain order).
    const Long64_t N_evts = t_evts->GetEntries();
    int      current_file_idx = 0;
    Long64_t prev_local_idx   = -1;
    Long64_t n_pass_cuts      = 0;
    Long64_t n_pass_trigbit   = 0;

    for (Long64_t i = 0; i < N_evts; ++i) {
        t_evts->GetEntry(i);

        const Long64_t li = static_cast<Long64_t>(f_local_idx);
        if (i > 0 && li < prev_local_idx) ++current_file_idx;
        prev_local_idx = li;

        // Quality cuts (mirror src/setup_cuts.h on pp45_epem).
        if (f_isBest != 1.0f)  continue;
        if (f_vz   <= -500.0f) continue;
        if (f_si   != 3.0f)    continue;
        // OA cut intentionally NOT applied — segmentation edges are
        // selection-independent (see header docstring).
        ++n_pass_cuts;

        // Only single-trigger PT3 or PT2 events feed the calibration.
        // Multi-trigger events (e.g. trigbit == 12288 = PT3 + PT2) are
        // ignored, matching the existing main.cc:687 convention.
        const int tb = static_cast<int>(f_trigbit);
        if (tb != trig_pt3 && tb != trig_pt2) continue;
        ++n_pass_trigbit;

        if (current_file_idx >= static_cast<int>(files.size())) continue;
        const std::string& path = files[current_file_idx].path;

        if (relErr(open_seg.n_pt3, open_seg.n_pt2) > max_rel_err) {
            // Phase 1: keep growing open until it's well-determined.
            // Open never needs per-event history (refinement scans only
            // the candidate's events), so track=false here.
            open_seg.addTrigger(path, li, tb, trig_pt3, trig_pt2,
                                current_file_idx, /*track=*/false);
            continue;
        }

        // Phase 2: open is stable. Build up candidate. Track each event
        // (file_idx + local_idx + PT3/PT2 flag) so we can refine the
        // boundary later by scanning within the candidate.
        cand.addTrigger(path, li, tb, trig_pt3, trig_pt2,
                        current_file_idx, /*track=*/true);
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
            // SPLIT with boundary refinement. The crude algorithm puts
            // the new segment start at cand.events[0]; the true change
            // point can be anywhere inside cand's span. findOptimalSplit
            // returns the index k* maximising the left/right z²; events
            // [0..k*] are merged into the closing segment, [k*+1..end]
            // become the new open segment with a refined start position.
            const size_t k = findOptimalSplit(open_seg, cand);

            long long abs_pt3 = 0, abs_pt2 = 0;
            for (size_t i = 0; i <= k; ++i) {
                if (cand.events[i].is_pt3) ++abs_pt3; else ++abs_pt2;
            }
            // Extend open with cand[0..k] then close it.
            open_seg.n_pt3      += abs_pt3;
            open_seg.n_pt2      += abs_pt2;
            open_seg.last_file   = files[cand.events[k].file_idx].path;
            open_seg.last_local  = cand.events[k].local_idx;
            closeOpen("SPLIT |z|=" + std::to_string(std::abs(z))
                      + " refined@" + std::to_string(k));

            // Build new open from cand[k+1..end].
            PoolState new_open;
            new_open.first_file  = files[cand.events[k+1].file_idx].path;
            new_open.first_local = cand.events[k+1].local_idx;
            new_open.last_file   = cand.last_file;
            new_open.last_local  = cand.last_local;
            for (size_t i = k + 1; i < cand.events.size(); ++i) {
                if (cand.events[i].is_pt3) ++new_open.n_pt3;
                else                       ++new_open.n_pt2;
            }
            // Don't carry events into the new open — the next candidate
            // (in Phase 2 again) will collect its own per-event history.
            open_seg = std::move(new_open);
            cand.reset();
        } else {
            // MERGE — absorb candidate into open (counts only; events
            // discarded since open never needs per-event history).
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

    std::cout << "  -> events: total=" << N_evts
              << "  pass_cuts=" << n_pass_cuts
              << "  PT3+PT2=" << n_pass_trigbit
              << "  (acceptance: "
              << std::fixed << std::setprecision(2)
              << (100.0 * n_pass_trigbit / std::max<Long64_t>(1, N_evts)) << "%)\n";
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
                         double max_rel_err = 0.01) {
    std::cout << "=== trigger_calibration"
              << "  z_threshold=" << z_threshold
              << "  max_rel_err=" << max_rel_err
              << " ===\n";

    // Paths are relative to the cwd from which root is invoked. Standard
    // usage is `cd research/calibration && root -l -b -q plots/...` so the
    // cal-mode ./ana outputs live in repo root (../../output_*_cal.root)
    // and the calibration ROOTs go to the repo root too.
    const std::vector<std::string> chans = {"epem", "epep", "emem"};
    for (const std::string& chan : chans) {
        const std::string in_path  = "../../output_" + chan + "_cal.root";
        const std::string out_path = "../../pt3_calibration_" + chan + ".root";
        processChannel(chan, in_path, out_path, z_threshold, max_rel_err);
    }
    std::cout << "\nDone.\n";
}

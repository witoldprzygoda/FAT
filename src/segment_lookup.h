/**
 * @file segment_lookup.h
 * @brief Per-segment PT3 trigger-bias weight lookup
 *
 * Loads pt3_calibration_<channel>.root (produced by research/calibration/)
 * and provides O(log n) lookup of (seg_idx, w_seg) for a (file_path, local
 * event index) pair. Used by main.cc to stamp each dilepton_nt row with the
 * trigger-bias weight that applies to its parent event.
 *
 * Schema of the input TTree "pt3_calibration":
 *   seg_idx    /I   segment index (global)
 *   file_path  /string   absolute path of the source ROOT file
 *   event_lo   /L   inclusive lower bound of local-event-index range
 *   event_hi   /L   inclusive upper bound of local-event-index range
 *   w          /D   per-segment PT3-bias weight
 *   sigma      /D   weight uncertainty
 *   n_pt3      /L   PT3 event count in segment (diagnostic)
 *   n_pt2      /L   PT2 event count in segment (diagnostic)
 *
 * File-path matching: keys are normalized to basename only (after the last
 * '/'). This is defensive against absolute-vs-relative path mismatches or
 * directory moves between the calibration run and the consumption run.
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2026
 */

#ifndef SEGMENT_LOOKUP_H
#define SEGMENT_LOOKUP_H

#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

class SegmentLookup {
public:
    struct SegRange {
        Long64_t lo;
        Long64_t hi;
        int      seg_idx;
        double   w;
    };

    // Load segments from pt3_calibration_<channel>.root.
    // Returns true on success; prints diagnostics to stdout/stderr.
    bool Load(const std::string& path) {
        per_file_.clear();
        n_loaded_ = 0;
        cal_path_ = path;

        if (path.empty()) return false;

        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie()) {
            std::cerr << "SegmentLookup: cannot open calibration file '"
                      << path << "'\n";
            return false;
        }

        auto* t = dynamic_cast<TTree*>(f->Get("pt3_calibration"));
        if (!t) {
            std::cerr << "SegmentLookup: TTree 'pt3_calibration' missing in '"
                      << path << "'\n";
            return false;
        }

        Int_t       seg_idx = 0;
        std::string file_path;
        std::string* file_path_ptr = &file_path;
        Long64_t    event_lo = 0;
        Long64_t    event_hi = 0;
        Double_t    w        = 1.0;

        t->SetBranchAddress("seg_idx",   &seg_idx);
        t->SetBranchAddress("file_path", &file_path_ptr);
        t->SetBranchAddress("event_lo",  &event_lo);
        t->SetBranchAddress("event_hi",  &event_hi);
        t->SetBranchAddress("w",         &w);

        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            const std::string key = basenameOf(file_path);
            per_file_[key].push_back({event_lo, event_hi, seg_idx, w});
            ++n_loaded_;
        }

        // Sort each file's segments by lo so we can binary-search.
        for (auto& kv : per_file_) {
            auto& v = kv.second;
            std::sort(v.begin(), v.end(),
                      [](const SegRange& a, const SegRange& b) {
                          return a.lo < b.lo;
                      });
        }

        std::cout << "SegmentLookup: loaded " << n_loaded_
                  << " segments from " << path
                  << " (covering " << per_file_.size() << " files)\n";
        return true;
    }

    // Look up segment for (file_path, local_idx).
    // Matching is on basename only (defensive). Returns {-1, 1.0} on miss.
    std::pair<int, double> Lookup(const std::string& file_path,
                                  Long64_t local_idx) const {
        if (per_file_.empty()) return {-1, 1.0};
        const std::string key = basenameOf(file_path);
        auto it = per_file_.find(key);
        if (it == per_file_.end()) {
            ++miss_file_;
            return {-1, 1.0};
        }
        const auto& v = it->second;
        // Linear scan is fine here: segments per file are O(1)–O(10).
        for (const auto& s : v) {
            if (local_idx >= s.lo && local_idx <= s.hi) {
                ++hit_;
                return {s.seg_idx, s.w};
            }
        }
        ++miss_range_;
        return {-1, 1.0};
    }

    // Diagnostics
    bool   isLoaded()      const { return !per_file_.empty(); }
    size_t segmentCount()  const { return n_loaded_; }
    size_t fileCount()     const { return per_file_.size(); }
    const std::string& path() const { return cal_path_; }

    Long64_t hits()        const { return hit_; }
    Long64_t missesFile()  const { return miss_file_; }
    Long64_t missesRange() const { return miss_range_; }

    void printSummary(std::ostream& os = std::cout) const {
        os << "SegmentLookup summary:\n";
        os << "  source           : " << cal_path_  << "\n";
        os << "  segments loaded  : " << n_loaded_  << "\n";
        os << "  files covered    : " << per_file_.size() << "\n";
        os << "  lookups (hits)   : " << hit_       << "\n";
        os << "  miss (no file)   : " << miss_file_ << "\n";
        os << "  miss (no range)  : " << miss_range_ << "\n";
    }

private:
    static std::string basenameOf(const std::string& p) {
        auto pos = p.find_last_of('/');
        return (pos == std::string::npos) ? p : p.substr(pos + 1);
    }

    std::map<std::string, std::vector<SegRange>> per_file_;
    size_t      n_loaded_   = 0;
    std::string cal_path_;
    mutable Long64_t hit_        = 0;
    mutable Long64_t miss_file_  = 0;
    mutable Long64_t miss_range_ = 0;
};

#endif // SEGMENT_LOOKUP_H

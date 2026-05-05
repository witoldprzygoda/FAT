// meson_research.cc — slice-based research on meson_dalitz_nt.
//
// Reads one channel of FAT output (output_*_exp.root) and produces TWO sets
// of 1D histograms of M(e+e-gamma) per opening-angle slice:
//   m_epemg_*       — RECONSTRUCTED (raw HADES tracking)
//   m_epemg_cor_*   — CORRECTED     (energy-loss corrected leptons)
//
// The same binary is run separately on the three CB channels (epem, epep,
// emem) via three different config files; downstream macros in plots/
// combine them into all/CB/signal triples per slice for either flavour.
//
// Slicing config — currently compile-time constants. To change, edit the
// SliceConfig namespace below and rebuild.
//
// Usage (from research/ — paths in JSON configs are relative to here):
//   ./meson_research config.json        # epem
//   ./meson_research config_epep.json   # ++
//   ./meson_research config_emem.json   # --
//
// JSON config keys (flat schema):
//   "input_file"   — FAT output ROOT file to read (e.g. output_epem_exp.root)
//   "ntuple_name"  — ntuple to read from (typically "meson_dalitz_nt")
//   "output_file"  — output ROOT file path (e.g. research_epem.root, relative to research/)
//   "extra_cut"    — optional TTree::Draw cut, ANDed with the slice cut
//                    (empty string = no extra cut)
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2025

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cctype>

// -----------------------------------------------------------------------------
// Slicing parameters — kept simple as compile-time constants.
// -----------------------------------------------------------------------------
namespace SliceConfig {
    // OA slice grids — two co-existing slicings on the same range:
    //   step 0.2° → fine (used by π⁰ Dalitz analysis, high statistics)
    //   step 0.5° → coarse (used by η Dalitz analysis, lower statistics)
    // Histogram names stay distinct because the slice edges differ
    // (e.g. m_epemg_oa_0p0_0p2 vs m_epemg_oa_0p0_0p5).
    constexpr double kSliceMin     = 0.0;
    constexpr double kSliceMax     = 15.0;
    constexpr double kSliceStep    = 0.2;
    constexpr double kSliceStepEta = 0.5;

    // Per-slice m_epemg histogram binning.
    constexpr int    kHistNBins = 160;      // 5 MeV/bin
    constexpr double kHistMin   = 0.0;
    constexpr double kHistMax   = 0.8;

    // Slice variable in the input ntuple.
    constexpr const char* kSliceVar = "oa_epem";

    // Plotted mass flavours: same hist binning, separate output histograms.
    // The 'name' is both the input branch name and the output hist prefix.
    struct MassFlavor {
        const char* name;   // ntuple branch + hist-name prefix
        const char* tag;    // human-readable label for log lines
    };
    constexpr MassFlavor kFlavors[] = {
        {"m_epemg",     "REC"},
        {"m_epemg_cor", "COR"},
    };
    constexpr int kNFlavors = sizeof(kFlavors) / sizeof(kFlavors[0]);
}

// -----------------------------------------------------------------------------
// Minimal JSON string-value reader for our flat schema. Looks for
//   "key" : "value"
// with whitespace tolerated. Returns empty string if key missing or value
// not a string. Sufficient for our 4-key configs; do not extend.
// -----------------------------------------------------------------------------
std::string getJsonString(const std::string& content, const std::string& key) {
    const std::string pat = "\"" + key + "\"";
    size_t pos = 0;
    while ((pos = content.find(pat, pos)) != std::string::npos) {
        size_t after = pos + pat.size();
        while (after < content.size() && std::isspace((unsigned char)content[after])) ++after;
        if (after < content.size() && content[after] == ':') {
            ++after;
            while (after < content.size() && std::isspace((unsigned char)content[after])) ++after;
            if (after < content.size() && content[after] == '"') {
                const size_t end = content.find('"', after + 1);
                if (end != std::string::npos)
                    return content.substr(after + 1, end - after - 1);
            }
            return "";
        }
        pos = after;
    }
    return "";
}

// Format slice edge "0.4" → "0p4", "15.0" → "15p0" — used to build histogram
// names that are valid C++ / ROOT identifiers (no dots).
std::string fmtEdge(double x) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%.1f", x);
    std::string s(buf);
    for (auto& c : s) if (c == '.') c = 'p';
    return s;
}

int main(int argc, char* argv[]) {
    std::string config_path = "config.json";
    if (argc > 1) config_path = argv[1];

    // -- Load config --------------------------------------------------------
    std::ifstream ifs(config_path);
    if (!ifs) {
        std::cerr << "ERROR: cannot open config file: " << config_path << "\n";
        return 1;
    }
    std::stringstream buf;
    buf << ifs.rdbuf();
    const std::string content = buf.str();

    const std::string input_file  = getJsonString(content, "input_file");
    const std::string ntuple_name = getJsonString(content, "ntuple_name");
    const std::string output_file = getJsonString(content, "output_file");
    const std::string extra_cut   = getJsonString(content, "extra_cut");

    if (input_file.empty() || ntuple_name.empty() || output_file.empty()) {
        std::cerr << "ERROR: config must define non-empty input_file, ntuple_name, output_file\n";
        return 1;
    }

    std::cout << "=== meson_research ===\n";
    std::cout << "  input:    " << input_file  << "\n";
    std::cout << "  ntuple:   " << ntuple_name << "\n";
    std::cout << "  output:   " << output_file << "\n";
    std::cout << "  cut:      '" << extra_cut << "'\n";

    using namespace SliceConfig;
    const int n_slices_pi0 = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStep));
    const int n_slices_eta = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStepEta));
    std::cout << "  slicing:  '" << kSliceVar << "' in [" << kSliceMin << ", "
              << kSliceMax << "] deg\n"
              << "      fine   step " << kSliceStep    << " deg → "
              << n_slices_pi0 << " slices  (π⁰)\n"
              << "      coarse step " << kSliceStepEta << " deg → "
              << n_slices_eta << " slices  (η)\n";
    std::cout << "  flavours:";
    for (const auto& fv : kFlavors) std::cout << "  " << fv.name << "(" << fv.tag << ")";
    std::cout << "\n";
    std::cout << "  per slice: each flavour in [" << kHistMin << ", " << kHistMax
              << "] GeV/c² with " << kHistNBins << " bins\n\n";

    // -- Open input ---------------------------------------------------------
    TFile* fin = TFile::Open(input_file.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open input file: " << input_file << "\n";
        return 1;
    }
    auto* t = (TTree*)fin->Get(ntuple_name.c_str());
    if (!t) {
        std::cerr << "ERROR: ntuple '" << ntuple_name << "' not in " << input_file << "\n";
        fin->Close();
        return 1;
    }
    std::cout << "  input ntuple entries: " << t->GetEntries() << "\n\n";

    // -- Open output --------------------------------------------------------
    TFile* fout = TFile::Open(output_file.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot create output file: " << output_file << "\n";
        fin->Close();
        return 1;
    }
    fout->cd();   // histograms below auto-attach here

    // -- Full-range histogram (per flavour) --------------------------------
    // m_epemg integrated over OA ∈ [kSliceMin, kSliceMax] — i.e. mathematically
    // the sum of all slice histograms below. Kept as a separate, named
    // histogram so plotting macros can read it directly without having to
    // accumulate slices. Built for each mass flavour (REC, COR).
    for (const auto& fv : kFlavors) {
        const std::string fname = std::string(fv.name) + "_full";
        const TString ftitle = TString::Format(
            "M(e^{+}e^{-}#gamma) %s, OA(e^{+}e^{-}) #in [%.1f, %.1f] deg (full);"
            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
            fv.tag, kSliceMin, kSliceMax);

        TH1D* h_full = new TH1D(fname.c_str(), ftitle, kHistNBins, kHistMin, kHistMax);
        h_full->Sumw2();

        std::stringstream cut;
        cut << kSliceVar << ">=" << kSliceMin << " && " << kSliceVar << "<" << kSliceMax;
        if (!extra_cut.empty()) cut << " && (" << extra_cut << ")";

        const std::string draw_expr = std::string(fv.name) + ">>" + fname;
        t->Draw(draw_expr.c_str(), cut.str().c_str(), "goff");

        std::cout << "  " << fname << ": " << h_full->GetEntries() << " entries (full)\n";
    }

    // -- Loop slices --------------------------------------------------------
    Long64_t total = 0;

    auto runSlicing = [&](double step, const char* tag) {
        const int n = static_cast<int>(std::round((kSliceMax - kSliceMin) / step));
        std::cout << "\n  [slicing " << tag << "  step=" << step
                  << " deg → " << n << " slices]\n";
        for (int i = 0; i < n; ++i) {
            const double oa_lo = kSliceMin + i * step;
            const double oa_hi = kSliceMin + (i + 1) * step;
            const std::string suff = "oa_" + fmtEdge(oa_lo) + "_" + fmtEdge(oa_hi);

            std::stringstream cut;
            cut << kSliceVar << ">=" << oa_lo << " && " << kSliceVar << "<" << oa_hi;
            if (!extra_cut.empty()) cut << " && (" << extra_cut << ")";

            for (const auto& fv : kFlavors) {
                const std::string hname  = std::string(fv.name) + "_" + suff;
                const TString     title  = TString::Format(
                    "M(e^{+}e^{-}#gamma) %s, OA(e^{+}e^{-}) #in [%.1f, %.1f] deg;"
                    "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts",
                    fv.tag, oa_lo, oa_hi);

                TH1D* h = new TH1D(hname.c_str(), title, kHistNBins, kHistMin, kHistMax);
                h->Sumw2();

                const std::string draw_expr = std::string(fv.name) + ">>" + hname;
                t->Draw(draw_expr.c_str(), cut.str().c_str(), "goff");

                total += static_cast<Long64_t>(h->GetEntries());
                std::cout << "  " << hname << ": " << h->GetEntries() << " entries\n";
            }
        }
    };

    runSlicing(kSliceStep,    "fine/π⁰");
    runSlicing(kSliceStepEta, "coarse/η");

    std::cout << "\n  Total fills across slices × flavours: " << total << "\n";

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Output: " << output_file << "\n";
    return 0;
}

// meson_research.cc — slice-based research on meson_dalitz_nt (SIMULATION).
//
// Reads a single FAT-sim output (output_epem_sim.root) and produces a set
// of 1D histograms of M(e+e-gamma), one per opening-angle slice.
//
// Three histograms per OA slice (REC / COR / TRU flavours of the e+e-γ mass):
//   m_epemg_oa_X_Y      — RECONSTRUCTED  (raw HADES tracking)
//   m_epemg_cor_oa_X_Y  — CORRECTED      (energy-loss corrected leptons)
//   m_epemg_tru_oa_X_Y  — SIMULATED-truth (sim_px/py/pz directly)
// All weighted by sim_genweight. Plus the integrated full-range counterparts:
//   m_epemg_full, m_epemg_cor_full, m_epemg_tru_full
//
// On simulation there is no like-sign CB — the binary produces a single
// channel only (one config file, one output ROOT).
//
// All histograms are filled in a SINGLE event loop with SetBranchAddress —
// 102 separate TTree::Draw calls on a 100M-row tree would be impractically
// slow. Optional extra_cut from the config is evaluated via TTreeFormula
// inside the loop.
//
// Usage (from research/ — paths in JSON config are relative to here):
//   ./meson_research config.json
//
// JSON config keys (flat schema):
//   "input_file"   — FAT sim output ROOT file (e.g. ../output_epem_sim.root)
//   "ntuple_name"  — ntuple to read from (typically "meson_dalitz_nt")
//   "output_file"  — output ROOT file path (e.g. outputs/research_sim.root)
//   "extra_cut"    — optional TTreeFormula expression, ANDed with the slice
//                    cut. sim_genweight is always applied as a fill weight.
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2025

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TTreeFormula.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <cctype>

// -----------------------------------------------------------------------------
// Slicing parameters — kept simple as compile-time constants.
// -----------------------------------------------------------------------------
namespace SliceConfig {
    constexpr double kSliceMin  = 0.0;
    constexpr double kSliceMax  = 10.0;
    constexpr double kSliceStep = 0.2;

    constexpr int    kHistNBins = 230;      // 2 MeV/bin (230 × 0.002 = 0.46)
    constexpr double kHistMin   = 0.0;
    constexpr double kHistMax   = 0.46;
}

// -----------------------------------------------------------------------------
// Minimal JSON string-value reader for our flat schema.
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

// Format slice edge "0.4" → "0p4", "10.0" → "10p0".
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

    using namespace SliceConfig;
    const int n_slices = static_cast<int>(std::round((kSliceMax - kSliceMin) / kSliceStep));

    std::cout << "=== meson_research (sim) ===\n";
    std::cout << "  input:    " << input_file  << "\n";
    std::cout << "  ntuple:   " << ntuple_name << "\n";
    std::cout << "  output:   " << output_file << "\n";
    std::cout << "  cut:      '" << extra_cut << "'   (× sim_genweight always)\n";
    std::cout << "  slicing:  oa_epem in [" << kSliceMin << ", " << kSliceMax
              << "] deg, step " << kSliceStep << " deg → "
              << n_slices << " slices\n";
    std::cout << "  per slice: REC m_epemg + COR m_epemg_cor + TRU m_epemg_sim, "
              << kHistNBins << " bins in ["
              << kHistMin << ", " << kHistMax << "] GeV/c²\n\n";

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
    const Long64_t n_entries = t->GetEntries();
    std::cout << "  input ntuple entries: " << n_entries << "\n\n";

    // -- Open output --------------------------------------------------------
    TFile* fout = TFile::Open(output_file.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot create output file: " << output_file << "\n";
        fin->Close();
        return 1;
    }
    fout->cd();

    // -- Pre-allocate histograms (REC + SIM-truth × full + per-slice) -------
    auto makeHist = [&](const std::string& name, const std::string& title) {
        auto* h = new TH1D(name.c_str(), title.c_str(),
                           kHistNBins, kHistMin, kHistMax);
        h->Sumw2();
        return h;
    };

    TH1D* h_full_rec = makeHist("m_epemg_full",
        TString::Format("M(e^{+}e^{-}#gamma) REC, OA #in [%.1f, %.1f] deg (full);"
                        "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                        kSliceMin, kSliceMax).Data());
    TH1D* h_full_cor = makeHist("m_epemg_cor_full",
        TString::Format("M(e^{+}e^{-}#gamma) COR, OA #in [%.1f, %.1f] deg (full);"
                        "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                        kSliceMin, kSliceMax).Data());
    TH1D* h_full_tru = makeHist("m_epemg_tru_full",
        TString::Format("M(e^{+}e^{-}#gamma) TRU, OA #in [%.1f, %.1f] deg (full);"
                        "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                        kSliceMin, kSliceMax).Data());

    std::vector<TH1D*> h_slice_rec(n_slices, nullptr);
    std::vector<TH1D*> h_slice_cor(n_slices, nullptr);
    std::vector<TH1D*> h_slice_tru(n_slices, nullptr);
    for (int i = 0; i < n_slices; ++i) {
        const double oa_lo = kSliceMin + i * kSliceStep;
        const double oa_hi = kSliceMin + (i + 1) * kSliceStep;
        const std::string suff = "oa_" + fmtEdge(oa_lo) + "_" + fmtEdge(oa_hi);

        h_slice_rec[i] = makeHist("m_epemg_" + suff,
            TString::Format("M(e^{+}e^{-}#gamma) REC, OA #in [%.1f, %.1f] deg;"
                            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                            oa_lo, oa_hi).Data());
        h_slice_cor[i] = makeHist("m_epemg_cor_" + suff,
            TString::Format("M(e^{+}e^{-}#gamma) COR, OA #in [%.1f, %.1f] deg;"
                            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                            oa_lo, oa_hi).Data());
        h_slice_tru[i] = makeHist("m_epemg_tru_" + suff,
            TString::Format("M(e^{+}e^{-}#gamma) TRU, OA #in [%.1f, %.1f] deg;"
                            "M_{e^{+}e^{-}#gamma} [GeV/c^{2}];Counts (sim_genweight)",
                            oa_lo, oa_hi).Data());
    }

    // -- Branch addresses ---------------------------------------------------
    // Note: the TRU flavour is stored in the input ntuple under the legacy
    // branch name 'm_epemg_sim'; we just rebrand it as TRU on the output.
    Float_t oa_epem = 0, m_epemg = 0, m_epemg_cor = 0, m_epemg_sim = 0,
            sim_genweight = 0;
    t->SetBranchAddress("oa_epem",       &oa_epem);
    t->SetBranchAddress("m_epemg",       &m_epemg);
    t->SetBranchAddress("m_epemg_cor",   &m_epemg_cor);
    t->SetBranchAddress("m_epemg_sim",   &m_epemg_sim);
    t->SetBranchAddress("sim_genweight", &sim_genweight);

    // Optional extra cut as a TTreeFormula. nullptr if no cut configured.
    std::unique_ptr<TTreeFormula> formula;
    if (!extra_cut.empty()) {
        formula.reset(new TTreeFormula("extra_cut", extra_cut.c_str(), t));
        if (!formula->GetTree()) {
            std::cerr << "ERROR: invalid extra_cut: " << extra_cut << "\n";
            fout->Close(); fin->Close();
            return 1;
        }
    }

    // -- Single-pass event loop --------------------------------------------
    Long64_t n_processed = 0, n_filled = 0;
    Long64_t print_every = std::max<Long64_t>(1, n_entries / 20);

    for (Long64_t ev = 0; ev < n_entries; ++ev) {
        t->GetEntry(ev);
        ++n_processed;

        if (formula && formula->EvalInstance() == 0) continue;

        // Determine slice index from OA. Reject if outside grid.
        if (oa_epem < kSliceMin || oa_epem >= kSliceMax) {
            // Not in slicing range → also skip the "full" histogram (which
            // is defined as the OA-window integration matching the slices).
            continue;
        }
        const int idx = static_cast<int>((oa_epem - kSliceMin) / kSliceStep);
        if (idx < 0 || idx >= n_slices) continue;

        h_full_rec ->Fill(m_epemg,     sim_genweight);
        h_full_cor ->Fill(m_epemg_cor, sim_genweight);
        h_full_tru ->Fill(m_epemg_sim, sim_genweight);
        h_slice_rec[idx]->Fill(m_epemg,     sim_genweight);
        h_slice_cor[idx]->Fill(m_epemg_cor, sim_genweight);
        h_slice_tru[idx]->Fill(m_epemg_sim, sim_genweight);
        ++n_filled;

        if (n_processed % print_every == 0) {
            std::cout << "  processed " << n_processed << " / " << n_entries
                      << "  (filled " << n_filled << ")\n";
        }
    }

    std::cout << "\n  total processed: " << n_processed
              << "   filled: " << n_filled
              << "   in-window: " << h_full_rec->GetEntries() << "\n";
    std::cout << "  full-range integrals: REC=" << h_full_rec->Integral()
              << "   COR=" << h_full_cor->Integral()
              << "   TRU=" << h_full_tru->Integral() << "\n";

    // -- Per-slice summary (entries + sim_genweight integral) -------------
    // The single-pass loop above interleaves slices, so we can't print a
    // per-slice line "as we fill" the way a Draw-per-slice version would.
    // Print the equivalent diagnostic at the end instead.
    std::cout << "\n  Per-slice summary:\n";
    for (int i = 0; i < n_slices; ++i) {
        const double oa_lo = kSliceMin + i * kSliceStep;
        const double oa_hi = kSliceMin + (i + 1) * kSliceStep;
        std::cout << "    OA [" << oa_lo << ", " << oa_hi << "] deg:"
                  << "  REC ∫=" << h_slice_rec[i]->Integral()
                  << "  |  COR ∫=" << h_slice_cor[i]->Integral()
                  << "  |  TRU ∫=" << h_slice_tru[i]->Integral() << "\n";
    }

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Output: " << output_file << "\n";
    return 0;
}

// research/calibration/trigger_scan.cc
//
// Standalone PT3 / PT2 trigger scanner. Reads the same .list file the FAT
// main analysis uses, opens each ROOT file independently, counts events with
// trigbit == PT3 (default 8192) and trigbit == PT2 (default 4096), and writes
// one entry per file to an output TTree:
//
//     pt3_perfile  (file_idx /I, n_pt3 /L, n_pt2 /L, file_path /string)
//
// Decoupled from ./ana so we can re-derive the trigger-bias calibration
// without re-running the heavy analysis. Designed to be parallelisable via
// run_parallel_scan.sh — that script splits the .list into chunks, launches
// one trigger_scan per chunk on a chunk-config, and hadd-merges the per-chunk
// outputs into the final scan ROOT.
//
// Usage (from research/calibration/):
//   make
//   ./trigger_scan config_epem.json
//
// Config schema (flat, JSON-string fields):
//   {
//     "input_source": "/path/to/source.list"   (or single .root)
//     "tree_name":    "EpEm_ID"                (TTree name in each file)
//     "trig_pt3":     "8192"                   (integer as string — uses our
//     "trig_pt2":     "4096"                    minimal getJsonString reader)
//     "output_file":  "trigger_scan_epem.root"
//   }
//
// @author Witold Przygoda (witold.przygoda@uj.edu.pl)
// @date 2026

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <cstdio>
#include <cstdlib>

namespace {

// Minimal JSON string-value reader, identical in spirit to meson_research.cc.
// Looks for  "key" : "value"  with whitespace tolerated. Sufficient for our
// flat schema; do not extend.
std::string getJsonString(const std::string& content, const std::string& key) {
    const std::string pat = "\"" + key + "\"";
    size_t pos = 0;
    while ((pos = content.find(pat, pos)) != std::string::npos) {
        size_t after = pos + pat.size();
        while (after < content.size() && std::isspace((unsigned char) content[after])) ++after;
        if (after < content.size() && content[after] == ':') {
            ++after;
            while (after < content.size() && std::isspace((unsigned char) content[after])) ++after;
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

// .list parser — accepts plain paths, chain->Add("...") format, blank lines,
// and # / // comments. Same logic as NTupleReader::openFromList.
std::vector<std::string> readListFile(const std::string& list_path) {
    std::vector<std::string> files;
    std::ifstream ifs(list_path);
    if (!ifs) return files;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        const size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        if (line[start] == '#' || line.compare(start, 2, "//") == 0) continue;

        std::string fp;
        const size_t q1 = line.find('"');
        if (q1 != std::string::npos) {
            const size_t q2 = line.find('"', q1 + 1);
            if (q2 != std::string::npos)
                fp = line.substr(q1 + 1, q2 - q1 - 1);
        } else {
            const size_t end = line.find_last_not_of(" \t;");
            if (end != std::string::npos)
                fp = line.substr(start, end - start + 1);
        }
        if (!fp.empty()) files.push_back(std::move(fp));
    }
    return files;
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
    std::string config_path = "config.json";
    if (argc > 1) config_path = argv[1];

    // -- Load config -------------------------------------------------------
    std::ifstream ifs(config_path);
    if (!ifs) {
        std::cerr << "ERROR: cannot open config file: " << config_path << "\n";
        return 1;
    }
    std::stringstream buf; buf << ifs.rdbuf();
    const std::string content = buf.str();

    const std::string input_source = getJsonString(content, "input_source");
    const std::string tree_name    = getJsonString(content, "tree_name");
    const std::string output_file  = getJsonString(content, "output_file");
    const std::string s_pt3        = getJsonString(content, "trig_pt3");
    const std::string s_pt2        = getJsonString(content, "trig_pt2");
    const int trig_pt3 = s_pt3.empty() ? 8192 : std::atoi(s_pt3.c_str());
    const int trig_pt2 = s_pt2.empty() ? 4096 : std::atoi(s_pt2.c_str());

    if (input_source.empty() || tree_name.empty() || output_file.empty()) {
        std::cerr << "ERROR: config must define non-empty "
                  << "input_source, tree_name, output_file\n";
        return 1;
    }

    std::cout << "=== trigger_scan ===\n"
              << "  input_source: " << input_source << "\n"
              << "  tree_name:    " << tree_name    << "\n"
              << "  trig_pt3:     " << trig_pt3     << "\n"
              << "  trig_pt2:     " << trig_pt2     << "\n"
              << "  output_file:  " << output_file  << "\n";

    // -- Resolve input list -----------------------------------------------
    std::vector<std::string> files;
    const bool is_list = (input_source.size() >= 5 &&
                          input_source.substr(input_source.size() - 5) == ".list");
    if (is_list) {
        files = readListFile(input_source);
        if (files.empty()) {
            std::cerr << "ERROR: input list is empty or unreadable: "
                      << input_source << "\n";
            return 1;
        }
    } else {
        files.push_back(input_source);
    }
    std::cout << "  files to scan: " << files.size() << "\n\n";

    // -- Open output ROOT --------------------------------------------------
    TFile* fout = TFile::Open(output_file.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot open output: " << output_file << "\n";
        if (fout) delete fout;
        return 1;
    }
    fout->cd();
    TTree* tout = new TTree("pt3_perfile", "Per-file PT3/PT2 trigger counts");
    Int_t       b_idx = 0;
    Long64_t    b_pt3 = 0, b_pt2 = 0;
    std::string b_path;
    tout->Branch("file_idx",  &b_idx,  "file_idx/I");
    tout->Branch("n_pt3",     &b_pt3,  "n_pt3/L");
    tout->Branch("n_pt2",     &b_pt2,  "n_pt2/L");
    tout->Branch("file_path", &b_path);

    // -- Scan files --------------------------------------------------------
    Long64_t total_pt3 = 0, total_pt2 = 0, total_evt = 0;
    int n_skipped = 0;

    for (size_t i = 0; i < files.size(); ++i) {
        const std::string& fpath = files[i];

        TFile* fin = TFile::Open(fpath.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            std::cerr << "  WARNING: cannot open " << fpath << " — skipping\n";
            if (fin) delete fin;
            ++n_skipped;
            continue;
        }
        TTree* tin = dynamic_cast<TTree*>(fin->Get(tree_name.c_str()));
        if (!tin) {
            std::cerr << "  WARNING: tree '" << tree_name << "' missing in "
                      << fpath << " — skipping\n";
            fin->Close(); delete fin;
            ++n_skipped;
            continue;
        }

        // The FAT input ntuples store trigbit as Float_t — match that.
        Float_t trigbit_f = 0.0f;
        tin->SetBranchStatus("*", 0);
        tin->SetBranchStatus("trigbit", 1);
        tin->SetBranchAddress("trigbit", &trigbit_f);

        Long64_t n_pt3 = 0, n_pt2 = 0;
        const Long64_t N = tin->GetEntries();
        for (Long64_t j = 0; j < N; ++j) {
            tin->GetEntry(j);
            const int trigbit = static_cast<int>(trigbit_f);
            if      (trigbit == trig_pt3) ++n_pt3;
            else if (trigbit == trig_pt2) ++n_pt2;
        }

        // Fill output entry. Buffers are bound to fout via the construction
        // context above; Fill() goes to fout regardless of current gFile.
        b_idx  = static_cast<Int_t>(i);
        b_pt3  = n_pt3;
        b_pt2  = n_pt2;
        b_path = fpath;
        tout->Fill();

        total_pt3 += n_pt3; total_pt2 += n_pt2; total_evt += N;

        // Print every 25th file (and the first / last) so logs don't explode
        // for big lists.
        if (i == 0 || i == files.size() - 1 || (i + 1) % 25 == 0) {
            std::printf("  [%4zu/%4zu] N=%-10lld PT3=%-8lld PT2=%-8lld  %s\n",
                        i + 1, files.size(),
                        (long long) N, (long long) n_pt3, (long long) n_pt2,
                        fpath.c_str());
        }

        fin->Close(); delete fin;
    }

    // -- Finalise output ---------------------------------------------------
    fout->cd();
    tout->Write();
    fout->Close();
    delete fout;

    std::cout << "\n  scanned " << files.size() - n_skipped
              << " / " << files.size() << " files"
              << "  (skipped " << n_skipped << ")\n"
              << "  totals: events=" << total_evt
              << "  N_PT3=" << total_pt3
              << "  N_PT2=" << total_pt2 << "\n";

    // Completion marker — used by monitor_parallel_scan.sh to detect "done".
    std::cout << "Scan Complete!\n";
    return 0;
}

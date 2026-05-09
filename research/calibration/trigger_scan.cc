// research/calibration/trigger_scan.cc
//
// Standalone PT3 / PT2 trigger scanner — phase (a) of the trigger-bias
// calibration pipeline. Reads the same .list file the FAT main analysis
// uses, opens each ROOT file independently, and emits ONE TTree entry per
// trigger event (PT3 or PT2). The downstream segmenter (phase (b)) decides
// where weight changes — there is NO a-priori chunking here.
//
// Output ROOT contains two TTrees:
//
//     trigger_events  (1 entry per PT3/PT2 event in the chain)
//        - file_path        /string   physical ROOT file the event lives in
//        - local_event_idx  /L        index inside that file
//        - trigbit          /I        the actual trigger code (PT3 or PT2)
//
//     files           (1 entry per file scanned)
//        - file_path        /string   absolute path of the file
//        - n_events_total   /L        total events in the file's tree
//        - n_pt3            /L        PT3 events found in this file
//        - n_pt2            /L        PT2 events found in this file
//
// The scanner is parallelisable — run_parallel_scan.sh splits the .list
// into chunks and runs one trigger_scan per chunk. After hadd, both TTrees
// concatenate cleanly: entries are in chunk-then-line order, which equals
// the original .list order (split is line-aligned, hadd merges lexically).
// The segmenter therefore sees trigger events in proper chain order and
// can detect file boundaries by file_path changes between consecutive
// entries.
//
// Usage (from research/calibration/):
//   make
//   ./trigger_scan config_epem.json     # or via run_parallel_scan.sh
//
// Config schema (flat):
//   {
//     "input_source": "/path/to/source.list",
//     "tree_name":    "EpEm_ID",
//     "trig_pt3":     "8192",
//     "trig_pt2":     "4096",
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

    // -- Resolve input list ------------------------------------------------
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

    TTree* t_events = new TTree("trigger_events",
                                "Per-PT3/PT2-event records");
    std::string b_path;
    Long64_t    b_local_idx = 0;
    Int_t       b_trigbit   = 0;
    t_events->Branch("file_path",       &b_path);
    t_events->Branch("local_event_idx", &b_local_idx, "local_event_idx/L");
    t_events->Branch("trigbit",         &b_trigbit,   "trigbit/I");

    TTree* t_files = new TTree("files",
                               "Per-file metadata for the scan");
    std::string  f_path;
    Long64_t     f_n_events = 0;
    Long64_t     f_n_pt3    = 0;
    Long64_t     f_n_pt2    = 0;
    t_files->Branch("file_path",      &f_path);
    t_files->Branch("n_events_total", &f_n_events, "n_events_total/L");
    t_files->Branch("n_pt3",          &f_n_pt3,    "n_pt3/L");
    t_files->Branch("n_pt2",          &f_n_pt2,    "n_pt2/L");

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

        Float_t trigbit_f = 0.0f;
        tin->SetBranchStatus("*", 0);
        tin->SetBranchStatus("trigbit", 1);
        tin->SetBranchAddress("trigbit", &trigbit_f);

        const Long64_t N = tin->GetEntries();
        Long64_t n_pt3_in_file = 0, n_pt2_in_file = 0;

        for (Long64_t j = 0; j < N; ++j) {
            tin->GetEntry(j);
            const int trigbit = static_cast<int>(trigbit_f);
            if (trigbit == trig_pt3 || trigbit == trig_pt2) {
                b_path     = fpath;
                b_local_idx = j;
                b_trigbit  = trigbit;
                fout->cd();
                t_events->Fill();
                if      (trigbit == trig_pt3) ++n_pt3_in_file;
                else if (trigbit == trig_pt2) ++n_pt2_in_file;
            }
        }

        // Per-file metadata entry — written even if no triggers, so the
        // segmenter knows about every file that was scanned (file boundary
        // accounting works correctly).
        f_path     = fpath;
        f_n_events = N;
        f_n_pt3    = n_pt3_in_file;
        f_n_pt2    = n_pt2_in_file;
        fout->cd();
        t_files->Fill();

        total_pt3 += n_pt3_in_file;
        total_pt2 += n_pt2_in_file;
        total_evt += N;

        // Sparse progress log (first, last, every 25th). For run_parallel_scan
        // the per-chunk logs go to TMPDIR/cfg_NNN.log; the monitor mines
        // these lines for the live progress.
        if (i == 0 || i == files.size() - 1 || (i + 1) % 25 == 0) {
            std::printf("  [%4zu/%4zu] N=%-10lld PT3=%-8lld PT2=%-8lld  %s\n",
                        i + 1, files.size(),
                        (long long) N,
                        (long long) n_pt3_in_file,
                        (long long) n_pt2_in_file,
                        fpath.c_str());
        }

        fin->Close(); delete fin;
    }

    // -- Finalise output ---------------------------------------------------
    fout->cd();
    t_events->Write();
    t_files->Write();
    fout->Close();
    delete fout;

    std::cout << "\n  scanned " << files.size() - n_skipped
              << " / " << files.size() << " files"
              << "  (skipped " << n_skipped << ")\n"
              << "  totals: events=" << total_evt
              << "  N_PT3=" << total_pt3
              << "  N_PT2=" << total_pt2
              << "  trigger_events_emitted=" << t_events->GetEntries()
              << "\n";

    std::cout << "Scan Complete!\n";
    return 0;
}

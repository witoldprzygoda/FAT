// ============================================================================
// s5_blind_prep_check.C  --  the pseudo-data of the blind test contain no truth and no weight
// ============================================================================
// Reads ONLY /hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/pseudo (never the sealed directory) and checks:
//   - every object of every ROOT file (class, name, title, axis titles), every branch and leaf of every tree
//   - whitelists: eta_toys.root tree 'eta_toys' = {toy, evt, m4, mee, mm4, dz_pi, pivtx}, tree 'toys' = {id, index};
//     histograms only 'toy_<6 digits>_m4mee' (TH2D) and 'toy_<6 digits>_mee' (TH1D); TNamed 'info'
//   - no forbidden token (truth, weight, class, set, BR, ...) in any name, title, axis title or branch / leaf name
//   - X17 toys: ids unique, both histograms present, every bin a non-negative integer, no Sumw2 array (unweighted),
//     the M(ee) histogram equals the projection of the 2D histogram over the X17 M4 window bins
//   - eta toys: entries per toy, finite values, pivtx in {0, 1}
// Output: plots/eta_sim/tables/s5_blind_prep_pseudo_branches.txt
// Run from /home/przygoda/HADES/FAT/FAT:
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_s5_blind_prep", kTRUE)' \
//        'plots/eta_sim/s5_blind_prep_check.C+'
// ============================================================================
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TClass.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TTree.h"
#include "TList.h"

namespace s5bc {
const char* const kPseudo = "/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/blind/pseudo";
const char* const kOut = "/home/przygoda/HADES/FAT/FAT/plots/eta_sim/tables/s5_blind_prep_pseudo_branches.txt";
FILE* gF = nullptr;
void P(const char* fmt, ...) {
    char buf[16384];
    va_list ap; va_start(ap, fmt); vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    fputs(buf, stdout);
    if (gF) fputs(buf, gF);
}
// tokens: lower-case alphanumeric runs ('_' and every other character split)
std::vector<std::string> Tokens(const std::string& s) {
    std::vector<std::string> t;
    std::string cur;
    for (char ch : s) {
        if (std::isalnum((unsigned char)ch)) cur += (char)std::tolower((unsigned char)ch);
        else { if (!cur.empty()) t.push_back(cur); cur.clear(); }
    }
    if (!cur.empty()) t.push_back(cur);
    return t;
}
const std::set<std::string> kForbidden = {
    "w", "wdb", "weight", "weights", "sumw", "cls", "lcls", "pcls", "tmpl", "template", "templates", "own", "owna", "ownb",
    "flag", "flags", "ghost", "set", "sets", "grp", "fold", "sim", "true", "truth", "br", "class", "classes", "inject",
    "injected", "signal", "bkg", "background", "key", "mu", "expect", "expected", "half", "parent", "parentid", "geninfo",
    "std", "lep", "gsig2", "gsig3", "gsig4", "gsig5", "geta", "gpi0d", "gconv", "gvec", "goth", "direct", "conversion",
    "dalitz", "bump", "peak45", "x17perbr", "seed"};
int gViol = 0;
void CheckText(const std::string& where, const std::string& s) {
    for (const std::string& t : Tokens(s))
        if (kForbidden.count(t)) { ++gViol; P("  VIOLATION: forbidden token '%s' in %s: \"%s\"\n", t.c_str(), where.c_str(), s.c_str()); }
}
bool IsToyHist(const std::string& n, std::string& kind, int& id) {
    // toy_<6 digits>_m4mee or toy_<6 digits>_mee
    if (n.size() < 14 || n.compare(0, 4, "toy_") != 0) return false;
    for (int i = 4; i < 10; ++i) if (!std::isdigit((unsigned char)n[i])) return false;
    if (n[10] != '_') return false;
    kind = n.substr(11);
    if (kind != "m4mee" && kind != "mee") return false;
    id = std::atoi(n.substr(4, 6).c_str());
    return true;
}

void CheckTree(TTree* t, const std::set<std::string>& allowed) {
    P("  tree '%s' (title \"%s\"), %lld entries; branches:\n", t->GetName(), t->GetTitle(), t->GetEntries());
    CheckText(std::string("tree name"), t->GetName());
    CheckText(std::string("tree title"), t->GetTitle());
    TObjArray* br = t->GetListOfBranches();
    for (int i = 0; i < br->GetEntriesFast(); ++i) {
        TBranch* b = (TBranch*)br->At(i);
        TObjArray* lv = b->GetListOfLeaves();
        std::string leaves;
        for (int j = 0; j < lv->GetEntriesFast(); ++j) {
            TLeaf* l = (TLeaf*)lv->At(j);
            leaves += std::string(j ? ", " : "") + l->GetName() + "/" + l->GetTypeName();
            CheckText("leaf name", l->GetName());
        }
        P("    %-8s  leaves: %-24s  title \"%s\"%s\n", b->GetName(), leaves.c_str(), b->GetTitle(),
          allowed.count(b->GetName()) ? "" : "   <-- NOT IN THE WHITELIST");
        if (!allowed.count(b->GetName())) ++gViol;
        CheckText("branch name", b->GetName());
        CheckText("branch title", b->GetTitle());
    }
}

void CheckEtaFile(const std::string& fn) {
    TFile* f = TFile::Open(fn.c_str(), "READ");
    if (!f || f->IsZombie()) { P("cannot open %s\n", fn.c_str()); ++gViol; return; }
    P("\n## %s\n", fn.c_str());
    TIter nx(f->GetListOfKeys());
    TKey* k;
    while ((k = (TKey*)nx())) {
        P("  key: %-10s class %-8s title \"%s\"\n", k->GetName(), k->GetClassName(), k->GetTitle());
        CheckText("key name", k->GetName());
        const std::string cn = k->GetClassName();
        if (cn == "TTree") {
            if (std::string(k->GetName()) != "eta_toys") { ++gViol; P("  VIOLATION: unexpected tree %s\n", k->GetName()); continue; }
        } else if (cn == "TNamed") {
            if (std::string(k->GetName()) != "info") { ++gViol; P("  VIOLATION: unexpected TNamed %s\n", k->GetName()); }
            TNamed* nm = (TNamed*)k->ReadObj();
            CheckText("info text", nm->GetTitle());
            P("  info: \"%s\"\n", nm->GetTitle());
        } else { ++gViol; P("  VIOLATION: unexpected object class %s\n", cn.c_str()); }
    }
    TTree* t = (TTree*)f->Get("eta_toys");
    if (!t) { P("  tree eta_toys missing\n"); ++gViol; f->Close(); return; }
    CheckTree(t, {"toy", "evt", "m4", "mee", "mm4", "dz_pi", "pivtx"});
    Int_t toy, evt, pivtx;
    Float_t m4, mee, mm4, dz;
    t->SetBranchAddress("toy", &toy); t->SetBranchAddress("evt", &evt); t->SetBranchAddress("pivtx", &pivtx);
    t->SetBranchAddress("m4", &m4); t->SetBranchAddress("mee", &mee); t->SetBranchAddress("mm4", &mm4); t->SetBranchAddress("dz_pi", &dz);
    std::map<int, long long> rows;
    std::map<int, int> maxevt;
    long long bad = 0, badpv = 0, outM4 = 0;
    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        rows[toy]++;
        maxevt[toy] = std::max(maxevt[toy], evt);
        if (!std::isfinite(m4) || !std::isfinite(mee) || !std::isfinite(mm4) || !std::isfinite(dz)) ++bad;
        if (pivtx != 0 && pivtx != 1) ++badpv;
        if (!(m4 >= 0.30f && m4 < 1.20f)) ++outM4;
    }
    double mn = 1e30, mx = 0, mean = 0;
    for (auto& x : rows) { mn = std::min(mn, (double)x.second); mx = std::max(mx, (double)x.second); mean += x.second; }
    if (!rows.empty()) mean /= rows.size();
    P("  toys: %zu (index %d..%d); rows per toy: mean %.1f, min %.0f, max %.0f; non-finite values %lld; pivtx outside {0,1} %lld; "
      "M4 outside 0.30-1.20: %lld\n", rows.size(), rows.empty() ? -1 : rows.begin()->first, rows.empty() ? -1 : rows.rbegin()->first,
      mean, mn, mx, bad, badpv, outM4);
    if (bad || badpv) ++gViol;
    f->Close();
    delete f;
}

void CheckX17File(const std::string& fn) {
    TFile* f = TFile::Open(fn.c_str(), "READ");
    if (!f || f->IsZombie()) { P("cannot open %s\n", fn.c_str()); ++gViol; return; }
    P("\n## %s\n", fn.c_str());
    std::map<int, int> seen;   // id -> bit 1 m4mee, 2 mee
    std::map<std::string, int> classes;
    int nHist = 0, nNonInt = 0, nNeg = 0, nSumw2 = 0, nProjBad = 0;
    TIter nx(f->GetListOfKeys());
    TKey* k;
    int shown = 0;
    while ((k = (TKey*)nx())) {
        const std::string cn = k->GetClassName(), nm = k->GetName();
        classes[cn]++;
        CheckText("key title", k->GetTitle());
        if (cn == "TTree") {
            P("  key: %-10s class %-8s title \"%s\"\n", nm.c_str(), cn.c_str(), k->GetTitle());
            if (nm != "toys") { ++gViol; P("  VIOLATION: unexpected tree %s\n", nm.c_str()); }
            continue;
        }
        if (cn == "TNamed") {
            TNamed* x = (TNamed*)k->ReadObj();
            P("  key: %-10s class %-8s\n  info: \"%s\"\n", nm.c_str(), cn.c_str(), x->GetTitle());
            if (nm != "info") { ++gViol; P("  VIOLATION: unexpected TNamed %s\n", nm.c_str()); }
            CheckText("info text", x->GetTitle());
            continue;
        }
        std::string kind;
        int id = 0;
        if (!IsToyHist(nm, kind, id)) { ++gViol; P("  VIOLATION: unexpected object %s (%s)\n", nm.c_str(), cn.c_str()); continue; }
        if ((kind == "m4mee" && cn != "TH2D") || (kind == "mee" && cn != "TH1D")) { ++gViol; P("  VIOLATION: %s has class %s\n", nm.c_str(), cn.c_str()); }
        seen[id] |= kind == "m4mee" ? 1 : 2;
        ++nHist;
        if (shown < 2) { P("  example key: %-18s class %-5s title \"%s\"\n", nm.c_str(), cn.c_str(), k->GetTitle()); ++shown; }
    }
    // content checks (all toys)
    int nToy = 0;
    for (auto& s : seen) {
        ++nToy;
        if (s.second != 3) { ++gViol; P("  VIOLATION: toy %06d lacks a histogram\n", s.first); continue; }
        TH2D* h = (TH2D*)f->Get(Form("toy_%06d_m4mee", s.first));
        TH1D* g = (TH1D*)f->Get(Form("toy_%06d_mee", s.first));
        for (TH1* x : {(TH1*)h, (TH1*)g}) {
            CheckText("histogram title", x->GetTitle());
            CheckText("x axis title", x->GetXaxis()->GetTitle());
            CheckText("y axis title", x->GetYaxis()->GetTitle());
            CheckText("z axis title", x->GetZaxis()->GetTitle());
            if (x->GetSumw2N() > 0) ++nSumw2;
        }
        // X17 M4 window bins from the README convention: 0.5330-0.5555 = x bins 21-29
        const int b0 = h->GetXaxis()->FindBin(0.5330 + 1e-6), b1 = h->GetXaxis()->FindBin(0.5555 - 1e-6);
        for (int iy = 0; iy <= h->GetNbinsY() + 1; ++iy) {
            double proj = 0;
            for (int ix = 0; ix <= h->GetNbinsX() + 1; ++ix) {
                const double c = h->GetBinContent(ix, iy);
                if (c < 0) ++nNeg;
                if (c != std::floor(c)) ++nNonInt;
                if (ix >= b0 && ix <= b1) proj += c;
            }
            if (iy >= 1 && iy <= g->GetNbinsX() && std::fabs(proj - g->GetBinContent(iy)) > 1e-9) ++nProjBad;
        }
    }
    P("  objects per class:");
    for (auto& c : classes) P(" %s %d;", c.first.c_str(), c.second);
    P("\n  toys %d (ids unique by construction of the map; histograms %d); non-integer bins %d, negative bins %d, histograms with "
      "Sumw2 %d, M(ee) projection mismatches %d\n", nToy, nHist, nNonInt, nNeg, nSumw2, nProjBad);
    if (nNonInt || nNeg || nSumw2 || nProjBad) ++gViol;
    TTree* t = (TTree*)f->Get("toys");
    if (t) {
        CheckTree(t, {"id", "index"});
        Int_t id;
        t->SetBranchAddress("id", &id);
        std::set<int> ids;
        for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); ids.insert(id); }
        P("  tree 'toys': %lld entries, %zu distinct ids, all with histograms: %s\n", t->GetEntries(), ids.size(),
          (ids.size() == seen.size() && (Long64_t)ids.size() == t->GetEntries()) ? "yes" : "NO");
        if (!(ids.size() == seen.size() && (Long64_t)ids.size() == t->GetEntries())) ++gViol;
    } else { ++gViol; P("  tree toys missing\n"); }
    f->Close();
    delete f;
}

void Run() {
    gF = fopen(kOut, "w");
    P("# s5_blind_prep_pseudo_branches: content of the blind-test pseudo-data directory (s5_blind_prep_check.C)\n");
    P("# directory %s (the sealed directory is NOT read)\n", kPseudo);
    P("# whitelist: eta_toys {toy, evt, m4, mee, mm4, dz_pi, pivtx}; toys {id, index}; histograms toy_<id>_m4mee (TH2D), toy_<id>_mee (TH1D);\n");
    P("# forbidden tokens in any name / title / branch / leaf:");
    for (const std::string& t : kForbidden) P(" %s", t.c_str());
    P("\n\n## files in the directory\n");
    TSystemDirectory dir("pseudo", kPseudo);
    TList* fl = dir.GetListOfFiles();
    std::vector<std::string> files;
    if (fl) {
        TIter it(fl);
        TObject* o;
        while ((o = it())) {
            const std::string n = o->GetName();
            if (n == "." || n == "..") continue;
            files.push_back(n);
        }
    }
    std::sort(files.begin(), files.end());
    for (const std::string& n : files) {
        Long_t id, flags, mt; Long64_t size = 0;
        gSystem->GetPathInfo((std::string(kPseudo) + "/" + n).c_str(), &id, &size, &flags, &mt);
        P("  %-16s %10.2f MB\n", n.c_str(), size / 1048576.);
        if (n != "eta_toys.root" && n != "x17_toys.root" && n != "README.md") { ++gViol; P("  VIOLATION: unexpected file %s\n", n.c_str()); }
    }
    CheckEtaFile(std::string(kPseudo) + "/eta_toys.root");
    CheckX17File(std::string(kPseudo) + "/x17_toys.root");
    // README text
    {
        const std::string rn = std::string(kPseudo) + "/README.md";
        FILE* r = fopen(rn.c_str(), "r");
        int nl = 0;
        if (r) {
            char line[8192];
            while (fgets(line, sizeof line, r)) {
                ++nl;
                std::string s(line);
                // the README describes the content by design; only numbers of the injected BR values or any truth count are forbidden
                for (const char* bad : {"1e-6", "3e-6", "1e-5", "3e-5", "1e-06", "3e-06", "1e-05", "3e-05", "seed"})
                    if (s.find(bad) != std::string::npos) { ++gViol; P("  VIOLATION: README line %d contains '%s'\n", nl, bad); }
            }
            fclose(r);
        }
        P("\n## README.md: %d lines, checked for injected BR values and seeds\n", nl);
    }
    P("\nRESULT: %s (%d violations)\n", gViol == 0 ? "PASS: no truth, weight, class or set column, no forbidden token" : "FAIL", gViol);
    fclose(gF);
    gF = nullptr;
}
}  // namespace s5bc

void s5_blind_prep_check() { s5bc::Run(); }

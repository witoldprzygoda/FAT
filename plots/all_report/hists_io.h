// ============================================================================
// hists_io.h  --  access to hists_v1.root of SPEC report_all (written by WP0)
// ============================================================================
// Usage (ACLiC macro in plots/all_report/):
//   #include "hists_io.h"
//   TFile* f = TFile::Open(alldefs::kHistFile);
//   TH1* n  = allio::Data(f, "S0", "all", "N",  allio::F10(2, 1));        // data rows, slice S3, class E2
//   TH1* cb = allio::Data(f, "S0", "all", "CB", allio::F10(2, 1));
//   TH1* c5 = allio::Sim (f, "S0", "nom", alldefs::kC5, allio::F10(2, 1)); // nb per bin
//   c5->Scale(alldefs::kNbToRaw);                                         // -> expected raw data rows (L x 63/64)
// Every function returns a clone owned by the caller (SetDirectory(nullptr)), or nullptr with a message.
// Data N of X2D (M(ee) 12-28 MeV) and of B_mee_R2..R4 (12-28 MeV) is -1 (blinded, D19): use MaskedX2D / MaskedBmee.
// ============================================================================
#ifndef ALLREPORT_HISTS_IO_H
#define ALLREPORT_HISTS_IO_H

#include <cstdio>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "defs.h"

namespace allio {

inline std::string F10(int slice, int ecl) { return alldefs::F10Name(slice, ecl); }   // slice 0..4, ecl 0 (E1) / 1 (E2)
inline std::string F5g(int slice) { return alldefs::F5gName(slice); }
inline std::string BVar(int var, int reg) { return alldefs::BVarName(var, reg); }

inline TH1* Get(TFile* f, const std::string& path) {
    TH1* h = f ? (TH1*)f->Get(path.c_str()) : nullptr;
    if (!h) { printf("allio: %s not found\n", path.c_str()); return nullptr; }
    TH1* c = (TH1*)h->Clone((path + "_c").c_str());
    c->SetDirectory(nullptr);
    return c;
}
// data: sel = S0, pre, dz40, PEAK, vz150, vz170, none, S0_best, pre_best, S0_corfb, pre_corfb;
//       sub = all, fold0, fold1, perI, perIIA, perIIB, perIII; what = N, Npp, Nmm, CB, CB20 (all only), CBrep (all only)
//       layout = F10_S<j>E<e>, F5g_S<j>, X2D, fine_M4, B_M4, win, tot; sub all also G3, B_map, B_<var>_R<k>
inline TH1* Data(TFile* f, const char* sel, const char* sub, const char* what, const std::string& layout) {
    return Get(f, std::string("data/") + sel + "/" + sub + "/" + what + "_" + layout);
}
// N - CB (raw rows; errors: N event-group and CB in quadrature); nullptr if a part is missing
inline TH1* DataSig(TFile* f, const char* sel, const char* sub, const std::string& layout, const char* cbName = "CB") {
    TH1* n = Data(f, sel, sub, "N", layout);
    TH1* cb = Data(f, sel, sub, cbName, layout);
    if (!n || !cb) { delete n; delete cb; return nullptr; }
    n->Add(cb, -1.);
    delete cb;
    return n;
}
// simulation: sel = S0, pre (all variants) or dz40, PEAK, vz150, vz170 (nom only); var = alldefs::SimVarTag(v);
//             comp = alldefs::Comp; suffix "" or "_hA", "_hB", "_rep" (S0__nom, pre__nom only)
inline TH1* Sim(TFile* f, const char* sel, const char* var, int comp, const std::string& layout, const char* suffix = "") {
    return Get(f, std::string("sim/") + sel + "__" + var + "/" + alldefs::CompTag(comp) + "_" + layout + suffix);
}
// window content of a "win" histogram by window id (alldefs::WinId)
inline double Win(const TH1* w, int id) { return w ? w->GetBinContent(id + 1) : 0.; }
inline double WinErr(const TH1* w, int id) { return w ? w->GetBinError(id + 1) : 0.; }
// blinding masks (1 = masked)
inline bool MaskedX2D(int binY) { return binY >= 6 && binY <= 13; }      // M(ee) 12-28 MeV, all M4 bins
inline bool MaskedBmee(int bin, int reg) { return reg >= 2 && reg <= 4 && bin >= 7 && bin <= 14; }   // 12-28 MeV in R2-R4

}  // namespace allio

#endif  // ALLREPORT_HISTS_IO_H

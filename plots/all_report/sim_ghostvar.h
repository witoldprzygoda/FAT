// ============================================================================
// sim_ghostvar.h  --  treatment of ghost-flagged true e+e- pairs (review fix)
// ============================================================================
// src/sim_truth.h TemplateClass sends the four-particle classes 2-6 to the
// signal (T1) BEFORE it looks at the ghost-lepton flag, but sends every other
// row with a ghost-flagged lepton track to T3 "fake lepton". A true one-source
// e+e- pair (lepton-pair class L1-L6) with a ghost-flagged track is therefore
// signal when its pions fit an eta channel and "fake lepton" when it is
// physical background. About 90 % of T3 are such pairs. The like-sign samples
// do not contain them, so the data CB cannot remove them: data N - CB contains
// them. The default comparison (T1 + T2) treats the ghost flag differently in
// signal and background. The two consistent variants:
//   (a) "ghosts in both":     T1 + T2 + ghost-flagged true pairs (GhostPair)
//   (b) "no ghosts anywhere": T1 without ghost-lepton rows + T2
// (T2 never has a ghost lepton: the flag sends such rows to T3.)
// The binding rule 1 of the study (fake lepton not in the compared spectrum)
// was written for misidentified leptons; the user must confirm the treatment.
// ============================================================================
#ifndef SIM_GHOSTVAR_H
#define SIM_GHOSTVAR_H

#include "SimEngine.h"

namespace simeng {

// ghost-flagged true one-source pair in T3 (not the same particle twice); the eta set adds the
// one-source classes L10 (eta -> pi+pi-e+e- leptons) and L11 (X17 leptons, W = 0 in the production)
inline bool GhostPair(const Row& r) {
    return r.tmpl == 3 && (r.flags & simtruth::kFlagGhostLepton) != 0 &&
           (r.flags & simtruth::kFlagSameParticle) == 0 && ((r.lcls >= 1 && r.lcls <= 6) || r.lcls == 10 || r.lcls == 11);
}
// ghost-variant of the compared simulated spectrum: 0 default T1 + T2, 1 (a), 2 (b)
inline bool InSimVar(const Row& r, int v) {
    if (X17Row(r) && r.tmpl == 1) return false;   // X17 signal: information only (W = 0 in the production)
    if (v == 1) return r.tmpl == 1 || r.tmpl == 2 || GhostPair(r);
    if (v == 2) return (r.tmpl == 1 && !GhostLepton(r)) || r.tmpl == 2;
    return r.tmpl == 1 || r.tmpl == 2;
}
inline const char* SimVarTag(int v) { return v == 1 ? "a_ghostsBoth" : v == 2 ? "b_noGhosts" : "T12"; }
inline const char* SimVarName(int v) {
    return v == 1 ? "(a) T1 + T2 + ghost-flagged true pairs" : v == 2 ? "(b) T1 without ghost rows + T2" : "default T1 + T2";
}
inline const char* SimVarLatex(int v) {
    return v == 1 ? "(a) T1 + T2 + ghost-flagged true pairs" : v == 2 ? "(b) T1 without ghost rows + T2" : "T1 + T2 (default)";
}
inline Comp CompSimVar(int v) { return {SimVarName(v), [v](const Row& r) { return InSimVar(r, v); }}; }
inline Comp CompGhostPair() { return {"ghost-flagged true pairs (T3)", [](const Row& r) { return GhostPair(r); }}; }
// background of a variant: T2 (+ ghost pairs in (a))
inline Comp CompBkgVar(int v) {
    return {std::string("background ") + SimVarTag(v), [v](const Row& r) { return r.tmpl == 2 || (v == 1 && GhostPair(r)); }};
}
// signal classes of a variant: (b) drops ghost rows
inline Comp CompSigVar(int v, int clsLo, int clsHi) {
    return {Form("signal cls %d-%d %s", clsLo, clsHi, SimVarTag(v)), [v, clsLo, clsHi](const Row& r) {
                return r.tmpl == 1 && r.cls >= clsLo && r.cls <= clsHi && !(v == 2 && GhostLepton(r));
            }};
}

// Relabel the template column in place so that every engine function (CompSim,
// CompPhysBkg, the template fit) follows a variant: (a) ghost pairs -> T2,
// (b) ghost signal rows -> T3. Restore with RestoreTmpl. Keep the saved copy.
struct TmplSave { std::vector<signed char> t[kNTreeAll]; };
inline void RelabelTmpl(SimData& d, int v, TmplSave& s) {
    for (int k = 0; k < kNTreeAll; ++k) {
        s.t[k].resize(d.t[k].size());
        for (size_t i = 0; i < d.t[k].size(); ++i) {
            Row& r = d.t[k][i];
            s.t[k][i] = (signed char)r.tmpl;
            if (v == 1 && GhostPair(r)) r.tmpl = 2;
            if (v == 2 && r.tmpl == 1 && GhostLepton(r)) r.tmpl = 3;
        }
    }
}
inline void RestoreTmpl(SimData& d, const TmplSave& s) {
    for (int k = 0; k < kNTreeAll; ++k)
        for (size_t i = 0; i < d.t[k].size() && i < s.t[k].size(); ++i) d.t[k][i].tmpl = s.t[k][i];
}

// vertex-z slices (mm) used by the review; base = "pre" without the vz > -160 cut
struct VzSlice { const char* tag; const char* name; float lo, hi; };
const VzSlice kVzSl[5] = {{"pre", "v_{z} > -160 (pre)", -160.f, 1e9f},
                          {"central", "-140 < v_{z} < -100", -140.f, -100.f},
                          {"edge", "-160 < v_{z} < -140", -160.f, -140.f},
                          {"down", "v_{z} > -100", -100.f, 1e9f},
                          {"novz", "no v_{z} cut (z > -500)", -500.f, 1e9f}};
template <class R> inline bool PassBaseNoVz(const R& r) {
    return r.isBest >= 0 && r.vz > -500.f && r.cor_ok == 1 && PT3(r);
}
template <class R> inline bool PassVz(const R& r, int iv) {
    return PassBaseNoVz(r) && r.vz > kVzSl[iv].lo && r.vz < kVzSl[iv].hi;
}

}  // namespace simeng

#endif  // SIM_GHOSTVAR_H

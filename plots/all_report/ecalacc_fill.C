// ============================================================================
// ecalacc_fill.C  --  feasibility of the FULL reconstruction of the eta bump with the ECAL photon
// ============================================================================
// SPEC report_all, package ecalacc (2026-11-15). New file; nothing else is changed.
//
// Question: eta -> pi+pi-pi0, pi0 -> e+e-gamma (class C2) and pi0 -> gamma gamma with one photon converted (C3).
// The bump is the part of these decays in which the photon is NOT seen. How many of the accepted pi+pi-e+e- rows
// have their photon inside the HADES ECAL (which covers only part of the polar angles and 5 of the 6 sectors), how
// many are found by the ECAL simulation, and what does this predict for the data?
//
//  A. DATA (scan_v2 compact, cor_ok == 1): the ECAL tag of WP6 / report_exp (good photon, 0.119 <= M(e+e-gamma) < 0.173,
//     M(ee) < 0.10; m5c with the photon energy rescaled to M(e+e-gamma) = m_pi0; frozen fit 0.5463 / 6.2 MeV, pol2,
//     0.49-0.62) for "pre" and S0, all and per M(pi+pi-) slice; photon (theta, phi, E) maps of the tagged rows;
//     M4 distance to the bump edge of the tagged eta (m5c peak window minus side bands).
//  B. SIMULATION (SimEngine copy: mode A, w x V3 = code 108, owned rows, std + lep + eta; expected raw data rows =
//     sum w x 27.0718): the same emulated tag on T1 + T2 (UNVERIFIED, as WP6) for pre and S0, all and per slice.
//  C. TRUTH of the bump rows (classes C2, C3; the standard set owns them): the true momenta of the four charged tracks
//     and the ECAL clusters with their GEANT truth are read from the PAT file (entry pat_entry).
//     - Photon truth for EVERY row ("conditional completion"): given the four true charged momenta, the missing
//       photon k must satisfy (Q + k)^2 = m_pi0^2 (Q = e+e-, or the converted photon for C3) and (P + k)^2 = m_eta^2
//       (P = pi+pi-e+e-). In the rest frame of Q this fixes |k| and the angle between k and P, leaving one free
//       azimuth phi with a FLAT measure (d^3k/k delta delta = const dphi). Its density along the circle is
//       f(p_eta = P + k) x |M|^2, f = the SMASH eta invariant density E d^3N/dp^3 (from the converted SMASH files of the standard set, about 60 files, 140 572 etas),
//       |M|^2 = 1 for the simulation (HGeant: phase space, isotropic gamma*). 720 phi points per row.
//       Variants: KLOE-2 Dalitz plot f(X, Y) along the circle; QED gamma* angle 1 + cos^2 + (4 m_e^2/q^2) sin^2
//       (normalised to mean 1), which also gives the C2 angle factor of the WHOLE sample (cross-check of 0.8855).
//     - ECAL geometry and cell status of the simulation (EmcGeomPar, EmcCellStatusPar of
//       allParam_feb22_sim_run_17000_gen4_10032026.root, read on vae25 by ecalacc_emcgeo_vae25.C): per sector a plane
//       (rms 0.003 mm) with 163 modules of 92 x 92 mm (pitch 94 mm); a ray from the true vertex hits a module if it
//       crosses its front face enlarged to the pitch. Status 0 modules (all of sector 0 and 8 single modules, the same in the data parameters) do not count.
//     - Closure: the photon found by the ECAL simulation (cluster truth: GEANT id 1, parent pi0, the codes and the
//       vertex of the pions) must lie on the circle, and the geometric test must say "inside" for it.
//     - Tag of the truth rows recomputed from the PAT clusters (reproduces the compact exactly for lambda = 1) and
//       with the cluster energies x lambda (lambda = data-like calibration), to see what the 100 MeV threshold and the
//       M(e+e-gamma) window do when the energy scale changes.
//
// Outputs (/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/): ecalacc_fill.txt (log + KEY lines),
//   rows_c23.bin + rows_c23_layout.txt (per-row table for ecalacc_ana.py), hist_*.txt (histograms), tag_fits.txt.
// Run (from /home/przygoda/HADES/FAT/FAT; about 10 GB, a few minutes):
//   root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/all_report/aclic_ecalacc", kTRUE)' \
//        'plots/all_report/ecalacc_fill.C+O'
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "sim_eta_physics.h"
#include "defs.h"

namespace ecalacc {
using namespace alldefs;

const char* const kOutDir  = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc";
const char* const kCells   = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/emc_cells_sim_v2.txt";
const char* const kEtaLab  = "/hdd3/przygoda/HADES/FAT_out/all_report/ecalacc/smash_eta_lab.txt";
const char* const kPatFile[4] = {"", "/hdd3/przygoda/HADES/files/sim/ppee_std_all.root",
                                 "/hdd3/przygoda/HADES/files/sim/ppee_lep_all.root",
                                 "/hdd3/przygoda/HADES/files/sim/ppee_eta_all.root"};
FILE* gOut = nullptr;
#define P(...) do { printf(__VA_ARGS__); if (gOut) fprintf(gOut, __VA_ARGS__); } while (0)

// ---- tag (frozen; copy of wp6_fill.C, which copies scan_bump.C) ------------------------------------------------
const double kPi0Lo = 0.1190, kPi0Hi = 0.1730;
const double kMeeMax = 0.100;
const double kTagMean = 0.5463, kTagSigma = 0.0062, kTagLo = 0.49, kTagHi = 0.62;
const int kTagOrder = 2, kTagNbin = 160, kTagReb = 8;
const double kTagHLo = 0.40, kTagHHi = 0.80;
template <class R> inline double M5c(const R& r) {
    const double a = (double)r.m_eeg * r.m_eeg - (double)r.mee * r.mee;
    if (!(a > 0)) return -1.0;
    const double f = (scan::kMpi0 * scan::kMpi0 - (double)r.mee * r.mee) / a;
    const double m2 = (double)r.m4 * r.m4 + f * ((double)r.m5 * r.m5 - (double)r.m4 * r.m4);
    return m2 > 0 ? std::sqrt(m2) : -1.0;
}
inline double M5cOf(double m4, double mee, double meeg, double m5) {
    const double a = meeg * meeg - mee * mee;
    if (!(a > 0)) return -1.0;
    const double f = (scan::kMpi0 * scan::kMpi0 - mee * mee) / a;
    const double m2 = m4 * m4 + f * (m5 * m5 - m4 * m4);
    return m2 > 0 ? std::sqrt(m2) : -1.0;
}
template <class R> inline bool Tagged(const R& r) {
    return r.g_idx > 0 && r.mee < kMeeMax && r.m_eeg >= kPi0Lo && r.m_eeg < kPi0Hi;
}
inline scan::PeakCfg TagCfg() {
    scan::PeakCfg c;
    c.mean = kTagMean; c.sigma = kTagSigma; c.lo = kTagLo; c.hi = kTagHi; c.order = kTagOrder;
    c.winLo = c.mean - scan::kWinNsig * c.sigma; c.winHi = c.mean + scan::kWinNsig * c.sigma;
    return c;
}
scan::Spectra SimSpectra(const double* y, int nbin, double lo, double hi, int reb) {   // copy of wp6_fill.C
    scan::Spectra s;
    s.nbin = nbin; s.reb = reb; s.nw = nbin / reb; s.lo = lo; s.hi = hi; s.bw = (hi - lo) / nbin;
    for (int k = 0; k < 3; ++k) { s.n[k].assign(nbin, 0.); s.w[k].assign(s.nw, 0.); }
    for (int i = 0; i < nbin; ++i) s.n[0][i] = y[i];
    for (int i = 0; i < nbin; ++i) s.w[0][i / reb] += y[i];
    s.cbW.assign(s.nw, 0.); s.cbWVar.assign(s.nw, 0.);
    s.cb.assign(nbin, 0.); s.cbVar.assign(nbin, 0.);
    s.sig.assign(y, y + nbin); s.sigVar.assign(y, y + nbin);
    for (int k = 0; k < 3; ++k) s.nsel[k] = 0;
    for (int i = 0; i < nbin; ++i) s.nsel[0] += y[i];
    return s;
}
inline int TagBin(double m5c) {
    if (!(m5c >= kTagHLo && m5c < kTagHHi)) return -1;
    const int i = (int)((m5c - kTagHLo) / ((kTagHHi - kTagHLo) / kTagNbin));
    return (i >= 0 && i < kTagNbin) ? i : -1;
}

// M(pi+pi-) slices of the tag study (DATA scale): 0 all, 1 < 0.34, 2 0.34-0.37, 3 0.37-0.40, 4 0.40-0.43, 5 >= 0.43
const int kNTs = 6;
const double kTsLo[kNTs] = {0.0, 0.0, 0.34, 0.37, 0.40, 0.43};
const double kTsHi[kNTs] = {9.0, 0.34, 0.37, 0.40, 0.43, 9.0};
const char* const kTsName[kNTs] = {"all", "<0.34", "0.34-0.37", "0.37-0.40", "0.40-0.43", ">=0.43"};
inline bool InTs(int s, double mpp) { return mpp >= kTsLo[s] && mpp < kTsHi[s]; }

// distance of M4 to the bump edge E(M(pipi), M(ee)) (data scale), GeV/c2; -9 above the kinematic limit
inline double DistEdge(double m4, double mpp, double mee) {
    const double e = Edge(mpp, mee);
    return e > 0 ? m4 - e : -9.;
}

// ---- simple text histograms -------------------------------------------------------------------------------------
struct H1 {
    std::string name; int n = 0; double lo = 0, hi = 1; std::vector<double> s, v;
    void Init(const std::string& nm, int nb, double a, double b) { name = nm; n = nb; lo = a; hi = b; s.assign(nb, 0.); v.assign(nb, 0.); }
    void Fill(double x, double w) { if (!(x >= lo && x < hi)) return; const int i = (int)((x - lo) / (hi - lo) * n); if (i >= 0 && i < n) { s[i] += w; v[i] += w * w; } }
    void Write(FILE* f) const {
        fprintf(f, "HIST %s %d %.6f %.6f\n", name.c_str(), n, lo, hi);
        for (int i = 0; i < n; ++i) fprintf(f, "%.8g %.8g\n", s[i], v[i]);
    }
};

// ============================================================================
// ECAL geometry
// ============================================================================
struct Cell { int sec, cell, status; double u[4], v[4], uc, vc; double ue[4], ve[4]; };
struct Sector { bool have = false; TVector3 c, n, a, b; std::vector<int> cells; double umin = 1e9, umax = -1e9, vmin = 1e9, vmax = -1e9; };
struct Geo {
    std::vector<Cell> cells; Sector sec[6];
    bool Load(const char* fn) {
        std::ifstream in(fn);
        if (!in) { P("cannot open %s\n", fn); return false; }
        std::string line;
        std::vector<std::array<TVector3, 4>> fronts; std::vector<std::array<int, 3>> ids;
        while (std::getline(in, line)) {
            std::istringstream is(line); std::string tag; is >> tag; if (tag != "CELL") continue;
            int s, c, st; double th, ph, x, y, z; is >> s >> c >> st >> th >> ph >> x >> y >> z;
            std::vector<TVector3> p(8);
            for (int i = 0; i < 8; ++i) { double a, b, cc; is >> a >> b >> cc; p[i].SetXYZ(a, b, cc); }
            std::vector<int> o(8); for (int i = 0; i < 8; ++i) o[i] = i;
            std::sort(o.begin(), o.end(), [&](int i, int j) { return p[i].Mag() < p[j].Mag(); });
            std::array<TVector3, 4> f = {p[o[0]], p[o[1]], p[o[2]], p[o[3]]};
            fronts.push_back(f); ids.push_back({s, c, st});
        }
        // plane per sector from all front corners (centroid + normal = the direction of the centroid, refined by
        // the cross product of two in-plane directions; the planes are flat to 0.003 mm)
        for (int s = 0; s < 6; ++s) {
            std::vector<TVector3> pts;
            for (size_t k = 0; k < ids.size(); ++k) if (ids[k][0] == s) for (int i = 0; i < 4; ++i) pts.push_back(fronts[k][i]);
            if (pts.empty()) continue;
            TVector3 c(0, 0, 0); for (auto& q : pts) c += q; c *= 1.0 / pts.size();
            // normal: from one cell front face (all coplanar)
            size_t k0 = 0; for (size_t k = 0; k < ids.size(); ++k) if (ids[k][0] == s) { k0 = k; break; }
            TVector3 e1 = fronts[k0][1] - fronts[k0][0], e2 = fronts[k0][2] - fronts[k0][0];
            TVector3 n = e1.Cross(e2).Unit(); if (n.Dot(c) < 0) n = -n;
            TVector3 a = (TVector3(0, 0, 1) - n * n.Z()); if (a.Mag() < 1e-6) a = TVector3(1, 0, 0) - n * n.X(); a = a.Unit();
            TVector3 b = n.Cross(a).Unit();
            sec[s].have = true; sec[s].c = c; sec[s].n = n; sec[s].a = a; sec[s].b = b;
        }
        for (size_t k = 0; k < ids.size(); ++k) {
            Cell q; q.sec = ids[k][0]; q.cell = ids[k][1]; q.status = ids[k][2];
            Sector& S = sec[q.sec];
            double uu[4], vv[4];
            for (int i = 0; i < 4; ++i) { TVector3 d = fronts[k][i] - S.c; uu[i] = d.Dot(S.a); vv[i] = d.Dot(S.b); }
            q.uc = 0.25 * (uu[0] + uu[1] + uu[2] + uu[3]); q.vc = 0.25 * (vv[0] + vv[1] + vv[2] + vv[3]);
            std::vector<int> o = {0, 1, 2, 3};
            std::sort(o.begin(), o.end(), [&](int i, int j) { return std::atan2(vv[i] - q.vc, uu[i] - q.uc) < std::atan2(vv[j] - q.vc, uu[j] - q.uc); });
            for (int i = 0; i < 4; ++i) {
                q.u[i] = uu[o[i]]; q.v[i] = vv[o[i]];
                q.ue[i] = q.uc + (q.u[i] - q.uc) * (94.0 / 92.0); q.ve[i] = q.vc + (q.v[i] - q.vc) * (94.0 / 92.0);
                S.umin = std::min(S.umin, q.ue[i]); S.umax = std::max(S.umax, q.ue[i]);
                S.vmin = std::min(S.vmin, q.ve[i]); S.vmax = std::max(S.vmax, q.ve[i]);
            }
            S.cells.push_back((int)cells.size());
            cells.push_back(q);
        }
        int nact = 0; for (auto& q : cells) nact += q.status == 1;
        P("geometry: %zu cells, %d active (status 1)\n", cells.size(), nact);
        return true;
    }
    static bool InQuad(const double* u, const double* v, double x, double y) {
        int pos = 0, neg = 0;
        for (int i = 0; i < 4; ++i) {
            const int j = (i + 1) & 3;
            const double cr = (u[j] - u[i]) * (y - v[i]) - (v[j] - v[i]) * (x - u[i]);
            if (cr > 0) ++pos; else if (cr < 0) ++neg;
        }
        return pos == 0 || neg == 0;
    }
    // ray from vertex V (mm) along direction d: returns the cell index (or -1); strict = the 92 mm faces
    int Hit(const TVector3& V, const TVector3& d, bool strict = false, int* secOut = nullptr) const {
        double ph = d.Phi() * TMath::RadToDeg(); if (ph < 0) ph += 360.;
        double x = ph - 60.; if (x < 0) x += 360.;
        const int s = std::min(5, (int)(x / 60.));
        if (secOut) *secOut = s;
        const Sector& S = sec[s];
        if (!S.have) return -1;
        const double den = S.n.Dot(d);
        if (den <= 1e-9) return -1;
        const double t = S.n.Dot(S.c - V) / den;
        if (t <= 0) return -1;
        const TVector3 X = V + t * d - S.c;
        const double u = X.Dot(S.a), v = X.Dot(S.b);
        if (u < S.umin || u > S.umax || v < S.vmin || v > S.vmax) return -1;
        for (int k : S.cells) {
            const Cell& q = cells[k];
            if (strict ? InQuad(q.u, q.v, u, v) : InQuad(q.ue, q.ve, u, v)) return k;
        }
        return -1;
    }
};

// ============================================================================
// SMASH eta invariant density f(y_lab, mT - m) ~ E d3N/dp3
// ============================================================================
struct EtaDens {
    TH2D* h = nullptr; long n = 0; double meanM = 0;
    bool Load(const char* fn) {
        std::ifstream in(fn);
        if (!in) { P("cannot open %s\n", fn); return false; }
        h = new TH2D("hEtaDens", "", 70, -0.5, 3.0, 45, 0.0, 1.8);
        h->SetDirectory(nullptr);
        double E, px, py, pz;
        while (in >> E >> px >> py >> pz) {
            const double pt = std::hypot(px, py), m2 = E * E - px * px - py * py - pz * pz;
            if (!(m2 > 0) || E <= std::fabs(pz)) continue;
            const double m = std::sqrt(m2), mt = std::sqrt(pt * pt + m2);
            const double y = 0.5 * std::log((E + pz) / (E - pz));
            h->Fill(y, mt - m, 1.0 / mt);      // dN/(dy dmT) / mT  ~ E d3N/dp3
            meanM += m; ++n;
        }
        if (n) meanM /= n;
        P("SMASH eta density: %ld etas from %s, mean mass %.5f GeV/c2\n", n, fn, meanM);
        return n > 1000;
    }
    double operator()(const TLorentzVector& pEta /* GeV */) const {
        const double E = pEta.E(), pz = pEta.Pz();
        if (E <= std::fabs(pz)) return 0.;
        const double y = 0.5 * std::log((E + pz) / (E - pz));
        const double mt = std::sqrt(pEta.Perp2() + meanM * meanM);
        const double x = mt - meanM;
        if (y <= h->GetXaxis()->GetBinCenter(1) || y >= h->GetXaxis()->GetBinCenter(70)) return 0.;
        if (x >= h->GetYaxis()->GetBinCenter(45)) return 0.;
        if (x <= h->GetYaxis()->GetBinCenter(1)) {   // below the first centre: value of the first bin row (interpolated in y)
            const int iy = h->GetYaxis()->FindBin(x);
            (void)iy;
            return std::max(0., h->Interpolate(y, h->GetYaxis()->GetBinCenter(1)));
        }
        return std::max(0., h->Interpolate(y, x));
    }
};

// KLOE-2 2016 eta -> pi+pi-pi0 Dalitz plot (the values of sim_eta_physics.h)
inline double FKloe(double X, double Y) {
    const double a = -1.095, b = 0.145, d = 0.081, f = 0.141, g = -0.044;
    return 1 + a * Y + b * Y * Y + d * X * X + f * Y * Y * Y + g * X * X * Y;
}

// ============================================================================
// PAT truth of one row
// ============================================================================
const int kNN = 5;
struct PatRow {
    float ep_px, ep_py, ep_pz, em_px, em_py, em_pz, pip_px, pip_py, pip_pz, pim_px, pim_py, pim_pz;
    float ep_p, pip_p, pip_vx, pip_vy, pip_vz, ep_vx, ep_vy, ep_vz;
    float pip_g0, pip_g1, pip_g2;
    float nmult;
    float pid[kNN], beta[kNN], ecl[kNN], thc[kNN], phc[kNN], secc[kNN], cellc[kNN];
    float sid[kNN], sppid[kNN], sg0[kNN], sg1[kNN], sg2[kNN], spx[kNN], spy[kNN], spz[kNN], svx[kNN], svy[kNN], svz[kNN];
};
struct PatReader {
    TFile* f = nullptr; TTree* t = nullptr; PatRow r;
    bool Open(const char* fn) {
        f = TFile::Open(fn, "READ");
        if (!f || f->IsZombie()) { P("cannot open %s\n", fn); return false; }
        t = (TTree*)f->Get("PipPimEpEm_ID");
        if (!t) { P("no PipPimEpEm_ID in %s\n", fn); return false; }
        t->SetBranchStatus("*", 0);
        auto B = [&](const char* n, float* a) { t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
        B("ep_sim_px", &r.ep_px); B("ep_sim_py", &r.ep_py); B("ep_sim_pz", &r.ep_pz);
        B("em_sim_px", &r.em_px); B("em_sim_py", &r.em_py); B("em_sim_pz", &r.em_pz);
        B("pip_sim_px", &r.pip_px); B("pip_sim_py", &r.pip_py); B("pip_sim_pz", &r.pip_pz);
        B("pim_sim_px", &r.pim_px); B("pim_sim_py", &r.pim_py); B("pim_sim_pz", &r.pim_pz);
        B("ep_sim_p", &r.ep_p); B("pip_sim_p", &r.pip_p);
        B("pip_sim_vertexx", &r.pip_vx); B("pip_sim_vertexy", &r.pip_vy); B("pip_sim_vertexz", &r.pip_vz);
        B("ep_sim_vertexx", &r.ep_vx); B("ep_sim_vertexy", &r.ep_vy); B("ep_sim_vertexz", &r.ep_vz);
        B("pip_sim_geninfo", &r.pip_g0); B("pip_sim_geninfo1", &r.pip_g1); B("pip_sim_geninfo2", &r.pip_g2);
        B("neutr_mult", &r.nmult);
        for (int k = 0; k < kNN; ++k) {
            const int i = k + 1;
            B(Form("neutr_pid_%d", i), &r.pid[k]); B(Form("neutr_beta_%d", i), &r.beta[k]);
            B(Form("neutr_cluster_energy_%d", i), &r.ecl[k]); B(Form("neutr_cluster_theta_%d", i), &r.thc[k]);
            B(Form("neutr_cluster_phi_%d", i), &r.phc[k]); B(Form("neutr_cluster_sector_%d", i), &r.secc[k]);
            B(Form("neutr_cluster_cell_%d", i), &r.cellc[k]);
            B(Form("neutr_sim_id_%d", i), &r.sid[k]); B(Form("neutr_sim_parentpid_%d", i), &r.sppid[k]);
            B(Form("neutr_sim_geninfo_%d", i), &r.sg0[k]); B(Form("neutr_sim_geninfo1_%d", i), &r.sg1[k]);
            B(Form("neutr_sim_geninfo2_%d", i), &r.sg2[k]);
            B(Form("neutr_sim_px_%d", i), &r.spx[k]); B(Form("neutr_sim_py_%d", i), &r.spy[k]); B(Form("neutr_sim_pz_%d", i), &r.spz[k]);
            B(Form("neutr_sim_vx_%d", i), &r.svx[k]); B(Form("neutr_sim_vy_%d", i), &r.svy[k]); B(Form("neutr_sim_vz_%d", i), &r.svz[k]);
        }
        return true;
    }
    bool Get(Long64_t e) { return t && e >= 0 && e < t->GetEntries() && t->GetEntry(e) > 0; }
};

// ============================================================================
// the per-row record written for ecalacc_ana.py
// ============================================================================
const int kNEb = 30;            // photon energy bins 0-600 MeV, 20 MeV
const double kEbW = 20.;
const int kNThb = 45;           // photon polar angle 0-90 deg, 2 deg
const int kNThr = 6;            // thresholds on the TRUE photon energy for the geometric fraction
const double kThr[kNThr] = {0., 60., 80., 100., 110., 130.};
const int kNLam = 4;            // calibration factors lambda for the recomputed tag
const double kLam[kNLam] = {1.0, 1.10, 1.20, 1.285};   // 1.285 = (0.1459 / 0.1287)^2: the data pi0 position
struct Rec {
    // identification and weights
    double set, cls, pre, s0, w, m4p, mpp, mee, m4true, dist, grp;
    // conditional completion (simulation |M|^2 = 1): fractions of the circle weight
    double ok, sumf, fGeo, fGeoStrict, fGeo6, fGeoThr[kNThr], fGeoKloe, fGeoQed, fGeoFlat, meanQed, eMean, eMeanAcc, thMean;
    double eHist[kNEb], eHistAcc[kNEb], thHist[kNThb];
    // truth photon found by the ECAL simulation
    double tMatch, tPass, tE, tTh, tPhi, tEcl, tGeo, tRes1, tRes2, tCosRes, tSec, tQed, tConv;
    // compact tag and recomputed tag
    double tagC, tagTrueC, m5cC, tagL[kNLam], tagTrueL[kNLam], m5cL[kNLam], meegL[kNLam];
};
const int kRecN = sizeof(Rec) / sizeof(double);

void Run(const char* etaFile, const char* outTag) {
    TStopwatch swAll;
    gSystem->mkdir(kOutDir, kTRUE);
    gOut = fopen(Form("%s/ecalacc_fill%s.txt", kOutDir, outTag), "w");
    P("ecalacc_fill.C -- ECAL photon of the eta bump (package ecalacc), %s\n", TDatime().AsString());
    FILE* fh = fopen(Form("%s/hist_data.txt", kOutDir), "w");

    // =========================================================================================================
    // A. DATA
    // =========================================================================================================
    const int nSelT = 2;                         // 0 pre, 1 S0
    const int selId[nSelT] = {kSelPre, kSelS0};
    std::vector<scan::Spectra> dTag;             // [sel * kNTs + slice]
    {
        scan::Data d;
        if (!scan::Load(d)) return;
        std::vector<scan::Sel> vt;
        // data selections: 0 pre, 1 S0, 2 pre with evtStart == 1, 3 pre with evtStart != 1 (the ECAL beta of the
        // neutral candidate is measured against the START time in the data)
        for (int s = 0; s < 4; ++s)
            for (int j = 0; j < kNTs; ++j) {
                const int sid = s == 1 ? kSelS0 : kSelPre;
                vt.push_back([=](const scan::Row& r) {
                    if (!Tagged(r) || !PassSel(sid, r, r.cor_ok == 1, r.dz_pi, r.mm4)) return false;
                    if (s == 2 && r.evtStart != 1) return false;
                    if (s == 3 && r.evtStart == 1) return false;
                    return InTs(j, r.mpipi);
                });
            }
        {
            double nAll[2] = {0, 0}, nB[2] = {0, 0};
            for (const auto& r : d.s[0]) {
                if (r.cor_ok != 1 || !PassSel(kSelPre, r, true, r.dz_pi, r.mm4)) continue;
                const int k = r.evtStart == 1 ? 0 : 1;
                nAll[k] += 1; if (r.m4 >= kBlo && r.m4 < kBhi) nB[k] += 1;
            }
            P("A. data pre e+e- rows: evtStart == 1 %.0f, != 1 %.0f (share %.4f); in the bump window B: %.0f / %.0f (share %.4f)\n",
              nAll[0], nAll[1], nAll[0] / (nAll[0] + nAll[1]), nB[0], nB[1], nB[0] / (nB[0] + nB[1]));
        }
        scan::FillOpt o;
        o.nbin = kTagNbin; o.lo = kTagHLo; o.hi = kTagHHi; o.reb = kTagReb;
        o.x = [](const scan::Row& r) { return M5c(r); };
        dTag = scan::FillMany(d, vt, o);
        // photon maps and distance-to-edge of the tagged rows ("pre"): e+e- and like-sign separately
        const char* hs[3] = {"us", "pp", "mm"};
        const double wLo = kTagMean - 2 * kTagSigma, wHi = kTagMean + 2 * kTagSigma;     // m5c peak window
        for (int k = 0; k < 3; ++k) {
            H1 hth, hph, hE, hEw, hEsb, hD, hDsb, hMeeg, hthw, hthsb;
            hth.Init(Form("data_%s_gtheta", hs[k]), 45, 0, 90); hph.Init(Form("data_%s_gphi", hs[k]), 72, 0, 360);
            hE.Init(Form("data_%s_gE", hs[k]), 40, 0, 800); hEw.Init(Form("data_%s_gE_win", hs[k]), 40, 0, 800);
            hEsb.Init(Form("data_%s_gE_sb", hs[k]), 40, 0, 800);
            hthw.Init(Form("data_%s_gtheta_win", hs[k]), 45, 0, 90); hthsb.Init(Form("data_%s_gtheta_sb", hs[k]), 45, 0, 90);
            hD.Init(Form("data_%s_dist_win", hs[k]), 40, -0.20, 0.04); hDsb.Init(Form("data_%s_dist_sb", hs[k]), 40, -0.20, 0.04);
            hMeeg.Init(Form("data_%s_meeg", hs[k]), 100, 0.0, 0.30);
            for (const auto& r : d.s[k]) {
                if (r.cor_ok != 1 || !PassSel(kSelPre, r, true, r.dz_pi, r.mm4)) continue;
                if (r.g_idx > 0 && r.mee < kMeeMax) hMeeg.Fill(r.m_eeg, 1.);
                if (!Tagged(r)) continue;
                hth.Fill(r.g_theta, 1.); hph.Fill(r.g_phi, 1.); hE.Fill(r.g_e, 1.);
                const double m5c = M5c(r);
                const bool inW = m5c >= wLo && m5c < wHi;
                const bool inSB = (m5c >= wLo - 0.030 && m5c < wLo - 0.0124) || (m5c >= wHi + 0.0124 && m5c < wHi + 0.030);
                const double de = DistEdge(r.m4, r.mpipi, r.mee);
                if (inW) { hD.Fill(de, 1.); hEw.Fill(r.g_e, 1.); hthw.Fill(r.g_theta, 1.); }
                if (inSB) { hDsb.Fill(de, 1.); hEsb.Fill(r.g_e, 1.); hthsb.Fill(r.g_theta, 1.); }
            }
            for (H1* h : {&hth, &hph, &hE, &hEw, &hEsb, &hD, &hDsb, &hMeeg, &hthw, &hthsb}) h->Write(fh);
        }
        P("A. data loaded and tag spectra filled (pre, S0; %d M(pipi) slices)\n", kNTs);
    }
    fclose(fh);

    // =========================================================================================================
    // B + C. SIMULATION
    // =========================================================================================================
    Geo geo;
    if (!geo.Load(kCells)) return;
    EtaDens dens;
    if (!dens.Load(etaFile && *etaFile ? etaFile : kEtaLab)) return;

    simphys::Init();
    simphys::Install();
    simphys::LoadAngles();
    simeng::SimData* D = new simeng::SimData();
    if (!simeng::Load(*D, simeng::kCompact, simeng::kMaskUS, -1e9, 1e9, 1)) return;
    const int trees[3] = {simeng::kStdUS, simeng::kLepUS, simeng::kEtaUS};

    // ---- B. emulated tag spectra of T1 + T2 and of the truth C2 + C3 (pre, S0, slices)
    std::vector<std::vector<double>> sT12(nSelT * kNTs, std::vector<double>(kTagNbin, 0.)), sC23 = sT12, sC2 = sT12, sC3 = sT12;
    std::vector<double> c23All(nSelT * kNTs, 0.), c2All = c23All, c3All = c23All, c23B = c23All;   // all C2/C3 rows (any photon)
    long nC23set[4] = {0, 0, 0, 0};
    struct Want { int tree; size_t idx; };
    std::vector<Want> want;
    for (int t = 0; t < 3; ++t) {
        const std::vector<simeng::Row>& v = D->t[trees[t]];
        for (size_t i = 0; i < v.size(); ++i) {
            const simeng::Row& r = v[i];
            if (r.ownA != 1) continue;
            const bool pre = PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4);
            if (!pre) continue;
            const bool s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4);
            int cn[6]; const int ncn = CompsOf(r, false, cn);
            bool isT12 = false; for (int q = 0; q < ncn; ++q) if (cn[q] == kT12) isT12 = true;
            if (!isT12) continue;
            const double w = r.w * simphys::RowFactor(r, kPhysNom) * kSimToRaw;
            if (w == 0.) continue;
            const bool isC23 = r.tmpl == 1 && (r.cls == 2 || r.cls == 3);
            const double mpp = kSPiPi * r.mpipi;
            const double m4p = M4sim(r.m4, r.m4_true, kS, kR);
            if (isC23) { want.push_back({trees[t], i}); ++nC23set[r.set]; }
            for (int s = 0; s < nSelT; ++s) {
                if (s == 1 && !s0) continue;
                for (int j = 0; j < kNTs; ++j) {
                    if (!InTs(j, mpp)) continue;
                    const int k = s * kNTs + j;
                    if (isC23) {
                        c23All[k] += w; if (r.cls == 2) c2All[k] += w; else c3All[k] += w;
                        if (m4p >= kBlo && m4p < kBhi) c23B[k] += w;
                    }
                    if (!Tagged(r)) continue;
                    const int it = TagBin(M5c(r));
                    if (it < 0) continue;
                    sT12[k][it] += w;
                    if (isC23) { sC23[k][it] += w; (r.cls == 2 ? sC2 : sC3)[k][it] += w; }
                }
            }
        }
    }
    P("B. C2 + C3 rows (owned, pre): std %ld, lep %ld, eta %ld\n", nC23set[1], nC23set[2], nC23set[3]);

    // ---- tag fits: data and simulation
    FILE* ft = fopen(Form("%s/tag_fits.txt", kOutDir), "w");
    fprintf(ft, "# sel slice | data T eT chi2 ndf free_mean free_sigma free_T | sim T12 frozen T, free T, free mean, free sigma | truth C23 C2 C3 in 0.49-0.62 and all m5c | all C23 rows C2 C3 (any photon), C23 in B\n");
    P("\nB. ECAL TAG: frozen fit %.4f / %.1f MeV, pol%d, %.2f-%.2f (data: N - CB; simulation: T1 + T2 expected raw rows, UNVERIFIED)\n",
      kTagMean, 1e3 * kTagSigma, kTagOrder, kTagLo, kTagHi);
    for (int s = 0; s < nSelT; ++s)
        for (int j = 0; j < kNTs; ++j) {
            const int k = s * kNTs + j;
            const scan::PeakFit fd = scan::FitPeakMain(dTag[k], TagCfg());
            const scan::FreeFit ffd = scan::FitPeakFree(dTag[k], TagCfg());
            const scan::Spectra ss = SimSpectra(sT12[k].data(), kTagNbin, kTagHLo, kTagHHi, kTagReb);
            const scan::PeakFit fs = scan::FitPeakMain(ss, TagCfg());
            const scan::FreeFit ffs = scan::FitPeakFree(ss, TagCfg());
            double trW = 0, trA = 0, tr2 = 0, tr3 = 0;
            for (int i = 0; i < kTagNbin; ++i) {
                const double x = kTagHLo + (i + 0.5) * (kTagHHi - kTagHLo) / kTagNbin;
                trA += sC23[k][i];
                if (x >= kTagLo && x < kTagHi) { trW += sC23[k][i]; tr2 += sC2[k][i]; tr3 += sC3[k][i]; }
            }
            P("  %-3s M(pipi) %-9s | data T = %7.1f +- %5.1f (chi2 %5.1f/%d; free %.4f / %.1f MeV, T %.1f) | sim T = %7.1f (frozen) %7.1f (free %.4f / %.1f MeV) | truth C2+C3 %7.1f (C2 %.1f, C3 %.1f) in 0.49-0.62, %7.1f all | all C2+C3 rows %8.1f (C2 %.1f, C3 %.1f), in B %.1f\n",
              SelTag(selId[s]), kTsName[j], fd.S, fd.eS, fd.chi2, fd.ndf, ffd.mean, 1e3 * ffd.sigma, ffd.S, fs.S, ffs.S, ffs.mean, 1e3 * ffs.sigma,
              trW, tr2, tr3, trA, c23All[k], c2All[k], c3All[k], c23B[k]);
            fprintf(ft, "TAG %d %d %.4f %.4f %.3f %d %.6f %.6f %.4f %.4f %.4f %.6f %.6f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", s, j, fd.S, fd.eS, fd.chi2,
                    fd.ndf, ffd.mean, ffd.sigma, ffd.S, fs.S, ffs.S, ffs.mean, ffs.sigma, trW, tr2, tr3, trA, c23All[k], c2All[k], c3All[k], c23B[k]);
            if (j == 0) {   // spectra for the figure
                fprintf(ft, "SPEC %d data_N", s); for (int i = 0; i < kTagNbin; ++i) fprintf(ft, " %.4f", dTag[k].n[0][i]); fprintf(ft, "\n");
                fprintf(ft, "SPEC %d data_CB", s); for (int i = 0; i < kTagNbin; ++i) fprintf(ft, " %.4f", dTag[k].cb[i]); fprintf(ft, "\n");
                fprintf(ft, "SPEC %d data_sigVar", s); for (int i = 0; i < kTagNbin; ++i) fprintf(ft, " %.4f", dTag[k].sigVar[i]); fprintf(ft, "\n");
                fprintf(ft, "SPEC %d sim_T12", s); for (int i = 0; i < kTagNbin; ++i) fprintf(ft, " %.5f", sT12[k][i]); fprintf(ft, "\n");
                fprintf(ft, "SPEC %d sim_C2", s); for (int i = 0; i < kTagNbin; ++i) fprintf(ft, " %.5f", sC2[k][i]); fprintf(ft, "\n");
                fprintf(ft, "SPEC %d sim_C3", s); for (int i = 0; i < kTagNbin; ++i) fprintf(ft, " %.5f", sC3[k][i]); fprintf(ft, "\n");
                fprintf(ft, "FITPAR %d data", s); for (int q = 0; q < fd.np; ++q) fprintf(ft, " %.8g", fd.par[q]); fprintf(ft, " x0 %.8g dx %.8g bw %.8g\n", fd.x0, fd.dx, fd.bw);
                fprintf(ft, "FITPAR %d sim", s); for (int q = 0; q < fs.np; ++q) fprintf(ft, " %.8g", fs.par[q]); fprintf(ft, " x0 %.8g dx %.8g bw %.8g\n", fs.x0, fs.dx, fs.bw);
            }
        }
    for (int s = 2; s < 4; ++s)
        for (int j = 0; j < kNTs; ++j) {
            const int k = s * kNTs + j;
            const scan::PeakFit fd = scan::FitPeakMain(dTag[k], TagCfg());
            P("  data pre, evtStart %s, M(pipi) %-9s | T = %7.1f +- %5.1f (chi2 %5.1f/%d)\n", s == 2 ? "== 1" : "!= 1", kTsName[j], fd.S, fd.eS, fd.chi2, fd.ndf);
            fprintf(ft, "TAGSTART %d %d %.4f %.4f %.3f %d\n", s, j, fd.S, fd.eS, fd.chi2, fd.ndf);
        }
    fclose(ft);

    // ---- C. truth of the C2 / C3 rows
    std::sort(want.begin(), want.end(), [&](const Want& a, const Want& b) {
        const simeng::Row& ra = D->t[a.tree][a.idx]; const simeng::Row& rb = D->t[b.tree][b.idx];
        if (ra.set != rb.set) return ra.set < rb.set;
        return ra.pat_entry < rb.pat_entry; });
    PatReader pr[4];
    std::vector<PatRow> pat(want.size());
    std::vector<char> patOk(want.size(), 0);
    long nBadMatch = 0;
    for (size_t q = 0; q < want.size(); ++q) {
        const simeng::Row& r = D->t[want[q].tree][want[q].idx];
        if (!pr[r.set].t && !pr[r.set].Open(kPatFile[r.set])) return;
        if (!pr[r.set].Get(r.pat_entry)) continue;
        pat[q] = pr[r.set].r;
        // identity check: true momenta of the compact and of the PAT entry
        if (std::fabs(pat[q].ep_p - r.ep_sim_p) > 1e-3f * std::max(1.f, std::fabs(r.ep_sim_p)) ||
            std::fabs(pat[q].pip_p - r.pip_sim_p) > 1e-3f * std::max(1.f, std::fabs(r.pip_sim_p))) { ++nBadMatch; continue; }
        patOk[q] = 1;
    }
    P("C. PAT truth read for %zu rows; identity check failed for %ld\n", want.size(), nBadMatch);

    const double mE = scan::kMe * 1000., mPi = 139.56995;
    auto vec4 = [](double px, double py, double pz, double m) { TLorentzVector v; v.SetXYZM(px, py, pz, m); return v; };
    // pass 1: the masses of the generator (HGeant) from the matched photons
    auto matchPhoton = [&](const PatRow& p) -> int {
        int best = -1;
        const int n = std::min(kNN, (int)std::lround(p.nmult));
        for (int k = 0; k < n; ++k) {
            if (std::lround(p.sid[k]) != 1 || std::lround(p.sppid[k]) != 7) continue;
            if (p.sg0[k] != p.pip_g0 || p.sg1[k] != p.pip_g1 || p.sg2[k] != p.pip_g2) continue;
            const double dv = std::sqrt(std::pow(p.svx[k] - p.pip_vx, 2) + std::pow(p.svy[k] - p.pip_vy, 2) + std::pow(p.svz[k] - p.pip_vz, 2));
            if (dv > 0.01) continue;
            if (best < 0 || p.ecl[k] > p.ecl[best]) best = k;
        }
        return best;
    };
    std::vector<double> mEtaV, mPi0V;
    for (size_t q = 0; q < want.size(); ++q) {
        if (!patOk[q]) continue;
        const simeng::Row& r = D->t[want[q].tree][want[q].idx];
        const PatRow& p = pat[q];
        const int k = matchPhoton(p);
        if (k < 0) continue;
        const TLorentzVector ep = vec4(p.ep_px, p.ep_py, p.ep_pz, mE), em = vec4(p.em_px, p.em_py, p.em_pz, mE);
        const TLorentzVector pp = vec4(p.pip_px, p.pip_py, p.pip_pz, mPi), pm = vec4(p.pim_px, p.pim_py, p.pim_pz, mPi);
        const TLorentzVector g = vec4(p.spx[k], p.spy[k], p.spz[k], 0.);
        mEtaV.push_back((ep + em + pp + pm + g).M());
        if (r.cls == 2) mPi0V.push_back((ep + em + g).M());
    }
    auto median = [](std::vector<double> v) { if (v.empty()) return -1.; std::sort(v.begin(), v.end()); return v[v.size() / 2]; };
    double mEtaG = median(mEtaV), mPi0G = median(mPi0V);
    P("C. generator masses from %zu (%zu) matched photons: m_eta %.4f, m_pi0 %.4f MeV/c2 (PDG 547.862, 134.977)\n", mEtaV.size(), mPi0V.size(), mEtaG, mPi0G);
    if (mEtaV.size() < 20) mEtaG = 547.862;
    if (mPi0V.size() < 20) mPi0G = 134.9768;
    {
        std::vector<double> a = mEtaV; std::sort(a.begin(), a.end());
        if (a.size() > 10) P("   m_eta quantiles 5/50/95 %%: %.3f %.3f %.3f\n", a[a.size() / 20], a[a.size() / 2], a[a.size() * 19 / 20]);
        std::vector<double> b = mPi0V; std::sort(b.begin(), b.end());
        if (b.size() > 10) P("   m_pi0 quantiles 5/50/95 %%: %.3f %.3f %.3f\n", b[b.size() / 20], b[b.size() / 2], b[b.size() * 19 / 20]);
    }

    // pass 2: the conditional completion and the records
    const int kNPhi = 720;
    std::mt19937_64 rng(0xEC41ACCULL);
    std::uniform_real_distribution<double> uni(0., 1.);
    std::vector<Rec> recs;
    long nNoSol = 0, nZeroF = 0;
    H1 hRes1, hRes2, hCosRes, hEclTrue;
    hRes1.Init("closure_res_pi0", 100, -5, 5); hRes2.Init("closure_res_eta", 100, -5, 5); hCosRes.Init("closure_cosalpha", 100, -0.01, 0.01);
    hEclTrue.Init("ecl_over_etrue", 100, 0., 2.);
    double qedAll[3] = {0, 0, 0}, qedW[3] = {0, 0, 0};    // C2 angle factor of all rows: sum w <W_qed>, sum w; [all, E1, E2]
    for (size_t q = 0; q < want.size(); ++q) {
        if (!patOk[q]) continue;
        const simeng::Row& r = D->t[want[q].tree][want[q].idx];
        const PatRow& p = pat[q];
        Rec R; std::fill((double*)&R, (double*)&R + kRecN, 0.);
        R.set = r.set; R.cls = r.cls; R.grp = r.grp;
        R.pre = 1; R.s0 = PassSel(kSelS0, r, r.cor_ok == 1, 0.f, r.mm4) ? 1 : 0;
        R.w = r.w * simphys::RowFactor(r, kPhysNom) * kSimToRaw;
        R.m4p = M4sim(r.m4, r.m4_true, kS, kR); R.mpp = kSPiPi * r.mpipi; R.mee = r.mee; R.m4true = r.m4_true;
        R.dist = DistEdge(R.m4p, R.mpp, R.mee);
        const TLorentzVector ep = vec4(p.ep_px, p.ep_py, p.ep_pz, mE), em = vec4(p.em_px, p.em_py, p.em_pz, mE);
        const TLorentzVector pp = vec4(p.pip_px, p.pip_py, p.pip_pz, mPi), pm = vec4(p.pim_px, p.pim_py, p.pim_pz, mPi);
        const TLorentzVector Q = ep + em, Pv = Q + pp + pm;
        const TVector3 V(p.pip_vx, p.pip_vy, p.pip_vz);
        const double Mq = Q.M();
        // ---- conditional completion in the rest frame of Q
        const TVector3 bQ = Q.BoostVector();
        TLorentzVector Ps = Pv; Ps.Boost(-bQ);
        TLorentzVector eps = ep; eps.Boost(-bQ);
        const double w0 = (mPi0G * mPi0G - Mq * Mq) / (2. * Mq);
        const double c2 = 0.5 * (mEtaG * mEtaG - Pv.M2());
        const double cosA = (Ps.E() * w0 - c2) / (Ps.Vect().Mag() * w0);
        const double aQ = 4. * mE * mE / (Mq * Mq);
        if (!(std::fabs(cosA) <= 1.)) { ++nNoSol; R.ok = 0; }
        else {
            R.ok = 1;
            const double sinA = std::sqrt(std::max(0., 1. - cosA * cosA));
            const TVector3 e1 = Ps.Vect().Unit();
            TVector3 tmp = std::fabs(e1.Z()) < 0.9 ? TVector3(0, 0, 1) : TVector3(1, 0, 0);
            const TVector3 e2 = tmp.Cross(e1).Unit(), e3 = e1.Cross(e2);
            const double off = uni(rng);
            double sf = 0, sGeo = 0, sGeoS = 0, sGeo6 = 0, sThr[kNThr] = {0}, sK = 0, sKGeo = 0, sQ = 0, sQGeo = 0, sE = 0, sEacc = 0, sTh = 0;
            double nFlat = 0, nFlatGeo = 0;   // variant: flat weight along the circle (no SMASH density)
            for (int i = 0; i < kNPhi; ++i) {
                const double ph = 2. * TMath::Pi() * (i + off) / kNPhi;
                const TVector3 ks = w0 * (cosA * e1 + sinA * std::cos(ph) * e2 + sinA * std::sin(ph) * e3);
                TLorentzVector k(ks, w0);
                k.Boost(bQ);
                const TLorentzVector pEta = Pv + k;
                const double f = dens(pEta * 1e-3);
                const double Eg = k.E();
                const TVector3 dir = k.Vect().Unit();
                int sec6 = -1;
                const int hit = geo.Hit(V, dir, false, &sec6);
                const bool in6 = hit >= 0;
                const bool in = in6 && geo.cells[hit].status == 1;
                nFlat += 1; if (in) nFlatGeo += 1;
                if (!(f > 0)) continue;
                const int hitS = geo.Hit(V, dir, true);
                const bool inS = hitS >= 0 && geo.cells[hitS].status == 1;
                // KLOE X dependence: energies in the eta rest frame
                double wK = 1.;
                {
                    const TVector3 bE = pEta.BoostVector();
                    TLorentzVector a = pp, b = pm, c = Q + k; a.Boost(-bE); b.Boost(-bE); c.Boost(-bE);
                    const double Qe = mEtaG - 2 * mPi - mPi0G;
                    const double X = std::sqrt(3.) * ((a.E() - mPi) - (b.E() - mPi)) / Qe, Y = 3. * (c.E() - mPi0G) / Qe - 1.;
                    wK = FKloe(X, Y);
                }
                double wQ = 1.;
                if (r.cls == 2) {
                    const double ct = eps.Vect().Unit().Dot(ks.Unit());
                    wQ = (1. + ct * ct + aQ * (1. - ct * ct)) * 3. / (4. + 2. * aQ);
                }
                sf += f; sK += f * wK; sQ += f * wQ;
                sE += f * Eg; sTh += f * dir.Theta() * TMath::RadToDeg();
                const int ie = (int)(Eg / kEbW); if (ie >= 0 && ie < kNEb) R.eHist[ie] += f;
                const int ith = (int)(dir.Theta() * TMath::RadToDeg() / 2.); if (ith >= 0 && ith < kNThb) R.thHist[ith] += f;
                if (in6) sGeo6 += f;
                if (inS) sGeoS += f;
                if (in) {
                    sGeo += f; sKGeo += f * wK; sQGeo += f * wQ; sEacc += f * Eg;
                    for (int t = 0; t < kNThr; ++t) if (Eg >= kThr[t]) sThr[t] += f;
                    if (ie >= 0 && ie < kNEb) R.eHistAcc[ie] += f;
                }
            }
            if (!(sf > 0)) { ++nZeroF; R.ok = 0; }
            else {
                R.sumf = sf; R.fGeo = sGeo / sf; R.fGeoStrict = sGeoS / sf; R.fGeo6 = sGeo6 / sf;
                for (int t = 0; t < kNThr; ++t) R.fGeoThr[t] = sThr[t] / sf;
                R.fGeoKloe = sK > 0 ? sKGeo / sK : 0.; R.fGeoQed = sQ > 0 ? sQGeo / sQ : 0.; R.meanQed = sQ / sf;
                R.fGeoFlat = nFlat > 0 ? nFlatGeo / nFlat : 0.;
                R.eMean = sE / sf; R.eMeanAcc = sGeo > 0 ? sEacc / sGeo : 0.; R.thMean = sTh / sf;
                for (int i = 0; i < kNEb; ++i) { R.eHist[i] /= sf; R.eHistAcc[i] /= sf; }
                for (int i = 0; i < kNThb; ++i) R.thHist[i] /= sf;
                if (r.cls == 2) {
                    qedAll[0] += R.w * R.meanQed; qedW[0] += R.w;
                    const int e = R.mee < 0.010 ? 1 : 2; qedAll[e] += R.w * R.meanQed; qedW[e] += R.w;
                }
            }
        }
        // ---- the photon found by the ECAL simulation
        const int km = matchPhoton(p);
        if (km >= 0) {
            R.tMatch = 1;
            R.tPass = (std::lround(p.pid[km]) == 1 && p.beta[km] > 0.8 && p.beta[km] < 1.2 && p.ecl[km] > 100.) ? 1 : 0;
            const TLorentzVector g = vec4(p.spx[km], p.spy[km], p.spz[km], 0.);
            R.tE = g.E(); R.tTh = g.Theta() * TMath::RadToDeg(); R.tPhi = g.Phi() * TMath::RadToDeg(); if (R.tPhi < 0) R.tPhi += 360.;
            R.tEcl = p.ecl[km]; R.tSec = p.secc[km];
            const int hit = geo.Hit(V, g.Vect().Unit());
            R.tGeo = (hit >= 0 && geo.cells[hit].status == 1) ? 1 : (hit >= 0 ? 0.5 : 0);
            R.tRes1 = (Q + g).M() - mPi0G; R.tRes2 = (Pv + g).M() - mEtaG;
            if (R.ok) {
                TLorentzVector gs = g; gs.Boost(-bQ);
                R.tCosRes = gs.Vect().Unit().Dot(Ps.Vect().Unit()) - cosA;
            }
            if (r.cls == 2) {
                TLorentzVector gs = g; gs.Boost(-bQ);
                const double ct = eps.Vect().Unit().Dot(gs.Vect().Unit());
                R.tQed = (1. + ct * ct + aQ * (1. - ct * ct)) * 3. / (4. + 2. * aQ);
            }
            hRes1.Fill(R.tRes1, 1.); hRes2.Fill(R.tRes2, 1.); hCosRes.Fill(R.tCosRes, 1.);
            if (R.tE > 0) hEclTrue.Fill(R.tEcl / R.tE, 1.);
        }
        // a cluster made by an e+- of the converted photon (GEANT id 2/3, parent photon, the codes of the pions)
        {
            const int n = std::min(kNN, (int)std::lround(p.nmult));
            for (int k = 0; k < n; ++k)
                if ((std::lround(p.sid[k]) == 2 || std::lround(p.sid[k]) == 3) && std::lround(p.sppid[k]) == 1 &&
                    p.sg0[k] == p.pip_g0 && p.sg1[k] == p.pip_g1 && p.sg2[k] == p.pip_g2) R.tConv = 1;
        }
        // ---- tag of the compact row and recomputed tags with the cluster energies x lambda
        R.tagC = (Tagged(r) && R.pre) ? 1 : 0;
        R.m5cC = R.tagC ? M5c(r) : -1.;
        if (R.tagC && km >= 0) R.tagTrueC = std::fabs(r.g_e - p.ecl[km]) < 1e-3f * std::max(1.f, p.ecl[km]) ? 1 : 0;
        {
            const TLorentzVector vEp = simeng::P4(r.ep_p, r.ep_theta, r.ep_phi, mE), vEm = simeng::P4(r.em_p, r.em_theta, r.em_phi, mE);
            const TLorentzVector vPip = simeng::P4(r.pip_p, r.pip_theta, r.pip_phi, mPi), vPim = simeng::P4(r.pim_p, r.pim_theta, r.pim_phi, mPi);
            const TLorentzVector vEE = vEp + vEm, v4 = vEE + vPip + vPim;
            const int n = std::min(kNN, (int)std::lround(p.nmult));
            for (int l = 0; l < kNLam; ++l) {
                double best = 1e9; int kb = -1; double meeg = -1, m5 = -1;
                // the clusters passing ecal_quality with the scaled energy; main.cc keeps the three most energetic
                std::vector<int> pass;
                for (int k = 0; k < n; ++k)
                    if (std::lround(p.pid[k]) == 1 && p.beta[k] > 0.8 && p.beta[k] < 1.2 && p.ecl[k] * kLam[l] > 100.) pass.push_back(k);
                std::stable_sort(pass.begin(), pass.end(), [&](int a, int b) { return p.ecl[a] > p.ecl[b]; });
                if (pass.size() > 3) pass.resize(3);
                for (int k : pass) {
                    const double E = p.ecl[k] * kLam[l];
                    const TLorentzVector vg = simeng::P4(E, p.thc[k], p.phc[k], 0.);
                    const double m = (vEE + vg).M() / 1000.;
                    if (std::fabs(m - 0.140) < best) { best = std::fabs(m - 0.140); kb = k; meeg = m; m5 = (v4 + vg).M() / 1000.; }
                }
                R.meegL[l] = meeg;
                const bool tg = kb >= 0 && r.mee < kMeeMax && meeg >= kPi0Lo && meeg < kPi0Hi;
                R.tagL[l] = tg ? 1 : 0;
                R.m5cL[l] = tg ? M5cOf(r.m4, r.mee, meeg, m5) : -1.;
                R.tagTrueL[l] = (tg && kb == km) ? 1 : 0;
            }
        }
        recs.push_back(R);
    }
    P("C. records: %zu; no solution on the circle %ld; zero SMASH density on the whole circle %ld\n", recs.size(), nNoSol, nZeroF);
    P("C. C2 angle factor of ALL C2 rows (pre) from the circle average of the QED gamma* weight: all %.4f, E1 %.4f, E2 %.4f\n",
      qedW[0] > 0 ? qedAll[0] / qedW[0] : 0., qedW[1] > 0 ? qedAll[1] / qedW[1] : 0., qedW[2] > 0 ? qedAll[2] / qedW[2] : 0.);
    // check of the recomputed tag (lambda = 1) against the compact
    {
        long nT = 0, nAgree = 0;
        for (const auto& R : recs) { nT += (R.tagC > 0 || R.tagL[0] > 0); nAgree += (R.tagC > 0) == (R.tagL[0] > 0) && (R.tagC == 0 || std::fabs(R.m5cC - R.m5cL[0]) < 2e-5); }
        P("C. recomputed tag (lambda = 1) vs compact: %ld rows tagged by either, %ld of all %zu records agree\n", nT, nAgree, recs.size());
    }
    // ---- write the records
    {
        FILE* fb = fopen(Form("%s/rows_c23.bin", kOutDir), "wb");
        for (const auto& R : recs) fwrite(&R, sizeof(double), kRecN, fb);
        fclose(fb);
        FILE* fl = fopen(Form("%s/rows_c23_layout.txt", kOutDir), "w");
        fprintf(fl, "nrec %zu ndouble %d\n", recs.size(), kRecN);
        fprintf(fl, "fields set cls pre s0 w m4p mpp mee m4true dist grp ok sumf fGeo fGeoStrict fGeo6");
        for (int t = 0; t < kNThr; ++t) fprintf(fl, " fGeoThr%d", (int)kThr[t]);
        fprintf(fl, " fGeoKloe fGeoQed fGeoFlat meanQed eMean eMeanAcc thMean");
        for (int i = 0; i < kNEb; ++i) fprintf(fl, " eHist%d", i);
        for (int i = 0; i < kNEb; ++i) fprintf(fl, " eHistAcc%d", i);
        for (int i = 0; i < kNThb; ++i) fprintf(fl, " thHist%d", i);
        fprintf(fl, " tMatch tPass tE tTh tPhi tEcl tGeo tRes1 tRes2 tCosRes tSec tQed tConv tagC tagTrueC m5cC");
        for (int l = 0; l < kNLam; ++l) fprintf(fl, " tagL%d", l);
        for (int l = 0; l < kNLam; ++l) fprintf(fl, " tagTrueL%d", l);
        for (int l = 0; l < kNLam; ++l) fprintf(fl, " m5cL%d", l);
        for (int l = 0; l < kNLam; ++l) fprintf(fl, " meegL%d", l);
        fprintf(fl, "\nlambda"); for (int l = 0; l < kNLam; ++l) fprintf(fl, " %.5f", kLam[l]);
        fprintf(fl, "\nthresholds"); for (int t = 0; t < kNThr; ++t) fprintf(fl, " %.1f", kThr[t]);
        fprintf(fl, "\nebins %d %.1f thbins %d 2.0\nmasses %.5f %.5f\n", kNEb, kEbW, kNThb, mEtaG, mPi0G);
        fclose(fl);
    }
    // ---- simulated photon maps of the tagged T12 rows (pre) for the comparison with the data
    {
        FILE* fs = fopen(Form("%s/hist_sim.txt", kOutDir), "w");
        H1 hth, hph, hE, hMeeg, hthC, hEC, hD;
        hth.Init("sim_T12_gtheta", 45, 0, 90); hph.Init("sim_T12_gphi", 72, 0, 360); hE.Init("sim_T12_gE", 40, 0, 800);
        hMeeg.Init("sim_C23_meeg", 100, 0., 0.30); hthC.Init("sim_C23tag_gtheta", 45, 0, 90); hEC.Init("sim_C23tag_gE", 40, 0, 800);
        hD.Init("sim_C23tag_dist", 40, -0.20, 0.04);
        for (int t = 0; t < 3; ++t)
            for (const auto& r : D->t[trees[t]]) {
                if (r.ownA != 1 || !PassSel(kSelPre, r, r.cor_ok == 1, 0.f, r.mm4)) continue;
                int cn[6]; const int ncn = CompsOf(r, false, cn);
                bool isT12 = false; for (int q = 0; q < ncn; ++q) if (cn[q] == kT12) isT12 = true;
                if (!isT12) continue;
                const double w = r.w * simphys::RowFactor(r, kPhysNom) * kSimToRaw;
                const bool isC23 = r.tmpl == 1 && (r.cls == 2 || r.cls == 3);
                if (isC23 && r.g_idx > 0 && r.mee < kMeeMax) hMeeg.Fill(r.m_eeg, w);
                if (!Tagged(r)) continue;
                hth.Fill(r.g_theta, w); hph.Fill(r.g_phi, w); hE.Fill(r.g_e, w);
                if (isC23) {
                    hthC.Fill(r.g_theta, w); hEC.Fill(r.g_e, w);
                    const double m5c = M5c(r);
                    if (m5c >= kTagLo && m5c < kTagHi) hD.Fill(DistEdge(M4sim(r.m4, r.m4_true, kS, kR), kSPiPi * r.mpipi, r.mee), w);
                }
            }
        for (H1* h : {&hth, &hph, &hE, &hMeeg, &hthC, &hEC, &hD, &hRes1, &hRes2, &hCosRes, &hEclTrue}) h->Write(fs);
        fclose(fs);
    }
    P("done in %.1f s\n", swAll.RealTime());
    fclose(gOut); gOut = nullptr;
}

}  // namespace ecalacc

// etaFile: SMASH eta momenta (default kEtaLab); outTag: suffix of ecalacc_fill<tag>.txt (test runs)
void ecalacc_fill(const char* etaFile = "", const char* outTag = "") { ecalacc::Run(etaFile, outTag); }

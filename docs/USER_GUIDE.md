# FAT Framework - User Guide

A practical guide for running and customizing physics analyses with the FAT framework.

---

## Table of Contents

1. [Quick Start](#1-quick-start)
2. [Configuration](#2-configuration)
3. [The Cut System](#3-the-cut-system)
4. [Adding Histograms](#4-adding-histograms)
5. [Working with Ntuples](#5-working-with-ntuples)
6. [Common Analysis Patterns](#6-common-analysis-patterns)

---

## 1. Quick Start

### 1.1 First Run

```bash
# Build the analysis
make

# Run with default config
./ana

# Run with custom config
./ana my_config.json
```

### 1.2 Files You Will Edit

| File | Purpose | When to Edit |
|------|---------|--------------|
| `config.json` | Input/output paths, beam energy | Every analysis |
| `src/setup_cuts.h` | Cut definitions | When changing selections |
| `src/setup_histograms.h` | Histogram definitions | When adding observables |
| `main.cc` | Physics logic in `processEvent()` | When changing analysis |

### 1.3 Minimal Workflow

1. Edit `config.json` to point to your input files
2. Define cuts in `src/setup_cuts.h`
3. Add histograms in `src/setup_histograms.h`
4. Implement physics in `main.cc:processEvent()`
5. Build and run: `make && ./ana`

---

## 2. Configuration

### 2.1 Essential config.json Parameters

```json
{
  "input": {
    "source": "my_files.list",      // Input: .root file or .list of files
    "tree_name": "EpEm_ID",          // TTree name in ROOT file
    "start_event": 0,                // First event (0 = beginning)
    "max_events": -1                 // -1 = process all events
  },
  "output": {
    "filename": "output.root",       // Output ROOT file
    "option": "RECREATE"             // Overwrite existing file
  },
  "beam": {
    "kinetic_energy": 4500.0,        // Beam kinetic energy [MeV]
    "luminosity": 1.0                // Weight scaling factor
  },
  "kinematics": {
    "reconstructed": true,           // Load measured momenta
    "corrected": true,               // Load energy-loss corrected
    "simulated": false,              // Load MC truth (if available)
    "analysis_type": "reconstructed" // Which to use for physics
  },
  "ecal": {
    "ecal_1": true,                  // Enable ECAL photon slot 1
    "ecal_2": true,                  // Enable ECAL photon slot 2
    "ecal_3": true                   // Enable ECAL photon slot 3
  }
}
```

### 2.2 Input File Formats

**Single ROOT file:**
```json
"source": "data/run001.root"
```

**Multiple ROOT files (comma-separated):**
```json
"source": "/path/to/file1.root, /path/to/file2.root, /path/to/file3.root"
```

**Multiple ROOT files (JSON array):**
```json
"source": [
    "/path/to/file1.root",
    "/path/to/file2.root",
    "/path/to/file3.root"
]
```

**File list (.list file):**
```json
"source": "my_files.list"
```

Contents of `my_files.list` (one path per line):
```
/path/to/file1.root
/path/to/file2.root
/path/to/file3.root
```

All formats create a TChain internally when multiple files are specified.

---

## 3. The Cut System

The cut system allows you to define selection criteria that filter events. Cuts are defined once and applied throughout the analysis with automatic efficiency tracking.

### 3.1 Cut Types Overview

| Type | Purpose | Example Use Case |
|------|---------|------------------|
| **Range Cut** | 1D value window | Mass windows, momentum limits |
| **Trigger Cut** | Bitmask selection | Hardware trigger requirements |
| **Graphical Cut** | 2D polygon | PID bands in beta vs p |

### 3.2 Defining Range Cuts

Range cuts select events where a value falls within `[min, max]`.

**Location:** `src/setup_cuts.h`

```cpp
inline void setupCuts(CutManager& cuts) {

    // ========================================================================
    // RANGE CUTS: defineRangeCut(name, min, max, description)
    // ========================================================================

    // Invariant mass windows
    cuts.defineRangeCut("dilepton_mass", 0.0, 1.5,
                        "Dilepton mass window [GeV]");

    cuts.defineRangeCut("eegamma_mass_pi0", 0.10, 0.17,
                        "e+e-gamma mass for pi0 Dalitz [GeV]");

    cuts.defineRangeCut("eegamma_mass_eta", 0.50, 0.60,
                        "e+e-gamma mass for eta Dalitz [GeV]");

    // Angular cuts
    cuts.defineRangeCut("opening_angle", 9.0, 180.0,
                        "Opening angle cut [deg]");

    cuts.defineRangeCut("ep_theta", 18.0, 85.0,
                        "Positron theta acceptance [deg]");

    cuts.defineRangeCut("em_theta", 18.0, 85.0,
                        "Electron theta acceptance [deg]");

    // Momentum cuts
    cuts.defineRangeCut("ep_momentum", 50.0, 2000.0,
                        "Positron momentum range [MeV/c]");

    cuts.defineRangeCut("em_momentum", 50.0, 2000.0,
                        "Electron momentum range [MeV/c]");

    // ECAL photon quality cuts
    cuts.defineRangeCut("photon_energy", 50.0, 10000.0,
                        "Photon energy window [MeV]");

    cuts.defineRangeCut("photon_theta", 12.0, 45.0,
                        "Photon theta acceptance [deg]");

    cuts.defineRangeCut("photon_beta", 0.9, 1.1,
                        "Photon beta (should be ~1)");

    // Quality cuts
    cuts.defineRangeCut("vertex_z", -60.0, 0.0,
                        "Vertex Z position [mm]");

    cuts.defineRangeCut("chi2_ndf", 0.0, 5.0,
                        "Track chi2/ndf quality");

    cuts.printDefinedCuts();  // Show all defined cuts
}
```

### 3.3 Applying Range Cuts in processEvent()

**Location:** `main.cc` in `processEvent()` function

```cpp
void processEvent(NTupleReader& reader, Manager& mgr, CutManager& cuts, ...) {

    // ... create particles and calculate observables ...

    double m_ee = dilepton.massGeV(analysis_type);
    double oa_epem = positron.openingAngle(electron, analysis_type);
    double ep_theta = positron.theta(analysis_type);
    double em_theta = electron.theta(analysis_type);

    // ========================================================================
    // PATTERN 1: Check if cut exists, then apply (recommended for optional cuts)
    // ========================================================================

    if (cuts.hasRangeCut("dilepton_mass")) {
        if (!cuts.passRangeCut("dilepton_mass", m_ee)) return;
    }

    // ========================================================================
    // PATTERN 2: Direct application (if cut always exists)
    // ========================================================================

    if (!cuts.passRangeCut("opening_angle", oa_epem)) return;

    // ========================================================================
    // PATTERN 3: Multiple cuts in sequence
    // ========================================================================

    // Apply angular acceptance for both leptons
    if (cuts.hasRangeCut("ep_theta")) {
        if (!cuts.passRangeCut("ep_theta", ep_theta)) return;
    }
    if (cuts.hasRangeCut("em_theta")) {
        if (!cuts.passRangeCut("em_theta", em_theta)) return;
    }

    // ========================================================================
    // PATTERN 4: Conditional cuts (only apply under certain conditions)
    // ========================================================================

    // Apply photon cuts only if we have a valid photon
    if (gamma1.isValid()) {
        if (cuts.hasRangeCut("photon_energy")) {
            if (!cuts.passRangeCut("photon_energy", gamma1.getEnergy())) {
                // Photon fails quality cut - mark as bad but don't reject event
                // (photon-specific logic here)
            }
        }
    }

    // ========================================================================
    // If we reach here, event passed all cuts - fill histograms
    // ========================================================================

    mgr.fillw("mass_ee_cut", m_ee, weight);
    // ... more histogram filling ...
}
```

### 3.4 Defining Trigger Cuts

Trigger cuts check bits in a trigger word (integer variable from ntuple).

```cpp
inline void setupCuts(CutManager& cuts) {

    // ========================================================================
    // TRIGGER CUTS: defineTriggerCut(name, mask, require_all, description)
    // ========================================================================

    // PT3 trigger (bit 2 = value 4)
    // require_all = false means OR logic (any matching bit passes)
    cuts.defineTriggerCut("PT3", 4, false, "PT3 physics trigger");

    // Lepton trigger (bit 4 = value 16)
    cuts.defineTriggerCut("lepton", 16, false, "Lepton trigger");

    // Combined trigger requirement (bits 2 AND 4 must both be set)
    // require_all = true means AND logic
    cuts.defineTriggerCut("PT3_and_lepton", 20, true, "PT3 AND lepton trigger");

    // Multiple OR bits (any of bits 0, 1, 2 set)
    cuts.defineTriggerCut("any_physics", 7, false, "Any physics trigger");
}
```

**Applying trigger cuts:**

```cpp
void processEvent(...) {
    // Read trigger word from ntuple
    int trigger_word = static_cast<int>(reader["trigger"]);

    // Apply trigger requirement
    if (cuts.hasTriggerCut("PT3")) {
        if (!cuts.passTriggerCut("PT3", trigger_word)) return;
    }
}
```

### 3.5 Defining Graphical Cuts (2D Polygon)

Graphical cuts use ROOT's TCutG to define 2D polygon regions.

```cpp
inline void setupCuts(CutManager& cuts) {

    // ========================================================================
    // GRAPHICAL CUTS: loadGraphicalCut(name, file, cutg_name, description)
    // ========================================================================

    // Load from ROOT file containing TCutG object
    try {
        cuts.loadGraphicalCut("ep_pid", "cuts/lepton_pid.root", "cutg_positron",
                              "Positron PID band");
        cuts.loadGraphicalCut("em_pid", "cuts/lepton_pid.root", "cutg_electron",
                              "Electron PID band");
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not load PID cuts: " << e.what() << "\n";
        // Analysis continues without these cuts
    }

    // ========================================================================
    // CREATE GRAPHICAL CUT PROGRAMMATICALLY
    // ========================================================================

    // Define polygon points (beta vs momentum band)
    TCutG* electron_band = new TCutG("electron_band", 5);
    electron_band->SetPoint(0, 100, 0.95);   // (p, beta)
    electron_band->SetPoint(1, 2000, 0.95);
    electron_band->SetPoint(2, 2000, 1.05);
    electron_band->SetPoint(3, 100, 1.05);
    electron_band->SetPoint(4, 100, 0.95);   // Close polygon

    cuts.addGraphicalCut("electron_beta_band", electron_band,
                         "Electron beta vs p band");
}
```

**Applying graphical cuts:**

```cpp
void processEvent(...) {
    double ep_p = positron.momentum(analysis_type);
    double ep_beta = reader["ep_beta"];  // or calculated

    // Check if point (x, y) is inside polygon
    // passGraphicalCut(name, x_value, y_value)
    if (cuts.hasGraphicalCut("ep_pid")) {
        if (!cuts.passGraphicalCut("ep_pid", ep_p, ep_beta)) return;
    }
}
```

### 3.6 Creating TCutG Files

**Interactive method in ROOT:**

```cpp
// In ROOT session with a 2D histogram displayed:
TH2F* h = (TH2F*)gDirectory->Get("h_beta_vs_p");
h->Draw("colz");

// Then in the canvas, use Edit -> Create TCutG
// Draw points around your selection region
// Save with: cutg->Write() in a ROOT file
```

**Script method:**

```cpp
// create_pid_cut.C
void create_pid_cut() {
    TFile* f = new TFile("cuts/lepton_pid.root", "RECREATE");

    // Positron band (beta vs p)
    TCutG* cut_ep = new TCutG("cutg_positron", 6);
    cut_ep->SetPoint(0, 50, 0.90);
    cut_ep->SetPoint(1, 1500, 0.95);
    cut_ep->SetPoint(2, 1500, 1.05);
    cut_ep->SetPoint(3, 50, 1.10);
    cut_ep->SetPoint(4, 50, 0.90);
    cut_ep->Write();

    // Electron band
    TCutG* cut_em = new TCutG("cutg_electron", 6);
    cut_em->SetPoint(0, 50, 0.90);
    cut_em->SetPoint(1, 1500, 0.95);
    cut_em->SetPoint(2, 1500, 1.05);
    cut_em->SetPoint(3, 50, 1.10);
    cut_em->SetPoint(4, 50, 0.90);
    cut_em->Write();

    f->Close();
}
```

### 3.7 Cut Flow Statistics

The CutManager automatically tracks how many events pass each cut.

**At the end of the analysis:**

```cpp
int main() {
    // ... event loop ...

    // Print cut efficiency report
    cuts.printCutFlow();
}
```

**Output:**

```
================================================================================
                              CUT FLOW SUMMARY
================================================================================
Cut Name                    | Tested     | Passed     | Efficiency
--------------------------------------------------------------------------------
dilepton_mass               |    100000  |     95000  |    95.00%
opening_angle               |     95000  |     82000  |    86.32%
ep_theta                    |     82000  |     78500  |    95.73%
em_theta                    |     78500  |     75200  |    95.80%
================================================================================
```

### 3.8 Advanced Cut Operations

```cpp
// Temporarily disable a cut (for systematic studies)
cuts.setRangeCutActive("opening_angle", false);

// Re-enable
cuts.setRangeCutActive("opening_angle", true);

// Disable all cuts
cuts.setAllCutsActive(false);

// Reset statistics (e.g., between analysis passes)
cuts.resetStatistics();

// Check if cut is active
if (cuts.isRangeCutActive("dilepton_mass")) {
    // ...
}

// Get cut parameters
double min, max;
cuts.getRangeCutLimits("dilepton_mass", min, max);
std::cout << "Mass cut: [" << min << ", " << max << "] GeV\n";
```

### 3.9 Complete Cut Example: e+e- Analysis

```cpp
// src/setup_cuts.h
inline void setupCuts(CutManager& cuts) {
    std::cout << "Setting up e+e- analysis cuts...\n";

    // === Dilepton Selection ===
    cuts.defineRangeCut("mass_ee", 0.0, 1.5, "Dilepton mass [GeV]");
    cuts.defineRangeCut("opening_angle", 9.0, 180.0, "e+e- opening angle [deg]");

    // === Lepton Acceptance ===
    cuts.defineRangeCut("lepton_theta", 18.0, 85.0, "Lepton polar angle [deg]");
    cuts.defineRangeCut("lepton_p", 50.0, 2000.0, "Lepton momentum [MeV/c]");

    // === ECAL Photon Quality ===
    cuts.defineRangeCut("photon_energy", 50.0, 10000.0, "Photon energy [MeV]");
    cuts.defineRangeCut("photon_theta", 12.0, 45.0, "Photon theta [deg]");

    // === Pi0 Selection ===
    cuts.defineRangeCut("mass_eegamma_pi0", 0.10, 0.17, "e+e-g mass for pi0 [GeV]");
    cuts.defineRangeCut("mass_gg_pi0", 0.10, 0.17, "gg mass for pi0 [GeV]");

    // === Eta Selection ===
    cuts.defineRangeCut("mass_eegamma_eta", 0.50, 0.60, "e+e-g mass for eta [GeV]");

    // === Vertex Quality ===
    cuts.defineRangeCut("vertex_z", -60.0, 0.0, "Vertex Z [mm]");

    cuts.printDefinedCuts();
}
```

```cpp
// main.cc - processEvent() excerpt
void processEvent(...) {
    // Create particles and calculate observables
    PParticle dilepton = positron + electron;
    double m_ee = dilepton.massGeV(analysis_type);
    double oa = positron.openingAngle(electron, analysis_type);

    // === Apply cuts sequentially ===

    // 1. Dilepton mass window
    if (!cuts.passRangeCut("mass_ee", m_ee)) return;

    // 2. Opening angle (reject conversions)
    if (!cuts.passRangeCut("opening_angle", oa)) return;

    // 3. Lepton acceptance
    if (!cuts.passRangeCut("lepton_theta", positron.theta(analysis_type))) return;
    if (!cuts.passRangeCut("lepton_theta", electron.theta(analysis_type))) return;
    if (!cuts.passRangeCut("lepton_p", positron.momentum(analysis_type))) return;
    if (!cuts.passRangeCut("lepton_p", electron.momentum(analysis_type))) return;

    // 4. Vertex quality (if available)
    if (reader.hasVariable("eVertZ")) {
        if (!cuts.passRangeCut("vertex_z", reader["eVertZ"])) return;
    }

    // === Event passed all cuts - fill physics histograms ===
    mgr.fillw("mass_dilepton", m_ee, weight);

    // === Check photon quality for e+e-gamma analysis ===
    for (auto* gamma : good_gammas) {
        if (!cuts.passRangeCut("photon_energy", gamma->getEnergy())) continue;
        if (!cuts.passRangeCut("photon_theta", gamma->ecal_theta)) continue;

        // Good photon - calculate e+e-gamma mass
        PParticle eegamma = dilepton + *gamma;
        double m_eeg = eegamma.massGeV(analysis_type);

        // Fill before pi0/eta selection
        mgr.fillw("mass_eegamma", m_eeg, weight);

        // Pi0 Dalitz selection
        if (cuts.passRangeCut("mass_eegamma_pi0", m_eeg)) {
            mgr.fillw("mass_eegamma_pi0_selected", m_eeg, weight);
        }

        // Eta Dalitz selection
        if (cuts.passRangeCut("mass_eegamma_eta", m_eeg)) {
            mgr.fillw("mass_eegamma_eta_selected", m_eeg, weight);
        }
    }
}
```

---

## 4. Adding Histograms

### 4.1 Basic Histogram Creation

**Location:** `src/setup_histograms.h`

```cpp
inline void setupHistograms(Manager& mgr) {

    // ========================================================================
    // 1D HISTOGRAM: create1D(name, title, nbins, xmin, xmax, folder)
    // ========================================================================

    // Title format: "Display Title;X-axis label;Y-axis label"
    mgr.create1D("mass_ee", "Dilepton mass;M_{e^{+}e^{-}} [GeV/c^{2}];Counts",
                 500, 0.0, 1.5, "dilepton");

    // With ROOT-style formatting in labels
    mgr.create1D("cos_theta_cms", "cos#theta (CMS);cos#theta;Counts",
                 100, -1, 1, "cms/angular");

    // ========================================================================
    // 2D HISTOGRAM: create2D(name, title, nbx, xmin, xmax, nby, ymin, ymax, folder)
    // ========================================================================

    mgr.create2D("mass_vs_pt", "M vs p_{T};M [GeV];p_{T} [MeV/c]",
                 100, 0.0, 1.0, 100, 0, 1500, "correlations");

    mgr.create2D("theta_ep_vs_em", "#theta correlation;#theta_{e^{-}};#theta_{e^{+}}",
                 90, 0, 90, 90, 0, 90, "correlations");
}
```

### 4.2 Folder Organization

Histograms are organized in folders within the ROOT file:

```cpp
// Creates folder structure: quality/
mgr.create1D("vertex_z", "...", 100, -100, 0, "quality");

// Creates nested folders: lab/angular/
mgr.create1D("ep_theta_lab", "...", 90, 0, 90, "lab/angular");

// Creates nested folders: cms/momentum/
mgr.create1D("ep_p_cms", "...", 100, 0, 2000, "cms/momentum");
```

### 4.3 Filling Histograms

**In processEvent():**

```cpp
// Unweighted fill
mgr.fill("histogram_name", value);

// Weighted fill (recommended - use event weight)
mgr.fillw("histogram_name", value, weight);

// 2D histogram
mgr.fillw("hist_2d", x_value, y_value, weight);
```

---

## 5. Working with Ntuples

### 5.1 Output Ntuple Variables

**Location:** Fill variables in `main.cc:processEvent()`

```cpp
void processEvent(...) {
    // ... physics calculations ...

    // Get ntuple reference
    DynamicHNtuple& nt = mgr.getDynamicNtuple("nt_particles");

    // Fill variables (can add new ones anytime!)
    nt["ep_p"] = positron.momentum(analysis_type);
    nt["ep_theta"] = positron.theta(analysis_type);
    nt["ep_phi"] = positron.phi(analysis_type);

    nt["ee_mass"] = dilepton.massGeV(analysis_type);
    nt["opening_angle"] = oa_epem;

    // Conditional variables (missing = -1.0 sentinel)
    if (gamma1.isValid()) {
        nt["g1_energy"] = gamma1.getEnergy();
    }

    nt["weight"] = weight;

    // Write entry
    nt.fill();
}
```

### 5.2 Reading Output in ROOT

```cpp
// In ROOT session
TFile* f = TFile::Open("output.root");
TNtuple* nt = (TNtuple*)f->Get("nt_particles");

// Draw histogram
nt->Draw("ee_mass");
nt->Draw("ee_mass", "opening_angle > 9");
nt->Draw("ee_mass", "weight");  // Weighted

// 2D plot
nt->Draw("ep_theta:em_theta", "", "colz");
```

---

## 6. Common Analysis Patterns

### 6.1 Before/After Cut Comparison

```cpp
// Fill BEFORE cuts
mgr.fillw("mass_ee_all", m_ee, weight);

// Apply cut
if (!cuts.passRangeCut("opening_angle", oa)) return;

// Fill AFTER cuts
mgr.fillw("mass_ee_oa_cut", m_ee, weight);
```

### 6.2 Multiple Mass Hypotheses

```cpp
// Same particle, different mass hypothesis
PParticleFwd particle_p(Physics::MASS_PROTON, "as_proton");
PParticleFwd particle_d(Physics::MASS_DEUTERON, "as_deuteron");

particle_p.setFromReader(reader, 1, KinematicType::RECONSTRUCTED);
particle_d.setFromReader(reader, 1, KinematicType::RECONSTRUCTED);

// Compare reconstructed masses
mgr.fillw("mass_hypo_proton", particle_p.mass(analysis_type), weight);
mgr.fillw("mass_hypo_deuteron", particle_d.mass(analysis_type), weight);
```

### 6.3 Systematic Variations

```cpp
// Nominal analysis
cuts.setRangeCutActive("opening_angle", true);
// ... run analysis ...

// Systematic: no opening angle cut
cuts.setRangeCutActive("opening_angle", false);
cuts.resetStatistics();
// ... run analysis again ...
```

### 6.4 Signal Region Selection

```cpp
// Define signal and sideband regions
cuts.defineRangeCut("pi0_signal", 0.120, 0.150, "Pi0 signal region [GeV]");
cuts.defineRangeCut("pi0_sideband_low", 0.050, 0.100, "Pi0 low sideband [GeV]");
cuts.defineRangeCut("pi0_sideband_high", 0.170, 0.220, "Pi0 high sideband [GeV]");

// In processEvent:
if (cuts.passRangeCut("pi0_signal", m_eeg)) {
    mgr.fillw("observable_signal", observable, weight);
}
if (cuts.passRangeCut("pi0_sideband_low", m_eeg) ||
    cuts.passRangeCut("pi0_sideband_high", m_eeg)) {
    mgr.fillw("observable_sideband", observable, weight);
}
```

---

## Quick Reference Card

### Cut Definition Syntax

```cpp
// Range cut
cuts.defineRangeCut("name", min, max, "description");

// Trigger cut
cuts.defineTriggerCut("name", bitmask, require_all, "description");

// Graphical cut (from file)
cuts.loadGraphicalCut("name", "file.root", "TCutG_name", "description");
```

### Cut Application Syntax

```cpp
// Range cut
if (!cuts.passRangeCut("name", value)) return;

// Trigger cut
if (!cuts.passTriggerCut("name", trigger_word)) return;

// Graphical cut
if (!cuts.passGraphicalCut("name", x, y)) return;

// Check existence first (optional cuts)
if (cuts.hasRangeCut("name")) {
    if (!cuts.passRangeCut("name", value)) return;
}
```

### Histogram Syntax

```cpp
// Create
mgr.create1D("name", "title;X;Y", nbins, xmin, xmax, "folder");
mgr.create2D("name", "title;X;Y", nbx, x0, x1, nby, y0, y1, "folder");

// Fill
mgr.fillw("name", value, weight);
mgr.fillw("name_2d", x, y, weight);
```

### Ntuple Syntax

```cpp
DynamicHNtuple& nt = mgr.getDynamicNtuple("name");
nt["variable"] = value;
nt.fill();
```

---

*Last updated: February 2026*

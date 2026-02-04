# FAT Analysis Framework - Pedagogical Tutorial

A comprehensive step-by-step guide to understanding the FAT (Fast Analysis Tool) framework for HADES physics analysis.

---

## Table of Contents

1. [Project Directory Structure](#1-project-directory-structure)
2. [Analysis Flow - Step by Step](#2-analysis-flow---step-by-step)
3. [JSON Configuration Parameters](#3-json-configuration-parameters)
4. [Particle Objects](#4-particle-objects)
5. [Cut System](#5-cut-system)
6. [Output Handling](#6-output-handling)
7. [Complete Data Flow Diagram](#7-complete-data-flow-diagram)
8. [Key Code Locations](#8-key-code-locations)
9. [Summary of Architecture](#9-summary-of-architecture)

---

## 1. Project Directory Structure

```
FAT/
├── main.cc                          # Entry point, event loop, physics analysis
├── config.json                      # Analysis configuration (JSON)
├── Makefile                         # Build system
├── MyLinkDef.h                      # ROOT dictionary linkage
├── MyDict.cc                        # Generated ROOT dictionary
├── h68_10.list                      # Input file list
│
├── src/                             # Core framework headers
│   ├── analysis_config.h            # JSON parser & configuration loader
│   ├── manager.h                    # Output orchestration (histograms, ntuples)
│   │
│   ├── pparticle.h                  # Base particle class
│   ├── pparticle_ecal.h             # ECAL detector particle extension
│   ├── pparticle_fwd.h              # Forward Tracker particle extension
│   │
│   ├── cut_manager.h                # Cut definition & evaluation system
│   ├── setup_cuts.h                 # Predefined cut configurations
│   │
│   ├── histogram_registry.h         # Central histogram storage
│   ├── histogram_factory.h          # Histogram creation utilities
│   ├── histogram_builder.h          # Fluent histogram builder
│   ├── setup_histograms.h           # Predefined histogram configurations
│   │
│   ├── hntuple.h / hntuple.cc       # Classic ntuple (fixed structure)
│   ├── dynamic_hntuple.h            # Flexible ntuple (add vars anytime)
│   ├── setup_ntuples.h              # Predefined ntuple configurations
│   │
│   ├── ntuple_reader.h              # Input ROOT file reader
│   ├── boost_frame.h                # Lorentz boost utilities
│   ├── reactionvertexfind.h         # Vertex reconstruction
│   ├── progressbar.h                # Progress display
│   └── console_box.h                # Formatted console output
│
├── examples/                        # Usage examples and tests
│   ├── PParticle_Usage_Examples.cc  # PParticle class tutorial
│   ├── PPip_ID_Refactored.cc        # Refactored analysis example
│   ├── test_boost_sign_convention.cc
│   ├── test_hntuple_improved_errors.cc
│   ├── test_improved_manager.cc
│   └── Makefile
│
├── doc/                             # Documentation
│   └── ANALYSIS_TUTORIAL.md         # This file
│
└── Documentation (*.md files in root)
    ├── DEVELOPER_TUTORIAL.md
    ├── HISTOGRAM_MANAGEMENT_GUIDE.md
    ├── HISTOGRAM_SYSTEM_SUMMARY.md
    ├── HNTUPLE_IMPROVEMENTS.md
    ├── BOOST_SIGN_CONVENTION.md
    ├── REFACTORING_GUIDE.md
    ├── REFACTORING_PROPOSAL.md
    ├── MANAGER_ANALYSIS.md
    └── NTUPLE_STATUS.md
```

---

## 2. Analysis Flow - Step by Step

### Phase 1: Configuration Loading (main.cc)

```
config.json  →  JsonParser::parseFile()  →  AnalysisConfig object
```

The `AnalysisConfig` class (src/analysis_config.h) contains a built-in JSON parser and provides typed accessors for all parameters.

```cpp
AnalysisConfig config;
config.load("config.json");
config.print();  // Display summary
```

### Phase 2: Beam System Setup

```cpp
double beam_kinetic_energy = config.getBeamKineticEnergy();
PParticle projectile = ParticleFactory::createBeamProton(beam_kinetic_energy);
PParticle target = ParticleFactory::createTargetProton();
PParticle beam = projectile + target;
EventFrames frames;
frames.addFrame("beam", BoostFrame(beam));
```

Creates the beam system and initializes reference frames for CMS boosts.

### Phase 3: Input Data Setup

```cpp
NTupleReader reader;
if (config.isInputRootFile()) {
    reader.open(source, tree_name);
} else {
    reader.openFromList(source, tree_name);  // .list file → TChain
}
```

### Phase 4: Output Setup

```cpp
Manager manager;
manager.openFile(output_filename, "RECREATE");
setupHistograms(manager);      // Creates ~40 histograms
setupNtuples(manager, config); // Creates 2 dynamic ntuples

CutManager cuts;
setupCuts(cuts);  // Defines mass window cuts
```

### Phase 5: Event Loop

```cpp
for (Long64_t i = start_event; i < end_event; ++i) {
    reader.getEntry(i);
    processEvent(reader, manager, cuts, beam, projectile, frames, config);
    progress.update(processed);
}
```

### Phase 6: Finalization

```cpp
cuts.printCutFlow();     // Print cut statistics
manager.printSummary();  // Print histogram/ntuple summary
manager.closeFile();     // Write all output to ROOT file
```

---

## 3. JSON Configuration Parameters

### Complete config.json Structure

```json
{
  "input": {
    "source": "h68_10.list",
    "tree_name": "PPip_ID",
    "start_event": 0,
    "max_events": -1
  },
  "output": {
    "filename": "output_ppip.root",
    "option": "RECREATE",
    "keep_intermediate_tree": false,
    "missing_value": -1.0
  },
  "beam": {
    "kinetic_energy": 4500.0,
    "luminosity": 1.0
  },
  "kinematics": {
    "reconstructed": true,
    "corrected": true,
    "simulated": false,
    "analysis_type": "reconstructed"
  },
  "forward_tracker": {
    "fwdet_1": false,
    "fwdet_2": false,
    "fwdet_3": false
  },
  "ecal": {
    "ecal_1": false,
    "ecal_2": false,
    "ecal_3": false
  }
}
```

### Parameter Reference Table

| Section | Parameter | Type | Default | Purpose |
|---------|-----------|------|---------|---------|
| **input** | `source` | string | "h68_10.list" | Input file (.root or .list) |
| | `tree_name` | string | "PPip_ID" | ROOT TTree name |
| | `start_event` | int | 0 | First event index |
| | `max_events` | int | -1 | Max events (-1 = all) |
| **output** | `filename` | string | "output.root" | Output ROOT file |
| | `option` | string | "RECREATE" | ROOT file mode |
| | `keep_intermediate_tree` | bool | false | Keep TTree before conversion |
| | `missing_value` | float | -1.0 | Sentinel for missing variables |
| **beam** | `kinetic_energy` | double | 4500.0 | Beam KE [MeV] |
| | `luminosity` | double | 1.0 | Weight scaling factor |
| **kinematics** | `reconstructed` | bool | true | Load reconstructed momenta |
| | `corrected` | bool | true | Load corrected momenta |
| | `simulated` | bool | false | Load MC truth momenta |
| | `analysis_type` | string | "reconstructed" | Which to use for analysis |
| **forward_tracker** | `fwdet_1/2/3` | bool | false | Enable FT detector slots |
| **ecal** | `ecal_1/2/3` | bool | false | Enable ECAL slots |

### Key Parameter Effects

- **`analysis_type`** selects which momentum representation is used throughout the analysis
- **`luminosity`** scales all histogram weights: `weight = event_weight * luminosity`
- **`missing_value`** fills gaps in DynamicHNtuple for events missing certain variables

---

## 4. Particle Objects

### 4.1 PParticle - Base Class (src/pparticle.h)

#### Core Design: Multiple Kinematic Representations

Each particle stores **3 momentum representations simultaneously**:

| KinematicType | Purpose | Source Variables |
|---------------|---------|------------------|
| `RECONSTRUCTED` | Raw detector measurement | `p_p`, `p_theta`, `p_phi` |
| `CORRECTED` | Energy-loss corrected | `p_p_corr_p`, `p_theta`, `p_phi` |
| `SIMULATED` | MC truth | `p_sim_px`, `p_sim_py`, `p_sim_pz` |

#### Data Members

```cpp
double mass_;                      // Rest mass [MeV/c²]
std::string name_;                 // Particle name

TLorentzVector p4_reconstructed_;  // Current 4-momentum (may be boosted)
TLorentzVector p4_corrected_;
TLorentzVector p4_simulated_;

TLorentzVector lab_frame_reconstructed_;  // Original LAB frame (preserved)
TLorentzVector lab_frame_corrected_;
TLorentzVector lab_frame_simulated_;
```

#### Construction Methods

```cpp
// From spherical coordinates (p [MeV/c], theta [deg], phi [deg])
PParticle proton(Physics::MASS_PROTON, "p");
proton.setFromSpherical(1500, 45.0, 30.0, KinematicType::RECONSTRUCTED);

// From Cartesian (px, py, pz [MeV/c])
proton.setFromCartesian(px, py, pz, KinematicType::SIMULATED);

// Factory functions
PParticle proton = ParticleFactory::createProton(p, theta, phi);
PParticle beam = ParticleFactory::createBeamProton(T_kin);
PParticle target = ParticleFactory::createTargetProton();
```

#### Physics Accessor Methods

```cpp
double mass(KinematicType type);      // Invariant mass [MeV/c²]
double massGeV(KinematicType type);   // Invariant mass [GeV/c²]
double momentum(KinematicType type);  // |p| [MeV/c]
double energy(KinematicType type);    // E [MeV]
double theta(KinematicType type);     // Polar angle [deg]
double phi(KinematicType type);       // Azimuthal angle [deg]
double cosTheta(KinematicType type);  // cos(theta)
double rapidity(KinematicType type);  // y
double beta(KinematicType type);      // v/c
double px(KinematicType type);        // Momentum x-component [MeV/c]
double py(KinematicType type);        // Momentum y-component [MeV/c]
double pz(KinematicType type);        // Momentum z-component [MeV/c]

double openingAngle(const PParticle& other, KinematicType type);  // [deg]
double deltaPhi(const PParticle& other, KinematicType type);      // [deg]
```

#### Composite Particle Operations

```cpp
// Missing mass technique
PParticle neutron = beam - proton - pion;

// Resonance reconstruction
PParticle deltaPP = proton + pion;

// Multi-particle composites
PParticle system = proton + pion1 + pion2;
```

#### Reference Frame Transformations

```cpp
// Apply Lorentz boost
proton.boost(beta_vector);

// Boost to another particle's rest frame
proton.boostToRestFrame(resonance, KinematicType::RECONSTRUCTED);

// Restore original LAB momentum
proton.resetToLAB();

// Access LAB frame (even after boosts)
TLorentzVector lab_p4 = proton.labFrame(KinematicType::RECONSTRUCTED);
```

---

### 4.2 PParticleEcal - ECAL Extension (src/pparticle_ecal.h)

**Inherits from**: PParticle

**Purpose**: Electromagnetic Calorimeter particles (photons, electrons, neutral pions)

#### Additional Data Members

```cpp
// Identification
int    ecal_pid;           // Particle ID
int    ecal_clusterId;     // EMC cluster index
double ecal_charge;        // Charge

// Energy/momentum
double ecal_energy;        // Calorimeter energy [MeV]
double ecal_p_pid;         // PID momentum [MeV/c]

// Mass
double ecal_mass;          // Reconstructed mass [MeV/c²]
double ecal_mass2;         // Mass squared [MeV²/c⁴]

// Timing
double ecal_beta;          // Measured β = v/c
double ecal_tof;           // Time of flight
double ecal_tofRec;        // Reconstructed TOF

// Position (primary)
double ecal_theta;         // Polar angle [deg]
double ecal_phi;           // Azimuthal angle [deg]
double ecal_r;             // Radial position [mm]
double ecal_z;             // Z position [mm]

// Position (secondary reconstruction)
double ecal_theta2, ecal_phi2, ecal_r2, ecal_z2;

// Quality
double ecal_chi2;          // Fit χ²
double ecal_dist;          // Distance to EMC [mm]
int    ecal_mult;          // Multiplicity
int    ecal_counter;       // Counter
bool   ecal_valid;         // Successfully filled?
```

#### Loading from NTuple

```cpp
PParticleEcal photon(0.0, "gamma");  // Mass=0 for photon
photon.setFromReader(reader, 1, KinematicType::RECONSTRUCTED);  // Reads neutr_*_1 variables
```

**Input Variable Pattern** (N = 1, 2, or 3):
- `neutr_energy_N`, `neutr_p_N`, `neutr_p_pid_N`
- `neutr_theta_N`, `neutr_phi_N`, `neutr_r_N`, `neutr_z_N`
- `neutr_beta_N`, `neutr_tof_N`, `neutr_chi2_N`
- `neutr_mult`, `neutr_counter`

#### ECAL-Specific Methods

```cpp
bool isValid();                   // Check if successfully filled
double getEnergy();               // Calorimeter energy [MeV]
double getChi2();                 // Fit quality
double getDistToEmc();            // Distance to EMC cluster [mm]
double getMeasuredBeta();         // TOF-based velocity
int getClusterId();               // EMC cluster index

// Quality cuts
bool passesQuality(double minEnergy, double maxChi2);

// Debug output
void printEcal();
```

---

### 4.3 PParticleFwd - Forward Tracker Extension (src/pparticle_fwd.h)

**Inherits from**: PParticle

**Purpose**: Forward Tracker particles (protons, deuterons, pions at forward angles)

#### Additional Data Members

```cpp
// Measured quantities
double fwd_beta;           // Measured β
double fwd_mass;           // Reconstructed mass [MeV/c²]
double fwd_mass2;          // Mass squared [MeV²/c⁴]
double fwd_charge;         // Charge
double fwd_tof;            // Time of flight
double fwd_tofRec;         // Reconstructed TOF

// Position (primary)
double fwd_theta;          // Polar angle [deg]
double fwd_phi;            // Azimuthal angle [deg]
double fwd_r;              // Radial position [mm]
double fwd_z;              // Z position [mm]

// Position (secondary reconstruction)
double fwd_theta2, fwd_phi2, fwd_r2, fwd_z2;

// Quality
double fwd_chi2;           // Fit χ²
int    fwd_ndf;            // Degrees of freedom
double fwd_chi2ndf;        // χ²/ndf (computed)
double fwd_distToRpc;      // Distance to RPC [mm]
int    fwd_mult;           // Multiplicity
int    fwd_counter;        // Counter
bool   fwd_valid;          // Successfully filled?
```

#### Loading from NTuple

```cpp
PParticleFwd fw_proton(Physics::MASS_PROTON, "fw_p1");
fw_proton.setFromReader(reader, 1, KinematicType::RECONSTRUCTED);  // Reads fwdet_*_1 variables
```

**Input Variable Pattern** (N = 1, 2, or 3):
- `fwdet_p_N`, `fwdet_theta_N`, `fwdet_phi_N` (required for momentum)
- `fwdet_beta_N`, `fwdet_mass_N`, `fwdet_mass2_N`
- `fwdet_chi2_N`, `fwdet_ndf_N`, `fwdet_dist_tofhit_N`
- `fwdet_mult`, `fwdet_counter`

#### FT-Specific Methods

```cpp
bool isValid();              // Check if successfully filled
double getChi2NDF();         // χ²/ndf quality metric
double getDistToRpc();       // Distance to RPC hit [mm]
double getMeasuredBeta();    // TOF-based velocity

// Quality cuts
bool passesQuality(double maxChi2ndf, double maxDistRpc);

// Debug output
void printFwd();
```

---

### 4.4 Physics Constants Namespace

```cpp
namespace Physics {
    // Nucleons
    constexpr double MASS_PROTON = 938.27208816;    // MeV/c²
    constexpr double MASS_NEUTRON = 939.56542052;

    // Pions
    constexpr double MASS_PION_PLUS = 139.57039;
    constexpr double MASS_PION_MINUS = 139.57039;
    constexpr double MASS_PION_ZERO = 134.9768;

    // Leptons
    constexpr double MASS_ELECTRON = 0.51099895;
    constexpr double MASS_MUON_PLUS = 105.6583755;
    constexpr double MASS_MUON_MINUS = 105.6583755;

    // Light nuclei
    constexpr double MASS_DEUTERON = 1875.61294257;
    constexpr double MASS_TRITON = 2808.92113298;

    // Conversions
    constexpr double D2R = 0.0174532925199;  // degrees to radians
    constexpr double R2D = 57.2957795131;    // radians to degrees
}
```

---

### 4.5 ParticleFactory Namespace

```cpp
namespace ParticleFactory {
    // Standard particles
    PParticle createProton(double p, double theta, double phi, KinematicType type);
    PParticle createPiPlus(double p, double theta, double phi, KinematicType type);
    PParticle createPiMinus(double p, double theta, double phi, KinematicType type);
    PParticle createEPlus(double p, double theta, double phi, KinematicType type);
    PParticle createEMinus(double p, double theta, double phi, KinematicType type);
    PParticle createMuPlus(double p, double theta, double phi, KinematicType type);
    PParticle createMuMinus(double p, double theta, double phi, KinematicType type);
    PParticle createDeuteron(double p, double theta, double phi, KinematicType type);
    PParticle createTriton(double p, double theta, double phi, KinematicType type);

    // Beam/target
    PParticle createBeamProton(double T_kin);    // From kinetic energy
    PParticle createTargetProton();              // At rest
}
```

---

### 4.6 Class Hierarchy Diagram

```
PParticle (base)
├── mass_, name_
├── 3× TLorentzVector (reconstructed, corrected, simulated)
├── 3× LAB frame copies (immune to boosts)
├── Arithmetic operators (+, -)
├── Frame transformations (boost, resetToLAB)
│
├── PParticleEcal (inherits PParticle)
│   ├── +20 ECAL fields (energy, chi2, dist, beta, etc.)
│   ├── setFromReader(reader, index, type)
│   └── passesQuality(minEnergy, maxChi2)
│
└── PParticleFwd (inherits PParticle)
    ├── +16 FT fields (chi2ndf, distRpc, beta, etc.)
    ├── setFromReader(reader, index, type)
    └── passesQuality(maxChi2ndf, maxDistRpc)
```

---

## 5. Cut System

### 5.1 Cut Types (src/cut_manager.h)

| Cut Type | Purpose | Definition | Evaluation |
|----------|---------|------------|------------|
| **RangeCut** | 1D value window | `(name, min, max, description)` | `value ∈ [min, max]` |
| **TriggerCut** | Bitmask selection | `(name, mask, require_all, description)` | AND/OR of bits |
| **GraphicalCut** | 2D polygon (TCutG) | `(name, TCutG*, description)` | Point inside polygon |

### 5.2 Cut Definition (src/setup_cuts.h)

```cpp
inline void setupCuts(CutManager& cuts) {
    // Range cuts for mass windows
    cuts.defineRangeCut("neutron_mass", 0.899, 0.986, "Neutron mass window [GeV]");
    cuts.defineRangeCut("deltaPP_mass", 0.8, 1.8, "Delta++ mass window [GeV]");

    // Trigger cuts (example)
    cuts.defineTriggerCut("physics", 4, false, "PT3 trigger");  // OR logic

    // Graphical cuts from ROOT file (example)
    cuts.loadGraphicalCut("proton_pid", "cuts/proton_pid.root", "cut_proton");
}
```

### 5.3 Cut Execution in processEvent

```cpp
// Pattern 1: Check existence, evaluate, early return on failure
if (cuts.hasRangeCut("neutron_mass")) {
    if (!cuts.passRangeCut("neutron_mass", m_n)) return;
}

// Pattern 2: Direct evaluation (if cut always exists)
if (!cuts.passRangeCut("deltaPP_mass", m_deltaPP)) return;

// Pattern 3: Graphical cut (2D)
if (!cuts.passGraphicalCut("proton_pid", beta, momentum)) return;

// Only reached if all cuts pass
mgr.fillw("mass_n_cut", m_n, weight);
```

### 5.4 Cut Statistics

Each cut automatically tracks:
- **`tested`**: Number of times evaluated
- **`passed`**: Number of times passed
- **`efficiency()`**: passed/tested ratio

```cpp
// At end of analysis
cuts.printCutFlow();  // Prints formatted ASCII table
```

**Output Example:**
```
╔════════════════════════════════════════════════════════════════╗
║                       CUT FLOW SUMMARY                         ║
╠════════════════════════════════════════════════════════════════╣
║ Cut Name                   │ Tested   │ Passed   │ Efficiency  ║
╠────────────────────────────┼──────────┼──────────┼─────────────╣
║ neutron_mass               │  100000  │   85000  │    85.00%   ║
║ deltaPP_mass               │   85000  │   71400  │    84.00%   ║
╚════════════════════════════════════════════════════════════════╝
```

### 5.5 Convenience Macros

```cpp
CUT_FLOW_START(cuts);
CUT_FLOW_CHECK(cuts, "neutron_mass", mass_n);
CUT_FLOW_CHECK(cuts, "deltaPP_mass", mass_deltaPP);
CUT_FLOW_CHECK_2D(cuts, "proton_pid", beta, momentum);
if (_cutflow_passed) {
    // All cuts passed - fill histograms/ntuples
}
```

### 5.6 Cut Flow Pattern

```
Event enters processEvent()
    │
    ├─► Fill "before cut" histograms (quality control)
    │
    ├─► Cut 1: neutron_mass in [0.899, 0.986]?
    │       NO → return (skip event)
    │
    ├─► Cut 2: deltaPP_mass in [0.8, 1.8]?
    │       NO → return (skip event)
    │
    └─► All cuts pass → Fill "after cut" histograms & ntuples
```

### 5.7 Cut Activation Control

```cpp
// Disable specific cuts without removing them
cuts.setRangeCutActive("neutron_mass", false);

// Re-enable
cuts.setRangeCutActive("neutron_mass", true);

// Toggle all cuts
cuts.setAllCutsActive(false);  // Disable all
cuts.setAllCutsActive(true);   // Enable all

// Reset statistics for fresh run
cuts.resetStatistics();
```

---

## 6. Output Handling

### 6.1 Manager Class (src/manager.h)

**Central orchestrator** for all output operations:

```cpp
Manager manager;
manager.openFile("output.root", "RECREATE");

// === Histogram Creation ===
manager.create1D("name", "title;X;Y", nbins, xmin, xmax, "folder", "description");
manager.create2D("name", "title", nbx, x0, x1, nby, y0, y1, "folder");
manager.create3D("name", "title", nbx, x0, x1, nby, y0, y1, nbz, z0, z1, "folder");

// Arrays of histograms
manager.create1DArray("basename", "title", count, nbins, xmin, xmax, "folder");

// === Ntuple Creation ===
manager.createDynamicNtuple("name", "title", missing_value, keep_tree);

// === Filling (in event loop) ===
manager.fill("histogram_name", value);           // Unweighted
manager.fillw("histogram_name", value, weight);  // Weighted
manager.fill("histogram_2d", x, y);
manager.fillw("histogram_2d", x, y, weight);

// === Finalization ===
manager.closeFile();  // Writes everything to ROOT file
```

### 6.2 Histogram System

#### Three Complementary APIs

**A. Manager (High-Level Convenience)**
```cpp
manager.create1D("h_theta", "Theta;#theta [deg];Counts",
                 100, 0, 180, "angular", "Polar angle distribution");
```

**B. HistogramBuilder (Fluent Interface)**
```cpp
auto h = HistogramBuilder()
    .name("h_theta_p")
    .title("Proton #theta")
    .bins(100, 0, 180)
    .folder("lab/angular")
    .description("Proton polar angle in LAB")
    .tag("proton")
    .tag("angular")
    .build1D();
```

**C. HistogramFactory (Static Methods)**
```cpp
// Create array of histograms
auto hists = HistogramFactory::create1DArray("h_mass", "Mass", 10, 100, 0, 2);
// Creates: h_mass_0, h_mass_1, ..., h_mass_9
```

#### Folder Organization

```cpp
mgr.create1D("p_momentum", "...", 100, 0, 3000, "lab/momentum");
mgr.create1D("p_theta_cms", "...", 100, 0, 180, "cms/angular");
mgr.create2D("dalitz", "...", 100, 1, 4, 100, 1, 4, "correlations");
```

Creates folder hierarchy in ROOT file:
```
output.root
├── quality/
│   ├── mass_n
│   └── mass_n_cut
├── lab/
│   ├── angular/
│   │   └── p_theta
│   └── momentum/
│       └── p_momentum
├── cms/
│   └── angular/
│       └── p_theta_cms
└── correlations/
    └── dalitz
```

#### Histogram Registry

Central storage with metadata support:
```cpp
// Type-safe access
TH1F* h1 = manager.getAs<TH1F>("h_theta");
TH2F* h2 = manager.getAs<TH2F>("dalitz");

// Check existence
if (manager.has("h_theta")) { ... }

// List all histograms
manager.printSummary();
```

---

### 6.3 Ntuple System

#### Two Ntuple Types

| Feature | HNtuple | DynamicHNtuple |
|---------|---------|----------------|
| Variable addition | Before first fill only | Anytime |
| Internal storage | TNtuple directly | TTree → TNtuple |
| Missing values | Error | Filled with sentinel |
| Use case | Known, fixed structure | Evolving analysis |

#### HNtuple (Classic)

```cpp
HNtuple nt("nt_basic", "Basic observables");
nt["momentum"] = p;
nt["theta"] = theta;
nt["phi"] = phi;
nt.fill();  // Variables locked after first fill
```

#### DynamicHNtuple (Recommended)

```cpp
// Creation
manager.createDynamicNtuple("nt_physics", "Physics variables",
                            config.getMissingValue(),           // -1.0
                            config.getKeepIntermediateTree());  // false

// In event loop
DynamicHNtuple& nt = *manager.getDynamicNtuple("nt_physics");

// Add variables freely - works even on first occurrence!
nt["p_p"] = proton.momentum();
nt["p_theta"] = proton.theta();
nt["weight"] = weight;

// Conditional variables are fine
if (has_ecal_hit) {
    nt["ecal_energy"] = ecal.getEnergy();  // Missing in other events → -1.0
}

nt.fill();  // Write entry
```

#### Finalization Process

When `manager.closeFile()` is called:
1. TTree entries are read
2. All variables sorted alphabetically
3. TNtuple created with final variable list
4. Missing values filled with sentinel (-1.0)
5. Intermediate TTree file deleted (unless `keep_intermediate_tree: true`)

**Console Output:**
```
╔════════════════════════════════════════════════════════════════╗
║         Converting TTree → TNtuple                             ║
╠════════════════════════════════════════════════════════════════╣
║ NTuple:     nt_physics                                         ║
║ Variables:  47                                                 ║
║ Entries:    50000                                              ║
╚════════════════════════════════════════════════════════════════╝

Variables (alphabetical order):
  [ 0] ecal_energy
  [ 1] p_p
  [ 2] p_theta
  ...

Converting: [████████████████████████████████████] 100%  Done in 0:05
✓ TNtuple 'nt_physics' created with 47 variables, 50000 entries
```

---

### 6.4 Setup Files

#### setup_histograms.h

```cpp
inline void setupHistograms(Manager& mgr) {
    // Quality control
    mgr.create1D("mass_n", "Neutron mass;m [GeV];Counts", 100, 0.8, 1.1, "quality");
    mgr.create1D("mass_p", "Proton mass;m [GeV];Counts", 100, 0.85, 1.05, "quality");

    // Composite particles
    mgr.create1D("mass_deltaPP", "#Delta^{++} mass;m [GeV];Counts", 100, 1.0, 2.5, "composites");

    // LAB frame
    mgr.create1D("p_theta_lab", "p #theta LAB;#theta [deg];Counts", 90, 0, 90, "lab/angular");

    // CMS frame
    mgr.create1D("cos_theta_deltaPP_cms", "#Delta^{++} cos#theta CMS;cos#theta;Counts",
                 100, -1, 1, "cms/angular");

    // Correlations
    mgr.create2D("dalitz_ppip_npip", "Dalitz;m^{2}(p#pi^{+});m^{2}(n#pi^{+})",
                 100, 1, 4, 100, 1, 4, "correlations");
}
```

#### setup_ntuples.h

```cpp
inline void setupNtuples(Manager& manager, const AnalysisConfig& config) {
    manager.createDynamicNtuple(
        "nt_particles",
        "Particle observables",
        config.getMissingValue(),
        config.getKeepIntermediateTree()
    );

    manager.createDynamicNtuple(
        "nt_compound",
        "Compound particle observables",
        config.getMissingValue(),
        config.getKeepIntermediateTree()
    );
}
```

---

## 7. Complete Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                     CONFIGURATION                                │
│  config.json → AnalysisConfig → parameters for all components   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        SETUP PHASE                               │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐             │
│  │ Beam System │  │ Input Chain │  │ Output File │             │
│  │ (PParticle) │  │(NTupleReader)│ │  (Manager)  │             │
│  └─────────────┘  └─────────────┘  └─────────────┘             │
│         │                │                │                      │
│         │                │         ┌──────┴──────┐               │
│         │                │         │             │               │
│         │                │    Histograms    Ntuples              │
│         │                │  (setupHistograms) (setupNtuples)     │
│         │                │                                       │
│         │                │         CutManager (setupCuts)        │
└─────────┴────────────────┴──────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      EVENT LOOP                                  │
│  for each event:                                                 │
│    1. reader.getEntry(i) → load branches                        │
│    2. Create PParticle objects from branches                    │
│    3. Set kinematics (RECONSTRUCTED/CORRECTED/SIMULATED)        │
│    4. Calculate derived particles (missing mass, composites)    │
│    5. Fill quality histograms (before cuts)                     │
│    6. Apply cuts → early return if fail                         │
│    7. Boost to CMS frame                                        │
│    8. Fill physics histograms (after cuts)                      │
│    9. Fill ntuples                                              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                     FINALIZATION                                 │
│  1. cuts.printCutFlow() → efficiency statistics                 │
│  2. manager.printSummary() → histogram/ntuple counts            │
│  3. manager.closeFile() → write all output to ROOT file         │
│     - Finalize DynamicHNtuples (TTree → TNtuple)                │
│     - Write histograms with folder hierarchy                    │
│     - Close TFile                                               │
└─────────────────────────────────────────────────────────────────┘
```

---

## 8. Key Code Locations

| Component | File | Description |
|-----------|------|-------------|
| Main entry | `main.cc` | Entry point, setup, event loop |
| Event processing | `main.cc:processEvent()` | Per-event physics logic |
| Configuration | `src/analysis_config.h` | JSON parser & accessors |
| PParticle | `src/pparticle.h` | Base particle class |
| PParticleEcal | `src/pparticle_ecal.h` | ECAL extension |
| PParticleFwd | `src/pparticle_fwd.h` | Forward Tracker extension |
| CutManager | `src/cut_manager.h` | Cut system |
| Manager | `src/manager.h` | Output orchestration |
| DynamicHNtuple | `src/dynamic_hntuple.h` | Flexible ntuples |
| HNtuple | `src/hntuple.h` | Classic ntuples |
| HistogramRegistry | `src/histogram_registry.h` | Histogram storage |
| HistogramFactory | `src/histogram_factory.h` | Histogram creation |
| HistogramBuilder | `src/histogram_builder.h` | Fluent histogram API |
| NTupleReader | `src/ntuple_reader.h` | Input reading |
| Setup: Histograms | `src/setup_histograms.h` | Histogram definitions |
| Setup: Cuts | `src/setup_cuts.h` | Cut definitions |
| Setup: Ntuples | `src/setup_ntuples.h` | Ntuple definitions |

---

## 9. Summary of Architecture

### Design Principles

1. **Configuration-driven**: All parameters from JSON, no hardcoding
2. **Separation of concerns**: Physics in processEvent(), infrastructure in classes
3. **Multiple kinematics**: RECONSTRUCTED/CORRECTED/SIMULATED stored simultaneously
4. **Memory safety**: RAII with `unique_ptr` throughout
5. **Flexible output**: DynamicHNtuple allows evolving analysis
6. **Weighted histograms**: Cross-section scaling via luminosity
7. **Organized output**: Folder hierarchy in ROOT files
8. **Statistics tracking**: Built-in cut efficiency monitoring

### Workflow Summary

```
1. Load config.json → AnalysisConfig
2. Create beam system (projectile + target)
3. Open input (NTupleReader)
4. Open output (Manager)
5. Setup: histograms, cuts, ntuples
6. Event loop:
   a. Create particles (PParticle, PParticleEcal, PParticleFwd)
   b. Fill from branches
   c. Apply cuts (CutManager)
   d. Build composites (operator+/-)
   e. Boost to CMS
   f. Fill histograms (weighted)
   g. Fill ntuples
7. Finalize: print statistics, write output
```

### Key Features

- **PParticle arithmetic**: `proton + pion` creates composite, `beam - proton - pion` gives missing mass
- **Frame preservation**: LAB frame always available via `labFrame()` even after boosts
- **Dynamic ntuples**: Add variables anytime, missing values handled automatically
- **Cut flow tracking**: Automatic efficiency calculation and formatted reports
- **Folder organization**: Clean ROOT file structure with folders

---

*Last updated: February 2026*

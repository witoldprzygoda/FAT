# FAT Framework - Developer Tutorial

A complete, step-by-step guide to working with the FAT (Final Analysis Tool) framework.

---

## Table of Contents

1. [Overview](#1-overview)
2. [File Structure](#2-file-structure)
3. [Configuration: config.json](#3-configuration-configjson)
4. [The Event Loop](#4-the-event-loop)
5. [Reading Variables (Reflection)](#5-reading-variables-reflection)
6. [Creating Particles](#6-creating-particles)
7. [Reference Frame Transformations](#7-reference-frame-transformations)
8. [Histogram Management](#8-histogram-management)
9. [Cut Management](#9-cut-management)
10. [Complete Workflow](#10-complete-workflow)
11. [Quick Reference](#11-quick-reference)

---

## 1. Overview

The FAT framework is a modern C++/ROOT toolkit for particle physics analysis. 

### Key Design Principles

| Feature | Approach |
|---------|----------|
| Configuration | External JSON file (`config.json`) |
| Variable reading | **Reflection** - no pre-declaration needed |
| Histograms | Named access via Manager class |
| Particles | PParticle class with 4-momentum algebra |
| Reference frames | BoostFrame class for Lorentz transformations |

---

## 2. File Structure

```
FAT/
├── config.json              # Configuration (input, output, beam)
├── main.cc                  # Physics analysis (processEvent + main)
├── Makefile                 # Build system
│
├── src/                     # Framework source headers
│   ├── analysis_config.h    # JSON configuration loader
│   ├── ntuple_reader.h      # Input data reader (reflection-based)
│   ├── pparticle.h          # Particle class + physics constants
│   ├── boost_frame.h        # Reference frame transformations
│   ├── manager.h            # Output histograms/ntuples
│   ├── cut_manager.h        # Cut definitions + statistics
│   ├── hntuple.h/.cc        # Output ntuple with named columns
│   ├── setup_histograms.h   # Histogram definitions (EDIT THIS)
│   └── setup_cuts.h         # Cut definitions (EDIT THIS)
│
└── examples/                # Example analyses
```

### Files You Edit for a New Analysis

| File | What to Edit |
|------|-------------|
| `config.json` | Input source, tree name, beam energy |
| `src/setup_histograms.h` | Add/modify histograms |
| `src/setup_cuts.h` | Add/modify cuts |
| `main.cc` → `processEvent()` | Physics analysis logic |

---

## 3. Configuration: config.json

The configuration file specifies **WHERE** to read data from and **WHERE** to write output:

```json
{
    "input": {
        "source": "h68_10.list",      // .root file OR .list file (auto-detected)
        "tree_name": "PPip_ID",        // TTree name in the ROOT file
        "start_event": 0,              // First event to process (default: 0)
        "max_events": -1               // -1 = all events
    },
    "output": {
        "filename": "output.root",
        "option": "RECREATE"
    },
    "beam": {
        "kinetic_energy": 1580.0       // MeV
    }
}
```

### Input Source Auto-Detection

- If `source` ends with `.root` → opens as single ROOT file
- If `source` ends with `.list` → reads file list, creates TChain

### File List Format (h68_10.list)

```
# Comments start with #
/path/to/file1.root
/path/to/file2.root
/path/to/file3.root
```

---

## 4. The Event Loop

**Location in main.cc: `main()` function, around lines 490-530**

```cpp
int main(int argc, char* argv[]) {
    // ... setup code ...
    
    // EVENT LOOP
    for (Long64_t i = start_event; i < end_event; ++i) {
        reader.getEntry(i);    // Load event i from tree
        
        processEvent(reader, manager, cuts, beam, projectile, frames, use_corrected);
    }
    
    // ... finalize code ...
}
```

### What Happens Each Iteration

1. `reader.getEntry(i)` - Loads all branch values for event `i`
2. `processEvent()` - Your analysis logic for one event

---

## 5. Reading Variables (Reflection)

**Location in main.cc: `processEvent()` function, top section**

### The Key Concept: No Pre-Declaration Needed!

**Old approach (traditional ROOT):**
```cpp
// Required 3 steps per variable:
Float_t pip_p;                              // 1. Declare
tree->SetBranchAddress("pip_p", &pip_p);    // 2. Bind
// ... then use pip_p                       // 3. Use
```

**New approach (FAT reflection):**
```cpp
// Just ONE line:
double pip_p = reader["pip_p"];    // Automatically binds & reads!
```

### How It Works in processEvent()

```cpp
void processEvent(NTupleReader& reader, Manager& mgr, ...) {
    
    // ============================================================
    // STEP 1: Read all needed variables at the TOP of the function
    // These are the ONLY places you use reader["..."]
    // ============================================================
    
    // Proton kinematics
    double p_p     = reader["p_p"];       // Momentum [MeV/c]
    double p_theta = reader["p_theta"];   // Theta [degrees]
    double p_phi   = reader["p_phi"];     // Phi [degrees]
    
    // Pion kinematics  
    double pip_p     = reader["pip_p"];
    double pip_theta = reader["pip_theta"];
    double pip_phi   = reader["pip_phi"];
    
    // Optional: check if variable exists
    double weight = 1.0;
    if (reader.hasVariable("weight")) {
        weight = reader["weight"];
    }
    
    // ============================================================
    // STEP 2: From here on, use regular C++ variables
    // No more reader["..."] needed!
    // ============================================================
    
    PParticle proton = ParticleFactory::createProton(p_p, p_theta, p_phi);
    // ... rest of analysis ...
}
```

### Important Notes

- The variable names (`"p_p"`, `"pip_theta"`, etc.) must **exactly match** the branch names in your ROOT tree
- To see what branches exist: `root -l file.root` then `tree->Print()`
- Use `reader.hasVariable("name")` to check if a branch exists

---

## 6. Creating Particles

**Location in main.cc: `processEvent()` function, after reading variables**

### Physics Constants (in `src/pparticle.h`)

```cpp
namespace Physics {
    constexpr double MASS_PROTON     = 938.27231;   // MeV/c²
    constexpr double MASS_NEUTRON    = 939.56542;   // MeV/c²
    constexpr double MASS_PION_PLUS  = 139.56995;   // MeV/c²
    constexpr double MASS_PION_MINUS = 139.56995;   // MeV/c²
    constexpr double MASS_ELECTRON   = 0.51099895;  // MeV/c²
    constexpr double MASS_MUON_PLUS  = 105.6583745; // MeV/c²
    // ... more in pparticle.h
}
```

### Creating Particles with ParticleFactory

```cpp
// Detected particles (momentum in MeV/c, angles in degrees):
PParticle proton = ParticleFactory::createProton(p_p, p_theta, p_phi);
PParticle pion   = ParticleFactory::createPiPlus(pip_p, pip_theta, pip_phi);

// Beam particles:
PParticle projectile = ParticleFactory::createBeamProton(1580.0);  // T_kin in MeV
PParticle target     = ParticleFactory::createTargetProton();       // At rest
PParticle beam       = projectile + target;                         // Total system
```

### Four-Momentum Algebra

```cpp
// Addition: composite particle
PParticle deltaPP = proton + pion;           // Δ++ = p + π+

// Subtraction: missing mass technique
PParticle neutron = beam - proton - pion;    // n = beam - p - π+

// Access kinematics:
double mass    = deltaPP.mass();       // MeV/c²
double massGeV = deltaPP.massGeV();    // GeV/c²
double mom     = proton.momentum();    // |p| in MeV/c
double theta   = proton.theta();       // degrees
double phi     = proton.phi();         // degrees
double costh   = proton.cosTheta();    // cos(theta)
double beta    = proton.beta();        // v/c
```

---

## 7. Reference Frame Transformations

**Location in main.cc: `processEvent()` function, after particle creation**

### Setting Up Frames (in main())

```cpp
// Create frame manager
EventFrames frames;

// Define beam CMS frame
frames.setBeamFrame(projectile, target);
```

### Boosting Particles (in processEvent())

```cpp
// Get the beam CMS frame
const BoostFrame& cms = frames.getFrame("beam");

// Boost particles - returns NEW copies (originals unchanged)
PParticle proton_cms  = cms.boost(proton);
PParticle deltaPP_cms = cms.boost(deltaPP);

// Access boosted kinematics
double cos_theta_cms = deltaPP_cms.cosTheta();
```

### Creating Composite Rest Frames

```cpp
// Create rest frame of p+π+ system
BoostFrame ppip_frame(proton + pion);

// Boost pion into p+π+ rest frame (helicity angle)
PParticle pion_in_ppip = ppip_frame.boost(pion);
double helicity_angle = pion_in_ppip.cosTheta();
```

---

## 8. Histogram Management

**Locations in main.cc:**
- Definition: `setupHistograms()` function
- Filling: `processEvent()` function

### Defining Histograms

```cpp
void setupHistograms(Manager& mgr) {
    // 1D histogram: (name, title, bins, xlow, xhigh, folder)
    mgr.create1D("mass_n", "Missing mass;M [GeV];Counts", 
                 200, 0.5, 1.5, "quality");
    
    // 2D histogram: (name, title, xbins, xlow, xhigh, ybins, ylow, yhigh, folder)
    mgr.create2D("dalitz", "Dalitz;M²(pπ+);M²(nπ+)",
                 100, 1.0, 4.0, 100, 1.0, 4.0, "correlations");
}
```

### Filling Histograms

```cpp
void processEvent(...) {
    // ... create particles ...
    
    // Fill 1D:
    mgr.fill("mass_n", neutron.massGeV());
    
    // Fill 2D:
    mgr.fill("dalitz", m2_ppip, m2_npip);
}
```

### Folder Organization

Histograms are organized into ROOT folders in the output file:

```
output.root
├── quality/
│   └── mass_n
├── correlations/
│   └── dalitz
└── cms/
    └── cos_theta
```

---

## 9. Cut Management

**Location in main.cc: cuts defined in `main()`, applied in `processEvent()`**

### Defining Cuts

```cpp
// In main():
CutManager cuts;
cuts.defineRangeCut("neutron_mass", 0.899, 0.986, "Neutron mass window");
cuts.defineRangeCut("deltaPP_mass", 0.8, 1.8, "Delta++ mass window");
```

### Applying Cuts

```cpp
void processEvent(...) {
    // ... create particles, fill quality histograms ...
    
    // Apply cut with early return
    if (!cuts.passRangeCut("neutron_mass", neutron.massGeV())) return;
    
    // Event passed - continue analysis
    mgr.fill("mass_n_cut", neutron.massGeV());
}
```

### Cut Statistics

At the end of analysis:

```cpp
cuts.printCutFlow();

// Output:
// ╔═══════════════════════════════════════════════════════════════╗
// ║                       CUT FLOW SUMMARY                         ║
// ╠═══════════════════════════════════════════════════════════════╣
// ║ neutron_mass    │ Tested: 100000 │ Passed: 85000 │ Eff: 85.0% ║
// ╚═══════════════════════════════════════════════════════════════╝
```

---

## 10. Complete Workflow

### Visual Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                         config.json                              │
│  - input.source: "data.root" or "files.list"                    │
│  - input.tree_name: "PPip_ID"                                   │
│  - beam.kinetic_energy: 1580.0                                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        main.cc: main()                           │
│                                                                  │
│  1. Load config.json                                            │
│  2. Create beam particles                                       │
│  3. Open input via NTupleReader                                 │
│  4. Call setupHistograms()                                      │
│  5. Define cuts                                                 │
│  6. EVENT LOOP: for each event call processEvent()              │
│  7. Print cut statistics                                        │
│  8. Close output file                                           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                   main.cc: processEvent()                        │
│                                                                  │
│  1. READ VARIABLES (top of function):                           │
│     double p_p = reader["p_p"];                                 │
│     double p_theta = reader["p_theta"];                         │
│     ...                                                         │
│                                                                  │
│  2. CREATE PARTICLES:                                           │
│     PParticle proton = ParticleFactory::createProton(p_p, ...); │
│     PParticle neutron = beam - proton - pion;                   │
│                                                                  │
│  3. FILL QUALITY HISTOGRAMS (before cuts)                       │
│                                                                  │
│  4. APPLY CUTS:                                                 │
│     if (!cuts.passRangeCut("name", value)) return;              │
│                                                                  │
│  5. BOOST TO CMS:                                               │
│     PParticle p_cms = frames.getFrame("beam").boost(proton);    │
│                                                                  │
│  6. FILL PHYSICS HISTOGRAMS (after cuts)                        │
└─────────────────────────────────────────────────────────────────┘
```

### Step-by-Step Instructions

1. **Edit config.json** - Set your input file, tree name, beam energy
2. **Edit setupHistograms()** - Define your histograms
3. **Edit processEvent() top** - Read the variables you need via `reader["name"]`
4. **Edit processEvent()** - Create particles, apply cuts, fill histograms
5. **Build**: `make`
6. **Run**: `./ana config.json`

---

## 11. Quick Reference

### Creating Particles
```cpp
PParticle p = ParticleFactory::createProton(mom, theta, phi);
PParticle composite = p1 + p2;
PParticle missing = beam - p1 - p2;
```

### Reading Variables
```cpp
double var = reader["branch_name"];   // Exactly matches tree branch name
if (reader.hasVariable("name")) { }   // Check if exists
```

### Boosting
```cpp
EventFrames frames;
frames.setBeamFrame(projectile, target);
PParticle p_cms = frames.getFrame("beam").boost(particle);
```

### Histograms
```cpp
mgr.create1D("name", "title;X;Y", bins, lo, hi, "folder");
mgr.fill("name", value);
```

### Cuts
```cpp
cuts.defineRangeCut("name", min, max);
if (!cuts.passRangeCut("name", value)) return;
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Variable 'xxx' not found` | Branch name must exactly match tree. Check with `tree->Print()` |
| `Histogram already exists` | Each histogram name must be unique |
| `File not found` | Check paths in config.json |
| Build errors | Run `make clean && make` |

### Debug Tips

```cpp
// See available variables in tree
std::vector<std::string> vars = reader.listVariables();
for (const auto& v : vars) std::cout << v << "\n";

// Print particle info
proton.print();

// Print config
config.print();
```

---

*Last updated: November 2025*
*Author: Witold Przygoda (witold.przygoda@uj.edu.pl)*

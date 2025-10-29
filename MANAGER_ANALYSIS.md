# Manager Class Analysis and Refactoring Proposal

## Executive Summary

The current `Manager` class provides basic histogram and ntuple management but has several architectural issues that limit maintainability, testability, and scalability. This document proposes a comprehensive refactoring using modern C++ patterns (Factory, Registry, Builder) and smart pointers.

---

## Current Architecture Analysis

### File Structure
```
src/manager.h          - Base Manager class (43 lines)
src/manager.cc         - Implementation (143 lines)
src/datamanager.h      - Derived class (26 lines)
datamanager.cc         - Initialization (149 lines)
data.h                 - Global histogram pointers (120+ lines)
```

---

## Critical Issues Identified

### 1. **Global State Pollution** (CRITICAL)

**Location**: `data.h` lines 23-74

```cpp
// 30+ global histogram pointers!
inline TH1F *mass_deltaPP {nullptr};
inline TH1F *mass_deltaP {nullptr};
inline TH1F *cos_theta_deltaPP {nullptr};
inline TH1F *mass_n_tab[40] {nullptr};      // Array of pointers!
inline TH1F *mass_p_tab[100] {nullptr};     // 100 pointers!
inline TH1F *miss_mass_tab[25][20] {nullptr}; // 500 pointers!
inline HNtuple *nt {nullptr};
```

**Problems**:
- ❌ Global state makes testing impossible
- ❌ Unclear initialization order
- ❌ No encapsulation
- ❌ Can't have multiple analyses simultaneously
- ❌ Name collisions across files

**Impact**: **HIGH** - Affects every part of the codebase

---

### 2. **Manual Memory Management** (CRITICAL)

**Location**: `manager.cc` lines 92-112

```cpp
Manager::~Manager()
{
    writeData();
    for (auto &histogram : histograms) {
        delete histogram;  // Manual deletion
    }
    for (auto &histogram : histograms2) {
        delete histogram;
    }
    for (auto &histogram : histograms3) {
        delete histogram;
    }
    // for (auto& ntuple : ntuples) {
    //     delete ntuple;  // COMMENTED OUT - MEMORY LEAK!
    // }
    outfile->Close();
    delete outfile;
}
```

**Problems**:
- ❌ Ntuple deletion commented out (line 107-109) → **Memory leak!**
- ❌ No RAII - exception unsafe
- ❌ Raw pointers everywhere
- ❌ Manual new/delete error-prone

**Impact**: **CRITICAL** - Active memory leak

---

### 3. **Three Separate Containers** (MEDIUM)

**Location**: `manager.h` lines 38-41

```cpp
std::vector<TH1F*> histograms;   // 1D histograms
std::vector<TH2F*> histograms2;  // 2D histograms
std::vector<TH3F*> histograms3;  // 3D histograms
std::vector<HNtuple*> ntuples;
```

**Problems**:
- ❌ Code duplication (3× write loop, 3× delete loop)
- ❌ Can't iterate over all histograms uniformly
- ❌ Adding TProfile, THnSparse requires new vector

**Impact**: **MEDIUM** - Maintenance burden

---

### 4. **Dangerous const_cast** (HIGH)

**Location**: `manager.cc` line 10

```cpp
TH1F *Manager::createHistogram(const TH1F* ptr)
{
    TH1F *histogram = const_cast<TH1F*>(ptr);  // ⚠️ DANGER!
    histograms.push_back(histogram);
    return histogram;
}
```

**Problems**:
- ❌ Violates const correctness
- ❌ Undefined behavior if ptr is actually const
- ❌ Unclear why this function exists

**Impact**: **HIGH** - Potential undefined behavior

---

### 5. **Repetitive Histogram Creation** (MEDIUM)

**Location**: `datamanager.cc` lines 107-131

```cpp
// Creating 40 histograms manually!
for (int i = 0; i < 40; ++i)
{
    sprintf(name, "mass_n_tab_%d", i);
    mass_n_tab[i] = createHistogram(name, name, 100, 0.5, 1.5);
    sprintf(name, "mass_nn_tab_%d", i);
    mass_nn_tab[i] = createHistogram(name, name, 100, 0.5, 1.5);
}

// Creating 100 more histograms!
for (int i = 0; i < 100; ++i)
{
    sprintf(namep, "mass_p_tab_%d", i);
    mass_p_tab[i] = createHistogram(namep, namep, 100, 0.5, 1.5);
    sprintf(namep, "mass_pp_tab_%d", i);
    mass_pp_tab[i] = createHistogram(namep, namep, 100, 0.5, 1.5);
}

// Creating 500 more histograms (25×20 matrix)!
for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 25; ++j) {
        sprintf(name2dim, "miss_mass_invmass_costh_%d_%d", j + 1, i + 1);
        miss_mass_tab[j][i] = createHistogram(name2dim, name2dim, 100, 0.5, 1.5);
    }
}
```

**Total**: **640+ histograms** created with repetitive code!

**Problems**:
- ❌ Repetitive sprintf calls
- ❌ Hard to change binning for all at once
- ❌ Can't share common configuration
- ❌ No metadata about what each histogram represents

**Impact**: **MEDIUM** - Maintenance burden

---

### 6. **No Organizational Structure** (MEDIUM)

**Problem**: All histograms flat in one file

```
output.root
├── mass_deltaPP          ← What analysis step?
├── mass_deltaP           ← What cut level?
├── pwa_pip_costh         ← PWA analysis
├── exp_sig               ← Experimental?
├── mass_n_tab_0          ← Part of array
├── mass_n_tab_1
├── ...
└── mass_n_tab_39         ← Which belongs where?
```

**Impact**: **MEDIUM** - Difficult to navigate output

---

### 7. **No Histogram Metadata** (LOW)

Current histograms only have:
- Name
- Title

Missing:
- ❌ Which analysis step created it?
- ❌ What cuts were applied?
- ❌ Physics significance?
- ❌ Units?
- ❌ Tags for categorization?

**Impact**: **LOW** - Quality of life

---

### 8. **Code Duplication in createHistogram** (LOW)

**Location**: `manager.cc` lines 15-49

Four nearly identical functions:
```cpp
TH1F *Manager::createHistogram(...) {
    TH1F *histogram = new TH1F(...);
    if (sumw2) histogram->Sumw2();  // Repeated 4 times
    histograms.push_back(histogram); // Different vector each time
    return histogram;
}
```

**Impact**: **LOW** - Maintenance burden

---

## Proposed Architecture: Modern Manager

### Design Principles

1. **No global state** - Everything encapsulated
2. **Smart pointers** - RAII, no manual delete
3. **Unified container** - Polymorphism via TH1* base
4. **Registry pattern** - Named histogram access
5. **Factory pattern** - Flexible creation
6. **Builder pattern** - Complex configurations
7. **Folder support** - Logical organization
8. **Metadata** - Physics context attached

---

## File Organization

### New Structure
```
src/HistogramFactory.h     - Factory for creating histograms
src/HistogramRegistry.h    - Registry for managing histograms
src/HistogramBuilder.h     - Builder for complex configurations
src/ImprovedManager.h      - Modernized Manager class
src/ImprovedManager.cc     - Implementation
```

---

## Implementation Preview

### Key Components

#### 1. Histogram Metadata
```cpp
struct HistogramMetadata {
    std::string name;
    std::string title;
    std::string folder;           // Logical grouping
    std::string description;      // Physics context
    std::vector<std::string> tags; // Categorization
    std::string units;
    std::string cutLevel;         // Which cuts applied
};
```

#### 2. Histogram Registry
```cpp
class HistogramRegistry {
public:
    // Add histogram with metadata
    void add(std::unique_ptr<TH1> hist, HistogramMetadata meta);

    // Access by name
    TH1* get(const std::string& name);
    TH1* operator[](const std::string& name);

    // Query
    std::vector<std::string> listByFolder(const std::string& folder);
    std::vector<std::string> listByTag(const std::string& tag);

    // Write to file
    void writeToFile(TFile* file);

private:
    std::map<std::string, std::unique_ptr<TH1>> histograms_;
    std::map<std::string, HistogramMetadata> metadata_;
};
```

#### 3. Histogram Factory
```cpp
class HistogramFactory {
public:
    // Create 1D histogram
    static std::unique_ptr<TH1F> create1D(
        const std::string& name,
        const std::string& title,
        int nbins, double xlow, double xup,
        bool sumw2 = true
    );

    // Create array of histograms
    static std::vector<std::unique_ptr<TH1F>> create1DArray(
        const std::string& basename,
        int count,
        int nbins, double xlow, double xup
    );

    // Create 2D matrix of histograms
    static std::vector<std::vector<std::unique_ptr<TH1F>>> create1DMatrix(
        const std::string& basename,
        int rows, int cols,
        int nbins, double xlow, double xup
    );
};
```

#### 4. Histogram Builder
```cpp
class HistogramBuilder {
public:
    HistogramBuilder& name(const std::string& n);
    HistogramBuilder& title(const std::string& t);
    HistogramBuilder& bins(int n, double low, double up);
    HistogramBuilder& folder(const std::string& f);
    HistogramBuilder& tag(const std::string& t);
    HistogramBuilder& description(const std::string& d);
    HistogramBuilder& sumw2(bool enable = true);

    std::unique_ptr<TH1F> build1D();
    std::unique_ptr<TH2F> build2D(int nbinsy, double ylow, double yup);
};
```

---

## Usage Comparison

### Current (Old) Way

```cpp
// In data.h - declare global
inline TH1F *mass_deltaPP {nullptr};

// In datamanager.cc initData() - create manually
mass_deltaPP = createHistogram("mass_deltaPP", "mass_deltaPP", 50, 0.8, 1.8);

// In PPip_ID.cc - use global
mass_deltaPP->Fill(deltaPP->M()/1000., weight);

// Memory leak risk, no organization, no metadata
```

### Proposed (New) Way - Option 1: Registry

```cpp
// In ImprovedManager - no globals needed!
class ImprovedManager {
protected:
    HistogramRegistry histograms_;

    void initHistograms() {
        histograms_.add(
            HistogramFactory::create1D("mass_deltaPP", "Delta++ Mass", 50, 0.8, 1.8),
            {.folder = "Resonances",
             .description = "Invariant mass of p+pi+ system",
             .tags = {"Delta++", "resonance"},
             .units = "GeV/c^2"}
        );
    }
};

// In analysis code - access via registry
histograms_["mass_deltaPP"]->Fill(deltaPP.massGeV(), weight);

// Or with error checking
if (auto* h = histograms_.get("mass_deltaPP")) {
    h->Fill(deltaPP.massGeV(), weight);
}
```

### Proposed (New) Way - Option 2: Builder

```cpp
// Fluent interface for complex configurations
histograms_.add(
    HistogramBuilder()
        .name("mass_deltaPP")
        .title("Delta++ Invariant Mass")
        .bins(50, 0.8, 1.8)
        .folder("Resonances")
        .tag("Delta++")
        .tag("p+pi+")
        .description("Invariant mass of proton-pion system")
        .sumw2(true)
        .build1D()
);
```

### Proposed (New) Way - Option 3: Array Creation

```cpp
// Instead of 40 lines of sprintf + createHistogram:
auto mass_n_histograms = HistogramFactory::create1DArray(
    "mass_n_tab",  // basename
    40,            // count
    100, 0.5, 1.5  // binning
);

// Register them
for (int i = 0; i < 40; i++) {
    histograms_.add(
        std::move(mass_n_histograms[i]),
        {.folder = "MissingMass",
         .tag = "neutron",
         .cutLevel = "bin_" + std::to_string(i)}
    );
}
```

---

## Benefits Summary

### 1. **No Global State**
```cpp
// Before: 30+ global pointers
inline TH1F *mass_deltaPP {nullptr};

// After: Encapsulated in Manager
class ImprovedManager {
    HistogramRegistry histograms_;  // All contained here
};
```

### 2. **RAII & Smart Pointers**
```cpp
// Before: Manual delete, memory leaks
for (auto &histogram : histograms) {
    delete histogram;  // Easy to forget
}
// delete ntuple;  // COMMENTED OUT - LEAK!

// After: Automatic cleanup
std::unique_ptr<TH1> hist;  // Deletes itself
```

### 3. **Unified Container**
```cpp
// Before: 3 separate vectors
std::vector<TH1F*> histograms;
std::vector<TH2F*> histograms2;
std::vector<TH3F*> histograms3;

// After: One container via polymorphism
std::map<std::string, std::unique_ptr<TH1>> histograms_;
// Works for TH1F, TH2F, TH3F, TProfile, THnSparse, etc.
```

### 4. **Named Access**
```cpp
// Before: Must declare global variable
inline TH1F *mass_deltaPP {nullptr};
mass_deltaPP->Fill(...);

// After: Access by name
histograms_["mass_deltaPP"]->Fill(...);
```

### 5. **Organization**
```cpp
// Before: Flat list of 640+ histograms
output.root
├── mass_deltaPP
├── mass_n_tab_0
├── ...

// After: Organized folders
output.root
├── Resonances/
│   ├── mass_deltaPP
│   └── mass_deltaP
├── Kinematics/
│   ├── p_momentum
│   └── pip_theta
└── MissingMass/
    ├── mass_n_tab_0
    └── ...
```

### 6. **Metadata**
```cpp
// Before: Just name and title
TH1F("mass_deltaPP", "mass_deltaPP", ...)

// After: Rich context
{
  .name = "mass_deltaPP",
  .title = "Delta++ Invariant Mass",
  .folder = "Resonances",
  .description = "Invariant mass of p+pi+ system in GeV/c^2",
  .tags = ["Delta++", "p+pi+", "resonance"],
  .units = "GeV/c^2",
  .cutLevel = "after_vertex_cut"
}
```

### 7. **Reduced Code Duplication**
```cpp
// Before: 40 lines to create 40 histograms
for (int i = 0; i < 40; ++i) {
    sprintf(name, "mass_n_tab_%d", i);
    mass_n_tab[i] = createHistogram(name, name, 100, 0.5, 1.5);
}

// After: 3 lines
auto hists = HistogramFactory::create1DArray("mass_n_tab", 40, 100, 0.5, 1.5);
for (auto& h : hists) histograms_.add(std::move(h));
```

---

## Migration Strategy

### Phase 1: Add New Infrastructure (Backward Compatible)
1. Add `HistogramRegistry.h`
2. Add `HistogramFactory.h`
3. Add `HistogramBuilder.h`
4. Keep old Manager working

### Phase 2: Create ImprovedManager (Parallel)
1. Create `ImprovedManager` inheriting from `Manager`
2. Implement registry internally
3. Forward `createHistogram()` to registry
4. Old code still works!

### Phase 3: Migrate Gradually
1. Move one histogram at a time to registry
2. Test each migration
3. Remove global pointers as you go

### Phase 4: Clean Up
1. Remove old `Manager` implementation
2. Rename `ImprovedManager` → `Manager`
3. Remove `data.h` global pointers

---

## Priorities

### HIGH Priority (Fix Now)
1. ✅ Fix ntuple memory leak (uncomment delete)
2. ✅ Replace raw pointers with smart pointers
3. ✅ Remove dangerous const_cast

### MEDIUM Priority (Next Iteration)
4. ✅ Implement HistogramRegistry
5. ✅ Implement HistogramFactory
6. ✅ Unify containers (use TH1* base)
7. ✅ Add folder support

### LOW Priority (Nice to Have)
8. ✅ Add metadata support
9. ✅ Implement Builder pattern
10. ✅ Add configuration file support

---

## Estimated Impact

### Code Reduction
- **datamanager.cc**: 149 lines → ~80 lines (47% reduction)
- **data.h**: 120 lines → ~20 lines (83% reduction)
- **manager.cc**: 143 lines → ~200 lines (initial increase, but more features)

### Maintainability
- **Before**: Must track 30+ global variables
- **After**: One registry, named access

### Safety
- **Before**: Manual memory management, leaks
- **After**: RAII, automatic cleanup

### Flexibility
- **Before**: Hard to add TProfile, THnSparse
- **After**: Works with any TH1-derived class

---

## Next Steps

Would you like me to:

1. **Implement HistogramRegistry** with full code?
2. **Implement HistogramFactory** for array creation?
3. **Implement HistogramBuilder** with fluent interface?
4. **Create ImprovedManager** that's backward compatible?
5. **Provide migration examples** for your specific histograms?

I can create production-ready code for any or all of these!

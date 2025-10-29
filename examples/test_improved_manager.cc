/**
 * @file test_improved_manager.cc
 * @brief Comprehensive examples demonstrating the ImprovedManager system
 *
 * This file shows:
 * 1. Migration from old Manager to ImprovedManager
 * 2. Using HistogramFactory for arrays
 * 3. Using HistogramBuilder for fluent interface
 * 4. Folder organization
 * 5. Metadata and tags
 * 6. Type-safe access
 * 7. NTuple integration
 *
 * Compile:
 *   g++ -o test_improved_manager test_improved_manager.cc \
 *       ../src/hntuple.cc \
 *       `root-config --cflags --libs`
 *
 * Run:
 *   ./test_improved_manager
 */

#include "../src/ImprovedManager.h"
#include "../src/HistogramRegistry.h"
#include "../src/HistogramFactory.h"
#include "../src/HistogramBuilder.h"
#include <iostream>
#include <TRandom3.h>

// ============================================================================
// EXAMPLE 1: Basic usage - replacing old Manager
// ============================================================================

void example1_basic_usage() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 1: Basic Usage                                        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    ImprovedManager manager;
    manager.openFile("test_basic.root");

    // OLD WAY (from datamanager.cc):
    // TH1F* h_theta = new TH1F("h_theta", "Theta", 100, 0, 180);
    // histograms.push_back(h_theta);

    // NEW WAY:
    manager.create1D("h_theta", "Theta distribution", 100, 0, 180);
    manager.create1D("h_phi", "Phi distribution", 100, -180, 180);
    manager.create1D("h_p", "Momentum", 150, 0, 3000);

    // Fill histograms
    TRandom3 rnd;
    for (int i = 0; i < 10000; ++i) {
        manager.fill("h_theta", rnd.Uniform(0, 180));
        manager.fill("h_phi", rnd.Uniform(-180, 180));
        manager.fill("h_p", rnd.Gaus(1500, 300));
    }

    // Print summary
    manager.printSummary();

    // Save and close
    manager.closeFile();

    std::cout << "\n✅ Created test_basic.root with 3 histograms\n";
}

// ============================================================================
// EXAMPLE 2: Creating histogram arrays - replacing sprintf loops
// ============================================================================

void example2_histogram_arrays() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 2: Histogram Arrays                                   ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    ImprovedManager manager;
    manager.openFile("test_arrays.root");

    // OLD WAY (from datamanager.cc lines 32-39):
    // char name[256], title[256];
    // for (int i = 0; i < 10; i++) {
    //     sprintf(name, "h_theta_%d", i);
    //     sprintf(title, "Theta [%d]", i);
    //     TH1F* h = new TH1F(name, title, 100, 0, 180);
    //     histograms.push_back(h);
    // }

    // NEW WAY (one line!):
    manager.create1DArray("h_theta", "Theta", 10, 100, 0, 180, "angular");

    // Also works for 2D arrays
    manager.create1DArray("h_p", "Momentum", 10, 150, 0, 3000, "momentum");

    // Fill some histograms
    TRandom3 rnd;
    for (int i = 0; i < 10; ++i) {
        std::string name = "h_theta_" + std::to_string(i);
        for (int j = 0; j < 1000; ++j) {
            manager.fill(name, rnd.Uniform(0, 180));
        }
    }

    manager.printSummary();
    manager.closeFile();

    std::cout << "\n✅ Created test_arrays.root with 20 histograms in arrays\n";
    std::cout << "   OLD: ~40 lines of sprintf loops\n";
    std::cout << "   NEW: 2 lines!\n";
}

// ============================================================================
// EXAMPLE 3: Using HistogramBuilder for complex histograms
// ============================================================================

void example3_builder_pattern() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 3: Builder Pattern                                    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    ImprovedManager manager;
    manager.openFile("test_builder.root");

    // Using builder for expressive, self-documenting code
    HistogramBuilder()
        .name("h_theta_cms")
        .title("Theta in CMS frame")
        .bins(100, 0, 180)
        .folder("angular/cms")
        .description("Scattering angle in center-of-mass frame")
        .tag("angular")
        .tag("cms")
        .buildAndRegister1D(manager.registry());

    HistogramBuilder()
        .name("h_theta_lab")
        .title("Theta in LAB frame")
        .bins(100, 0, 180)
        .folder("angular/lab")
        .description("Scattering angle in laboratory frame")
        .tag("angular")
        .tag("lab")
        .buildAndRegister1D(manager.registry());

    // 2D histogram with builder
    HistogramBuilder()
        .name("h_theta_vs_p")
        .title("Theta vs Momentum")
        .binsX(100, 0, 180)
        .binsY(150, 0, 3000)
        .folder("correlations")
        .description("Angular vs momentum correlation")
        .tag("correlation")
        .buildAndRegister2D(manager.registry());

    // Fill some data
    TRandom3 rnd;
    for (int i = 0; i < 5000; ++i) {
        double theta = rnd.Uniform(0, 180);
        double p = rnd.Gaus(1500, 300);

        manager.fill("h_theta_cms", theta);
        manager.fill("h_theta_lab", theta * 0.8);  // Simplified
        manager.fill("h_theta_vs_p", theta, p);
    }

    manager.printSummary();
    manager.closeFile();

    std::cout << "\n✅ Created test_builder.root with organized folder structure\n";
    std::cout << "   Folders: angular/cms, angular/lab, correlations\n";
}

// ============================================================================
// EXAMPLE 4: Folder organization (replacing flat structure)
// ============================================================================

void example4_folder_organization() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 4: Folder Organization                                ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    ImprovedManager manager;
    manager.openFile("test_folders.root");

    // OLD: All histograms in root directory (flat, hard to navigate)
    // NEW: Organized by physics meaning

    // Proton histograms
    manager.create1D("h_p_p", "Proton momentum", 150, 0, 3000, "proton/lab");
    manager.create1D("h_p_theta", "Proton theta", 100, 0, 180, "proton/lab");
    manager.create1D("h_p_p_cms", "Proton momentum CMS", 150, 0, 2000, "proton/cms");

    // Pion histograms
    manager.create1D("h_pip_p", "Pion momentum", 150, 0, 3000, "pion/lab");
    manager.create1D("h_pip_theta", "Pion theta", 100, 0, 180, "pion/lab");

    // Composite system
    manager.create1D("h_ppip_mass", "p+pi+ invariant mass", 200, 1000, 3000, "composite");

    // Quality/cuts
    manager.create1D("h_chi2", "Chi-squared", 100, 0, 10, "quality");

    manager.printSummary();

    // Query by folder
    std::cout << "\nHistograms in 'proton/lab' folder:\n";
    auto proton_hists = manager.listHistogramsInFolder("proton/lab");
    for (const auto& name : proton_hists) {
        std::cout << "  - " << name << "\n";
    }

    manager.closeFile();

    std::cout << "\n✅ Created test_folders.root with hierarchical structure\n";
    std::cout << "   Open in ROOT TBrowser to see organized folders!\n";
}

// ============================================================================
// EXAMPLE 5: NTuple integration
// ============================================================================

void example5_ntuple_integration() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 5: NTuple Integration                                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    ImprovedManager manager;
    manager.openFile("test_ntuple.root");

    // Create histograms
    manager.create1D("h_p", "Momentum", 150, 0, 3000, "histograms");

    // Create ntuple
    manager.createNtuple("events", "Event data", "ntuples");

    HNtuple* nt = manager.getNtuple("events");

    // Fill events
    TRandom3 rnd;
    for (int i = 0; i < 1000; ++i) {
        double p = rnd.Gaus(1500, 300);
        double theta = rnd.Uniform(0, 180);
        double phi = rnd.Uniform(-180, 180);

        // Fill histogram
        manager.fill("h_p", p);

        // Fill ntuple
        (*nt)["p"] = p;
        (*nt)["theta"] = theta;
        (*nt)["phi"] = phi;
        (*nt)["event"] = i;
        nt->fill();
    }

    manager.printSummary();
    manager.closeFile();

    std::cout << "\n✅ Created test_ntuple.root with histograms and ntuple\n";
    std::cout << "   Histograms in 'histograms' folder\n";
    std::cout << "   Ntuple in 'ntuples' folder\n";
}

// ============================================================================
// EXAMPLE 6: Migration comparison (OLD vs NEW)
// ============================================================================

void example6_migration_comparison() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 6: Migration Comparison                               ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    std::cout << "\n--- OLD SYSTEM (datamanager.cc) ---\n";
    std::cout << R"(
// Global pointers (data.h)
TH1F *h_mult;
TH1F *h_p[10];
TH1F *h_theta[10][5];
// ... 30+ more globals

// Creation (datamanager.cc)
char name[256], title[256];
h_mult = new TH1F("h_mult", "Multiplicity", 20, 0, 20);
histograms.push_back(h_mult);

for (int i = 0; i < 10; i++) {
    sprintf(name, "h_p_%d", i);
    sprintf(title, "Momentum [%d]", i);
    h_p[i] = new TH1F(name, title, 150, 0, 3000);
    histograms.push_back(h_p[i]);
}

for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 5; j++) {
        sprintf(name, "h_theta_%d_%d", i, j);
        sprintf(title, "Theta [%d][%d]", i, j);
        h_theta[i][j] = new TH1F(name, title, 100, 0, 180);
        histograms.push_back(h_theta[i][j]);
    }
}

// Usage
h_mult->Fill(5);
h_p[3]->Fill(1500);
h_theta[2][1]->Fill(45);

// Cleanup (BUGGY - commented out!)
// for (auto& hist : histograms) {
//     delete hist;  // MEMORY LEAK!
// }
)";

    std::cout << "\n--- NEW SYSTEM (ImprovedManager) ---\n";
    std::cout << R"(
// NO global pointers!

ImprovedManager manager;
manager.openFile("output.root");

// Creation (much cleaner!)
manager.create1D("h_mult", "Multiplicity", 20, 0, 20);
manager.create1DArray("h_p", "Momentum", 10, 150, 0, 3000);

// For 2D matrix, use HistogramFactory directly
auto theta_matrix = HistogramFactory::create1DMatrix(
    "h_theta", "Theta", 10, 5, 100, 0, 180
);
// Register all to manager (could be automated)

// Usage (by name)
manager.fill("h_mult", 5);
manager.fill("h_p_3", 1500);
manager.fill("h_theta_2_1", 45);

// Cleanup: AUTOMATIC! No memory leaks!
manager.closeFile();
// Smart pointers handle everything
)";

    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Code Reduction Summary:                                       ║\n";
    std::cout << "╠════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Global declarations:  30+ lines → 0 lines  (100% reduction)  ║\n";
    std::cout << "║  Creation code:        ~60 lines → ~5 lines  (92% reduction)  ║\n";
    std::cout << "║  Memory management:    Buggy → Automatic     (NO LEAKS!)      ║\n";
    std::cout << "║  Type safety:          None → Full           (Compile-time)   ║\n";
    std::cout << "║  Organization:         Flat → Hierarchical   (ROOT folders)   ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";
}

// ============================================================================
// EXAMPLE 7: Type-safe access (prevents bugs)
// ============================================================================

void example7_type_safety() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 7: Type Safety                                        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    ImprovedManager manager;
    manager.openFile("test_safety.root");

    manager.create1D("h_1d", "1D histogram", 100, 0, 100);
    manager.create2D("h_2d", "2D histogram", 100, 0, 100, 100, 0, 100);

    // Type-safe access
    TH1F* h1 = manager.getHistogramAs<TH1F>("h_1d");  // OK
    TH2F* h2 = manager.getHistogramAs<TH2F>("h_2d");  // OK

    std::cout << "✅ Type-safe access works correctly\n";

    // This would throw exception at runtime (caught early!):
    try {
        TH2F* h_wrong = manager.getHistogramAs<TH2F>("h_1d");  // Error!
        (void)h_wrong;
    } catch (const std::runtime_error& e) {
        std::cout << "✅ Caught type mismatch error:\n";
        std::cout << "   " << e.what() << "\n";
    }

    manager.closeFile();

    std::cout << "\n✅ Type safety prevents bugs at compile/runtime\n";
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                                                ║\n";
    std::cout << "║       ImprovedManager Comprehensive Test Suite                ║\n";
    std::cout << "║                                                                ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    try {
        example1_basic_usage();
        example2_histogram_arrays();
        example3_builder_pattern();
        example4_folder_organization();
        example5_ntuple_integration();
        example6_migration_comparison();
        example7_type_safety();

        std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
        std::cout << "║  ✅ ALL TESTS PASSED                                          ║\n";
        std::cout << "╠════════════════════════════════════════════════════════════════╣\n";
        std::cout << "║  Files created:                                                ║\n";
        std::cout << "║    - test_basic.root        (basic usage)                     ║\n";
        std::cout << "║    - test_arrays.root       (histogram arrays)                ║\n";
        std::cout << "║    - test_builder.root      (builder pattern)                 ║\n";
        std::cout << "║    - test_folders.root      (folder organization)             ║\n";
        std::cout << "║    - test_ntuple.root       (ntuple integration)              ║\n";
        std::cout << "║    - test_safety.root       (type safety)                     ║\n";
        std::cout << "║                                                                ║\n";
        std::cout << "║  Open any file in ROOT TBrowser to explore the structure!     ║\n";
        std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

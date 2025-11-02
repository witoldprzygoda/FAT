/**
 * @file test_hntuple_improved_errors.cc
 * @brief Demonstrates improved error handling in HNtuple
 *
 * This example shows:
 * 1. Correct usage of HNtuple (lazy construction)
 * 2. Error messages when trying to add variables after freeze
 * 3. New diagnostic methods
 */

#include "../src/hntuple.h"
#include <TFile.h>
#include <iostream>

void example1_correct_usage() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 1: Correct Usage (Lazy Construction)                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    TFile* outfile = new TFile("test_output.root", "RECREATE");

    HNtuple ntuple("my_ntuple", "Test NTuple");
    ntuple.setFile(outfile);

    // Define variables BEFORE first fill()
    std::cout << "Setting variables before first fill()...\n";
    ntuple["energy"] = 100.5;
    ntuple["momentum"] = 50.2;
    ntuple["theta"] = 45.0;

    // Check status before freeze
    std::cout << "\nBefore fill() - Status: "
              << (ntuple.isFrozen() ? "FROZEN" : "UNFROZEN") << "\n";
    std::cout << "Number of variables: " << ntuple.getNVariables() << "\n\n";

    // First fill() - this FREEZES the structure
    std::cout << "Calling first fill()...\n";
    ntuple.fill();

    // Check status after freeze
    std::cout << "\nAfter fill() - Status: "
              << (ntuple.isFrozen() ? "FROZEN" : "UNFROZEN") << "\n";
    std::cout << "Number of variables: " << ntuple.getNVariables() << "\n";

    // Fill more events (this is OK - same structure)
    std::cout << "\nFilling more events with same structure...\n";
    for (int i = 0; i < 5; i++) {
        ntuple["energy"] = 100.0 + i;
        ntuple["momentum"] = 50.0 + i * 2;
        ntuple["theta"] = 45.0 + i * 5;
        ntuple.fill();
    }

    std::cout << "Filled 5 more events successfully!\n";

    // Print structure
    std::cout << "\nFinal structure:\n";
    ntuple.printStructure();

    outfile->Write();
    outfile->Close();
    delete outfile;
}

void example2_error_add_after_freeze() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 2: ERROR - Adding Variable After Freeze              ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    TFile* outfile = new TFile("test_output2.root", "RECREATE");

    HNtuple ntuple("my_ntuple2", "Test NTuple 2");
    ntuple.setFile(outfile);

    // Define initial variables
    ntuple["energy"] = 100.5;
    ntuple["momentum"] = 50.2;

    // First fill() - FREEZES structure
    std::cout << "Freezing structure with first fill()...\n";
    ntuple.fill();

    // Try to add NEW variable after freeze (THIS WILL FAIL!)
    std::cout << "\nAttempting to add 'phi' variable AFTER freeze...\n";
    try {
        ntuple["phi"] = 30.0;  // ERROR!
    }
    catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
    }

    // But we CAN modify existing variables
    std::cout << "\nModifying existing 'energy' variable (this is OK)...\n";
    ntuple["energy"] = 200.0;  // This is fine
    ntuple.fill();
    std::cout << "Success! Existing variables can be modified.\n";

    outfile->Close();
    delete outfile;
}

void example3_query_methods() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 3: New Query/Diagnostic Methods                      ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    TFile* outfile = new TFile("test_output3.root", "RECREATE");

    HNtuple ntuple("physics_data", "Physics Analysis NTuple");
    ntuple.setFile(outfile);

    // Add several variables
    ntuple["p_p"] = 1580.0;
    ntuple["p_theta"] = 45.0;
    ntuple["p_phi"] = 30.0;
    ntuple["pip_p"] = 850.0;
    ntuple["pip_theta"] = 60.0;
    ntuple["pip_phi"] = 120.0;
    ntuple["ppip_m"] = 2200.0;

    std::cout << "Before freeze:\n";
    std::cout << "  isFrozen():      " << (ntuple.isFrozen() ? "YES" : "NO") << "\n";
    std::cout << "  getNVariables(): " << ntuple.getNVariables() << "\n";

    // Check if variables exist
    std::cout << "\nChecking variables:\n";
    std::cout << "  hasVariable(\"p_p\"):    " << (ntuple.hasVariable("p_p") ? "YES" : "NO") << "\n";
    std::cout << "  hasVariable(\"mass\"):   " << (ntuple.hasVariable("mass") ? "YES" : "NO") << "\n";

    // Get variable names
    std::cout << "\nVariable names:\n";
    auto names = ntuple.getVariableNames();
    for (size_t i = 0; i < names.size(); i++) {
        std::cout << "  [" << i << "] " << names[i] << "\n";
    }

    // Freeze
    ntuple.fill();

    std::cout << "\nAfter freeze:\n";
    std::cout << "  isFrozen():      " << (ntuple.isFrozen() ? "YES" : "NO") << "\n";

    // Print full structure
    std::cout << "\n";
    ntuple.printStructure();

    // Get structure as string (useful for logging)
    std::string structure = ntuple.getStructureString();
    std::cout << "\nStructure string can be saved to logs or files.\n";

    outfile->Close();
    delete outfile;
}

void example4_error_no_variables() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 4: ERROR - Filling Without Variables                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    TFile* outfile = new TFile("test_output4.root", "RECREATE");

    HNtuple ntuple("empty_ntuple", "Empty NTuple");
    ntuple.setFile(outfile);

    // Try to fill without adding any variables (THIS WILL FAIL!)
    std::cout << "Attempting to fill() without setting any variables...\n";
    try {
        ntuple.fill();  // ERROR!
    }
    catch (const std::exception& e) {
        std::cout << "\n" << e.what() << std::endl;
    }

    outfile->Close();
    delete outfile;
}

void example5_typo_detection() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  EXAMPLE 5: Detecting Typos in Variable Names                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n\n";

    TFile* outfile = new TFile("test_output5.root", "RECREATE");

    HNtuple ntuple("data", "Data NTuple");
    ntuple.setFile(outfile);

    // Define variables
    ntuple["energy"] = 100.0;
    ntuple["momentum"] = 50.0;
    ntuple["theta"] = 45.0;

    // Freeze
    ntuple.fill();

    // Simulate common typo: "enrgy" instead of "energy"
    std::cout << "Attempting to set 'enrgy' (typo for 'energy')...\n";
    try {
        ntuple["enrgy"] = 200.0;  // Typo!
    }
    catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
        std::cout << "\nThe error message shows all valid variables,\n";
        std::cout << "making it easy to spot typos!\n";
    }

    outfile->Close();
    delete outfile;
}

int main() {
    std::cout << "═══════════════════════════════════════════════════════════════════\n";
    std::cout << "  HNtuple Improved Error Handling Examples\n";
    std::cout << "═══════════════════════════════════════════════════════════════════\n";

    try {
        example1_correct_usage();
        example2_error_add_after_freeze();
        example3_query_methods();
        example4_error_no_variables();
        example5_typo_detection();
    }
    catch (const std::exception& e) {
        std::cerr << "\nUnexpected error: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\n═══════════════════════════════════════════════════════════════════\n";
    std::cout << "  Summary of Improvements\n";
    std::cout << "═══════════════════════════════════════════════════════════════════\n";
    std::cout << "  ✓ Detailed error messages showing current structure\n";
    std::cout << "  ✓ Lists all existing variables when error occurs\n";
    std::cout << "  ✓ Clear indication of ntuple name and fill count\n";
    std::cout << "  ✓ Helpful suggestions for fixing the error\n";
    std::cout << "  ✓ New query methods: isFrozen(), getNVariables(), etc.\n";
    std::cout << "  ✓ Diagnostic method: printStructure()\n";
    std::cout << "  ✓ Validation: prevents filling without variables\n";
    std::cout << "  ✓ User-friendly boxed formatting for clarity\n";
    std::cout << "═══════════════════════════════════════════════════════════════════\n";

    return 0;
}

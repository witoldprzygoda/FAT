# research/calibration/Makefile
#
# Standalone "trigger_scan" binary — scans the analysis input list,
# counts PT3/PT2 trigbits per ROOT file, writes a small TTree.
#
# Build:  make                              (run from research/calibration/)
# Run:    ./trigger_scan <config.json>      (run from research/calibration/)
# Parallel: ./run_parallel_scan.sh <config.json> [N_parts] [N_concurrent]
#
# Outputs land per the config — typically trigger_scan_<channel>.root next
# to this Makefile (gitignored).

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

CXX      = g++
CXXFLAGS = -std=c++17 -O3 -g -Wall

EXE = trigger_scan

$(EXE): trigger_scan.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $(EXE) $< $(ROOTLIBS)

clean:
	rm -f $(EXE)

.PHONY: clean

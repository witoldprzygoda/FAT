# Standalone "research" binary — separate from the FAT main analysis.
# Reads output_*_exp.root files produced by FAT and runs slice-based
# studies (currently: per-OA-slice m_epemg histograms).
#
# Build:  make                              (run from research/ directory)
# Run:    ./meson_research <config.json>    (run from research/ directory)
#
# All paths in the JSON configs are relative to research/ (e.g. input
# "../output_epem_exp.root", output "outputs/research_epem.root").

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

CXX      = g++
CXXFLAGS = -std=c++17 -O3 -g -Wall

EXE = meson_research

$(EXE): meson_research.cc
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $(EXE) $< $(ROOTLIBS)

clean:
	rm -f $(EXE)

.PHONY: clean

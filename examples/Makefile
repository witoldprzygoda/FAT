# Makefile for PParticle Examples
# Demonstrates how to compile code using the new PParticle system

# Detect ROOT installation
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTGLIBS  := $(shell root-config --glibs)

# Compiler settings
CXX        = g++
CXXFLAGS   = -std=c++11 -Wall -Wextra -O2 $(ROOTCFLAGS) -I../src
LDFLAGS    = $(ROOTLIBS) -lPhysics

# Targets
TARGETS    = PParticle_Usage_Examples

.PHONY: all clean test

all: $(TARGETS)

# Build the usage examples
PParticle_Usage_Examples: PParticle_Usage_Examples.cc ../src/PParticle.h ../src/BoostFrame.h
	@echo "Compiling PParticle usage examples..."
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Build successful! Run with: ./PParticle_Usage_Examples"

# Run the examples
test: PParticle_Usage_Examples
	@echo ""
	@echo "======================================================"
	@echo "Running PParticle Usage Examples"
	@echo "======================================================"
	@./PParticle_Usage_Examples
	@echo ""
	@echo "======================================================"
	@echo "All examples completed successfully!"
	@echo "======================================================"

clean:
	@echo "Cleaning up..."
	rm -f $(TARGETS) *.o

# Help target
help:
	@echo "PParticle Examples Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all   - Build all examples (default)"
	@echo "  test  - Build and run examples"
	@echo "  clean - Remove built files"
	@echo "  help  - Show this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make         # Build examples"
	@echo "  make test    # Build and run"
	@echo "  make clean   # Clean up"

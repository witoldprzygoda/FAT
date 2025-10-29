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
TARGETS    = PParticle_Usage_Examples test_boost_sign_convention

.PHONY: all clean test test-boost

all: $(TARGETS)

# Build the usage examples
PParticle_Usage_Examples: PParticle_Usage_Examples.cc ../src/PParticle.h ../src/BoostFrame.h
	@echo "Compiling PParticle usage examples..."
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Build successful! Run with: ./PParticle_Usage_Examples"

# Build the boost sign convention test
test_boost_sign_convention: test_boost_sign_convention.cc ../src/PParticle.h ../src/BoostFrame.h
	@echo "Compiling boost sign convention test..."
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Build successful! Run with: ./test_boost_sign_convention"

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

# Run boost sign convention test
test-boost: test_boost_sign_convention
	@echo ""
	@echo "======================================================"
	@echo "Running Boost Sign Convention Test"
	@echo "======================================================"
	@./test_boost_sign_convention
	@echo ""

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

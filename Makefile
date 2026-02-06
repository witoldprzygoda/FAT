# Compiler and Flags
CXX    = g++ 
CFLAGS = $(shell root-config --cflags) -std=c++17 -g -Wall -fPIC -I./src
LIBS   = $(shell root-config --libs) -lProof -lEG

# Source and header files
SOURCES = main.cc src/hntuple.cc

# All headers for dependency tracking (changes trigger recompilation)
ALL_HEADERS = $(wildcard src/*.h)

# Headers for ROOT dictionary generation (subset needed by rootcling)
DICT_HEADERS = src/hntuple.h src/manager.h src/histogram_registry.h \
               src/histogram_factory.h src/histogram_builder.h \
               src/pparticle.h src/boost_frame.h \
               src/ntuple_reader.h src/cut_manager.h src/analysis_config.h \
               src/progressbar.h

OBJECTS = $(SOURCES:.cc=.o)
DICT    = MyDict.cc
DICTOBJ = $(DICT:.cc=.o)
LINKDEF = MyLinkDef.h

# Output executable
EXECUTABLE = ana

# Default action
all: $(DICT) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(DICTOBJ)
	$(CXX) -o $@ $^ $(LIBS)

# Pattern rule for object files (depends on ALL headers)
%.o: %.cc $(ALL_HEADERS) $(LINKDEF)
	$(CXX) $(CFLAGS) -c $< -o $@

# Dictionary generation (uses only DICT_HEADERS)
$(DICT): $(DICT_HEADERS) $(LINKDEF)
	rootcling -f $@ $(DICT_HEADERS) $(LINKDEF)

# Cleanup
clean:
	rm -f $(EXECUTABLE) $(OBJECTS) $(DICT) $(DICTOBJ) *.pcm

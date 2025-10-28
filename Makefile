# Compiler and Flags
CXX    = g++ 
CFLAGS = $(shell root-config --cflags) -O3 -Wall -fPIC
LIBS   = $(shell root-config --libs) -lProof -lEG

# Source and header files
SOURCES = main.cc src/manager.cc datamanager.cc src/Base_ID.cc PPip_ID.cc src/hntuple.cc
HEADERS = data.h src/manager.h src/datamanager.h src/Base_ID.h PPip_ID.h src/hntuple.h
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

# Pattern rule for object files
%.o: %.cc $(HEADERS) $(LINKDEF)
	$(CXX) $(CFLAGS) -c $< -o $@

# Dictionary generation
$(DICT): $(HEADERS) $(LINKDEF)
	rootcling -f $@ $(HEADERS) $(LINKDEF)

# Cleanup
clean:
	rm -f $(EXECUTABLE) $(OBJECTS) $(DICT) $(DICTOBJ) *.pcm


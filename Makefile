BAMTOOLS_DIR := bamtools
ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
CXX := g++
CXXFLAGS := -Wno-deprecated -Wall -g -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(BAMTOOLS_DIR)/include/ -L$(BAMTOOLS_DIR)/lib/
OBJS := $(FASTAHACK)

all: Microassembler


libbamtools.a:
	cd $(BAMTOOLS_DIR) && mkdir -p build && cd build && cmake .. && $(MAKE)
	cp $(BAMTOOLS_DIR)/lib/libbamtools.a ./

$(FASTAHACK):
	cd fastahack && $(MAKE)

$(VCFLIB): libbamtools.a
	cd vcflib && $(MAKE)


Microassembler: Microassembler.cc Microassembler.hh align.cc util.cc Mer.hh Ref.hh ReadInfo.hh ReadStart.hh Transcript.hh Edge.cc Edge.hh ContigLink.hh Node.cc Node.hh Path.cc Path.hh ContigLink.cc Graph.cc Graph.hh libbamtools.a $(FASTAHACK)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJS) Microassembler.cc Edge.cc Node.cc Graph.cc Path.cc ContigLink.cc align.cc util.cc -o Microassembler -lbamtools -lz

clean:
	rm -rf Microassembler libbamtools.a $(OBJS)
	cd $(BAMTOOLS_DIR)/build && make clean

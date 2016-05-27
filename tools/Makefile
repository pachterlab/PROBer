CXX = g++
CXXFLAGS = -std=c++11 -Wall -I. -O3

PROGRAMS = parseGTF annotateSite showGeneNames showSite

.PHONY : all clean

all : $(PROGRAMS)

$(PROGRAMS) :
	$(CXX) $(CXXFLAGS) -o $@ $<

parseGTF : parseGTF.cpp GTFItem.h genome_structs.h
annotateSite : annotateSite.cpp genome_structs.h
showGeneNames : showGeneNames.cpp genome_structs.h
showSite : showSite.cpp genome_structs.h

clean :
	rm -f $(PROGRAMS) *~



CXX = g++
CXXFLAGS= -Wall -std=c++11 -O3 -msse4.2 #-pg -g #-Wall #-O3
LINKPATH=
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = main.o 

#asan=1
ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address -g
	LDFLAGS+=-fsanitize=address -ldl -g
endif

all: centrifuger centrifuger-build centrifuger-inspect centrifuger-quant

centrifuger-build: CentrifugerBuild.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

centrifuger: CentrifugerClass.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

centrifuger-inspect: CentrifugerInspect.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

centrifuger-quant: CentrifugerQuant.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)


CentrifugerBuild.o: CentrifugerBuild.cpp Builder.hpp ReadFiles.hpp Taxonomy.hpp defs.h compactds/*.hpp 
CentrifugerClass.o: CentrifugerClass.cpp Classifier.hpp ReadFiles.hpp Taxonomy.hpp defs.h ResultWriter.hpp ReadPairMerger.hpp ReadFormatter.hpp BarcodeCorrector.hpp BarcodeTranslator.hpp compactds/*.hpp 
CentrifugerInspect.o: CentrifugerInspect.cpp Taxonomy.hpp defs.h compactds/*.hpp 
CentrifugerQuant.o: CentrifugerQuant.cpp Quantifier.hpp Taxonomy.hpp defs.h compactds/*.hpp

clean:
	rm -f *.o centrifuger-build centrifuger centrifuger-inspect centrifuger-quant

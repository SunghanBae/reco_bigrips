OBJ = bigrips reco
GO4SYS = /media/bae/data1/research/RIKEN/201605~06_EURICA/go4-5.2.0

all: $(OBJ)

include $(GO4SYS)/Makefile.config

CXX	= g++

#LIBS            = -L/usr/local/lib -lXMLParser
ROOTCFLAGS	= $(shell root-config --cflags)
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)
TARTLIBS        = -L$(TARTSYS)/lib -lanacore -lanaroot -lanabrips -lanawinds -lanaloop -lanacore -L/usr/local/lib -lXMLParser
INCLUDES	= -I$(TARTSYS)/include -I$(ROOTSYS)/include


ifdef GO4_WIN32
   GO4SYS = ../go4
endif

bigrips: MakeFullBigRIPSTree.o
	$(CXX) $< -o $@ $(INCLUDES) $(ROOTLIBS) $(ROOTFLAGS) $(ROOTGLIBS) $(TARTLIBS) -std=c11 -std=c++11

reco: reco_main.cpp reco.h
	$(CXX) $< -o $@ $(INCLUDES) $(ROOTLIBS) $(ROOTFLAGS) $(ROOTGLIBS) $(TARTLIBS) -std=c11 -std=c++11

clean: 
	rm -f $(TARGET) *.o

OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
#CXX = $(shell $(ROOTSYS)/bin/root-config --cxx)
CXX = g++
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags) \
	-O3 -W -Wall -Wextra -Wno-long-long \
	-fno-common \
	-I$(INCDIR) -I$(HEEDDIR) -I$(ROOTINC)

# Debug flags
CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --glibs) -lGeom -lgfortran -lm

make: trific gasTable

trific: trific.cc
	$(CXX) $(CFLAGS) -c trific.cc
	$(CXX) $(CFLAGS) -o trific trific.o $(LDFLAGS)

gasTable: gasTable.cc
	$(CXX) $(CFLAGS) -c gasTable.cc
	$(CXX) $(CFLAGS) -o gasTable gasTable.o $(LDFLAGS)

clean:
	rm -f *.o *.*~ *~ trific gasTable

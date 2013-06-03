



CXXFLAGS = -g -Wall -DDEBUG_HKM
LCXX_FLAGS = -g -Wall -DDEBUG_HKM

CXX=gcc

.cpp.o:
	$(CXX) $(CXXFLAGS) $(CCFLAGS) -c $<

all: metal_Li_LiCl_electrons


metal_Li_LiCl_electrons: metal_Li_LiCl_electrons.o
	$(CXX) $(LCXX_FLAGS) -o metal_Li_LiCl_electrons   metal_Li_LiCl_electrons.o -lm -lstdc++ 

clean:
	rm -f *.o metal_Li_LiCl_electrons out.txt polarizationData.csv diff*

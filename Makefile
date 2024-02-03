ROOTDIR=$(shell pwd)
OUTPUTDIR=$(ROOTDIR)/paraview-output/

.PHONY: all cleanall clean clean_paraview
all: step-1 step-2 step-3 step-4
step-%: step-%


# Target to be used with the GNU Compiler Collection.
# On Hamilton, the default version of this compiler is 8.5.0.
# If can load a more recent version with
#     $ module load gcc/12.2
step-% step-%.o NBodySimulation.o: CXX?=g++
step-% step-%.o NBodySimulation.o: CXXFLAGS?=-fopenmp -O3 -march=native -std=c++0x -fno-math-errno
NBodySimulation.o: NBodySimulation.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
step-%.o: step-%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
NBodySimulationVectorised.o: NBodySimulationVectorised.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
step-%: NBodySimulationVectorised.o NBodySimulation.o step-%.o
	$(CXX) $(CXXFLAGS) -o $@ $^

.silent: cleanall clean clean_paraview
cleanall: clean clean_paraview

clean:
	rm -rf $(ROOTDIR)/step-[1-4] $(ROOTDIR)/*.o

clean_paraview:
	if test -d "$(OUTPUTDIR)"; then \
		find $(OUTPUTDIR) -iname "result-*.vtp" -delete; \
		rm -rf $(OUTPUTDIR)/result.pvd; \
	fi

test:
	./validate.sh
	python3 validate.py

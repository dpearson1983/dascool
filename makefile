CXX = g++
CXXLIBS = -lgsl -lgslcblas
CXXOPTS = -march=native -mtune=native -O3

localinstall:
	$(MAKE) libdascool.so
	$(shell cp libdascool.so $(HOME)/lib)
	$(shell cp dascool.hpp $(HOME)/include)
	
libdascool.so: dascool.hpp dascool.cpp
	$(CXX) $(CXXLIBS) $(CXXOPTS) dascool.cpp -fPIC -shared -o libdascool.so 
	
clean:
	$(shell rm $(HOME)/lib/libdascool.so)
	$(shell rm $(HOME)/include/dascool.hpp)
	$(shell rm libdascool.so)

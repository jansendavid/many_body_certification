
INC+=-I$(MOSEK)/src/fusion_cxx/
INC+=-I$(MOSEK)/h/
INC+=-I./include/
INC+=-I./include/sos
INC+=-I$(EIGEN)
#########################

CCFLAGS=  -Wno-unused-variable -std=c++17 -O2 -std=gnu++1z
CCGFLAGS+=-I. $(DEBUGFLAGS) 

LIBFLAGS=-lfusion64 -lmosek64
LIBSPATH_1=-L$(MOSEK)/src/fusion_cxx/
LIBSPATH_2=-L$(MOSEK)/bin/
#-lboost_program_options -lboost_filesystem 




DB=-g
CXX=g++
ND=-DNDEBUG

#-fopenmp -DOMPPAR

#
#BINS=
#
#OBJECTS=holGS2d
#all: $(OBJECTS)
##############################
2dxxz: src/2dxxz.cpp
	g++ -std=c++11 src/2dxxz.cpp -o bin/$@ -std=c++20 -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)
1dxxz: src/1dxxz.cpp
	g++ -std=c++11 src/1dxxz.cpp -o bin/$@ -std=c++20 -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)
##tests##
#######################################
test_spins: tests/test_spins.cpp
	g++ -std=c++11 tests/test_spins.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_ham_spins: tests/test_ham_spins.cpp
	g++ -std=c++11 tests/test_ham_spins.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_ham_spins_eff: tests/test_ham_spins_eff.cpp
	g++ -std=c++11 tests/test_ham_spins_eff.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_ham_spins_eff_dual: tests/test_ham_spins_eff_dual.cpp
	g++ -std=c++11 tests/test_ham_spins_eff_dual.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)


test_momentum_spins: tests/test_momentum_spins.cpp
	g++ -std=c++11 tests/test_momentum_spins.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_momentum_spins_eff: tests/test_momentum_spins_eff.cpp
	g++ -std=c++11 tests/test_momentum_spins_eff.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_momentum_spins_eff_xy: tests/test_momentum_spins_eff_xy.cpp
	g++ -std=c++11 tests/test_momentum_spins_eff_xy.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_momentum_spins_sos_dual: tests/test_momentum_spins_sos_dual.cpp
	g++ -std=c++11 tests/test_momentum_spins_sos_dual.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_momentum_spins_dual_xy: tests/test_momentum_spins_dual_xy.cpp
	g++ -std=c++11 tests/test_momentum_spins_dual_xy.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)


test_momentum_spins_dual: tests/test_momentum_spins_dual.cpp
	g++ -std=c++11 tests/test_momentum_spins_dual.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)


test: src/test.cpp
	g++ -std=c++11 src/test.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)
test_hamiltonian: src/test_hamiltonian.cpp
	g++ -std=c++11 src/test_hamiltonian.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_hamiltonian_spin: src/test_hamiltonian_spin.cpp
	g++ -std=c++11 src/test_hamiltonian_spin.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

test_functions: tests/test_functions.cpp
	g++ -std=c++11 tests/test_functions.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)


momentum_test: src/momentum_test.cpp
	g++ -std=c++11 src/momentum_test.cpp -o bin/$@ -std=c++20 -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)
#sos
test_ham_spins_dual: sos/test_ham_spins_dual.cpp
	g++ -std=c++11 sos/test_ham_spins_dual.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)
test_ham_spins_sos: sos/test_ham_spins_sos.cpp
	g++ -std=c++11 sos/test_ham_spins_sos.cpp -o bin/$@ -I$(MOSEK)/h/ -L$(MOSEK)/bin/ -Wl,-rpath-link,$(MOSEK)/bin/ -Wl,-rpath=$(MOSEK)/bin/ -lmosek64 -lfusion64 $(CCFLAGS)  $(INC)

clean:
	rm bin/*

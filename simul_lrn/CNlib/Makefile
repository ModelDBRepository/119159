#--------------------------------------------------------------------------
# Author: Thomas Nowotny
#
# Institute: Institute for Nonlinear Dynamics
#            University of California San Diego
#            La Jolla, CA 92093-0402
#
# email to:  tnowotny@ucsd.edu
#
# initial version: 2005-08-17
#
#--------------------------------------------------------------------------

C++ = g++
C = gcc
RM= rm -f
MAKE= make
FLAGS= -Wall -g -O3 -funroll-all-loops -ffast-math -march=athlon
HOME= /autohome/quarz/nowotny
INCLUDE= -I$(HOME)/include/datalib/ -I$(HOME)/include/numlib/ -I$(HOME)/include/ISAAC_C++

all: CN_rk65n.o testCN

CN_rk65n.o: CN_rk65n.h CN_rk65n.cc
	$(C++) $(FLAGS) $(INCLUDE) -g -c CN_rk65n.cc

testCN.o: testCN.cc
	$(C++) $(FLAGS) $(INCLUDE) -g -c testCN.cc

testCN: testCN.o CN_rk65n.o
	$(C++) $(FLAGS) $(INCLUDE) -o testCN testCN.o CN_rk65n.o

clean:
	$(RM) *.o testCN

obj = 		VectorOperation.o main.o

run:		$(obj)
	g++ -g -Wall -O0 $^ -o $@ -L/home/wang/intel/mkl/8.0.1/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread	


VectorOperation.o main.o: VectorOperation.h
CXXFLAGS += 	-I/home/wang/intel/mkl/8.0.1/include   -g -Wall
vpath %.h 	/home/wang/intel/mkl/8.0.1/include    
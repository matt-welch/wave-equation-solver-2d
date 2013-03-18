#/*******************************************************************************
# * FILENAME:    Makefile
# * DESCRIPTION: build options for 2D wave equation solver (wavesolv.c) 
# * AUTHOR:      James Matthew Welch [JMW]
# * SCHOOL:      Arizona State University
# * CLASS:       CSE598: High Performance Computing
# * INSTRUCTOR:  Dr. Gil Speyer
# * SECTION:     20520
# * TERM:        Spring 2013
# *******************************************************************************/
#
all: mpi serial

mpi: wavesolv_mpi 

serial: wavesolv

wavesolv: wavesolv_serial.c 
	gcc wavesolv_serial.c -o wavesolv $(CONFIG)

wavesolv_mpi: wavesolv_mpi.c
	mpicc wavesolv_mpi.c -o wavesolv_mpi $(CONFIG)
	ctags --c-kinds=+defglmpstux *.c

parallel: wavesolv_serial.c wavesolv_mpi.c
	mpicc wavesolv_mpi.c -o wavesolv_mpi_omp $(CONFIG) $(PARCONFIG)
	gcc wavesolv_serial.c -o wavesolv_omp $(CONFIG) $(PARCONFIG) 

animate: wavesolv_serial.c wavesolv_mpi.c
	mpicc wavesolv_mpi.c -o wavesolv_mpi $(CONFIG) -DCREATEANIMATION
	gcc wavesolv_serial.c -o wavesolv $(CONFIG)  -DCREATEANIMATION

wavesolv_mpi.o: wavesolv_mpi.c
	mpicc.mpich2 wavesolv_mpi.c -o wavesolv_mpi $(CONFIG) 

# compile mpi program with logging to use in JumpShot
log: wavesolv_mpi.c
	mpicc  wavesolv_mpi.c -o wavesolv_mpi_log -lm -mpe=mpilog

profile: wavesolv_mpi.c
	mpicc wavesolv_mpi.c -o wavesolv_mpi_omp_prof -pg $(CONFIG) $(PARCONFIG) 
	mpicc wavesolv_mpi.c -o wavesolv_mpi_prof -pg $(CONFIG) 

clean:
	rm -f wavesolv wavesolv_mpi wavesolv_omp wavesolv_mpi_omp *.o core core.*

tidy: clean
	rm -f *.*~ *~ *.asv output*

DEBUG_FLAGS = -g3 -ggdb -O0 -DDEBUG -Wall -pedantic -lm
CONFIG		= -Wall -pedantic -lm 
PARCONFIG	= -fopenmp -DUSEOPENMP
EASY_FLAGS	= -lm

debug: CONFIG=$(DEBUG_FLAGS)
debug: wavesolv wavesolv_mpi

easy: CONFIG=$(EASY_FLAGS)
easy: wavesolv_serial wavesolv_mpi

quiet: CONFIG=$(EASY_FLAGS) -DQUIET
quiet: all

test: wavesolv wavesolv_mpi
	./wavesolv 100 100

runmpi: debug 
	mpiexec -n 4 -comm=pmi ./wavesolv_mpi 

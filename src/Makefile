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
all: wavesolv_mpi wavesolv_serial

wavesolv: wavesolv.c 
	gcc wavesolv.c -o wavesolv $(CONFIG)

wavesolv_mpi: wavesolv_mpi.c
	mpicc wavesolv_mpi.c -o wavesolv_mpi $(CONFIG) $(PARCONFIG)
	ctags --c-kinds=+defglmpstux *.c

# compile mpi program with logging to use in JumpShot
log: wavesolv_mpi.c
	mpicc.openmpi wavesolv_mpi.c -o wavesolv_mpi $(CONFIG) $(PARCONFIG) -mpilog 

wavesolv_serial: wavesolv

clean:
	rm -f wavesolv wavesolv_mpi *.o core core.*

tidy: clean
	rm -f *.*~ *~

DEBUG_FLAGS = -g3 -ggdb -O0 -DDEBUG -Wall -pedantic -lm
CONFIG		= -Wall -pedantic -lm
PARCONFIG	= -DPARALLEL 
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
	mpiexec -n 4 ./wavesolv_mpi 

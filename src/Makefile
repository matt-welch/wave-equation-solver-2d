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
all: wavesolv wavesolv_serial

wavesolv: 
	gcc wavesolv.c -o wavesolv $(CONFIG)

wavesolv_serial:

clean:
	rm -f wavesolv *.o core core.*

tidy: clean
	rm -f *.*~ *~

DEBUG_FLAGS = -g3 -ggdb -O0 -Wall -pedantic -DDEBUG -fopenmp -lm
CONFIG		= -Wall -pedantic -lm
PARCONFIG	= $(CONFIG) -DPARALLEL
EASY_FLAGS	= -lm

debug: CONFIG=$(DEBUG_FLAGS)
debug: wavesolv

easy: CONFIG=$(EASY_FLAGS)
easy: wavesolv

option1: CONFIG= -O0 -lm -DQUIET
option1: wavesolv

option2: CONFIG= -O3 -lm -DQUIET
option2: wavesolv


quiet: CONFIG=$(EASY_FLAGS) -DQUIET
quiet: all

test: quiet 
	time ./wavesolv 1000
	time ./wavesolv_parallel 1000
	

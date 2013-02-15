This project will involve writing a parallel 2D wave equation solver.
The wave equation can be succinctly represented as: 

    u(tt) = c^2 * ( u(xx) + u(yy) )

where the constant c has units of velocity (m/s)

This program will solve the 2D wave equation where u is a function of t, x, & y: 
u = f(t, x, y)

At time t(-1) and t(0), the magnitude of the wave across the entire domain 
is 0 as this is the base condition. At some time later, a pulse will be
introduced at one edge of the 2D domain.  This wave will propagate across the
domain until it reflects across the alternate edge fo the domain.  This
simulation will continue until the maximum magnitude of the wave is less than
20% of the initial magnitude.  At this point, another pulse will be introduced
at another edge, 90 degrees clockwise form the last pulse.  

Requirements::
Milestone 1:
    1. Write an MPI program that implements the wave equation on a 480X480 grid, propagatin         a single pulse from one side for 100 time steps. Output the data for the last time step to a text file
    with 3 columns (x, y, and u).
    2. Provide preliminary performance numbers for this program on 1, 4 and 16 processors.
    3. Give one paragraph describing your strategy for implementing the
       additional parts of the program (described below) using MPI.

Milestone 2:
    1. Now have your program find the maximum amplitude (u) across the entire domain. When this
    maximum is found to be 20% of the height of the original pulse, launch another pulse from the
    boundary edge 90 clockwise. Run for a 2000 time steps and provide an output file for the last
    time step.
    2. Insert OpenMP pragmas into the program (Make this another version so you can compare the
    two programs.) You can compile by using mpicc with the –openmp flag. Before you run the
    program, you must tell the compute nodes that you will be running MPI jobs between nodes
    with OpenMP running on the nodes. This is done by using the –pernode flag with mpiexec:
    mpiexec –pernode ./yourmpiprogram
    When –pernode is run, the compute nodes will automatically send all threads to one CPU per
    node. To override this you must type (before running mpiexec):
    export KMP_AFFINITY=norespect,scatter
    Now you will be running a hybrid code.
    3. Prepare a final report comparing the two implementations for two dataset sizes (1024x1024 and
    2048x2048) up to 256 processors. Present speedup plots (1, 4, 16, 64, 256), profiling
    information (describing where time is spent in the program) and jumpshot pictures. Discuss
    your communication patterns (data manipulation between neighbors, use of collectives),
    granularity, concurrency vs. communication overhead. Justify decisions you made in your code
    design. You do not need to include the output routine or run for all the time steps to get your
    performance data.


 Project 2 was posted yesterday.  I have made some revisions and posted these today.  Specifically:
1) Due date for Milestone 2 is now 3/17
2) You can simply use two zero matrices for timesteps l=0 and l=1
3) You can simply hard code delta x/y and delta t to fulfill the CFL condition
4) An output text file must be provided with your codes.

ideas:
use shorts for t, x, y
use integer for u like 32 bit A/D
output





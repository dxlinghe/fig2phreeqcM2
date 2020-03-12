EXECS = ex0 
GCC=g++
MPICC=mpicc

all:${EXECS}

ex0:add.c
	${MPICC} -g -O1 -std=gnu99  -I/home/sliu/softw/include -L/home/sliu/softw/lib -o ex0 add.c -liphreeqc #-limf


clean:
	rm ${EXECS} EL_* 

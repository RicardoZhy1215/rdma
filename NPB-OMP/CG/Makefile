SHELL=/bin/sh
BENCHMARK=cg
BENCHMARKU=CG

include ../config/make.def

OBJS = cg.o ${COMMON}/c_print_results.o  \
       ${COMMON}/c_${RAND}.o ${COMMON}/c_timers.o ${COMMON}/c_wtime.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB} -libverbs

cg.o:             cg.cpp npbparams.hpp
	${CCOMPILE} cg.cpp

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.hpp core

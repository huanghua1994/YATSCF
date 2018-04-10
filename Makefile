CC  = icc
CXX = icpc
EXE = YATSCF.exe

BLAS_LIBS      = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread  
LIBCINT_INCDIR = /home/huangh/libcint
LIBCINT_LIBDIR = /home/huangh/libcint
LIBCINT_LIB    = ${LIBCINT_LIBDIR}/libcint.a 
ERI_LIB        = /home/huangh/gtfock-simint/build-avx512/install/lib64/libsimint.a

INCS           = -I./ -I${LIBCINT_INCDIR} 
LIBS           = ${BLAS_LIBS} ${LIBCINT_LIB} ${ERI_LIB}

CFLAGS         = -Wall -g -O3 -qopenmp -std=gnu99 -xHost
LDFLAGS        = -L${LIBCINT_LIB} -lpthread -qopenmp

OBJS = utils.o TinySCF.o main.o 

$(EXE): Makefile $(OBJS) ${LIBCINT_LIB} ${ERI_LIB}
	$(CC) ${CFLAGS} ${LDFLAGS} $(OBJS) -o $(EXE) ${LIBS}

utils.o: Makefile utils.c utils.h
	$(CC) ${CFLAGS} ${INCS} -c utils.c   -o $@ 
	
TinySCF.o: Makefile TinySCF.c TinySCF.h utils.h
	$(CC) ${CFLAGS} ${INCS} ${BLAS_LIBS} -c TinySCF.c -o $@ 
	
main.o: Makefile main.c TinySCF.h
	$(CC) ${CFLAGS} ${INCS} -c main.c    -o $@ 

clean:
	rm -f *.o $(EXE)

CC  = icc
CXX = icpc
EXE = TinySCF.exe

BLAS_LIBS      = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread  
LIBCINT_INCDIR = /home/huangh/libcint
LIBCINT_LIBDIR = /home/huangh/libcint
LIBCINT_LIB    = ${LIBCINT_LIBDIR}/libcint.a 
ERI_LIB        = /home/huangh/gtfock-simint/build-avx512/install/lib64/libsimint.a

INCS           = -I./ -I${LIBCINT_INCDIR} 
LIBS           = ${BLAS_LIBS} ${LIBCINT_LIB} ${ERI_LIB}

CFLAGS         = -Wall -g -O3 -qopenmp -std=gnu99 -xHost
LDFLAGS        = -L${LIBCINT_LIB} -lpthread -qopenmp

OBJS = utils.o build_density.o build_Fock.o DIIS.o TinySCF.o main.o 

$(EXE): Makefile $(OBJS) ${LIBCINT_LIB} ${ERI_LIB}
	$(CC) ${CFLAGS} ${LDFLAGS} $(OBJS) -o $(EXE) ${LIBS}

utils.o: Makefile utils.c utils.h
	$(CC) ${CFLAGS} ${INCS} -c utils.c -o $@ 

build_density.o: Makefile build_density.c build_density.h TinySCF.h
	$(CC) ${CFLAGS} ${INCS} ${BLAS_LIBS} -c build_density.c -o $@ 

build_Fock.o: Makefile build_Fock.c build_Fock.h TinySCF.h
	$(CC) ${CFLAGS} ${INCS} -c build_Fock.c -o $@ 

DIIS.o: Makefile DIIS.c DIIS.h TinySCF.h
	$(CC) ${CFLAGS} ${INCS} -c DIIS.c -o $@ 

TinySCF.o: Makefile TinySCF.c TinySCF.h utils.h
	$(CC) ${CFLAGS} ${INCS} ${BLAS_LIBS} -c TinySCF.c -o $@ 
	
main.o: Makefile main.c TinySCF.h
	$(CC) ${CFLAGS} ${INCS} -c main.c    -o $@ 

clean:
	rm -f *.o $(EXE)

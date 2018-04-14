alignlen = 64

CC  = icc

AR  = xiar rcs

SIMINT_DIR = /home/huangh/gtfock-simint/build-avx512/install

INC=-I. -I${SIMINT_DIR}/include

OPTFLAGS = -m64 

SRC = cint_basisset.c cint_simint.c
CFLAGS  = -O3 -Wall -qopenmp -g -std=gnu99 -qno-offload
CFLAGS += ${OPTFLAGS}
CFLAGS += -D__ALIGNLEN__=${alignlen}

LIBCINT = libcint.a
OBJS := $(addsuffix .o, $(basename $(SRC)))

all: ${LIBCINT} 

${LIBCINT}: ${OBJS}
	${AR} $@ $^

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@ 

clean:
	rm -f *.o *.s *.d *~ *.a


L=-pthread
HG_DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE} -DUSE_HIC

##########
#
# EDIT the lines below so that they point to your kent source
# and gsl header files and libraries 
#
HG_INC += -I/home/cl454/src/kent/src/hg/inc -I/home/cl454/src/kent/src/inc -I/home/cl454/src/kent/src/htslib
L += /home/cl454/src/kent/src/lib/${MACHTYPE}/jkweb.a /home/cl454/src/kent/src/htslib/libhts.a

HG_INC += -I/home/cl454/src/gsl/gsl-2.6_install/include
L += /home/cl454/src/gsl/gsl-2.6_install/lib/libgsl.a /home/cl454/src/gsl/gsl-2.6_install/lib/libgslcblas.a

HG_INC += -I/home/cl454/src/R/R-3.6.3/src/include
L += /home/cl454/src/R/R-3.6.3/src/nmath/standalone/libRmath.a

L += /usr/lib/x86_64-linux-gnu/libssl.a

#
# END of basic editing
#
##########

CC=gcc
ifeq (${COPT},)
    COPT=-O -g
endif
ifeq (${CFLAGS},)
    CFLAGS=
endif

SYS = $(shell uname -s)

ifeq (${HG_WARN},)
  ifeq (${SYS},Darwin)
      HG_WARN = -Wall -Wno-unused-variable
      HG_WARN_UNINIT=
  else
    ifeq (${SYS},SunOS)
      HG_WARN = -Wall -Wformat -Wimplicit -Wreturn-type
      HG_WARN_UNINIT=-Wuninitialized
    else
      HG_WARN = -Wall -Werror -Wformat -Wimplicit -Wreturn-type
      HG_WARN_UNINIT=-Wuninitialized
    endif
  endif
  # -Wuninitialized generates a warning without optimization
  ifeq ($(findstring -O,${COPT}),-O)
     HG_WARN += ${HG_WARN_UNINIT}
  endif
endif

CFLAGS += ${HG_WARN}

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${HG_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<

L += -lcrypto -lm -lz

A = wgHmm
H = hmm.h
O = hmm.o wgHmm.o

wgHmm: ${O} ${MYLIBS}
	${CC} ${COPT} -o ${A} $O ${MYLIBS} $L

wgHmm.o: wgHmm.c hmm.h

hmm.o: hmm.c

clean:
	rm -f ${A} ${O}


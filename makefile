
L=-pthread
HG_DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}

##########
#
# EDIT the lines below so that they point to your kent source
# header files and libraries 
#
HG_INC += -I/home/lowec/kent/src/hg/inc -I/home/lowec/kent/src/inc
L += /home/lowec/kent/src/lib/${MACHTYPE}/jkweb.a
#
# END of basic editing
#
##########

##########
#
# If you compiled your kent source with sam/bam and tabix support
# then you will need to edit the below lines as well.
# If none of that sounded familiar to you, then you probably don't
# need to edit this.
#
ifeq (${USE_SAMTABIX},1)
    KNETFILE_HOOKS=1
    USE_BAM=1
    USE_TABIX=1
    ifeq (${SAMTABIXINC},)
        SAMTABIXINC = ${SAMTABIXDIR}
    endif
    ifeq (${SAMTABIXLIB},)
        SAMTABIXLIB = ${SAMTABIXDIR}/libsamtabix.a
    endif
    HG_INC += -I${SAMTABIXINC}
    L+=${SAMTABIXLIB} -lz
    HG_DEFS+=-DUSE_SAMTABIX -DUSE_BAM -DUSE_TABIX -DKNETFILE_HOOKS
endif
#
# End of sam/bam and tabix editing
#
##########

##########
#
# If you compiled your kent source with SSL support
# then you will need to use these lines
#
ifeq (${USE_SSL},1)
	L+=-lssl -lcrypto
	HG_DEFS+=-DUSE_SSL
endif
#
###########


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

L += -lm -lz

A = wgHmm
H = hmm.h
O = hmm.o wgHmm.o

wgHmm: ${O} ${MYLIBS}
	${CC} ${COPT} -o ${A} $O ${MYLIBS} $L

wgHmm.o: wgHmm.c hmm.h

hmm.o: hmm.c

clean:
	rm -f ${A} ${O}


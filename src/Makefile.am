ZLIB = /usr
BZIP = /homes/ssinghal/Builds
LZO = /homes/ssinghal/Builds
ACOMP_INCLUDES = -I${ZLIB}/include -I${BZIP}/include -I${LZO}/include
ACOMP_LIBS = -lm ${ZLIB}/lib64 -lz ${BZIP}/lib -lbz ${LZO}/lib

AM_CFLAGS = -g -Wall -fPIC ${ACOMP_INCLUDES}
AM_LDFLAGS = ${ACOMP_LIBS}

#bin_PROGRAMS = acomp 

#acomp_SOURCES = acomp.c
#acomp_LDADD = ${LIBS}

lib_LTLIBRARIES = libacomp.la
libacomp_la_SOURCES = acomp.c 
libacomp_la_LDFLAGS = -version-info 1:0:0


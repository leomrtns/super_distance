AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib -I$(srcdir)/../@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@  
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@
LDADD = ../lib/libsuper_distance.la ../biomcmc-lib/lib/libbiomcmc_static.la \
				$(GTKDEPS_LIBS) @ZLIB_LIBS@ @LZMA_LIBS@ $(AM_LDFLAGS)

bin_PROGRAMS = super_distance 
super_distance_SOURCES = src_super_distance.c src_super_distance.h
super_distance_LDADD = $(LDADD) 

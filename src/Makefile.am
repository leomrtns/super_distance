AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib -I$(srcdir)/../@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../lib/libsuper_distance.la ../biomcmc-lib/lib/libbiomcmc_static.la  $(AM_LDFLAGS)

bin_PROGRAMS = super_distance 

super_distance_SOURCES = src_super_distance.c src_super_distance.h

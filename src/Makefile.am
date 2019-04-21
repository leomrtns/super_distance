AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib -I$(srcdir)/../biomcmc-lib/lib  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../lib/libsuper_sptree.la ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)

bin_PROGRAMS = control super_distance 

super_distance_SOURCES = super_distance.c super_distance.h
control_SOURCES = control_file.c control_file.h

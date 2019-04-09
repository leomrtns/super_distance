AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ -I$(srcdir)/../biomcmc-lib/lib
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS) ../biomcmc-lib/lib/libbiomcmc.la

common_headers = topology_space.h topology_build.h \
                 super_sptree.h

common_src     = topology_space.c topology_build.c \
                 super_sptree.c

lib_LTLIBRARIES = libsuper_sptree.la   ## change to noinst_LT... if should be linked statically (not installed globally)
otherincludedir = $(includedir)/super_sptree  ## make sure all headers go to directory 
otherinclude_HEADERS = config.h $(common_headers)
libsuper_sptree_la_SOURCES = $(common_src) 


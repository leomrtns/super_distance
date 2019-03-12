AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ -I$(srcdir)/biomcmc-lib

AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS) biomcmc-lib/libbiomcmc.la


common_headers = topology_common.h topology_space.h topology_mrca.h topology_build.h \
                 topology_splitset.h super_sptree.h

common_src     = topology_common.c topology_space.c topology_mrca.c topology_build.c \
                 topology_splitset.c super_sptree.c

lib_LTLIBRARIES = libsuper_sptree.la
otherincludedir = $(includedir)/super_sptree
otherinclude_HEADERS = config.h $(common_headers)
libsuper_sptree_la_SOURCES = $(common_src) $(common_headers)


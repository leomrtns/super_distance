dnl Autoconf init, autoheader output location, config and prefix directories
AC_INIT([super_distance], [0.1], [leomrtns@gmail.com])
AC_CONFIG_AUX_DIR([.])

dnl AC_CONFIG_SRCDIR([src/main_distance.c])
AC_CONFIG_HEADER([lib/config.h])
AC_PREFIX_DEFAULT(${HOME})
dnl libtoolize recomendation
AC_CONFIG_MACRO_DIR([m4])
dnl Override default O2
CFLAGS=${CFLAGS-""}
CXXFLAGS=${CXXFLAGS-""}

dnl automake initialization (completely unrelated to compiler arguments)
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_USE_SYSTEM_EXTENSIONS
AM_MAINTAINER_MODE
dnl must be called after AC_USE_SYSTEM_EXTENSIONS (new restriction in autoconf2.68) 
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl Basic compiler and related tools
AC_LANG_C
AC_PROG_CC
AC_PROG_INSTALL

dnl directs the libc header files to provide the standard GNU system interfaces including all GNU extensions
AC_GNU_SOURCE

dnl Libtool (for library): in lib/Makefile.am, we include an "abstract" libfoo.la
AC_PROG_LIBTOOL

dnl openMP: sets $OPENMP_CFLAGS which should be passed to CFLAGS, CPPFLAGS; creates preprocessor macro _OPENMP 
dnl         (checked with "ifdef _OPENMP"); user can disable it through "--disable-openmp"
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_CPPFLAGS)
dnl unit tests in for "make check"
PKG_CHECK_MODULES([CHECK], [check >= 0.9.10])

dnl library checks (math etc.) are done by biomcmc-lib

AC_MSG_RESULT([                ===    configuration options specific to super_distance])

AC_MSG_CHECKING([whether to build (slower) debug code])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable debugging with gdb and friends (default=no)])],
    [debugit="$enableval"], [debugit=no])
AC_MSG_RESULT([$debugit])
if test x"$debugit" = x"yes"; then
    AC_DEFINE([BIOMCMC_DEBUG],[],[Debug Mode, with assert()-like error checking])
    AM_CFLAGS="${AM_CFLAGS} -g -W -Wall -Werror -Wno-uninitialized -O0"
else
    AM_CFLAGS="${AM_CFLAGS} -funroll-loops -fomit-frame-pointer -finline-functions -O4"
fi

AC_MSG_CHECKING([whether you want static binaries (not the library, but the executables)])
AC_ARG_ENABLE(static-binary,
    [  --enable-static-binary  static binaries, that run on same arch without the libraries [[default=no]]],
    [ statbin_use="yes" ], [ statbin_use=""    ])
if test -n "${statbin_use}"; then
    AC_MSG_RESULT([yes])
    AM_LDFLAGS="-static ${AM_LDFLAGS}";
else
    AC_MSG_RESULT([no])
fi

AC_MSG_CHECKING([whether to use GNU89 standard (temporary option, debug only)])
AC_ARG_ENABLE(gnu89,
    [AS_HELP_STRING([--enable-gnu89],[gnu89 standard for the GCC library (default=no)])],
    [ gnu89_use="yes" ], [ gnu89_use=""    ])
if test -n "${gnu89_use}"; then
    AC_MSG_RESULT([yes])
    AM_CFLAGS="${AM_CFLAGS} -std=gnu89"
else
    AC_MSG_RESULT([no])
    AM_CFLAGS="${AM_CFLAGS} -std=gnu11"
fi

AC_MSG_RESULT([                ===    end of specific configuration options])

dnl propagate changed vars among final makefiles
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([MPI_CXXLIBS])
AC_SUBST([MPI_CXXFLAGS])

dnl generate makefiles (last recipes, after defining CFLAGS etc.)
AC_CONFIG_FILES([ Makefile lib/Makefile src/Makefile tests/Makefile])

# Call biomcmc-lib ./configure script recursively.
AC_CONFIG_SUBDIRS([biomcmc-lib])
AC_OUTPUT

echo \
"---- 
Configuration parameters for super_distance:

 Source code location:   ${srcdir}
 Compiler     (library): ${CC}
 Compiler flags: ${AM_CFLAGS} 
 Linker  flags:  ${AM_LDFLAGS}
 Install path:   ${prefix}

----"

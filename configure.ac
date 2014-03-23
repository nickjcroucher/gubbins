AC_INIT([Gubbins], [0.1])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])   

AC_PROG_LIBTOOL
AC_PROG_CC                                   
AC_PROG_CXX

# Define a new recursive target. In our case it produces a .deb installable.
AM_EXTRA_RECURSIVE_TARGETS([release])

AC_CONFIG_FILES([Makefile src/Makefile release/Makefile])
AC_OUTPUT

AC_CHECK_HEADERS([zlib.h math.h])

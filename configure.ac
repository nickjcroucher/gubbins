AC_INIT([Gubbins], [0.1])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])

# Add our own macro dir, so we can check for PTHREAD
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_LIBTOOL
AC_PROG_CC                                   
AC_PROG_CXX

# Checks for pthread
AX_PTHREAD

# Ensure zlib is available
PKG_CHECK_EXISTS([zlib])

# Define a new recursive target. In our case it produces a .deb installable.
AM_EXTRA_RECURSIVE_TARGETS([release])

AC_CONFIG_FILES([Makefile src/Makefile release/Makefile])
AC_OUTPUT

AC_CHECK_HEADERS([zlib.h math.h])

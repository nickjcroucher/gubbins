AC_INIT([Gubbins], [0.1])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])   
AC_PROG_CC                                   

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

PKG_CHECK_EXISTS([zlib], [have_zlib=yes], [have_zlib=no])
if test x$have_zlib = xyes
then
	PKG_CHECK_MODULES(ZLIB, zlib)
else
	dnl We don't always have the luxury of .pc files for zlib, it seems
	have_zlib=yes
	AC_CHECK_LIB([z], [gzopen], [true], [have_zlib=no])
	AC_CHECK_HEADER([zlib.h], [true], [have_zlib=no])
	ZLIB_CFLAGS=
	ZLIB_LIBS="-lz"
fi
AC_SUBST(ZLIB_CFLAGS)
AC_SUBST(ZLIB_LIBS)

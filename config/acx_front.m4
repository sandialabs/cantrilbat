dnl @synopsis ACX_FRONT([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the FRONT
dnl linear-algebra interface (see http://www.netlib.org/front/).
dnl On success, it sets the FRONT_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with FRONT, you should link with:
dnl
dnl 	$FRONT_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-front=<lib> in order to use some
dnl specific FRONT library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the FRONT and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a FRONT
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_FRONT.
dnl
dnl @version $Id: acx_front.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_FRONT], [
acx_front_ok=no

AC_ARG_WITH(front,
	[AC_HELP_STRING([--with-front=<lib>], [use FRONT library <lib>])])
case $with_front in
	yes | "") FRONT_LIBS="-lmpfront -lparlib" ;;
	no) acx_front_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) FRONT_LIBS="$with_front" ;;
	*) FRONT_LIBS="-l$with_front" ;;
esac

# Set linker name of FRONT function to check for.
front_func="mf_setup"

# First, check FRONT_LIBS environment variable
if test "x$FRONT_LIBS" != "x"; then
	save_LIBS="$LIBS"; LIBS="$FRONT_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $front_func in $FRONT_LIBS])
	AC_TRY_LINK_FUNC($front_func, [acx_front_ok=yes], [FRONT_LIBS=""])
	AC_MSG_RESULT($acx_front_ok)
	LIBS="$save_LIBS"
	if test "$acx_front_ok" = "no"; then
		FRONT_LIBS=""
	fi
fi

# FRONT linked to by default?
if test "$acx_front_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$LIBS $FLIBS"
	AC_CHECK_FUNC($front_func, [acx_front_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic FRONT library?
if test "$acx_front_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_LIB(front, $front_func,
		[acx_front_ok=yes; FRONT_LIBS="-lmpfront -lparlib"], [], [$FLIBS])
	LIBS="$save_LIBS"
fi

AC_SUBST(FRONT_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$acx_front_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_FRONT,1,[Define if you have FRONT library.]),[$1])
        :
else
        acx_front_ok=no
        $2
fi
])dnl ACX_FRONT

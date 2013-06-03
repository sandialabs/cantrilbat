dnl @synopsis ACX_Y12M([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the Y12M
dnl linear-algebra interface (see http://www.netlib.org/y12m/).
dnl On success, it sets the Y12M_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with Y12M, you should link with:
dnl
dnl 	$Y12M_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro ,
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-y12m=<lib> in order to use some
dnl specific Y12M library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the Y12M.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a Y12M
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_Y12M.
dnl
dnl @version $Id: acx_y12m.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_Y12M], [
acx_y12m_ok=no

AC_ARG_WITH(y12m,
	[AC_HELP_STRING([--with-y12m=<lib>], [use Y12M library <lib>])])
case $with_y12m in
	yes | "") ;;
	no) acx_y12m_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) Y12M_LIBS="$with_y12m" ;;
	*) Y12M_LIBS="-l$with_y12m" ;;
esac

# Get fortran linker name of Y12M function to check for.
AC_F77_FUNC(y12cck)

# First, check Y12M_LIBS environment variable
if test "x$Y12M_LIBS" != "x"; then
	save_LIBS="$LIBS"; LIBS="$Y12M_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $y12cck in $Y12M_LIBS])
	AC_TRY_LINK_FUNC($y12cck, [acx_y12m_ok=yes], [Y12M_LIBS=""])
	AC_MSG_RESULT($acx_y12m_ok)
	LIBS="$save_LIBS"
	if test "$acx_y12m_ok" = "no"; then
		Y12M_LIBS=""
	fi
fi

# Y12M linked to by default?
if test "$acx_y12m_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$LIBS $FLIBS"
	AC_CHECK_FUNC($y12cck, [acx_y12m_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic Y12M library?
if test "$acx_y12m_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_LIB(y12m, $y12cck,
		[acx_y12m_ok=yes; Y12M_LIBS="-ly12m"], [], [$FLIBS])
	LIBS="$save_LIBS"
fi

AC_SUBST(Y12M_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$acx_y12m_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_Y12M,1,[Define if you have Y12M library.]),[$1])
        :
else
        acx_y12m_ok=no
        $2
fi
])dnl ACX_Y12M

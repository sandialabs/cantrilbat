dnl @synopsis ACX_UMFPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the UMFPACK
dnl linear-algebra interface (see http://www.netlib.org/umfpack/).
dnl On success, it sets the UMFPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with UMFPACK, you should link with:
dnl
dnl 	$UMFPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-umfpack=<lib> in order to use some
dnl specific UMFPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the UMFPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a UMFPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_UMFPACK.
dnl
dnl @version $Id: acx_umfpack.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_UMFPACK], [
AC_REQUIRE([ACX_BLAS])

acx_umfpack_ok=no
acx_umfpack_inc_ok=no
acx_umfpack_lib_ok=no

AC_ARG_WITH(umfpack-lib,
	[AC_HELP_STRING([--with-umfpack-lib=<lib>], [use UMFPACK library <lib>])])
case $with_umfpack_lib in
	yes | "") ;;
	no) acx_umfpack_lib_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) UMFPACK_LIBS="$with_umfpack_lib" 
	;;
#	*) UMFPACK_LIBS="-l$with_umfpack_lib" ;;
	*) UMFPACK_LIBS="-L$with_umfpack_lib -lumfpack-linux -lamd-linux" ;;
esac


AC_ARG_WITH(umfpack-inc,
	[AC_HELP_STRING([--with-umfpack-inc=<inc>],[use UMFPACK include path <inc>])])
case $with_umfpack_inc in
	yes | "") ;;
	no) acx_umfpack_inc_ok=disable ;;
	-* ) UMFPACK_INCS="$with_umfpack_inc" ;;
	*) UMFPACK_INCS="-I$with_umfpack_inc" ;;
esac

# We cannot use UMFPACK if BLAS is not found
if test "x$acx_blas_ok" != "xyes"; then
	acx_umfpack_ok=noblas
	AC_MSG_WARN([The Umfpack solver package requires the BLAS library.])
fi

# Check for the main UMFPACK header file
umf_hdr="umfpack.h"

# Check the UMFPACK_INCS environment variable
if test "x$UMFPACK_INCS" != "x"; then
	CPPFLAGS_SAVE="$CPPFLAGS"
	CPPFLAGS="$CPPFLAGS $UMFPACK_INCS"
	AC_CHECK_HEADERS([$umf_hdr],
		[acx_umfpack_inc_ok=yes],
		[acx_umfpack_inc_ok=no])
	if test "$acx_umfpack_inc_ok" = "no" ; then
		UMFPACK_INCS=""
	fi
	CPPFLAGS="$CPPFLAGS_SAVE"
fi

# Next, check the standard CPP paths
if test "$acx_umfpack_inc_ok" = "no" ; then
	CPPFLAGS_SAVE="$CPPFLAGS"
	CPPFLAGS="$CPPFLAGS"
	UMFPACK_INCS=""
	AC_CHECK_HEADERS([$umf_hdr],
		[acx_umfpack_inc_ok=yes],
		[acx_umfpack_inc_ok=no])
	CPPFLAGS="$CPPFLAGS_SAVE"
fi


# Set linker name of UMFPACK function to check for.
#umf_func="UMF_set_proc_config"
umf_func="umfpack_di_free_symbolic"

EXTRA_LIBS="$BLAS_LIBS"

# First, check the UMFPACK_LIBS environment variable
if test "$acx_umfpack_lib_ok" = "no" ; then
    if test "x$UMFPACK_LIBS" != "x"; then
	save_LIBS="$LIBS"; LIBS="$UMFPACK_LIBS $EXTRA_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $umf_func in $UMFPACK_LIBS])
	AC_TRY_LINK_FUNC($umf_func, [acx_umfpack_lib_ok=yes], [UMFPACK_LIBS=""])
	AC_MSG_RESULT($acx_umfpack_lib_ok)
	LIBS="$save_LIBS"
	if test "$acx_umfpack_lib_ok" = "no"; then
		UMFPACK_LIBS=""
	fi
    fi
fi

# UMFPACK linked to by default?  (is sometimes included in BLAS lib)
if test "$acx_umfpack_lib_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$EXTRA_LIBS $LIBS $FLIBS"
	AC_CHECK_FUNC($umf_func, [acx_umfpack_lib_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic UMFPACK library?
if test "$acx_umfpack_lib_ok" = "no"; then
	AC_CHECK_LIB(umfpack, $umf_func,
		[acx_umfpack_lib_ok=yes; UMFPACK_LIBS="-lumfpack -lamd"], [], [$EXTRA_LIBS $FLIBS])
fi

# OK, we really need both the includes and the libs
if test "$acx_umfpack_inc_ok" = "yes" -a "$acx_umfpack_lib_ok" = "yes" ; then
	acx_umfpack_ok="yes"
else
	AC_MSG_WARN([The Umfpack solver package will not be included.])
	UMFPACK_LIBS=""
	UMFPACK_INCS=""
fi

AC_SUBST(UMFPACK_LIBS)
AC_SUBST(UMFPACK_INCS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_umfpack_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_UMFPACK,1,[Define if you have UMFPACK library.]),[$1])
        :
else
        acx_umfpack_ok=no
        $2
fi
])dnl ACX_UMFPACK

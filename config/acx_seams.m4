dnl @synopsis ACX_SEAMS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a SEAMS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_SEAMS.
dnl
dnl @version $Id: acx_seams.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_SEAMS], [dnl
acx_seams_ok=no
acx_seams_provided=no
# The user can specify the location of a SEAMS installation 
# which will be used later as places to look for EXODUS II and netCDF
# installations if those are not explicitly specified.
#
# The user may simply say "--with-seams" without an argument. In that
# case, likely places are searched in an attempt to find the installation.
#
seams_dir_list=""
seams_dir_list="$seams_dir_list /usr/local/eng_sci/struct/current"
seams_dir_list="$seams_dir_list /usr/local/eng_sci/struct/current-gcc"
seams_dir_list="$seams_dir_list /usr/local/current"
seams_dir_list="$seams_dir_list /projects/seacas/atlantis/current"
seams_dir_list="$seams_dir_list /projects/seacas/current"
seams_dir_list="$seams_dir_list /usr/local/SEAMS"
seams_dir_list="$seams_dir_list /usr/local/SEACAS"
seams_dir_list="$seams_dir_list ./.."
seams_dir_list="$seams_dir_list ./../.."
AC_ARG_WITH([seams],
	    [AC_HELP_STRING([--with-seams=PATH], 
			    [specify SEAMS distribution at PATH])],
	    [acx_seams_provided=yes
	     case x$withval in
		xyes | x) SEAMS_HOME="guess" ;;
		* ) SEAMS_HOME="$withval" ;;
	     esac],
	    [acx_seams_provided=no])
#
# User wants SEAMS and has specified a location...is it valid SEAMS?
#
if test "$acx_seams_provided" = "yes" ; then
	if test x"$SEAMS_HOME" != x"guess" ; then
		AC_MSG_CHECKING([for SEAMS in $SEAMS_HOME])
		rm -fr conftest.dir
		if mkdir conftest.dir ; then
			cd conftest.dir
			cat >Imakefile <<'_ACEOF'
acfindx:
	@echo 'acx_seams_incroot="${INCROOT}"; acx_seams_usrlibdir="${USRLIBDIR}"; acx_seams_libdir="${LIBDIR}"; acx_seams_accessdir="${ACCESSDIR}"'
_ACEOF
			if (${SEAMS_HOME}/etc/accmkmf) >/dev/null 2>/dev/null && test -f Makefile; then
				eval `${MAKE-make} acfindx 2>/dev/null | grep -v make`
			fi
			cd ..
			rm -fr conftest.dir
			if test x"$acx_seams_libdir" != x"" ; then
				AC_MSG_RESULT([yes])
				acx_seams_ok=yes
			else
				AC_MSG_RESULT([no])
			fi
		fi
	fi
fi


if test x"$acx_seams_ok" = x"no" ; then
	AC_MSG_CHECKING([for pre-existing SEAMS])
	rm -fr conftest.dir
	if mkdir conftest.dir; then
		cd conftest.dir
		cat >Imakefile <<'_ACEOF'
acfindx:
	@echo 'acx_seams_incroot="${INCROOT}"; acx_seams_usrlibdir="${USRLIBDIR}"; acx_seams_libdir="${LIBDIR}"; acx_seams_accessdir="${ACCESSDIR}"'
_ACEOF
		# First, try what the user sees in their PATH...
		if (accmkmf) >/dev/null 2>/dev/null && test -f Makefile; then
			eval `${MAKE-make} acfindx 2>/dev/null | grep -v make`
		else
			# Look in a bunch of places...
			for d in $seams_dir_list ; do
				if ($d/etc/accmkmf) > /dev/null 2>/dev/null && test -f Makefile ; then
					eval `${MAKE-make} acfindx 2>/dev/null | grep -v make`		
					break
				fi
			done
		fi
		cd ..
		rm -fr conftest.dir	
	fi
	if test x"$acx_seams_incroot" = x"" ; then
		AC_MSG_RESULT([none])
		AC_MSG_WARN([SEAMS option specified, but no valid SEAMS found])
	else
		acx_seams_ok=yes
		AC_MSG_RESULT([$acx_seams_accessdir])
	fi
fi

# Set some general paths
if test x"$acx_seams_provided" = x"yes" ; then
	SEAMS_LIB="$acx_seams_usrlibdir"
	SEAMS_INC="$acx_seams_incroot"
	AC_SUBST(SEAMS_LIB)
	AC_SUBST(SEAMS_INC)
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$acx_seams_ok" = "xyes"; then
        ifelse([$1],[],AC_DEFINE(HAVE_SEAMS,1,[Define if you have SEAMS library.]),[$1])
        :
else
        acx_seams_ok=no
        $2
fi
])dnl ACX_SEAMS

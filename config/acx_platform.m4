dnl @synopsis ACX_PLATFORM()
dnl
dnl This macro is used to set a CPP defines for GOMAs platform dependence.
dnl Optional.
dnl @version $Id: acx_platform.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Patrick Notz <pknotz@sandia.gov>
dnl
AC_DEFUN([ACX_PLATFORM], [dnl
acx_platform_provided="no"
acx_platform_value=""
# User defined platform
AC_ARG_WITH([platform],
	    [AC_HELP_STRING([--with-platform=PLATFORM],
	                    [use PLATFORM as hint])],
	    [acx_platform_provided=yes
	     AC_MSG_CHECKING([specified Goma platform $withval])
	     case "$withval" in
	     	*linux*)            acx_platform_value="linux" ;;
		*solaris* | *sun* ) acx_platform_value="sun";;
		*hp*)               acx_platform_value="hpux";;
		*dec* | *osf* )     acx_platform_value="dec_osf1";;
		*aix* | *ibm* )     acx_platform_value="_AIX";;
		*sgi* | *irix* )    acx_platform_value="sgi";;
		*)                  acx_platform_value="$withval";;
	     esac
	     AC_MSG_RESULT([$acx_platform_value])],
	     [acx_platform_provided="no"
              acx_platform_value="none"])
])dnl ACX_PLATFORM

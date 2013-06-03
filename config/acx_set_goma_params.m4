dnl @synopsis ACX_SET_GOMA_PARAMS()
dnl
dnl This macro is used to set CPP defines for GOMA.
dnl
dnl @version $Id: acx_set_goma_params.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Patrick Notz <pknotz@sandia.gov>
dnl
AC_DEFUN([ACX_SET_GOMA_PARAMS], [

# Set MDE
AC_ARG_WITH(mde,[AC_HELP_STRING([--with-mde=<int>],[set MDE to <int>])])
case $with_mde in
	yes) AC_MSG_ERROR([Usage is:  --with-mde=<int>]);;
	no | "") with_mde=no ;;
	*) MDE="$with_mde" ;;
esac

# Let's make sure this is an <int>
if test "$with_mde" != "no" ; then
	test $with_mde -gt 0 > /dev/null 2> /dev/null
    result="$?"
    case $result in
	    2) AC_MSG_ERROR([Usage is:  --with-mde=<int>]) ;;
	    1) AC_MSG_ERROR([You must supply a positive integer for MDE]) ;;
	    *) AC_DEFINE_UNQUOTED(MDE,$MDE,[Set the maximum number of nodes per element])
	       AC_MSG_NOTICE([Setting the maximum number of nodes per element (MDE) to $MDE]) ;;
    esac
fi

# Set MAX_PROB_VAR
AC_ARG_WITH(max_prob_var,[AC_HELP_STRING([--with-max_prob_var=<int>],[set MAX_PROB_VAR to <int>])])
case $with_max_prob_var in
	yes) AC_MSG_ERROR([Usage is:  --with-max_prob_var=<int>]);;
	no | "") with_max_prob_var=no ;;
	*) MAX_PROB_VAR="$with_max_prob_var" ;;
esac

# Let's make sure this is an <int>
if test "$with_max_prob_var" != "no" ; then
    test $with_max_prob_var -gt 0 > /dev/null 2> /dev/null
    result="$?"
    case $result in
	    2) AC_MSG_ERROR([Usage is:  --with-max_prob_var=<int>]);;
	    1) AC_MSG_ERROR([You must supply a positive integer for MAX_PROB_VAR]);;
	    *) AC_DEFINE_UNQUOTED(MAX_PROB_VAR,$MAX_PROB_VAR,[Set the maximum number of problem variables])
	       AC_MSG_NOTICE([Setting the maximum number of problem variables (MAX_PROB_VAR) to $MAX_PROB_VAR]) ;;
    esac
fi
# Set MAX_CONC
AC_ARG_WITH(max_conc,[AC_HELP_STRING([--with-max_conc=<int>],[set MAX_CONC to <int>])])
case $with_max_conc in
	yes) AC_MSG_ERROR([Usage is:  --with-max_conc=<int>]);;
	no | "") with_max_conc=no ;;
	*) MAX_CONC="$with_max_conc" ;;
esac

# Let's make sure this is an <int>
if test "$with_max_conc" != "no" ; then
    test $with_max_conc -gt 0 > /dev/null 2> /dev/null
    result="$?"
    case $result in
	    2) AC_MSG_ERROR([Usage is:  --with-max_conc=<int>]);;
	    1) AC_MSG_ERROR([You must supply a positive integer for MAX_CONC]);;
	    *) AC_DEFINE_UNQUOTED(MAX_CONC,$MAX_CONC,[Set the maximum number of species])
	       AC_MSG_NOTICE([Setting the maximum number of species (MAX_CONC) to $MAX_CONC]) ;;
    esac
fi

# Set COUPLED_FILL
AC_ARG_WITH(coupled_fill,[AC_HELP_STRING([--with-coupled_fill],[define COUPLED_FILL])])
case $with_coupled_fill in
	yes) AC_DEFINE(COUPLED_FILL,1,[Make the level set algorithm fully coupled])
             AC_MSG_NOTICE([Compiling for the fully coupled level set algorithm]) ;;
	*) AC_MSG_NOTICE([Compiling for the de-coupled (subcycling) level set algorithm]) ;;

esac

# Set DEBUG
AC_ARG_WITH(debug,[AC_HELP_STRING([--with-debug],[define DEBUG])])
case $with_debug in
	yes) AC_DEFINE(DEBUG,1,[Compile in extra debugging code])
             AC_MSG_NOTICE([Compiling in extra debugging code]) ;;
	*) : ;;
esac


])dnl ACX_SET_GOMA_PARAMS

/**
 * @file mdp_allo.h
 *   Declarations for Multi Dimensional Pointer (mdp) malloc routines, which
 *   allow for dimensioning of arbitrarily dimensioned pointer arrays.
 */
/*
 * $Author: hkmoffa $
 * $Revision: 508 $
 * $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef MDP_ALLO_H_UTILSRC
#define MDP_ALLO_H_UTILSRC

/*
 * Include the header here in order to pick up size_t definition
 */
#include <cstring>

/**
 *
 *  General Principles
 *     If an allocation size is set to 0, then the actual
 *     allocation size is set to 1 within the program.
 *     Something is always allocated.
 *
 *     All checks for allocations are always checked for success.
 *     If a failure is found, thene mdp_alloc_eh() is called
 *     for disposition of the error condition.
 *
 */

namespace mdpUtil {

/****************************************************************************/
/*
*  If we have array_alloc() from another Sandia program, we will not use
*  the one from this mdp_array_alloc. Instead we will redefine the names
*/

#ifdef HAVE_ARRAY_ALLOC
#  define mdp_array_alloc array_alloc
#  define mdp_safe_free   safe_free
#endif

/*!
 *  MDP_INT_NOINIT is a poor man's way of specifying whether a value should be
 *  initialized. These are seldom used numbers which can be used in place
 *  of real ints and dbls to indicate that initialization shouldn't take
 *  place.
 */
#define MDP_INT_NOINIT          -68361

/*!
 *  MDP_DBL_NOINIT is a poor man's way of specifying whether a value should be
 *  initialized. These are seldom used numbers which can be used in place
 *  of real ints and dbls to indicate that initialization shouldn't take
 *  place.
 */
#define MDP_DBL_NOINIT          -1.241E11


#ifndef _C16_NAME_DEF
#  define _C16_NAME_DEF
   typedef char    C16_NAME[16]; /* Character array used for fortran names */
   typedef char    C16_NAME_STR[17];
#endif
/****************************************************************************/
/*
 * Externals that should be set by the calling program.
 * These are only used for debugging purposes.
 */
#ifdef MDP_MPDEBUGIO
extern int MDP_MP_Nprocs;
extern int MDP_MP_myproc;
#endif

/****************************************************************************/
/*
 * MDP_SAFE_DELETE() 
 *       This useful define is great for delete single instances of
 *       mallocing using new.
 */
#define MDP_SAFE_DELETE(a) if (a) { delete (a); a = 0; }


/****************************************************************************/

#define mdp_alloc_struct(x, num) (x *) mdp_array_alloc(1, (num), sizeof(x))


/* function declarations for dynamic array allocation */

//! allocates multidimensional pointer arrays of arbitrary length
//! via a single malloc call
/*!
 *  The first dimension is the number of dimensions in the allocation
 */
extern double *mdp_array_alloc(int numdim, ...);

//! Free a vector and set its value to 0
/*!
 *  This function carries out the following operation
 *  @code
 *    free(*hndVec);
 *    *hndVec = 0;
 *  @endcode
 *
 * @param hndVec This is the address of the pointer, expressed as
 *               a void **
 *
 * Note, a key idea behind the mdp suite of routines is that
 * a pointer that can be malloced is either malloced or its value
 * is 0. This routine enforces this convention.
 */
extern void mdp_safe_free(void ** hndVec);

//! Allocate a vector of integers
/*!
 *  The vector is initialized, unless the default int value is set
 * to MDP_INT_NOINIT
 * 
 * @param len  Length of the vector
 * @param defval Default value for the int, defaults to MDP_INT_NOINIT
 *
 * @return returns a pointer to the vector
 */
extern int *mdp_alloc_int_1(int len, const int defval = MDP_INT_NOINIT);

//! Allocate a vector of integers, potentially freeing memory first
/*!
 *  The vector is initialized, unless the default int value is set to MDP_INT_NOINIT
 *
 * Input
 * --------------
 * @param array_hdl Previous value of pointer. If non-NULL will try
 *                to free the memory at this address before doing
 *                a new alloc
 * @param len  Length of the vector
 * @param defval Default value for the int, defaults to MDP_INT_NOINIT
 *
 * Output 
 * ---------
 * @return  *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 */
extern void mdp_safe_alloc_int_1(int **array_hdl, int len, 
				 const int defval = MDP_INT_NOINIT);

//!    Reallocates a one dimensional array of ints, copying old
//!    information to the new array
/*!
 *    Reallocates a one dimensional array of ints. This routine always allocates space for at least one int.
 *    Calls the smalloc() routine to ensure that all malloc
 *    calls go through one location. This routine will then copy
 *    the pertinent information from the old array to the
 *    new array.
 *
 *     NewArray[0:old_len-1] = OldArray[0:old_len-1];
 *     NewArray[old_len:new_len-1] = defval;
 *
 * Input
 * --------------
 * @param array_hdl Previous value of pointer. If non-NULL will try
 *                to free the memory at this address before doing
 *                a new alloc
 * @param new_len  New Length of the vector
 * @param old_len  New Length of the vector
 * @param defval Default value for the int, defaults to MDP_INT_NOINIT
 *
 * Output 
 * ---------
 * @return  *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 */
extern void mdp_realloc_int_1(int **array_hdl, int new_len, int old_len, 
			      const int defval = MDP_INT_NOINIT);


//! Allocates a 2D matrix of integers with column pointers located at the front of the allocation
/*!
 *    The matrix is initialized on allocation, unless the default int value is set
 *    to MDP_INT_NOINIT, which is the default.
 *
 *     matrix[len1][len2]
 *
 *  All int data entries are contiguous. Therefore, it may
 *  be used as input into BLAS matrix function calls.
 *  This can be considered to be in fortran order format with
 *  len2 as the number of rows, and len1 as the number of columns.
 *
 *  matrix[jcol] refers to the jcol column of the matrix.
 *  Therefore, matrix[0] is a pointer to the beginning of the
 *  data portion of the structure.
 *  The structure will have len1 pointers at the beginning
 *  that holds pointers into the top of the columns of the
 *  contiguous data.
 *
 *  The entire structure may be deallocated via one free call.
 * 
 *   @param[in]   len1        Outer Length of the vector , i.e, the row in fortran format
 *   @param[in]   len2        Inner length of the matrix , .e., the column in fortran format
 *   @param[in]   defval      Default value for the int, defaults to MDP_INT_NOINIT
 *
 *   @return                  Returns a pointer to the malloced matrix
 */
extern int **mdp_alloc_int_2(int len1, int len2, const int defval = MDP_INT_NOINIT);

//! Allocate a vector of doubles
/*!
 *  The vector is initialized, unless the default double value is set to MDP_DBL_NOINIT
 * 
 *   @param[in]    len1          Length of the vector
 *   @param[in]    defval        Default value for the double, defaults to MDP_DBL_NOINIT
 *
 *   @return                     Returns a pointer to the malloced vector
 */
extern double *mdp_alloc_dbl_1(int len1, const double defval = MDP_DBL_NOINIT);


extern void    mdp_safe_alloc_dbl_1(double **, int, const double); 

//! Reallocate a vector of doubles possibly retaining a subset of values
/*!
 *  Reallocates the array and sets:
 *
 *  (*hndVec)[0:oldLen-1]    = oldVec[0:oldLen-1]
 *  (*hndVec)[oldLen:newLen] = defVal
 *  
 *   @param[in,out]  hndVec     Handle to the old allocated vector on input. On return
 *                              this contains the handle to the new allocated vector.
 *   @param[in]      newLen     New Length of the vector
 *   @param[in]      oldLen     Old Length of the vector
 *   @param[in]      defVal     Default fill value to be used on the new members of the double
 */
extern void    mdp_realloc_dbl_1(double ** hndVec, int newLen, int oldLen, 
                                 const double defVal);

extern void    mdp_realloc_dbl_2(double ***, int, int, int, int,
				 const double);

extern char   *mdp_alloc_char_1(int, const char = '\0');
extern void    mdp_safe_alloc_char_1(char **, int, const char = '\0');    
extern char  **mdp_alloc_VecFixedStrings(int, int);
extern void    mdp_safe_alloc_VecFixedStrings(char ***, int, int);
extern void    mdp_realloc_VecFixedStrings(char ***, int,  int, int);

extern double **mdp_alloc_dbl_2(int, int, const double);
extern void    mdp_safe_alloc_dbl_2(double ***, int, int, 
				    const double = MDP_DBL_NOINIT);

extern C16_NAME *mdp_alloc_C16_NAME_1(int, const int);
extern void    mdp_safe_alloc_C16_NAME_1(C16_NAME **, int, const int);

extern void  **mdp_alloc_ptr_1(int);
extern void    mdp_safe_alloc_ptr_1(void ***, int);
extern void    mdp_realloc_ptr_1(void ***, int, int);
extern void    mdp_copy_ptr_1(void *const, const void *const, int);
extern void  **mdp_dupl_ptr_1(const void * const, int);

extern char   *mdp_copy_C16_NAME_to_string(const C16_NAME);
extern void    mdp_copy_VecFixedStrings(char ** const, const char ** const,
					int, size_t);

//! Copy a cstring 
/*!
 *  @param cstring cstring to be copied. It is not changed by the routine
 *
 *  @return Returns a malloced copy of the string
 */
extern char   *mdp_copy_string(const char * cstring);

extern void    mdp_safe_copy_string(char **, const char *);

//! Copy a double vector to a double vector
/*!
 *  copyTo[len] = copyFrom[len] 
 *
 * Input
 * -------
 * @param copyFrom Vector to  copy ( length >= len)
 * @param len      Length of the copy
 *
 * Output
 * -------
 * @param copyTo   Vector to receive the copy ( length >= len)
 */
extern void    mdp_copy_dbl_1(double * const copyTo, 
                              const double * const copyFrom, 
                              int len);

//! Copy a double array to a double array
/*!
 *  copyTo[len1][len2] = copyFrom[len1][len2]
 *
 * Input
 * --------
 * @param copyFrom Vector to  copy ( length >= len1 * len2)
 * @param len1      Length of the first dimension
 * @param len2      Length of the second dimension
 *
 * Output
 * ----------
 * @param copyTo   Array to receive the copy ( length >= len1 * len2)
 */
extern void    mdp_copy_dbl_2(double ** const copyTo, const double ** const copyFrom,
			      int len1, int len2);

extern void    mdp_copy_int_1(int * const,    const int * const, int);
extern void    mdp_copy_int_2(int ** const,   const int ** const, int, int);
extern void    mdp_init_dbl_1(double * const, double, int);
extern void    mdp_zero_dbl_1(double * const, int);
extern void    mdp_zero_int_1(int * const, int);
extern void    mdp_init_dbl_2(double ** const v, double value, 
			      int len1, int len2);
extern void    mdp_init_int_1(int * const,    int, int);
/*
 * subtractRD.cpp
 */
extern double subtractRD(double, double);
extern "C" double subtractrd_(double *, double *);
/* 
 * checkFinite
 */
//extern void checkZeroFinite(double);
//extern void checkFinite(double);
//extern "C" void checkfinite_(double *);
//extern void checkMagnitude(double tmp, double trigger = 1.0E20);
/****************************************************************************/
}
#endif
/****************************************************************************/


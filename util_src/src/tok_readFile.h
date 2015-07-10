/**
 * @file tok_readFile.h
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef TOK_READFILE_H
#define TOK_READFILE_H

#include "tok_input_util.h"
#include <vector>

//----------------------------------------------------------------------------------------------------------------------------------
namespace TKInput
{
//==================================================================================================================================
//! Read a list of doubles from a file.
/*!
 *   This routine reads a file consisting of lines consisting of two doubles per line  putting the contents into two arrays
 *   of doubles. There is an option of having a header consisting of one integer, containing the number of lines of doubles.
 *
 *     example file:
 * 
 *            3
 *            0.0 1.0
 *            1.0 6.54444
 *            2.0 12.5555555
 *
 *    @param[in]      file_name         Name of the file to read the doubles from
 *    @param[out]     a1                First column consisting of doubles
 *    @param[out]     a2                Second column consisting of doubles
 *
 *    @return                           Returns the number of lines of doubles read as an int.
 */
int tok_read2dbl(const char* const file_name, std::vector<double>& a1, std::vector<double>& a2);

//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
/**************************************************************************/

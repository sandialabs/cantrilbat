/**
 * @file tok_readFile.h
 *
 * $Author: hkmoffa $
 * $Revision: 5 $
 * $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
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
namespace TKInput {
  extern int tok_read2dbl(const char * const,
			  std::vector<double>&, std::vector<double>&);
}
/**************************************************************************/
#endif 
/**************************************************************************/

/**
 *  @file m1d_CanteraElectrodeGlobals.h
 *    Contains the global variables that are defined within this directory 
 */
/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_CANTERAELECTRODEGLOBALS_H
#define M1D_CANTERAELECTRODEGLOBALS_H

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d {
 
class ProblemStatementCell;

//==================================================================================================================================
//!  We define a global ProblemStatementCell object pointer
/*!
 *   It's a pointer because we may be using an inherited object in its stead. We only
 *   need one of them because we are only solving one problem at a time.
 *
 *  All object need information about the problem environment. This information is available from  this pointer.
 */
extern ProblemStatementCell *PSCinput_ptr;

//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

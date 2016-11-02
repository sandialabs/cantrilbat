/**
 *  @file importPL.h
 *  Declarations for utility routine to read files containing phase descriptions into PhaseList objects.
 *  (see \ref ExtendedPhaseGroups ).
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef IMPORTPL_H
#define IMPORTPL_H

#include "Electrode_defs.h"

#include <string>


#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

class PhaseList;

//! Import all phases found in a single file into a PhaseList object
/*!
 *  Import all phases found in a single file into a PhaseList object, in an additive fashion.
 *  This returns the number of phases found, processed, and added  to the PhaseList object.
 *
 * @param       pl              Pointer to the PhaseList object
 * @param       canteraFile     Cantera CTML file
 *
 * @return                      Returns the number of phases added to the PhaseList object.
 */
size_t importAllCTMLIntoPhaseList(PhaseList* const pl, const std::string& canteraFile);

}
#endif

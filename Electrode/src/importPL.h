/**
 *  @file importAllCTML.h
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/**
 *  $Id: importPL.h 571 2013-03-26 16:44:21Z hkmoffa $
 *
 */

#ifndef IMPORTPL_H
#define IMPORTPL_H


#include "cantera/kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"


#include <string>

namespace Cantera
{

class PhaseList;
class XML_Node;
class Transport;



//! Import all phases found in a single file into a PhaseList object
/*!
 *  Import all phases found in a single file into a PhaseList object,
 *  in an additive fashion.
 *  This returns the number of phases found, processed, and added
 *  to the PhaseList object.
 *
 * @param pl             Pointer to the PhaseList object
 * @param canteraFile    Cantera CTML file
 */
int importAllCTMLIntoPhaseList(PhaseList* pl, std::string canteraFile);


}
#endif

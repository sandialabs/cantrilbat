/**
 *  @file importAllCTML.h
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/**
 *  $Id: importPL.h 507 2013-01-07 22:48:29Z hkmoffa $
 *
 */

#ifndef IMPORTPL_H
#define IMPORTPL_H


#include "zuzax/kinetics.h"
#include "zuzax/kinetics/InterfaceKinetics.h"


#include <string>

namespace Zuzax
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
   * @param canteraFile    Zuzax CTML file
   */
  int importAllCTMLIntoPhaseList(PhaseList *pl, std::string canteraFile);


}
#endif

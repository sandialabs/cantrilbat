/**
 *  @file importAllCTML.h
 */

/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef IMPORTALLCTML_H
#define IMPORTALLCTML_H

#include "cantera/base/ct_defs.h"
#include "cantera/thermo/ThermoPhase.h"

#include <string>
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

class ThermoPhase;
class Kinetics;
class InterfaceKinetics;
class PhaseList;
class XML_Node;
class Transport;

//==================================================================================================================================
//! Create a single ThermoPhase object, not currently supported by the main Zuzax distribution
/*!
 *    Import a single phase found in a single file.
 *    Add new thermophase definitions, ones not currently supported
 *    by the main Zuzax distribution.
 *
 *       All phases are imported and instantiated by using the
 *       constructor call with a reference to the phase XML object.
 *       This starts a process whereby the ThermoPhase object
 *       completely initializes itself.
 *
 *       For models that we understand, we pop the error off
 *       of the Zuzax stack.
 *       Models that we don't understand, or that cause an error  we rethrow the CanteraError
 *
 *  @param[in]               xmlphase            XML Node containing the phase to be imported.
 *
 *  @return                                      Returns a malloced pointer to the new phase
 */
thermo_t_double* processExpandedThermoPhase(XML_Node* const xmlphase);

//==================================================================================================================================
//! Process a kinetics manager, which may or may not be currently supported by the main Zuzax distribution.
/*!
 *    Import a Kinetics manager from a CTML description.
 *    Add new Kinetics definitions, ones not currently supported
 *    by the main Zuzax distribution.
 *
 *       For models that we understand, we pop the error off
 *       of the Cantera stack. This includes the "NONE" model.
 *       It's ok to return with the null pointer from this routine.
 *
 *  @param[in]               xmlPhase            pointer to the XML node containing the phase information
 *
 *  @param[in]               tpList              STL vector of pointers to ThermoPhase objects. These must be the correct objects
 *                                               for the kinetics manager class.
 *
 *  @return                                      Returns a malloced pointer to the new kinetics manager
 */
Kinetics* processExpandedKinetics(XML_Node* const xmlPhase, std::vector<thermo_t_double*> tpList);

//==================================================================================================================================
//! Process an interface kinetics manager, which may or may not be currently supported by the main Cantera distribution.
/*!
 *    Import an interface Kinetics manager from a CTML description.
 *    Add new Kinetics definitions, ones not currently supported by the main Zuzax distribution.
 *
 *       For models that we understand, we pop the error off
 *       of the Cantera stack. This includes the "NONE" model.
 *       It's ok to return with the null pointer from this routine.
 *
 *  @param[in]               xmlPhase            pointer to the XML node containing the phase information
 *
 *  @param[in]               tpList              STL vector of pointers to ThermoPhase objects. These must be the correct objects
 *                                               for the kinetics manager class.
 *
 *  @return                                     Returns a malloced pointer to the new InterfaceKinetics manager
 */
InterfaceKinetics* processExpandedInterfaceKinetics(XML_Node* const xmlPhase, std::vector<thermo_t_double*> tpList);

//==================================================================================================================================
//! Process transport properties for a phase.
/*!
 *  Process an expanded set of transport properties for a phase.
 *  This calls the Cantera routine newtransportMgr to process
 *  Zuzax's transport first. Then, it processes any new transport properties.
 *
 *  @param[in]               xmlPhase            Pointer to the XML node containing the phase information
 * 
 *  @param[in]               th                  Pointer to the previously processed %ThermoPhase object for the phase in question
 *
 *  @return                                      Returns a malloced pointer to the new InterfaceKinetics manager
 */
Transport* processExpandedTransport(const XML_Node* const xmlPhase, thermo_t_double* const th);

//==================================================================================================================================
}
//--------------------------------------------------------------------------------------------------------------------------------
#endif

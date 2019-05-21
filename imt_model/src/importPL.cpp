/**
 *  @file importCTML.cpp
 *
 *     This file contains routines which are global routines, i.e.,
 *     not part of any object. These routine take as input, ctml
 *     pointers to data, and pointers to Zuzax objects. The purpose
 *     of these routines is to intialize the Zuzax objects with data
 *     from the ctml tree structures.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/**
 *  $Id: importPL.cpp 507 2013-01-07 22:48:29Z hkmoffa $
 *
 */


//   Zuzax includes

#include "importPL.h"

#include "importAllCTML.h"
#include "zuzax/thermo/mix_defs.h"

#include "zuzax/kinetics.h"


#include "zuzax/thermo/speciesThermoTypes.h"
#include "zuzax/thermo/ThermoPhase.h"
#include "zuzax/thermo/ThermoFactory.h"

#include "zuzax/transport.h"

#include "zuzax/multiphase/PhaseList.h"
#include "zuzax/thermo/IdealSolidSolnPhase.h"
#include "zuzax/thermo/IdealMolalSoln.h"
#include "zuzax/thermo/DebyeHuckel.h"
#include "zuzax/thermo/StoichSubstanceSSTP.h"
#include "zuzax/thermo/HMWSoln.h"

#include "zuzax/base/xml.h"
#include "zuzax/base/ctml.h"
#include "zuzax/kinetics/InterfaceKinetics.h"


#include <vector>
#include "stdio.h"

//using namespace ctml;
using namespace std;

namespace Zuzax
{


//==================================================================================================================================
//!  Given an XML_Node pointing to a phase, add the phase to a PhaseList object
/*!
 *   Import all phases found in a single file into a PhaseList object.
 *   All phases are imported and instantiated by using the constructor call with a reference to the phase XML object.
 *   This starts a process whereby the ThermoPhase object completely initializes itself. The phase is separated into separate lists
 *   based on the "dim" attribute.
 *
 *   @param[in]           xmlphase                  XML_Node pointer, pointing to a phase named "phase"
 *   @param[in]           pl                        Pointer to the PhaseList
 *   @param[in]           canteraFile               Name of the file, used for error printouts
 */
static void processPhasePL(XML_Node* const xmlphase, PhaseList* const pl, const std::string& canteraFile)
{
    ThermoPhase* tPhase = processExpandedThermoPhase(xmlphase);
    //ThermoPhase* tPhase = newPhase(*xmlphase);
    if (!tPhase) {
        throw ZuzaxError("processPhasePL()",
                           "ERROR: tPhase = 0 while processing phase in file, " + canteraFile);
    }
    std::string dimS = xmlphase->operator[]("dim");
    if (dimS == "3") {
        pl->addVolPhase(tPhase);
    } else if (dimS == "2") {
        pl->addSurPhase(tPhase);
    } else {
        throw ZuzaxError("processPhasePL",
                           "While processing file, " + canteraFile + ", unknown dim string: " + dimS);
    }
}
//================================================================================================================================

  /*
   *
   * static routine findXMLAllPhasePL
   *
   *    Import all phases found in a single file into a PhaseList object,
   *    in a recursive and additive fashion. pl may or may not contain
   *    phases going into this routine.
   *
   */
  static void findXMLAllPhasePL(XML_Node *root, PhaseList *pl, std::string canteraFile)
  {
    XML_Node *sc = 0;
    if (!root) return;
    string idattrib;
    string rname = root->name();
    if (rname == "phase") {
      processPhasePL(sc, pl, canteraFile);
    }                                              
    const vector<XML_Node*> &vsc = root->children();
    for (size_t n = 0; n < root->nChildren(); n++) {
      sc = vsc[n];
      if (sc->name() == "phase") {
	processPhasePL(sc, pl, canteraFile);
      } else {
	findXMLAllPhasePL(sc, pl, canteraFile);
      }
    }
  }

  /*
   *
   * importAllCTMLIntoPhaseList()
   *
   *  Import all phases found in a single file into a PhaseList object,
   *  in an additive fashion. 
   *  This returns the number of phases found, processed, and added
   *  to the PhaseList object.
   *
   *      pl -> Pointer to the PhaseList object
   *      canteraFile -> Zuzax CTML file
   *
   */
  int importAllCTMLIntoPhaseList(PhaseList *pl, std::string canteraFile) {
    XML_Node *xc = 0;
    try {
      xc = get_XML_File(canteraFile);
    }  catch (ZuzaxError) {
      showErrors();
      throw ZuzaxError("importAllCTMLIntoPhaseList",
			 string("Could not find/process file, ") +
			 canteraFile + string(" -> aborting"));
    }
    if (!xc) {
      throw ZuzaxError("importAllCTMLIntoPhaseList",
			 string("Could not find/process file, ") +
			 canteraFile + string(" -> aborting")); 
    }
    findXMLAllPhasePL(xc, pl, canteraFile);
    int nphases = pl->nPhases();
    return nphases;
  }


 
  /*************************************************************************/
}    
 

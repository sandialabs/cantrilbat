/**
 *  @file importCTML.cpp
 *
 *     This file contains routines which are global routines, i.e.,
 *     not part of any object. These routine take as input, ctml
 *     pointers to data, and pointers to Cantera objects. The purpose
 *     of these routines is to intialize the Cantera objects with data
 *     from the ctml tree structures.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/**
 *  $Id: importAllCTML.cpp 508 2013-01-07 22:54:04Z hkmoffa $
 *
 */


//   Cantera includes

#include "importAllCTML.h"
#include "cantera/thermo/mix_defs.h"

#include "cantera/kinetics.h"
#include "SolidKinetics.h"


#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"

#include "cantera/transport.h"

//#include "PhaseList.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/StoichSubstanceSSTP.h"
#include "cantera/thermo/HMWSoln.h"

#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"
#include "cantera/kinetics/InterfaceKinetics.h"

#include "ElectrodeKinetics.h"
#include "ElectrolyteKinetics.h"

#include <vector>
#include "stdio.h"

using namespace ctml;
using namespace std;

namespace Cantera {

 
  // Create a single ThermoPhase object, not currently supported
  // by the main Cantera distribution
  /*
   *    Import a single phase found in a single file.
   *    Add new thermophase definitions, ones not currently supported
   *    by the main cantera distribution.
   * 
   * Note: all phases are imported and instantiated by using the
   *       constructor call with a reference to the phase XML object.
   *       This starts a process whereby the ThermoPhase object
   *       completely initializes itself.
   *
   *       For models that we understand, we pop the error off
   *       of the Cantera stack.
   *
   *       Models that we don't understand, or that cause an error
   *       return with the null pointer. This is considered an
   *       error.
   *
   * @param xmlphase  XML Node containing the phase to be imported.
   */
  ThermoPhase *processExpandedThermoPhase(XML_Node *xmlphase)
  {
    ThermoPhase *tPhase = 0;
    try {
      tPhase = newPhase(*xmlphase);
    } catch(UnknownThermoPhaseModel &uName) {
      const XML_Node &th = xmlphase->child("thermo");
      string model = th["model"];
      printf("\t\t model = %s\n", model.c_str());
      /*
       *  IdealSolidSolnPhase: This is not yet in the Factory, but is in
       *                       the main Cantera distribution.
       */
      if (model == "IdealSolidSolution") {
	popError();
	tPhase = new IdealSolidSolnPhase(*xmlphase);
      }
      /*
       * IdealMolalSoln: This is not yet in the Factory, but is in
       *                 the main Cantera distribution.
       */
      else if (model == "IdealMolalSoln") {
	popError();
	tPhase = new IdealMolalSoln(*xmlphase);
      } 
      /*
       * DebyeHuckel:    This is not yet in the Factory, but is in
       *                 the main Cantera distribution.
       */
      else if (model == "DebyeHuckel") {
	popError();
	tPhase = new DebyeHuckel(*xmlphase);
      } 
      /*
       * HMWSoln:        This is not yet in the Factory, but is in
       *                 the main Cantera distribution.
       */
      else if (model == "HMW" || model == "HMWSoln") {
	popError();
	tPhase = new HMWSoln(*xmlphase);
      }
      /*
       * StoichSubstance: This is not yet in the Factory, but is in
       *                  the main Cantera distribution.
       */
      else if (model == "StoichSubstance") {
	popError();
	tPhase = new StoichSubstanceSSTP(*xmlphase);
      } else {
	throw;
      }
    }
    return tPhase;
  }

  //  Process a kinetics manager, which may or may not be
  //  currently supported by the main Cantera distribution.
  /*
   *    Import a Kinetics manager from a CTML description.
   *    Add new Kinetics definitions, ones not currently supported
   *    by the main cantera distribution.
   *
   *       For models that we understand, we pop the error off
   *       of the Cantera stack. This includes the "NONE" model.
   *       It's ok to return with the null pointer from this routine.
   *
   * @param xmlPhase   Pointer to the XML node containing the phase 
   *                   information
   *  @param tpList    STL vector of pointers to ThermoPhase
   *                   objects. These must be the correct objects
   *                   for the kinetics manager class.
   */
  Kinetics *processExpandedKinetics(XML_Node *xmlPhase,
				    vector<ThermoPhase*> tpList) 
  {
    Kinetics *kin = 0;
    if (tpList.size() == 0) {
      printf("processExpandedKinetics ERROR: expecting tpList to "
	     "be filled before entry\n");
      exit(-1);
    }
    ThermoPhase *tp = tpList[0];
    if (!tp) {
      printf("processExpandedKinetics ERROR: expecting tpList to "
	     "be filled before entry\n");
      exit(-1);
    }

    XML_Node &kinNode = xmlPhase->child("kinetics");
    string kModel = kinNode.attrib("model");
    if (kModel == "" || (!strcasecmp(kModel.c_str(), "None"))) {
      return (Kinetics *) 0;
    } 
    /*
     * SolidKinetics:       This is a New one.
     *                      It's for solid phases
     */
    else if (kModel == "SolidKinetics") {
      SolidKinetics *sk_ptr = new SolidKinetics(tp);
      sk_ptr->importMechanism(*xmlPhase, "");
      kin = sk_ptr;
    } else {
      try {
	kin = newKineticsMgr(*xmlPhase, tpList);
      } catch (UnknownKineticsModel) {
	XML_Node &kinNode = xmlPhase->child("kinetics");
	string kModel = kinNode.attrib("model");
	if (kModel == "SolidKinetics") {
	  SolidKinetics *sk_ptr = new SolidKinetics(tp);
	  sk_ptr->importMechanism(*xmlPhase, "");
	  kin = sk_ptr;
	  popError();
	} else if (kModel == "Electrolyte") {
	  ElectrolyteKinetics *ek_ptr = new ElectrolyteKinetics(tp);
	  ek_ptr->importMechanism(*xmlPhase, "");
	  kin = ek_ptr;
	} else if (kModel == "NONE") {
	  popError();
	} else { 
	  showErrors(cout);
	}
      }
    }
    return kin;
  }
  
  //   Process an interface kinetics manager, which may or may not be
  //   currently supported by the main Cantera distribution.
  /*
   *    Import an InterfaceKinetics manager from a CTML description.
   *    Add new Kinetics definitions, ones not currently supported
   *    by the main cantera distribution.
   *
   *       For models that we understand, we pop the error off
   *       of the Cantera stack. This includes the "NONE" model.
   *       It's ok to return with the null pointer from this routine.
   *
   * @param xmlPhase pointer to the XML node containing the phase 
   *                   information
   *  @param tpList   STL vector of pointers to ThermoPhase
   *                  objects. These must be the correct objects
   *                  for the kinetics manager class.
   */
  InterfaceKinetics *
  processExpandedInterfaceKinetics(XML_Node *xmlPhase,
				   std::vector<ThermoPhase*> tpList) 
  {
    InterfaceKinetics *kin = 0;
    if (tpList.size() == 0) {
      printf("processExpandedKinetics ERROR: expecting tpList to "
	     "be filled before entry\n");
      exit(-1);
    }
    ThermoPhase *tp = tpList[0];
    if (!tp) {
      printf("processExpandedKinetics ERROR: expecting tpList to "
	     "be filled before entry\n");
      exit(-1);
    }

    XML_Node &kinNode = xmlPhase->child("kinetics");
    string kModel = kinNode.attrib("model");
    if (kModel == "" || (!strcasecmp(kModel.c_str(), "None"))) {
      return (InterfaceKinetics *) 0;
    } else {
      try {
	Kinetics *kinBase = newKineticsMgr(*xmlPhase, tpList);
        kin = dynamic_cast<InterfaceKinetics *>(kinBase);
        if (kin == 0) {
          throw CanteraError(" processExpandedInterfaceKinetics",
                             " Dynamic cast to InterfaceKinetics failed"); 
        }
      } catch (UnknownKineticsModel &ee) {
	XML_Node &kinNode = xmlPhase->child("kinetics");
	string kModel = kinNode.attrib("model");
	if (kModel == "Electrode") {
	  popError();
	  ElectrodeKinetics *ek_ptr = new ElectrodeKinetics();
	  ek_ptr->importKinetics(*xmlPhase, tpList);
	  kin = ek_ptr;
	} else if (kModel == "NONE") {
	  kin = 0;
	  popError();
	} else { 
	  showErrors(cout);
	  throw ee;
	}
      }
    }
    return kin;
  }

  //! Process transport properties for a phase.
  /*
   *  Process an expanded set of transport properties for a phase.
   *  This calls the Cantera routine newtransportMgr to process
   *  cantera's transport first. Then, it processes any new
   *  transport properties.
   *
   *  Currently, there are no new ones.
   *
   *  @param xmlPhase  pointer to the XML node containing the phase 
   *                   information
   *  @param th        Pointer to the previously processed ThermoPhase
   *                   object.
   */
  Transport *processExpandedTransport(XML_Node *xmlPhase, ThermoPhase *tp) {
    Transport *tran = 0;
    if (xmlPhase->hasChild("transport")) {
      XML_Node &tranNode = xmlPhase->child("transport");
      string tModel = tranNode.attrib("model");
      if (tModel == "" || (!strcasecmp(tModel.c_str(), "None"))) {
	return tran;
      } else {
	tran = newTransportMgr(tModel, tp);
      }
    }
    return tran;
  }




 
  /*************************************************************************/
}    
 

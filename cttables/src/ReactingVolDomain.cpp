/**
 * @file ReactingVolDomain.cpp
 *
 * $Author: hkmoffa $
 * $Revision: 497 $
 * $Date: 2013-01-07 14:17:04 -0700 (Mon, 07 Jan 2013) $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "ReactingVolDomain.h"
#include "importAllCTML.h"
#include "PhaseList.h"
#include "cantera/base/ctml.h"

#include <cstdio>

extern int DebugPrinting;
using namespace std;
namespace Cantera {

  /**
   * Constructor #1 for the IdealReactingGas object.
   * 
   * infile = input file location
   * id = ?
   */
  ReactingVolDomain::ReactingVolDomain() :
    m_kinetics(0),
    m_InterfaceKinetics(0),
    m_transport(0),
    numPhases(0),
    iphaseKin(-1),
    m_DoSurfKinetics(false),
    m_DoHomogKinetics(false),
    m_ok(false),
    m_XMLPhaseTree(0),
    m_XMLTree_owned(true),
    m_I_Own_Thermo(false)
  {

  }

  /**
   * Constructor #2
   * 
   * root = Previously built xml tree structure (not owned
   *        by this object 
   * id = ?
   */
/*
  ReactingVolDomain::
  ReactingVolDomain(XML_Node& phaseRoot, string id) :
    m_kinetics(0),
    m_InterfaceKinetics(0),
    m_transport(0),
    numPhases(0),
    iphaseKin(-1),
    m_DoSurfKinetics(false),
    m_DoHomogKinetics(false),
    m_ok(false),
    m_XMLPhaseTree(0),
    m_XMLTree_owned(false),
    m_I_Own_Thermo(false)
  {
    m_ok = importFromXML(phaseRoot);
    if (!m_ok) {
      throw CanteraError("ReactingVolDomain",
			 "buildSolutionFromXML returned false");
    }
  }
*/ 
  /**
   * Destructor for the ReactingVolDomain object.
   * 
   * We must decide whether this object owns its own xml tree
   * structure.
   */       
  ReactingVolDomain::~ReactingVolDomain() { 
    if (m_kinetics) {
      delete m_kinetics; m_kinetics = 0;
    }
    if (m_InterfaceKinetics) {
      delete m_InterfaceKinetics; m_InterfaceKinetics = 0;
    }

    if (m_XMLTree_owned) {
      delete m_XMLPhaseTree;
    }
    if (m_I_Own_Thermo) {
      for (int i = 0; i < numPhases; i++) {
	delete tpList[i];
	tpList[i] = 0;
      }
    }
  
  }
   
  /** 
   *  This ostream function describes how to extend cout output
   *  functions to this object. The way this is done is to
   *  call the Cantera's report function, which takes a ThermoPhase
   *  object as its arugment.
   *  This is a "friend" function to the class IdealReactingGas.
   *  Both this function and report are in the Cantera namespace.
   *
   *  Note -> The output doesn't cover kinetics.
   */
  std::ostream& operator<<(std::ostream& s, 
			   ReactingVolDomain& mix) {
    ThermoPhase *th;
    for (int i = 0; i < mix.numPhases; i++) {
      th = mix.tpList[i];
      std::string r = th->report(true);
      s << r;
    }
    return s;
  }

  /**
   *
   *
   */
  bool ReactingVolDomain::importFromXML(XML_Node& phaseRoot)  {
    try {
      int iph;
      /*
       * Vector of pointers to all of the phases in the file
       */
      XML_Node *TopCTML = &phaseRoot;
      if (phaseRoot.name() != "ctml") {
	TopCTML = phaseRoot.findNameID("ctml", "");
      }
      if (!TopCTML) {
	throw CanteraError("ReactingVolDomain::importFromXML","ctml not found");
      }

      vector<XML_Node *> phaseChildren;
      string nm = "phase";
      TopCTML->getChildren(nm, phaseChildren);

      int nPhasesFound = phaseChildren.size();
      /*
       * Resize the internal list of pointers and
       * get a pointer to the vacant ThermoPhase pointer
       */
      tpList.resize(nPhasesFound, 0);
      tplRead.resize(nPhasesFound, 0);
      kinOrder.resize(nPhasesFound, -1);
      xmlList.resize(nPhasesFound, 0);

      iphaseKin = -1; 
      XML_Node *phaseArrayXML = 0;
      for (iph = 0; iph < nPhasesFound; iph++) {
	XML_Node *xmlPhase = phaseChildren[iph];
	XML_Node *kineticsXML = xmlPhase->findNameID("kinetics", "");
	phaseArrayXML = xmlPhase->findNameID("phaseArray", "");
	if (kineticsXML) {
	  if (iphaseKin >=0) {
	    if (phaseArrayXML) {
	      iphaseKin = iph;
	    }
	  } else {
	    iphaseKin = iph;
	  }
	}
      }
	  
    
      phaseArrayXML = 0;
      if (iphaseKin >= 0) {
	XML_Node *xmlPhase = phaseChildren[iphaseKin];
	phaseArrayXML = xmlPhase->findNameID("phaseArray", "");
	if (phaseArrayXML) {
	  vector<string> phase_ids;
	  ctml::getStringArray(*phaseArrayXML, phase_ids);
	  int np = phase_ids.size();
	  for (iph = 0; iph < np; iph++) {
	    string phaseID = phase_ids[iph];
	    bool found = false;
	    for (int jph = 0; jph < nPhasesFound; jph++) {
	      if (phaseChildren[jph]->id() == phaseID) {
		XML_Node *xmlPhase = phaseChildren[jph];
		xmlList[numPhases] = xmlPhase;
		tpList[numPhases] = processExpandedThermoPhase(xmlPhase);
		tplRead[jph] = 1;
		kinOrder[numPhases] = iph;
		numPhases++;
	      }
	    }
	    if (!found) {
	      throw CanteraError("import", "phase not found");
	    }
	  }
	}
	if (!tplRead[iphaseKin]) {
	  XML_Node *xmlPhase = phaseChildren[iphaseKin];
	  xmlList[numPhases] = xmlPhase;
	  tpList[numPhases] = processExpandedThermoPhase(xmlPhase);
	  tplRead[iphaseKin] = 1;
	  kinOrder[numPhases] = iphaseKin;
	  numPhases++;
	}
      }

      
      
      for (iph = 0; iph < nPhasesFound; iph++) {
	if (iphaseKin != iph) {
	  /*
	   * Fill in the ThermoPhase object by querying the
	   * const XML_Node tree located at x.
	   */
	  XML_Node *xmlPhase = phaseChildren[iph];
	  xmlList[numPhases] = xmlPhase;
	  tpList[numPhases] = processExpandedThermoPhase(xmlPhase);
	  kinOrder[numPhases] = iphaseKin;
	  tplRead[iphaseKin] = 1;
	  numPhases++;
	}
      }

      for (iph = 0; iph < nPhasesFound; iph++) {
	if (! tplRead[iphaseKin]) {
	  XML_Node *xmlPhase = phaseChildren[iph];
	  xmlList[numPhases] = xmlPhase;
	  tpList[numPhases] = processExpandedThermoPhase(xmlPhase);
	  kinOrder[numPhases] = iph;
	  tplRead[iphaseKin] = 1;
	  numPhases++;
	}
      }

      if (nPhasesFound > numPhases) {	
        for (iph = 0; iph < nPhasesFound - numPhases; iph++) {
	  int jph = numPhases + iph;
	  if (tpList[jph] != 0) {
	    throw CanteraError(" ReactingVolDomain::importFromXML" , "Confused tpList[]");
	  }
 	}
	tpList.resize(numPhases);
      }

      /*
       * Fill in the kinetics object k, by querying the
       * const XML_Node tree located by x. The source terms and
       * eventually the source term vector will be constructed
       * from the list of ThermoPhases in the vector, phases.
       */
      if (iphaseKin >= 0) {
	XML_Node *xmlPhase = phaseChildren[iphaseKin];
	if (!phaseArrayXML) {
	  m_kinetics = processExpandedKinetics(xmlPhase, tpList);
	  if (!m_kinetics) {
	    if (DebugPrinting) {
	      printf("No kinetics object was found - that's ok\n");
	    }
	  }
	} else {
	  m_InterfaceKinetics = 
	    processExpandedInterfaceKinetics(xmlPhase, tpList);
	  if (!m_InterfaceKinetics) {
	    if (DebugPrinting) {
	      printf("No interface kinetics object was found - that's ok\n");
	    }
	  }
	}
      }
      return true;

    }
    catch (CanteraError) {
      showErrors(cout);
      throw CanteraError("ReactingVolDomain::importFromXML",
			 "error encountered");
      return false;
    }
  }

  /*
   *
   *
   */
  bool ReactingVolDomain::importFromPL(Cantera::PhaseList *pl, int ivkin, int iskin)  {
    try {
      int iph;
      
      XML_Node *kinXMLPhase = 0;
      ThermoPhase *kinPhase = 0;

      //XML_Node *vPhase = pl->VolPhaseXMLNodes[0];
      if (iskin >= 0) {
	kinXMLPhase = pl->surPhaseXMLNode(iskin);
	kinPhase = &(pl->surPhase(iskin));
      } else    if (ivkin >= 0) {
	kinXMLPhase = pl->volPhaseXMLNode(ivkin);
	kinPhase = &(pl->volPhase(ivkin));
      }
      //AssertThrow(vPhase, "vPhase must be defined");


      int nPhasesFound = pl->nVolPhases() + pl->nSurPhases();
      /*
       * Resize the internal list of pointers and
       * get a pointer to the vacant ThermoPhase pointer
       */
      tpList.resize(nPhasesFound, 0);
      tplRead.resize(nPhasesFound, 0);
      kinOrder.resize(nPhasesFound, -1);
      xmlList.resize(nPhasesFound, 0);
      iphaseKin = -1; 
      if (iskin >= 0) {
	iphaseKin = iskin;
	m_DoSurfKinetics = true;
      } else if (ivkin >= 0) {
	iphaseKin = ivkin;
	m_DoHomogKinetics = true;
      }
  

      numPhases = 0;
      if (iphaseKin >= 0) {
	xmlList[numPhases] = kinXMLPhase;
	tpList[numPhases] = kinPhase;
	tplRead[numPhases] = 1;
	kinOrder[numPhases] = iphaseKin;
	numPhases++;
      }

     
      XML_Node *phaseArrayXML = 0;
      if (iphaseKin >= 0) {
	XML_Node *xmlPhase = kinXMLPhase;
	phaseArrayXML = xmlPhase->findNameID("phaseArray", "");
	if (phaseArrayXML) {
	  std::vector<std::string> phase_ids;
	  ctml::getStringArray(*phaseArrayXML, phase_ids);
	  int npToFind = phase_ids.size();
	  for (iph = 0; iph < npToFind; iph++) {
	    string phaseID = phase_ids[iph];
	    bool found = false;
	    for (int jph = 0; jph < pl->nVolPhases(); jph++) {
	      XML_Node *xmlPhase_j = pl->volPhaseXMLNode(jph);
	      string pname = xmlPhase_j->operator[]("id");
	      if (phaseID == pname) {
		found = true;
		xmlList[numPhases] = xmlPhase_j;
		tpList[numPhases] = &(pl->volPhase(jph));
		tplRead[jph] = 1;
		kinOrder[numPhases] = iph;
		numPhases++;
		break;
	      }
	    }

	    if (!found) {
	      throw CanteraError("ReactingVolDomain::importFromPL",
				 "Phase, requested in phaseArray, was not found: " + phaseID);
	    }
	  }
	}
      } else {
	for (iph = 0; iph < pl->nVolPhases(); iph++) {
	  xmlList[numPhases] = pl->volPhaseXMLNode(iph);
	  tpList[numPhases]  = &(pl->volPhase(iph));
	  tplRead[numPhases] = 1;
	  kinOrder[numPhases] = iph;
	  numPhases++;
	}
      }

      if (nPhasesFound > numPhases) {
	for (iph = 0; iph < nPhasesFound - numPhases; iph++) {
	  int jph = iph +  numPhases;
	  if (tpList[jph] != 0) {
	    throw CanteraError("ReactingVolDomain::importFromPL",
			       " Confused tpList");
	  }
	}
	tpList.resize(numPhases);
      }
     
      /*
       * Fill in the kinetics object k, by querying the
       * const XML_Node tree located by x. The source terms and
       * eventually the source term vector will be constructed
       * from the list of ThermoPhases in the vector, phases.
       */
      if (iphaseKin >= 0) {
	XML_Node *xmlPhase =  pl->volPhaseXMLNode(iphaseKin);
	if (!phaseArrayXML) {
	  m_kinetics = processExpandedKinetics(xmlPhase, tpList);
	  if (!m_kinetics) {
	    if (DebugPrinting) {
	      printf("No kinetics object was found - that's ok\n");
	    }
	  }
	} else {
	  xmlPhase =  pl->surPhaseXMLNode(iphaseKin);
	  m_InterfaceKinetics = 
	    processExpandedInterfaceKinetics(xmlPhase, tpList);
	  if (!m_InterfaceKinetics) {
	    if (DebugPrinting) {
	      printf("No interface kinetics object was found - that's ok\n");
	    }
	  }
	}
      }
      return true;

    }
    catch (CanteraError) {
      showErrors(cout);
      throw CanteraError("ReactingVolDomain::importFromXML",
			 "error encountered");
      return false;
    }
  }
}



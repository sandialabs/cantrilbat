/**
 *  @file importCTML.cpp
 *
 *     This file contains routines which are global routines, i.e.,
 *     not part of any object. These routine take as input, ctml
 *     pointers to data, and pointers to Zuzax objects. The purpose
 *     of these routines is to intialize the Zuzax objects with data
 *     from the ctml tree structures.
 */
/**
 *  $Id: importAllPhases.cpp 511 2013-01-07 23:32:30Z hkmoffa $
 *
 */
/*
 * Zuzax includes: The paths change if you are compiling it 
 * as part of the cantera build tree or as an application
 */
#ifdef COMPILE_IN_CANTERA_BUILDTREE
#include "mix_defs.h"
#include "speciesThermoTypes.h"
#include "ThermoPhase.h"
#include "kernelThermoFactory.h"
#include "xml.h"
#else
#include "zuzax/thermo/speciesThermoTypes.h"
#include "zuzax/thermo/ThermoPhase.h"
#include "zuzax/thermo/ThermoFactory.h"
#include "zuzax/base/xml.h"
#include "zuzax/base/ctml.h"
#endif

#include "importAllPhases.h"
#include "IdealSolidSolnPhase.h"
#include "mdp_allo.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
using namespace std;
using namespace mdpUtil;

/**************************************************************************
 *
 * importAllCTML() static routine processPhasePL
 *
 *
 *    Import all phases found in a single file into a PhaseList object.
 *
 */
static void processPhasePL(XML_Node *xmlphase,  EPEQUIL_INPUT *pi)
{
    ThermoPhase *tPhase = 0;
    try {
      tPhase = newPhase(*xmlphase);
    } catch(UnknownThermoPhaseModel &uName) {
	 
      printf("ERROR -> importAllCTML , caught unknown model\n");
      const XML_Node &th = xmlphase->child("thermo");
      string model = th["model"];
      printf("\t\t model = %s\n", model.c_str());
      if (model == "IdealSolidSolution") {
	tPhase = new IdealSolidSolnPhase(*xmlphase);
      } else {
	printf("unknown model\n");
	exit(-1);
      }
    }

    mdp_realloc_ptr_1((void ***)&(pi->tplist),
		      pi->nphase+1, pi->nphase);
    pi->tplist[pi->nphase] = tPhase;
    pi->nphase++;
    MultiPhase *mp = pi->m_mp;
    mp->addPhaseWithMoles(tPhase, 0.0);
}

/**************************************************************************
 *
 * importAllCTML() -> static routine findXMLAllPhasePL
 *
 *    Import all phases found in a single file into a PhaseList object,
 *    in a recursive and additive fashion. pl may or may not contain
 *    phases going into this routine.
 *
 */
static void findXMLAllPhasePL(XML_Node *root, EPEQUIL_INPUT *pi)
{
    XML_Node *sc = 0;
    if (!root) return;
    string idattrib;
    string rname = root->name();
    if (rname == "phase") {
      processPhasePL(sc, pi);
    }                                              
    const vector<XML_Node*> &vsc = root->children();
    for (int n = 0; n < (int) root->nChildren(); n++) {
      sc = vsc[n];
      if (sc->name() == "phase") {
	processPhasePL(sc, pi);
      } else {
	findXMLAllPhasePL(sc, pi);
      }
    }
}

/**************************************************************************
 *
 * importAllCTML()
 *
 *    Import all phases found in a single file into a PhaseList object.
 *  This returns the number of phases found and processed.
 */
int importAllCTML(EPEQUIL_INPUT *pi, string canteraFile) {
    XML_Node *xc = 0;
    try {
      xc = get_XML_File(canteraFile);
    }  catch (ZuzaxError) {
      throw ZuzaxError("importAllCTML",
			 string("Could not find/process file, ") +
			 canteraFile + string(" -> aborting"));
    }
    if (!xc) {
      throw ZuzaxError("importAllCTML",
			 string("Could not find/process file, ") +
			 canteraFile + string(" -> aborting")); 
    }
    findXMLAllPhasePL(xc, pi);
    int nphases = pi->nphase;
    return nphases;
}

/*************************************************************************/


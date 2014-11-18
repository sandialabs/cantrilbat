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
 *  $Id: importPL.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 *
 */


//   Cantera includes

#include "importPL.h"
#include "PhaseList.h"

#include "importAllCTML.h"

#include <vector>
#include "cstdio"

using namespace std;

namespace Cantera
{



/*
 *
 * importAllCTML() static routine processPhasePL
 *
 *
 *    Import all phases found in a single file into a PhaseList object.
 * Note: all phases are imported and instantiated by using the
 *       constructor call with a reference to the phase XML object.
 *       This starts a process whereby the thermophase object
 *       completely initializes itself.
 */
static void processPhasePL(XML_Node* xmlphase, PhaseList* pl, std::string canteraFile)
{
    ThermoPhase* tPhase = processExpandedThermoPhase(xmlphase);
    //ThermoPhase* tPhase = newPhase(*xmlphase);
    if (!tPhase) {
        printf("processPhasePL ERROR: tPhase = 0\n");
        exit(-1);
    }
    string dimS = xmlphase->operator[]("dim");
    if (dimS == "3") {
        pl->addVolPhase(tPhase, xmlphase, canteraFile);
    } else if (dimS == "2") {
        pl->addSurPhase(tPhase, xmlphase, canteraFile);
    } else {
        throw CanteraError("processPhasePL", "unknown dim string: " + dimS);
    }
}

/*
 *
 * static routine findXMLAllPhasePL
 *
 *    Import all phases found in a single file into a PhaseList object,
 *    in a recursive and additive fashion. pl may or may not contain
 *    phases going into this routine.
 *
 */
static void findXMLAllPhasePL(XML_Node* root, PhaseList* pl, std::string canteraFile)
{
    XML_Node* sc = 0;
    if (!root) {
        return;
    }
    string idattrib;
    string rname = root->name();
    if (rname == "phase") {
        processPhasePL(sc, pl, canteraFile);
    }
    const vector<XML_Node*>& vsc = root->children();
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
 *      canteraFile -> Cantera CTML file
 *
 */
int importAllCTMLIntoPhaseList(PhaseList* pl, std::string canteraFile)
{
    XML_Node* xc = 0;
    try {
        xc = get_XML_File(canteraFile);
    }  catch (CanteraError) {
        showErrors();
        throw CanteraError("importAllCTMLIntoPhaseList",
                           string("Could not find/process file, ") +
                           canteraFile + string(" -> aborting"));
    }
    if (!xc) {
        throw CanteraError("importAllCTMLIntoPhaseList",
                           string("Could not find/process file, ") +
                           canteraFile + string(" -> aborting"));
    }
    findXMLAllPhasePL(xc, pl, canteraFile);
    int nphases = pl->nPhases();
    return nphases;
}

/*************************************************************************/
}


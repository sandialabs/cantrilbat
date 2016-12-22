/**
 *  @file importPL.cpp
 *    Definitions for utility routine to read files containing phase descriptions into PhaseList objects.
 *    (see \ref ExtendedPhaseGroups ).
 */
/*
 *     This file contains routines which are global routines, i.e.,
 *     not part of any object. These routine take as input, ctml
 *     pointers to data, and pointers to Cantera objects. The purpose
 *     of these routines is to intialize the Cantera objects with data from the ctml tree structures.
 */

/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


#include "cantera/multiphase/PhaseList.h"
#include "cantera/base/xml.h"

#include "importPL.h"
#include "importAllCTML.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
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
    thermo_t_double* tPhase = processExpandedThermoPhase(xmlphase);
    if (!tPhase) { 
	throw CanteraError("processPhasePL()",
			   "ERROR: tPhase = 0 while processing phase in file, " + canteraFile);
    }
    std::string dimS = xmlphase->operator[]("dim");
    size_t iphGlob;
    if (dimS == "3") {
       iphGlob =  pl->addVolPhase(tPhase, xmlphase);
    } else if (dimS == "2") {
        iphGlob = pl->addSurPhase(tPhase, xmlphase);
    } else if (dimS == "1") {
        iphGlob = pl->addEdgePhase(tPhase, xmlphase);
    } else {
        throw CanteraError("processPhasePL()", "While processing file, " + canteraFile + ", unknown dim string: " + dimS);
    }
    pl->movePhaseOwnership(2, iphGlob);
}
//==================================================================================================================================
//!  Recursive search for XML_Node files named phase
/*!
 *   Import all phases found in a single file into a PhaseList object,
 *   in a recursive and additive fashion. pl may or may not contain phases going into this routine.
 *
 *   @param[in]           root                       Beginning XML_Node to check
 *   @param[in]           pl                         Pointer to the phase list object
 *   @param[in]           canteraFile                Name of the file, used for error output only
 *   @param[in]           nrecursive                 Number of levels to search recursively. If it is one, then only
 *                                                   the current node and its children will be search for the name "phase"
 */
static void findXMLAllPhasePL(XML_Node* const root, PhaseList* const pl, const std::string& canteraFile, int nrecursive)
{
    XML_Node* sc = 0;
    if (!root) {
        return;
    }
    std::string idattrib;
    std::string rname = root->name();
    if (rname == "phase") {
        processPhasePL(sc, pl, canteraFile);
    }
    const std::vector<XML_Node*>& vsc = root->children();
    for (size_t n = 0; n < root->nChildren(); n++) {
        sc = vsc[n];
        if (sc->name() == "phase") {
            processPhasePL(sc, pl, canteraFile);
        } else {
	    if (nrecursive > 0) {
		findXMLAllPhasePL(sc, pl, canteraFile, nrecursive - 1);
	    }
        }
    }
}
//==================================================================================================================================
size_t importAllCTMLIntoPhaseList(PhaseList* const pl, const std::string& canteraFile)
{
    XML_Node* xc = 0;
    try {
        xc = get_XML_File(canteraFile);
    }  catch (CanteraError) {
        showErrors();
        throw CanteraError("importAllCTMLIntoPhaseList", "Could not find/process file, " + canteraFile + " -> aborting");
    }
    if (!xc) {
        throw CanteraError("importAllCTMLIntoPhaseList", "Could not find/process file, " + canteraFile + " -> aborting");
    }
    // Search the first 3 levels of the XML tree for a phase node. -> I don't think there is any need to go further.
    findXMLAllPhasePL(xc, pl, canteraFile, 2);
    size_t nphases = pl->nPhases();
    return nphases;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------


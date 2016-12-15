/**
 * @file ReactingVolDomain.cpp
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "ReactingVolDomain.h"
#include "importAllCTML.h"
#include "cantera/multiphase/PhaseList.h"
#include "cantera/base/ctml.h"

#include <cstdio>

extern int DebugPrinting;

using namespace std;
//---------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
ReactingVolDomain::ReactingVolDomain() :
    m_kinetics(0),
    m_InterfaceKinetics(0),
    m_transport(0),
    m_NumKinPhases(0),
    m_NumPLPhases(0),
    m_NumKinSpecies(0),
    m_iphGlobKin(npos),
    m_DoSurfKinetics(false),
    m_DoHomogKinetics(false),
    m_ok(false),
    m_XMLPhaseTree(0),
    m_XMLTree_owned(true),
    m_I_Own_Thermo(false)
{
}
//==================================================================================================================================
ReactingVolDomain::~ReactingVolDomain()
{
    if (m_kinetics) {
        delete m_kinetics;
        m_kinetics = 0;
    }
    if (m_InterfaceKinetics) {
        delete m_InterfaceKinetics;
        m_InterfaceKinetics = 0;
    }

    if (m_XMLTree_owned) {
        delete m_XMLPhaseTree;
    }
    if (m_I_Own_Thermo) {
        for (size_t i = 0; i < tpList.size(); i++) {
            delete tpList[i];
            tpList[i] = 0;
        }
    }
}
//==================================================================================================================================
std::ostream& operator<<(std::ostream& s, ReactingVolDomain& mix)
{
    ThermoPhase* th;
    for (size_t i = 0; i < mix.tpList.size(); i++) {
        th = mix.tpList[i];
        std::string r = th->report(true);
        s << r;
    }
    return s;
}
//==================================================================================================================================
bool ReactingVolDomain::importFromXML(XML_Node& phaseRoot)
{
    size_t numPhases = 0;
    try {
        //int iph;
        /*
         * Vector of pointers to all of the phases in the file
         */
        XML_Node* TopCTML = &phaseRoot;
        if (phaseRoot.name() != "ctml") {
            TopCTML = phaseRoot.findNameID("ctml", "");
        }
        if (!TopCTML) {
            throw CanteraError("ReactingVolDomain::importFromXML","ctml not found");
        }

        std::vector<XML_Node*> phaseChildren;
        std::string nm = "phase";
        TopCTML->getChildren(nm, phaseChildren);

        size_t nPhasesFound = phaseChildren.size();
        /*
         * Resize the internal list of pointers and
         * get a pointer to the vacant ThermoPhase pointer
         */
        tpList.resize(nPhasesFound, 0);
        tplRead.resize(nPhasesFound, 0);
        kinOrder.resize(nPhasesFound, npos);
        xmlList.resize(nPhasesFound, 0);



        m_iphGlobKin = npos;
        XML_Node* phaseArrayXML = 0;
        if (nPhasesFound >= 0) {
            for (size_t iph = 0; iph < nPhasesFound; iph++) {
                XML_Node* xmlPhase = phaseChildren[iph];
                xmlList[iph] = xmlPhase;
                tpList[iph] = processExpandedThermoPhase(xmlPhase);
                XML_Node* kineticsXML = xmlPhase->findNameID("kinetics", "");
                phaseArrayXML = xmlPhase->findNameID("phaseArray", "");
                if (kineticsXML) {
                    if (m_iphGlobKin != npos) {
                        if (phaseArrayXML) {
                            m_iphGlobKin = iph;
                        }
                    } else {
                        m_iphGlobKin = iph;
                    }
                }
            }
        }



        numPhases = 0;
        phaseArrayXML = 0;
        if (m_iphGlobKin >= 0) {
            XML_Node* xmlPhase = phaseChildren[m_iphGlobKin];
            phaseArrayXML = xmlPhase->findNameID("phaseArray", "");
            if (phaseArrayXML) {
                std::vector<std::string> phase_ids;
                ZZctml::getStringArray(*phaseArrayXML, phase_ids);
                size_t np = phase_ids.size();
                for (size_t iph = 0; iph < np; iph++) {
                    std::string phaseID = phase_ids[iph];
                    bool found = false;
                    for (size_t jph = 0; jph < nPhasesFound; jph++) {
                        if (phaseChildren[jph]->id() == phaseID) {
                            XML_Node* xmlPhase = phaseChildren[jph];
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
            if (!tplRead[m_iphGlobKin]) {
                XML_Node* xmlPhase = phaseChildren[m_iphGlobKin];
                xmlList[numPhases] = xmlPhase;
                tpList[numPhases] = processExpandedThermoPhase(xmlPhase);
                tplRead[m_iphGlobKin] = 1;
                kinOrder[numPhases] = m_iphGlobKin;
                numPhases++;
            }
        }

        for (size_t iph = 0; iph < nPhasesFound; iph++) {
            if (m_iphGlobKin != iph) {
                /*
                 * Fill in the ThermoPhase object by querying the
                 * const XML_Node tree located at x.
                 */
                XML_Node* xmlPhase = phaseChildren[iph];
                xmlList[numPhases] = xmlPhase;
                tpList[numPhases] = processExpandedThermoPhase(xmlPhase);
                kinOrder[numPhases] = m_iphGlobKin;
                tplRead[m_iphGlobKin] = 1;
                numPhases++;
            }
        }

        for (size_t iph = 0; iph < nPhasesFound; iph++) {
            if (! tplRead[m_iphGlobKin]) {
                XML_Node* xmlPhase = phaseChildren[iph];
                xmlList[numPhases] = xmlPhase;
                tpList[numPhases] = processExpandedThermoPhase(xmlPhase);
                kinOrder[numPhases] = iph;
                tplRead[m_iphGlobKin] = 1;
                numPhases++;
            }
        }

        if (nPhasesFound > (size_t) numPhases) {
            for (size_t iph = 0; iph < nPhasesFound - numPhases; iph++) {
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
        if (m_iphGlobKin != npos) {
            XML_Node* xmlPhase = phaseChildren[m_iphGlobKin];
            if (!phaseArrayXML) {
                m_kinetics = processExpandedKinetics(xmlPhase, tpList);
                if (!m_kinetics) {
                    if (DebugPrinting) {
                        printf("No kinetics object was found - that's ok\n");
                    }
                }
            } else {
                m_InterfaceKinetics = processExpandedInterfaceKinetics(xmlPhase, tpList);
                if (!m_InterfaceKinetics) {
                    if (DebugPrinting) {
                        printf("No interface kinetics object was found - that's ok\n");
                    }
                }
            }
        }
        return true;

    } catch (CanteraError) {
        showErrors(cout);
        throw CanteraError("ReactingVolDomain::importFromXML", "error encountered");
        return false;
    }
}
//==================================================================================================================================
bool ReactingVolDomain::importSurKinFromPL(PhaseList* pl, size_t iphSurKin)
{
    m_iphGlobKin = pl->phaseIndexFromSurPhaseIndex(iphSurKin);
    m_DoSurfKinetics = true;
    
    m_NumPLPhases = pl->nPhases();
    tpList.resize(m_NumPLPhases);
    for (size_t iph = 0; iph < m_NumPLPhases; iph++) {
        tpList[iph] = &(pl->thermo(iph));
    }
    
    XML_Node* xmlPhase = &(pl->surPhaseXMLNode(iphSurKin));
    m_InterfaceKinetics = processExpandedInterfaceKinetics(xmlPhase, tpList);
    if (!m_InterfaceKinetics) {
        if (DebugPrinting) {
              printf("No interface kinetics object was found - that's ok\n");
        }
    } else { 

        m_NumKinPhases = m_InterfaceKinetics->nPhases();
        PLtoKinPhaseIndex_.resize(m_NumPLPhases, npos);
        PLtoKinSpeciesIndex_.resize(pl->nSpecies(), npos);
        KintoPLSpeciesIndex_.resize(m_InterfaceKinetics->nKinSpecies(), npos);

        for (size_t kph = 0; kph < m_NumKinPhases; kph++) {
            thermo_t_double& tt = thermo(kph);
            std::string kname = tt.id();
            phaseID kID = tt.phaseIdentifier();
            size_t jph = npos;
            for (size_t iph = 0; iph < m_NumPLPhases; iph++) {
                ThermoPhase& pp = pl->thermo(iph);
                std::string iname = pp.id();
                phaseID iID = pp.phaseIdentifier();
                if (iID == kID) {
                    jph = iph;
                    break;
                }
            }
            if (jph == npos) {
                throw CanteraError("ElectrodeKinetics_intoPL::importFromPL()", "phase not found");
            }
            kinOrder[kph] = jph;
            PLtoKinPhaseIndex_[jph] = kph;

            size_t PLkstart = pl->globalSpeciesIndex(jph, 0);
            size_t nspPhase = tt.nSpecies();
            size_t kstart = m_InterfaceKinetics->kineticsSpeciesIndex(kph, 0);
            for (size_t k = 0; k < nspPhase; k++) {
                if (PLtoKinSpeciesIndex_[k + PLkstart] != npos) {
                    throw CanteraError("ElectrodeKinetics_intoPL::importFromPL()",
                                       "Indexing error found while initializing  PLtoKinSpeciesIndex_");
                }
                
                PLtoKinSpeciesIndex_[k + PLkstart] = kstart + k;
                KintoPLSpeciesIndex_[kstart + k] = PLkstart + k;
            }
        }
    }
    return true;
}
//==================================================================================================================================
bool ReactingVolDomain::importVolKinFromPL(PhaseList* pl, size_t iphVolKin)
{
    try {
        m_NumPLPhases = pl->nPhases();
        tpList.resize(m_NumPLPhases);
        for (size_t iph = 0; iph < m_NumPLPhases; iph++) {
            tpList[iph] = &(pl->thermo(iph));
        }
        if (iphVolKin != npos) {
	    XML_Node* xmlPhase = &(pl->volPhaseXMLNode(iphVolKin));
            m_kinetics = processExpandedKinetics(xmlPhase, tpList);
            if (!m_kinetics) {
                if (DebugPrinting) {
                    printf("No kinetics object was found - that's ok\n");
                }
                return false;
            }
            /*
             *  Create a mapping between the ReactingSurfPhase to the PhaseList phase
             */
        }
        m_NumKinPhases = m_kinetics->nPhases();
        m_NumKinSpecies = m_kinetics->nKinSpecies();
        kinOrder.resize(m_NumKinPhases, npos);
        PLtoKinPhaseIndex_.resize(pl->nPhases(), npos);
        PLtoKinSpeciesIndex_.resize(pl->nSpecies(), npos);
        KintoPLSpeciesIndex_.resize(m_NumKinSpecies, npos);
        //size_t kKinIndex = 0;
        for (size_t kph = 0; kph < m_NumKinPhases; kph++) {
            thermo_t_double& tt = thermo(kph);
            std::string kname = tt.id();
            phaseID kID = tt.phaseIdentifier();
            size_t jph = npos;
            for (size_t iph = 0; iph < pl->nPhases(); iph++) {
                ThermoPhase& pp = pl->thermo(iph);
                std::string iname = pp.id();
                phaseID iID = pp.phaseIdentifier();
                if (iID == kID) {
                    jph = iph;
                    break;
                }
            }
            if (jph == npos) {
                throw CanteraError("ElectrodeKinetics_intoPL::importFromPL()", "phase not found");
            }
            kinOrder[kph] = jph;
            PLtoKinPhaseIndex_[jph] = kph;

            size_t PLkstart = pl->globalSpeciesIndex(jph, 0);
            size_t nspPhase = tt.nSpecies();
            size_t kstart = m_kinetics->kineticsSpeciesIndex(kph, 0);
            for (size_t k = 0; k < nspPhase; k++) {
                if (PLtoKinSpeciesIndex_[k + PLkstart] != npos) {
                    throw CanteraError("ElectrodeKinetics_intoPL::importFromPL()",
                                       "Indexing error found while initializing  PLtoKinSpeciesIndex_");
                }
                PLtoKinSpeciesIndex_[k + PLkstart] = kstart + k;
                KintoPLSpeciesIndex_[kstart + k] = k + PLkstart;
            }
        }



    } catch (CanteraError) {
        showErrors(cout);
        throw CanteraError("ReactingVolDomain::importFromXML", "error encountered");
        return false;
    }
    return true;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------



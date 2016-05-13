/*
 * $Id: PhaseList.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "PhaseList.h"
#include "Electrode_Exception.h"
#include "cantera/thermo.h"

#include <new>

using namespace std;

namespace Cantera
{

/**************************************************************************
 *
 * Constructors:
 *
 *   Allow for the initialization of one volumetric phase during
 *  the construction process. The addition of more volume phases
 *  and surface phases will be handled by subroutine calls.
 *
 *   vp = default is null, i.e., don't need a phase to initialize.
 *   vol = default is 0.0
 *   ownership = default is to not own the phases in the list.
 */
PhaseList::PhaseList(bool ownership) :
    m_NumTotPhases(0),
    m_NumTotSpecies(0),
    NumVolPhases_(0),
    m_totNumVolSpecies(0),
    VolPhaseList(0),
    VolPhaseXMLNodes(0),
    VolPhaseHasKinetics(0),
    m_NumSurPhases(0),
    m_NumEdgePhases(0),
    m_totNumSurSpecies(0),
    SurPhaseList(0),
    SurPhaseXMLNodes(0),
    SurPhaseHasKinetics(0),
    PhaseList_(0),
    PhaseNames_(0),
    m_numElements(0),
    m_PhaseSpeciesStartIndex(0),
    IOwnPhasePointers(ownership),
    m_GlobalElementObj(0)
{
    m_PhaseSpeciesStartIndex.resize(1, 0);
    m_GlobalElementObj = new Elements();
}

//================================================================================================
PhaseList::~PhaseList()
{
    if (IOwnPhasePointers) {
        size_t i;
        for (i = 0; i < NumVolPhases_; i++) {
            delete VolPhaseList[i];
        }
        for (i = 0; i < m_NumSurPhases; i++) {
            delete SurPhaseList[i];
        }
        for (i = 0; i < m_NumEdgePhases; i++) {
            delete EdgePhaseList[i];
        }
    }
    if (m_GlobalElementObj) {
        delete m_GlobalElementObj;
    }
}
//================================================================================================
PhaseList::PhaseList(const PhaseList& right) :
    m_NumTotPhases(0),
    m_NumTotSpecies(0),
    NumVolPhases_(0),
    m_totNumVolSpecies(0),
    VolPhaseList(0),
    VolPhaseXMLNodes(0),
    VolPhaseHasKinetics(0),
    m_NumSurPhases(0),
    m_NumEdgePhases(0),
    m_totNumSurSpecies(0),
    SurPhaseList(0),
    SurPhaseXMLNodes(0),
    SurPhaseHasKinetics(0),
    PhaseList_(0),
    PhaseNames_(0),
    m_numElements(0),
    IOwnPhasePointers(true),
    m_GlobalElementObj(0)
{
    /*
     * Call the assignment operator
     */
    PhaseList::operator=(right);
}
//================================================================================================
/*
 * operator=()
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working assignment operator
 */
PhaseList& PhaseList::operator=(const PhaseList& right)
{
    size_t i;
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    /*
     * First delete current contents of this object
     */

    if (IOwnPhasePointers) {
        size_t i;
        for (i = 0; i < NumVolPhases_; i++) {
            delete VolPhaseList[i];
            VolPhaseList[i] = 0;
        }
        for (i = 0; i < m_NumSurPhases; i++) {
            delete SurPhaseList[i];
            SurPhaseList[i] = 0;
        }
        for (i = 0; i < m_NumEdgePhases; i++) {
            delete EdgePhaseList[i];
            EdgePhaseList[i] = 0;
        }
    }
    if (m_GlobalElementObj) {
        delete m_GlobalElementObj;
        m_GlobalElementObj = 0;
    }

    /*
     * Next copy over the size data
     */
    m_NumTotPhases  =      right.m_NumTotPhases;
    m_NumTotSpecies =      right.m_NumTotSpecies;
    NumVolPhases_   =      right.NumVolPhases_;
    m_totNumVolSpecies =   right.m_totNumVolSpecies;
    IOwnPhasePointers  = right.IOwnPhasePointers;

    VolPhaseList = right.VolPhaseList;
    if (IOwnPhasePointers) {
        for (size_t i = 0; i < NumVolPhases_; i++) {
            VolPhaseList[i] = (right.VolPhaseList[i])->duplMyselfAsThermoPhase();
        }
    }

    //
    // These are shallow pointers, so they need to be reevaluated
    //
    VolPhaseXMLNodes=  right.VolPhaseXMLNodes;
    for (size_t i = 0; i < NumVolPhases_; i++) {
        VolPhaseXMLNodes[i] = &(VolPhaseList[i]->xml());
    }

    VolPhaseHasKinetics  = right.VolPhaseHasKinetics;

    m_NumSurPhases = right.m_NumSurPhases;
    m_NumEdgePhases = right.m_NumEdgePhases;
    m_totNumSurSpecies = right.m_totNumSurSpecies;


    SurPhaseList.resize(right.m_NumSurPhases, 0);
    if (IOwnPhasePointers) {
        for (i = 0; i < m_NumSurPhases; i++) {
            SurPhaseList[i] = (right.SurPhaseList[i])->duplMyselfAsThermoPhase();
        }
    } else {
        for (i = 0; i < m_NumSurPhases; i++) {
            SurPhaseList[i] = right.SurPhaseList[i];
        }
    }

    //
    // These are shallow pointers, so they need to be reevaluated
    //
    SurPhaseXMLNodes = right.SurPhaseXMLNodes;
    for (i = 0; i < m_NumSurPhases; i++) {
        SurPhaseXMLNodes[i] = &(SurPhaseList[i]->xml());
    }
    SurPhaseHasKinetics = right.SurPhaseHasKinetics;

    EdgePhaseList.resize(right.m_NumSurPhases, 0);
    if (IOwnPhasePointers) {
        for (i = 0; i < m_NumSurPhases; i++) {
            EdgePhaseList[i] = (right.EdgePhaseList[i])->duplMyselfAsThermoPhase();
        }
    } else {
        for (i = 0; i < m_NumEdgePhases; i++) {
            EdgePhaseList[i] = right.EdgePhaseList[i];
        }
    }

    EdgePhaseXMLNodes = right.EdgePhaseXMLNodes;
    for (i = 0; i < m_NumEdgePhases; i++) {
        EdgePhaseXMLNodes[i] = &(EdgePhaseList[i]->xml());
    }
    EdgePhaseHasKinetics = right.EdgePhaseHasKinetics;

    m_numElements = right.m_numElements;

    m_PhaseSpeciesStartIndex  = right.m_PhaseSpeciesStartIndex;

    IOwnPhasePointers = right.IOwnPhasePointers;

    if (m_GlobalElementObj) {
        delete m_GlobalElementObj;
    }
    m_GlobalElementObj = new Elements(*(right.m_GlobalElementObj));

    PhaseList_.resize(m_NumTotPhases);
    PhaseNames_.resize(m_NumTotPhases);
    for (size_t i = 0; i < NumVolPhases_; i++) {
        PhaseList_[i] =  VolPhaseList[i];
        PhaseNames_[i] =  VolPhaseList[i]->name();
    }
    for (i = 0; i < m_NumSurPhases; i++) {
        PhaseList_[i + NumVolPhases_] = SurPhaseList[i];
        PhaseNames_[i + NumVolPhases_] = SurPhaseList[i]->name();
    }
    size_t istart = NumVolPhases_ + m_NumSurPhases;
    for (i = 0; i < m_NumEdgePhases; i++) {
        PhaseList_[i + istart] = EdgePhaseList[i];
        PhaseNames_[i + istart] = EdgePhaseList[i]->name();
    }

    return *this;
}
//==================================================================================================================================
void PhaseList::
addVolPhase(std::string canteraFile, const std::string& phaseID)
{
    XML_Node* xroot = get_XML_File(canteraFile);
    XML_Node* vPhase = findXMLPhase(xroot, phaseID);
    Cantera::ThermoPhase *tp = Cantera::newPhase(canteraFile, phaseID);
    addVolPhase(tp, vPhase);
}
//===================================================================================================================================
void PhaseList::
addSurPhase(std::string canteraFile,  const std::string& phaseID)
{
    XML_Node* xroot = get_XML_File(canteraFile);
    XML_Node* vPhase = findXMLPhase(xroot, phaseID);
    Cantera::ThermoPhase *tp = Cantera::newPhase(canteraFile, phaseID);
    addSurPhase(tp, vPhase);
}
//==================================================================================================================================
/*
 *  addVolPhase:
 *
 *     This routine adds a volumetric phase to the list of phases
 *     in a PhaseList object. The phase must already have been
 *     specified. This routine then goes on to adding the
 *     additional information specified in the phase to the
 *     meta information lists kept by the PhaseList object
 */
void PhaseList::
addVolPhase(Cantera::ThermoPhase* const vp, Cantera::XML_Node* vPhase)
{
    AssertThrow(vp!=0, "Volume Phase Pointer must be nonzero");

    // Check for incompatibilities
    if (m_NumTotPhases > 0) {
        for (size_t k = 0; k < vp->nSpecies(); k++) {
            std::string sname = vp->speciesName(k);
            size_t gindex = globalSpeciesIndex(sname);
            if (gindex != npos) {
                throw CanteraError("addVolPhase()", "Species name, " + sname + " is a duplicated in different ThermoPhases\n");
            }
        }
        string tname = vp->name();
        for (size_t k = 0; k < m_NumTotPhases; k++) {
            string pname = PhaseList_[k]->name();
            if (tname == pname) {
                throw CanteraError("addVolPhase()", "phase name, " + tname + " is a duplicated for different ThermoPhases\n");
            }
        }
    }

    if (!vPhase) {
      vPhase = &(vp->xml());
    }

    if (vp->nDim() != 3) {
        throw CanteraError("PhaseList::addVolPhase()",
                           "number of dimensions isn't three");
    }

    NumVolPhases_++;
    VolPhaseList.resize(NumVolPhases_, 0);
    VolPhaseList[NumVolPhases_ - 1] = vp;

    size_t nSpecies = vp->nSpecies();
    m_NumTotPhases++;
    m_NumTotSpecies += nSpecies;
    m_totNumVolSpecies += nSpecies;

    /*
     * Store XML node
     */
    VolPhaseXMLNodes.resize(NumVolPhases_, 0);
    AssertThrow(vPhase != 0,"must be nonzero");
    VolPhaseXMLNodes[NumVolPhases_ - 1] = vPhase;

    /*
     * Push back indecises of surface phases
     */
    m_PhaseSpeciesStartIndex.resize(m_NumTotPhases + 1, 0);
    if (m_NumSurPhases > 0) {
        int indexP = NumVolPhases_ + m_NumSurPhases;
        for (size_t isp = 0; isp < m_NumSurPhases; isp++) {
            m_PhaseSpeciesStartIndex[indexP] = m_PhaseSpeciesStartIndex[indexP - 1] + nSpecies;
            indexP--;
        }
    }
    /*
     * Push back indecises of edge phases
     */
    if (m_NumEdgePhases > 0) {
        int indexP = NumVolPhases_ + m_NumSurPhases + m_NumEdgePhases;
        for (size_t iep = 0; iep < m_NumEdgePhases; iep++) {
            m_PhaseSpeciesStartIndex[indexP] = m_PhaseSpeciesStartIndex[indexP - 1] + nSpecies;
            indexP--;
        }
    }
    /*
     *  Add the new entry
     */
    m_PhaseSpeciesStartIndex[NumVolPhases_] = m_PhaseSpeciesStartIndex[NumVolPhases_ - 1] + nSpecies;

    /*
     * Check elements list -> enforce strict conformance.
     */
    int numE = vp->nElements();
    for (int e = 0; e < numE; e++) {
        string symb1 = vp->elementName(e);
        double weight1 = vp->atomicWeight(e);
        m_GlobalElementObj->addUniqueElement(symb1, weight1);
    }
    m_numElements = m_GlobalElementObj->nElements();

    /*
     * Check to see whether the phase has a kinetics object
     */
    VolPhaseHasKinetics.resize(NumVolPhases_, 0);
    if (vPhase->hasChild("kinetics")) {
        const XML_Node& kinNode = vPhase->child("kinetics");
        std::string smodel = kinNode["model"];
        if (smodel != "" && smodel != "none" && smodel != "None") {
            VolPhaseHasKinetics[NumVolPhases_-1] = 1;
        }
    }


    PhaseList_.resize(m_NumTotPhases);
    PhaseNames_.resize(m_NumTotPhases);

    for (size_t i = 0; i < NumVolPhases_; i++) {
        PhaseList_[i] = VolPhaseList[i];
        PhaseNames_[i] =  VolPhaseList[i]->name();
    }
    for (size_t i = 0; i < m_NumSurPhases; i++) {
        PhaseList_[i + NumVolPhases_] = SurPhaseList[i];
        PhaseNames_[i + NumVolPhases_] = SurPhaseList[i]->name();
    }
    size_t istart = m_NumSurPhases + NumVolPhases_;
    for (size_t i = 0; i < m_NumEdgePhases; i++) {
        PhaseList_[i + istart] = EdgePhaseList[i];
        PhaseNames_[i + istart] = EdgePhaseList[i]->name();
    }

    vp->realNumberRangeBehavior_ = DONOTHING_CTRB;
    vp->realNumberRangeBehavior_ = CHANGE_OVERFLOW_CTRB;

}
//==================================================================================================================================
/*
 *
 *  addSurfPhase:
 *
 *     This routine adds a surface phase to the list of phases
 *     in a PhaseList object. The phase must already have been
 *     specified. This routine then goes on to adding the
 *     additional information specified in the phase to the
 *     meta information lists kept by the PhaseList object
 */
void PhaseList::
addSurPhase(Cantera::ThermoPhase* const sp, Cantera::XML_Node* sPhase)
{

    AssertThrow(sp != 0, "must be nonzero");
    // Check for incompatibilities
    if (m_NumTotPhases > 0) {
        for (size_t k = 0; k < sp->nSpecies(); k++) {
            string sname = sp->speciesName(k);
            int gindex = globalSpeciesIndex(sname);
            if (gindex != -1) {
                throw CanteraError("addVolPhase()",
                                   "Species name, " + sname + " is a duplicated in different ThermoPhases\n");
            }
        }
        string tname = sp->name();
        for (size_t k = 0; k < m_NumTotPhases; k++) {
            string pname = PhaseList_[k]->name();
            if (tname == pname) {
                throw CanteraError("addVolPhase()",
                                   "phase name, " + tname + " is a duplicated for different ThermoPhases\n");
            }
        }
    }

    // Get the storred phase XML tree
    if (!sPhase) {
      sPhase = &(sp->xml());
    }

    if (sp->nDim() >=  3) {
        throw CanteraError("PhaseList::addSurPhase()",
                           "Number of dimensions is three or greater");
    }
    m_NumSurPhases++;
    SurPhaseList.resize(m_NumSurPhases, 0);
    SurPhaseList[m_NumSurPhases - 1] = sp;

    int nSpecies = sp->nSpecies();
    m_NumTotPhases++;
    m_NumTotSpecies += nSpecies;
    m_totNumSurSpecies += nSpecies;

    /*
     * Store XML node
     */
    SurPhaseXMLNodes.resize(m_NumSurPhases, 0);
    AssertThrow(sPhase != 0,"must be nonzero");
    SurPhaseXMLNodes[m_NumSurPhases - 1] = sPhase;

    /*
     * Push back indecises of surface phases
     */
    m_PhaseSpeciesStartIndex.resize(m_NumTotPhases + 1, 0);

    /*
     * Push back indecises of edge phases
     */
    if (m_NumEdgePhases > 0) {
        int indexP = NumVolPhases_ + m_NumSurPhases + m_NumEdgePhases;
        for (size_t iep = 0; iep < m_NumEdgePhases; iep++) {
            m_PhaseSpeciesStartIndex[indexP] = m_PhaseSpeciesStartIndex[indexP - 1] + nSpecies;
            indexP--;
        }
    }
    size_t istart = NumVolPhases_ + m_NumSurPhases;
    m_PhaseSpeciesStartIndex[istart] = m_PhaseSpeciesStartIndex[istart - 1] + nSpecies;

    /*
     * Check elements list -> enforce strict conformance.
     */
    int numE = sp->nElements();
    for (int e = 0; e < numE; e++) {
        string symb1 = sp->elementName(e);
        double weight1 = sp->atomicWeight(e);
        m_GlobalElementObj->addUniqueElement(symb1, weight1);
    }
    m_numElements = m_GlobalElementObj->nElements();

    /*
     * Check to see whether the phase has a kinetics object
     */
    SurPhaseHasKinetics.resize(m_NumSurPhases, 0);
    SurPhaseHasKinetics[m_NumSurPhases-1] = 0;
    if (sPhase->hasChild("kinetics")) {
        const XML_Node& kinNode = sPhase->child("kinetics");
        string smodel = kinNode["model"];
        if (smodel != "" && smodel != "none" && smodel != "None") {
            SurPhaseHasKinetics[m_NumSurPhases-1] = 1;
        }
    }

    PhaseList_.resize(m_NumTotPhases);
    PhaseNames_.resize(m_NumTotPhases);
    for (size_t i = 0; i < NumVolPhases_; i++) {
        PhaseList_[i] = VolPhaseList[i];
        PhaseNames_[i] =  VolPhaseList[i]->name();
    }
    for (size_t i = 0; i < m_NumSurPhases; i++) {
        PhaseList_[i + NumVolPhases_] = SurPhaseList[i];
        PhaseNames_[i + NumVolPhases_] = SurPhaseList[i]->name();
    }
    istart = NumVolPhases_ + m_NumSurPhases;
    for (size_t i = 0; i < m_NumEdgePhases; i++) {
        PhaseList_[i + istart] = SurPhaseList[i];
        PhaseNames_[i + istart] = SurPhaseList[i]->name();
    }
}
//==================================================================================================================================
/*
 *
 *  getVolPhaseIndex:
 *
 *     This routine returns the phase index of a phase. This
 *     number is the index value of the phase in the PhaseList
 *     object.
 *     NOTE: This function may have several different overloaded
 *           forms, depending on how the phase is to be looked
 *           up. The first form looks up the phase according to
 *           the pointer value of the object.
 */
size_t PhaseList::getVolPhaseIndex(const ThermoPhase* const vp) const
{
    for (size_t i = 0; i < NumVolPhases_; i++) {
        const ThermoPhase* const temp = VolPhaseList[i];
        if (temp == vp) {
            return i;
        }
    }
    return npos;
}
//==================================================================================================================================
size_t PhaseList::getSurPhaseIndex(const ThermoPhase* const sp) const
{
    for (size_t i = 0; i < m_NumSurPhases; i++) {
        const ThermoPhase* const temp = SurPhaseList[i];
        if (temp == sp) {
            return i;
        }
    }
    return npos;
}
//==================================================================================================================================
std::string PhaseList::phaseName(size_t globalPhaseIndex) const
{
    return PhaseNames_[globalPhaseIndex];
}
//==================================================================================================================================
size_t PhaseList::getGlobalPhaseIndex(const ThermoPhase* const tp) const
{
    size_t vi = getVolPhaseIndex(tp);
    if (vi != npos) {
        return vi;
    } else {
        size_t si = getSurPhaseIndex(tp);
        if (si != npos) {
            return NumVolPhases_ + si;
        }
    }
    return npos;
}
//==================================================================================================================================
size_t PhaseList::globalPhaseIndex(std::string phaseName, bool phaseIDAfter) const
{
    size_t ip;
    string pname;
    const ThermoPhase* tp_ptr;
    for (ip = 0; ip < m_NumTotPhases; ip++) {
        tp_ptr = PhaseList_[ip];
        pname = tp_ptr->name();
        if (phaseName == pname) {
            return ip;
        }
    }
    if (phaseIDAfter) {
        for (ip = 0; ip < m_NumTotPhases; ip++) {
            tp_ptr = PhaseList_[ip];
            pname = tp_ptr->id();
            if (phaseName == pname) {
                return ip;
            }
        }
    }
    return npos;
}
//==================================================================================================================================
size_t PhaseList::getGlobalSpeciesIndex(const ThermoPhase* const ttp, size_t k) const
{
    size_t iphase = getVolPhaseIndex(ttp);
    if (iphase != npos) {
        return getGlobalSpeciesIndex(iphase, k);
    }
    if (iphase == npos) {
        iphase = getSurPhaseIndex(ttp);
        if (iphase == npos) {
            throw CanteraError("PhaseList::getGlobalSpeciesIndex()", " error could not find phase index ");
        }
    }
    return getGlobalSpeciesIndex(iphase + NumVolPhases_, k);
}
//==================================================================================================================================
size_t PhaseList::globalSpeciesIndex(const std::string speciesName, const std::string phaseName) const
{
    size_t lindex = npos;
    if (phaseName == "") {
        for (size_t i = 0; i < m_NumTotPhases; i++) {
            ThermoPhase* tp = PhaseList_[i];
            lindex = tp->speciesIndex(speciesName);
            if (lindex != npos) {
                return ((m_PhaseSpeciesStartIndex[i]) + lindex);
            }
        }
    } else {
        size_t p = globalPhaseIndex(phaseName);
        if (p == npos) {
            //throw CanteraError("PhaseList::globalSpeciesIndex", "phase not found: " + phaseName);
            return npos;
        }
        size_t k = PhaseList_[p]->speciesIndex(speciesName);
        if (k == npos) {
            return npos;
            //throw CanteraError("PhaseList::globalSpeciesIndex", "species not found: "
            //                  + speciesName + " in phase " + phaseName);
        }
        return (int)( m_PhaseSpeciesStartIndex[p] + k);

    }
    return lindex;
}
//==================================================================================================================================
size_t PhaseList::getGlobalSpeciesIndex(size_t globPhaseIndex, size_t k) const
{
    AssertTrace(globPhaseIndex < m_NumTotPhases);
    int istart = m_PhaseSpeciesStartIndex[globPhaseIndex];
    return (istart + (int) k);
}
//==================================================================================================================================
size_t PhaseList::getGlobalSpeciesIndexVolPhaseIndex(size_t volPhaseIndex, size_t k) const
{
    AssertTrace(volPhaseIndex != npos);
    AssertTrace(k != npos);
    AssertTrace(volPhaseIndex < NumVolPhases_);
    const ThermoPhase* const tp = VolPhaseList[volPhaseIndex];
    AssertTrace(k < tp->nSpecies());
    size_t istart = m_PhaseSpeciesStartIndex[volPhaseIndex];
    return (istart + k);
}
//==================================================================================================================================
size_t PhaseList::getGlobalSpeciesIndexSurPhaseIndex(size_t surPhaseIndex, size_t k) const
{
    AssertTrace(surPhaseIndex != npos);
    AssertTrace(k != npos);
    AssertTrace(surPhaseIndex < m_NumSurPhases);
    const ThermoPhase* const tp = SurPhaseList[surPhaseIndex];
    AssertTrace(k < tp->nSpecies());
    size_t phaseIndex = NumVolPhases_ + surPhaseIndex;
    int istart = m_PhaseSpeciesStartIndex[phaseIndex];
    return (istart + (int) k);
}
//==================================================================================================================================
size_t PhaseList::getPhaseIndexFromGlobalSpeciesIndex(size_t globalSpeciesIndex) const
{
    AssertTrace(globalSpeciesIndex != npos);
    if (globalSpeciesIndex >= m_PhaseSpeciesStartIndex[m_NumTotPhases]) {
        throw CanteraError("PhaseList::getPhaseIndexFromGlobalSpeciesIndex ", " error");
    }
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        if (globalSpeciesIndex < m_PhaseSpeciesStartIndex[i+1]) {
            return (int) i;
        }
    }
    return npos;
}
//==================================================================================================================================
ThermoPhase* PhaseList::getPhase(const char* phaseName) const
{
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        //if phase names match, return this phase
        if (!PhaseNames_[i].compare(phaseName)) {
            return PhaseList_[i];
        }
    }
    // did not find matching phase
    throw CanteraError("PhaseList::getPhase()", "Phase name was not found\n");
    return 0;
}
//==================================================================================================================================
ThermoPhase* PhaseList::getPhase(const std::string& phaseName) const
{
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        //if phase names match, return this phase
        if (!PhaseNames_[i].compare(phaseName)) {
            return PhaseList_[i];
        }
    }
    throw CanteraError("PhaseList::getPhase()", "Phase name was not found\n");
    return 0;
}
//====================================================================================================================================
void PhaseList::
getLocalIndecisesFromGlobalSpeciesIndex(size_t globalSpeciesIndex, size_t& phaseIndex, size_t& localSpeciesIndex) const
{
    phaseIndex = getPhaseIndexFromGlobalSpeciesIndex(globalSpeciesIndex);
    localSpeciesIndex = (globalSpeciesIndex - m_PhaseSpeciesStartIndex[phaseIndex]);
}
//====================================================================================================================
//! Kernel function that checks consistency of the phase name and the the name and number/name/order of species in a phase 
/*!
 *   It returns true if the phases are the same and have the same species in it in the same order.
 *
 *   @param[in]      tpA                    Pointer to the first %ThermoPhase object
 *   @param[in]      tpB                    Pointer to the second ThermoPhase object
 *
 *   @return                                Returns true if the two ThermoPhase classes have the same eos and the same
 *                                          phase name and species names. Returns false otherwise.
 */
static bool ThermoPhasesTheSameNames(const ThermoPhase* const tpA, const ThermoPhase* const tpB)
{
    // Check the id() attribute of the phase to see that it is the same
    string sna = tpA->id();
    string snb = tpB->id();
    if (sna != snb) {
	return false;
    }
    // We also check the eosType() of the phase
    int ia = tpA->eosType();
    int ib = tpB->eosType();
    if (ia != ib) {
	return false;
    }
    size_t nspA = tpA->nSpecies();
    size_t nspB = tpB->nSpecies();
    if (nspA != nspB) {
        return false;
    }
    for (size_t k = 0; k <  nspA; ++k) {
	string spA = tpA->speciesName(k);
	string spB = tpB->speciesName(k);
	if (spA != spB) {
	    return false;
	}
    }
    return true;
}
//==================================================================================================================================
// Compare this Phaselist against another PhaseList, seeing if it contains the same lists
/*
 *        0  Phase lists are completely the same
 *        1  Phase lists are the same, but in a different phase order, and possibly different species
 *        2  Owning phase list is a superset 
 *        3  Guest phaseList is a superset
 *        4  apples and oranges.
 */
int PhaseList::compareOtherPL(const PhaseList* const plGuest) const
{
    int result = 0;
    //
    //  First compare for absolute agreement:
    //       We say that there is absolute agreement if 
    //          volume phases are same number and are the same ThermoPhase
    //          surface phases are same number and are the same ThermoPhase
    int compGood = 1;
    if (m_NumTotPhases == plGuest->m_NumTotPhases) {
	if (NumVolPhases_ == plGuest->NumVolPhases_) {
	    for (size_t i = 0; i < (size_t) NumVolPhases_; ++i) {
		const ThermoPhase* tpa = PhaseList_[i];
		const ThermoPhase* tpb = plGuest->PhaseList_[i];
		bool sameNames = ThermoPhasesTheSameNames(tpa, tpb);
		if (!sameNames) {
		    compGood = 0;
		    break;
		}
	    }
	} else {
	    compGood = 0;
	}
	if (m_NumSurPhases == plGuest->m_NumSurPhases) {
	    for (size_t i = 0; i < (size_t) m_NumSurPhases; ++i) {
		const ThermoPhase* tpa = SurPhaseList[i];
		const ThermoPhase* tpb = plGuest->SurPhaseList[i];
		bool sameNames = ThermoPhasesTheSameNames(tpa, tpb);
		if (!sameNames) {
		    compGood = 0;
		    break;
		}
	    }
	} else {
	    compGood = 0;
	}
    } else {
	compGood = 0;
    }
    if (compGood == 1) {
	return 0;
    } else {
	result = 4;
    }

    std::vector<size_t> volPhaseMatchForA(NumVolPhases_, npos);
    std::vector<size_t> volPhaseMatchForB(plGuest->NumVolPhases_, npos);
    //
    //  Loop through owning volume phases, looking for agreement. We find the best match if any
    //
    bool foundAllA = true;
    bool foundAllB = true;
    for (size_t i = 0; i < (size_t) NumVolPhases_; ++i) {
	bool iFound = false;
	int numBestAgreement = 0;
	int numCurrentAgreement = 0;
	const ThermoPhase* tpa = VolPhaseList[i];
	string ida = tpa->id();
	string pna = tpa->name();
	int iEOSa = tpa->eosType();
	const std::vector<std::string>&  sNa_list = tpa->speciesNames();
	for (size_t iGuest =  0; i < (size_t) plGuest->NumVolPhases_; ++i) {
	    const ThermoPhase* tpb = plGuest->VolPhaseList[i];
	    numCurrentAgreement = 0;
	    if (ida == tpb->id()) {
		numCurrentAgreement++;
	    }
	    if (pna == tpb->name()) {
		numCurrentAgreement++;
	    }
	    for (size_t kGuest =  0; kGuest < (size_t) tpb->nSpecies(); ++kGuest) {
		std::string spnb = tpb->speciesName(kGuest);
		std::vector<std::string>::const_iterator it = find(sNa_list.begin(), sNa_list.end(), spnb);
		if (it != sNa_list.end()) {
		    numCurrentAgreement++;
		}
	    }
	    if (numCurrentAgreement >= numBestAgreement && numCurrentAgreement > 0) {
	
		if (numCurrentAgreement == numBestAgreement) {
		    if (iEOSa == tpb->eosType()) {
			numCurrentAgreement++;
		    } else {
			continue;
		    }
		} else {
		    if (iEOSa == tpb->eosType()) {
			numCurrentAgreement++;
		    }
		}
		if (iFound) {
		    size_t iGuest_alt = volPhaseMatchForA[i];
		    volPhaseMatchForB[iGuest_alt] = npos;
		}
		iFound = true;
		numBestAgreement = numCurrentAgreement;
		volPhaseMatchForA[i] = iGuest;
		volPhaseMatchForB[iGuest] = i;
	    }
	}
	//  If we didn't find a matching A then indicate it
	if (!iFound) {
	    foundAllA = false;
	}
    }
    for (size_t iB = 0; iB <  volPhaseMatchForB.size(); ++iB) {
	if (volPhaseMatchForB[iB] == npos) {
	    foundAllB = false;
	}
    }

    if (foundAllA && foundAllB) {
	result = 1;
    } else {
	if (foundAllB) {
	    result = 2;
	} else if (foundAllA) {
	    result = 3;
	} else {
	    result = 4;
	}
    }

    int sresult = result;

    std::vector<size_t> surPhaseMatchForA(m_NumSurPhases, npos);
    std::vector<size_t> surPhaseMatchForB(plGuest->m_NumSurPhases, npos);
    //
    //  Loop through owning surface phases, looking for agreement. We find the best match if any
    //
    foundAllA = true;
    foundAllB = true;
    for (size_t i = 0; i < (size_t) m_NumSurPhases; ++i) {
	bool iFound = false;
	int numBestAgreement = 0;
	int numCurrentAgreement = 0;
	const ThermoPhase* tpa = SurPhaseList[i];
	string ida = tpa->id();
	string pna = tpa->name();
	int iEOSa = tpa->eosType();
	const std::vector<std::string>&  sNa_list = tpa->speciesNames();
	for (size_t iGuest =  0; i < (size_t) plGuest->m_NumSurPhases; ++i) {
	    const ThermoPhase* tpb = plGuest->SurPhaseList[i];
	    numCurrentAgreement = 0;
	    if (ida == tpb->id()) {
		numCurrentAgreement++;
	    }
	    if (pna == tpb->name()) {
		numCurrentAgreement++;
	    }
	    for (size_t kGuest =  0; kGuest < (size_t) tpb->nSpecies(); ++kGuest) {
		std::string spnb = tpb->speciesName(kGuest);
		std::vector<std::string>::const_iterator it = find(sNa_list.begin(), sNa_list.end(), spnb);
		if (it != sNa_list.end()) {
		    numCurrentAgreement++;
		}
	    }
	    if (numCurrentAgreement >= numBestAgreement && numCurrentAgreement > 0) {
	
		if (numCurrentAgreement == numBestAgreement) {
		    if (iEOSa == tpb->eosType()) {
			numCurrentAgreement++;
		    } else {
			continue;
		    }
		} else {
		    if (iEOSa == tpb->eosType()) {
			numCurrentAgreement++;
		    }
		}
		if (iFound) {
		    size_t iGuest_alt = volPhaseMatchForA[i];
		    volPhaseMatchForB[iGuest_alt] = npos;
		}
		iFound = true;
		numBestAgreement = numCurrentAgreement;
		surPhaseMatchForA[i] = iGuest;
		surPhaseMatchForB[iGuest] = i;
	    }
	}
	//  If we didn't find a matching A then indicate it
	if (!iFound) {
	    foundAllA = false;
	}
    }
    for (size_t iB = 0; iB <  surPhaseMatchForB.size(); ++iB) {
	if (surPhaseMatchForB[iB] == npos) {
	    foundAllB = false;
	}
    }

    if (foundAllA && foundAllB) {
	sresult = 1;
    } else {
	if (foundAllB) {
	    sresult = 2;
	} else if (foundAllA) {
	    sresult = 3;
	} else {
	    sresult = 4;
	}
    }

    if (result == 1) {
	if (sresult != 1) {
	    result = sresult;
	}
    } else if (result == 2) {
	if (sresult != 2 || sresult != 1) {
	    result = 4;
	}
    } else if (result == 3) {
	if (sresult != 3 || sresult != 1) {
	    result = 4;
	}
    }

    return result;
}
//==================================================================================================================================
ThermoPhase&
PhaseList::thermo(int globalPhaseIndex) const
{
    if (globalPhaseIndex < 0) {
        throw CanteraError("PhaseList::thermo()", "error");
    }
    return thermo((size_t) globalPhaseIndex);
}
//==================================================================================================================================
ThermoPhase&
PhaseList::thermo(size_t globalPhaseIndex) const
{
    if (globalPhaseIndex == npos) {
        throw CanteraError("PhaseList::thermo()", "error");
    }
    if (globalPhaseIndex < NumVolPhases_) {
        return *(VolPhaseList[globalPhaseIndex]);
    }
    if (globalPhaseIndex >= (NumVolPhases_ + m_NumSurPhases)) {
        throw CanteraError("PhaseList::thermo()", "error");
    }
    size_t isurphase = globalPhaseIndex - NumVolPhases_;
    return *(SurPhaseList[isurphase]);
}
//==================================================================================================================================
ThermoPhase&
PhaseList::phase(int globalPhaseIndex) const
{
    if (globalPhaseIndex < 0) {
        throw CanteraError("PhaseList::phase()", "error");
    }
    return phase((size_t) globalPhaseIndex);
}
//==================================================================================================================================
ThermoPhase&
PhaseList::phase(size_t globalPhaseIndex) const
{
    if (globalPhaseIndex == npos) {
        throw CanteraError("PhaseList::phase()", "error");
    }
    if (globalPhaseIndex < NumVolPhases_) {
        return *(VolPhaseList[globalPhaseIndex]);
    }
    if (globalPhaseIndex >= (NumVolPhases_ + m_NumSurPhases)) {
        throw CanteraError("PhaseList::phase()", "error");
    }
    size_t isurphase = globalPhaseIndex - NumVolPhases_;
    return *(SurPhaseList[isurphase]);
}
//==================================================================================================================================
ThermoPhase&
PhaseList::thermo(const char* phaseName) const
{
    std::string pp(phaseName);
    return thermo(pp);
}
//==================================================================================================================================
ThermoPhase&
PhaseList::thermo(const std::string& phaseName) const
{
    for (size_t ip = 0; ip <  m_NumTotPhases; ++ip) {
	if (phaseName == PhaseNames_[ip]) {
	    return *PhaseList_[ip];
	}
    }
    throw CanteraError("PhaseList::thermo()", "Can't find phase named, " + phaseName + ", in list of phases");
}
//==================================================================================================================================
const Elements* PhaseList ::getGlobalElements() const
{
    return m_GlobalElementObj;
}
//==================================================================================================================================
std::string
PhaseList::elementName(size_t e) const
{
    return m_GlobalElementObj->elementName(e);
}
//==================================================================================================================================
size_t
PhaseList::elementIndex(const std::string& elemName) const
{
    return m_GlobalElementObj->elementIndex(elemName);
}
//==================================================================================================================================
size_t PhaseList::nElements() const
{
    return m_numElements;
}
//==================================================================================================================================
XML_Node* PhaseList::surPhaseXMLNode(size_t iSurIndex) const
{
    return SurPhaseXMLNodes[iSurIndex];
}
//====================================================================================================================
XML_Node* PhaseList::volPhaseXMLNode(size_t iVolIndex) const
{
    return VolPhaseXMLNodes[iVolIndex];
}
//====================================================================================================================
ThermoPhase& PhaseList::volPhase(size_t iVolIndex)
{
    return *(VolPhaseList[iVolIndex]);
}
//====================================================================================================================
ThermoPhase& PhaseList::surPhase(size_t iSurIndex)
{
    return *(SurPhaseList[iSurIndex]);
}
//====================================================================================================================
size_t PhaseList::nSurPhases() const
{
    return m_NumSurPhases;
}
//====================================================================================================================
size_t PhaseList::nVolPhases() const
{
    return NumVolPhases_;
}
//====================================================================================================================
size_t PhaseList::nVolSpecies() const
{
    return m_totNumVolSpecies;
}
//====================================================================================================================
size_t PhaseList::nSpecies() const
{
    return m_NumTotSpecies;
}
//====================================================================================================================
size_t PhaseList::nPhases() const
{
    return m_NumTotPhases;
}
//====================================================================================================================
bool PhaseList::volPhaseHasKinetics(size_t iVolIndex) const
{
    return VolPhaseHasKinetics[iVolIndex];
}
//====================================================================================================================
bool PhaseList::surPhaseHasKinetics(size_t iSurIndex) const
{
    return SurPhaseHasKinetics[iSurIndex];
}
//====================================================================================================================
std::string PhaseList::speciesName(size_t iGlobSpeciesIndex) const
{
    int iPhase = getPhaseIndexFromGlobalSpeciesIndex(iGlobSpeciesIndex);
    int kStart = (int) m_PhaseSpeciesStartIndex[iPhase];
    int kLocal = iGlobSpeciesIndex - kStart;
    ThermoPhase& tp = thermo(iPhase);
    return tp.speciesName(kLocal);
}
//======================================================================================================================
void PhaseList::setState_TP(doublereal temperature, doublereal pressure)
{
    for (size_t i = 0; i <  NumVolPhases_; i++) {
        ThermoPhase *tp = VolPhaseList[i];
        tp->setState_TP(temperature, pressure);
    } 
    for (size_t i = 0; i < m_NumSurPhases; i++) {
        ThermoPhase *tp = SurPhaseList[i];
        tp->setState_TP(temperature, pressure);
    }
}
//======================================================================================================================
}
//======================================================================================================================

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
PhaseList::PhaseList(ThermoPhase* const vp, XML_Node* vPhase, double vol,
                     bool ownership) :
    m_NumTotPhases(0),
    m_NumTotSpecies(0),
    NumVolPhases_(0),
    m_totNumVolSpecies(0),
    VolPhaseList(0),
    VolPhaseXMLNodes(0),
    VolPhaseHasKinetics(0),
    m_NumSurPhases(0),
    m_totNumSurSpecies(0),
    SurPhaseList(0),
    SurPhaseXMLNodes(0),
    SurPhaseHasKinetics(0),
    PhaseList_(0),
    PhaseNames_(0),
    m_numElements(0),
    m_PhaseSpeciesStartIndex(0),
    CanteraFNSurface(),
    IOwnPhasePointers(ownership),
    m_GlobalElementObj(0)
{
    m_PhaseSpeciesStartIndex.resize(1, 0);

    if (vp) {
        AssertThrow(vPhase->name() == "phase", "PhaseList constructor with bad xml node");
        string dimS = vPhase->operator[]("dim");
        if (dimS == "3") {
            addVolPhase(vp, vPhase);
        } else if (dimS == "2") {
            addSurPhase(vp, vPhase);
        } else {
            throw CanteraError("processPhasePL", "unknown dim string: " + dimS);
        }
    }

    m_GlobalElementObj = new Elements();
}

//================================================================================================
PhaseList::~PhaseList()
{
    if (IOwnPhasePointers) {
        int i;
        for (i = 0; i < NumVolPhases_; i++) {
            delete VolPhaseList[i];
        }
        for (i = 0; i < m_NumSurPhases; i++) {
            delete SurPhaseList[i];
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
    operator=(right);
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
    int i;
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
        int i;
        for (i = 0; i < NumVolPhases_; i++) {
            delete VolPhaseList[i];
            VolPhaseList[i] = 0;
        }
        for (i = 0; i < m_NumSurPhases; i++) {
            delete SurPhaseList[i];
            SurPhaseList[i] = 0;
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
        for (i = 0; i < NumVolPhases_; i++) {
            VolPhaseList[i] = (right.VolPhaseList[i])->duplMyselfAsThermoPhase();
        }
    }

    VolPhaseXMLNodes=  right.VolPhaseXMLNodes;

    VolPhaseHasKinetics  = right.VolPhaseHasKinetics;

    m_NumSurPhases = right.m_NumSurPhases;
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

    SurPhaseXMLNodes = right.SurPhaseXMLNodes;
    SurPhaseHasKinetics = right.SurPhaseHasKinetics;

    m_numElements = right.m_numElements;

    m_PhaseSpeciesStartIndex  = right.m_PhaseSpeciesStartIndex;

    CanteraFNSurface = right.CanteraFNSurface;

    IOwnPhasePointers = right.IOwnPhasePointers;

    if (m_GlobalElementObj) {
        delete m_GlobalElementObj;
    }
    m_GlobalElementObj = new Elements(*(right.m_GlobalElementObj));

    PhaseList_.resize(m_NumTotPhases);
    PhaseNames_.resize(m_NumTotPhases);
    for (i = 0; i < NumVolPhases_; i++) {
        PhaseList_[i] =  VolPhaseList[i];
        PhaseNames_[i] =  VolPhaseList[i]->name();
    }
    for (i = 0; i < m_NumSurPhases; i++) {
        PhaseList_[i + NumVolPhases_] = SurPhaseList[i];
        PhaseNames_[i + NumVolPhases_] = SurPhaseList[i]->name();
    }

    return *this;
}
//================================================================================================
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
addVolPhase(Cantera::ThermoPhase* const vp, XML_Node* vPhase, std::string canteraFile)
{

    // Check for incompatibilities
    if (m_NumTotPhases > 0) {
        for (size_t k = 0; k < vp->nSpecies(); k++) {
            string sname = vp->speciesName(k);
            int gindex = globalSpeciesIndex(sname);
            if (gindex != -1) {
                throw CanteraError("addVolPhase()",
                                   "Species name, " + sname + " is a duplicated in different ThermoPhases\n");
            }
        }
        string tname = vp->name();
        for (int k = 0; k < m_NumTotPhases; k++) {
            string pname = PhaseList_[k]->name();
            if (tname == pname) {
                throw CanteraError("addVolPhase()",
                                   "phase name, " + tname + " is a duplicated for different ThermoPhases\n");
            }
        }
    }


    NumVolPhases_++;
    VolPhaseList.resize(NumVolPhases_, 0);
    VolPhaseList[NumVolPhases_ - 1] = vp;

    int nSpecies = vp->nSpecies();
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
        for (int isp = 0; isp < m_NumSurPhases; isp++) {
            m_PhaseSpeciesStartIndex[indexP] =
                m_PhaseSpeciesStartIndex[indexP - 1] + nSpecies;
            indexP--;
        }
    }
    m_PhaseSpeciesStartIndex[NumVolPhases_] =
        m_PhaseSpeciesStartIndex[NumVolPhases_ - 1] + nSpecies;

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
        string smodel = kinNode["model"];
        if (smodel != "" && smodel != "none" && smodel != "None") {
            VolPhaseHasKinetics[NumVolPhases_-1] = 1;
        }
    }


    PhaseList_.resize(m_NumTotPhases);
    PhaseNames_.resize(m_NumTotPhases);

    for (int i = 0; i < NumVolPhases_; i++) {
        PhaseList_[i] = VolPhaseList[i];
        PhaseNames_[i] =  VolPhaseList[i]->name();
    }
    for (int i = 0; i < m_NumSurPhases; i++) {
        PhaseList_[i + NumVolPhases_] = SurPhaseList[i];
        PhaseNames_[i + NumVolPhases_] = SurPhaseList[i]->name();
    }

}
//==================================================================================================================
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
addSurPhase(Cantera::ThermoPhase* const sp, Cantera::XML_Node* sPhase, std::string canteraFile)
{

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
        for (int k = 0; k < m_NumTotPhases; k++) {
            string pname = PhaseList_[k]->name();
            if (tname == pname) {
                throw CanteraError("addVolPhase()",
                                   "phase name, " + tname + " is a duplicated for different ThermoPhases\n");
            }
        }
    }

    if (m_NumSurPhases == 0) {
        CanteraFNSurface = canteraFile;
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

    m_PhaseSpeciesStartIndex[m_NumTotPhases] =
        m_PhaseSpeciesStartIndex[m_NumTotPhases - 1] + nSpecies;

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
    for (int i = 0; i < NumVolPhases_; i++) {
        PhaseList_[i] = VolPhaseList[i];
        PhaseNames_[i] =  VolPhaseList[i]->name();
    }
    for (int i = 0; i < m_NumSurPhases; i++) {
        PhaseList_[i + NumVolPhases_] = SurPhaseList[i];
        PhaseNames_[i + NumVolPhases_] = SurPhaseList[i]->name();
    }
}
//==================================================================================================================
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
int PhaseList::getVolPhaseIndex(const ThermoPhase* const vp) const
{
    for (int i = 0; i < NumVolPhases_; i++) {
        const ThermoPhase* const temp = VolPhaseList[i];
        if (temp == vp) {
            return i;
        }
    }
    return -1;
}
//==================================================================================================================
int PhaseList::getSurPhaseIndex(const ThermoPhase* const sp) const
{
    for (int i = 0; i < m_NumSurPhases; i++) {
        const ThermoPhase* const temp = SurPhaseList[i];
        if (temp == sp) {
            return i;
        }
    }
    return -1;
}
//==================================================================================================================
std::string PhaseList::phaseName(int globalPhaseIndex) const
{
    return PhaseNames_[globalPhaseIndex];
}
//===================================================================================================================
int PhaseList::getGlobalPhaseIndex(const ThermoPhase* const tp) const
{
    int vi = getVolPhaseIndex(tp);
    if (vi != -1) {
        return vi;
    } else {
        int si = getSurPhaseIndex(tp);
        if (si != -1) {
            return NumVolPhases_ + si;
        }
    }
    return -1;
}
//================================================================================================
int  PhaseList::globalPhaseIndex(std::string phaseName, bool phaseIDAfter) const
{
    int ip;
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
    return -1;
}
//================================================================================================
int PhaseList::
getGlobalSpeciesIndex(const ThermoPhase* const ttp, int k) const
{
    int iphase = getVolPhaseIndex(ttp);
    if (iphase != -1) {
        return getGlobalSpeciesIndex(iphase, k);
    }
    if (iphase == -1) {
        iphase = getSurPhaseIndex(ttp);
        if (iphase == -1) {
            throw CanteraError(" ", " error ");
        }
    }
    return getGlobalSpeciesIndex(iphase + NumVolPhases_, k);
}
//================================================================================================
int PhaseList::globalSpeciesIndex(const std::string speciesName, const std::string phaseName) const
{
    int lindex = -1;
    if (phaseName == "") {
        for (int i = 0; i < m_NumTotPhases; i++) {
            ThermoPhase* tp = PhaseList_[i];
            lindex = tp->speciesIndex(speciesName);
            if (lindex >= 0) {
                return ((m_PhaseSpeciesStartIndex[i]) + lindex);
            }
        }
    } else {
        int p = globalPhaseIndex(phaseName);
        if (p < 0) {
            return -2;
            //throw CanteraError("PhaseList::globalSpeciesIndex", "phase not found: " + phaseName);
        }
        int k = PhaseList_[p]->speciesIndex(speciesName);
        if (k < 0) {
            return k;
            //throw CanteraError("PhaseList::globalSpeciesIndex", "species not found: "
            //                  + speciesName + " in phase " + phaseName);
        }
        return m_PhaseSpeciesStartIndex[p] + k;

    }
    return lindex;
}
//================================================================================================
int PhaseList::
getGlobalSpeciesIndex(int globPhaseIndex, int k) const
{
    AssertTrace(globPhaseIndex < m_NumTotPhases);
    int istart = m_PhaseSpeciesStartIndex[globPhaseIndex];
    return (istart + k);
}
//================================================================================================
int PhaseList::
getGlobalSpeciesIndexVolPhaseIndex(int volPhaseIndex, int k) const
{
    AssertTrace(volPhaseIndex >= 0);
    AssertTrace(volPhaseIndex < NumVolPhases_);
    ThermoPhase* tp = VolPhaseList[volPhaseIndex];
    AssertTrace(k < (int) tp->nSpecies());
    int istart = m_PhaseSpeciesStartIndex[volPhaseIndex];
    return (istart + k);
}
//================================================================================================
int PhaseList::
getGlobalSpeciesIndexSurPhaseIndex(int surPhaseIndex, int k) const
{
    AssertTrace(surPhaseIndex >= 0);
    AssertTrace(surPhaseIndex < m_NumSurPhases);
    ThermoPhase* tp = SurPhaseList[surPhaseIndex];
    AssertTrace(k < (int) tp->nSpecies());
    int phaseIndex = NumVolPhases_ + surPhaseIndex;
    int istart = m_PhaseSpeciesStartIndex[phaseIndex];
    return (istart + k);
}
//================================================================================================
int PhaseList::
getPhaseIndexFromGlobalSpeciesIndex(int globalSpeciesIndex) const
{
    int i;
    if (globalSpeciesIndex < 0 ||
            globalSpeciesIndex >= m_PhaseSpeciesStartIndex[m_NumTotPhases]) {
        throw CanteraError(" ", " error");
    }
    for (i = 0; i < m_NumTotPhases; i++) {
        if (globalSpeciesIndex < m_PhaseSpeciesStartIndex[i+1]) {
            return i;
        }
    }
    return -1;
}
//================================================================================================
ThermoPhase* PhaseList::getPhase(const char* phaseName) const
{
    for (int i = 0; i < m_NumTotPhases; i++) {
        //if phase names match, return this phase
        if (!PhaseNames_[i].compare(phaseName)) {
            return PhaseList_[i];
        }
    }
    // did not find matching phase
    throw CanteraError("PhaseList::getPhase()",
                       "Phase name was not found\n");
    return 0;

}
//====================================================================================================================
void PhaseList::
getLocalIndecisesFromGlobalSpeciesIndex(int globalSpeciesIndex,
                                        int& phaseIndex,
                                        int& localSpeciesIndex) const
{
    phaseIndex =
        getPhaseIndexFromGlobalSpeciesIndex(globalSpeciesIndex);
    localSpeciesIndex = globalSpeciesIndex
                        - m_PhaseSpeciesStartIndex[phaseIndex];
}
//===================================================================================================================
ThermoPhase&
PhaseList::thermo(int globalPhaseIndex) const
{
    if (globalPhaseIndex < 0) {
        throw CanteraError("  ", "error");
    }
    if (globalPhaseIndex < NumVolPhases_) {
        return *(VolPhaseList[globalPhaseIndex]);
    }
    if (globalPhaseIndex >= (NumVolPhases_ + m_NumSurPhases)) {
        throw CanteraError("  ", "error");
    }
    int isurphase = globalPhaseIndex - NumVolPhases_;
    return *(SurPhaseList[isurphase]);
}
//===================================================================================================================
const Elements* PhaseList ::getGlobalElements() const
{
    return m_GlobalElementObj;
}
//====================================================================================================================
std::string
PhaseList::elementName(int e) const
{
    return m_GlobalElementObj->elementName(e);
}
//====================================================================================================================
//  Return the file name of the file containing the first surface phase
//  encountered when initializing this object
std::string  PhaseList::firstSurfaceFile() const
{
    return CanteraFNSurface;
}
//====================================================================================================================
int PhaseList::nElements() const
{
    return m_numElements;
}
//====================================================================================================================
XML_Node* PhaseList::surPhaseXMLNode(int iSurIndex) const
{
    return SurPhaseXMLNodes[iSurIndex];
}
//====================================================================================================================
XML_Node* PhaseList::volPhaseXMLNode(int iVolIndex) const
{
    return VolPhaseXMLNodes[iVolIndex];
}
//====================================================================================================================
ThermoPhase& PhaseList::volPhase(int iVolIndex)
{
    return *(VolPhaseList[iVolIndex]);
}
//====================================================================================================================
ThermoPhase& PhaseList::surPhase(int iSurIndex)
{
    return *(SurPhaseList[iSurIndex]);
}
//====================================================================================================================
int PhaseList::nSurPhases() const
{
    return m_NumSurPhases;
}
//====================================================================================================================
int PhaseList::nVolPhases() const
{
    return NumVolPhases_;
}
//====================================================================================================================
int PhaseList::nVolSpecies() const
{
    return m_totNumVolSpecies;
}
//====================================================================================================================
int PhaseList::nSpecies() const
{
    return m_NumTotSpecies;
}
//====================================================================================================================
int PhaseList::nPhases() const
{
    return m_NumTotPhases;
}
//====================================================================================================================
bool PhaseList::volPhaseHasKinetics(int iVolIndex) const
{
    return VolPhaseHasKinetics[iVolIndex];
}
//====================================================================================================================
bool PhaseList::surPhaseHasKinetics(int iSurIndex) const
{
    return SurPhaseHasKinetics[iSurIndex];
}
//====================================================================================================================
std::string PhaseList::speciesName(int iGlobSpeciesIndex) const
{
    int iPhase = getPhaseIndexFromGlobalSpeciesIndex(iGlobSpeciesIndex);
    int kStart = m_PhaseSpeciesStartIndex[iPhase];
    int kLocal = iGlobSpeciesIndex - kStart;
    ThermoPhase& tp = thermo(iPhase);
    return tp.speciesName(kLocal);
}
//======================================================================================================================
}
//======================================================================================================================

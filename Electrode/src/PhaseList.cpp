/**
 * @file PhaseList.h 
 *    Declarations for the  base class within Zuzax that handles indexing
 *    within multiphase applications (see \ref ExtendedPhaseGroups and class \link Zuzax::PhaseList\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "PhaseList.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
/*
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
    m_OrderedByDimensionality(true),
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
    m_totNumEdgeSpecies(0),
    SurPhaseList(0),
    SurPhaseXMLNodes(0),
    SurPhaseHasKinetics(0),
    PhaseList_(0),
    PhaseNames_(0),
    m_numElements(0),
    m_PhaseSpeciesStartIndex(0),
    IOwnPhasePointers(ownership),
    m_lastPhaseIndexAdded(npos),
    m_nSpLastPhaseAdded(npos),
    m_dimLastPhaseAdded(npos),
    m_lastPhaseIndexDeleted(npos),
    m_nSpLastPhaseDeleted(npos),
    m_dimLastPhaseDeleted(npos)

{
    m_PhaseSpeciesStartIndex.resize(1, 0);
}
//==================================================================================================================================
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
}
//==================================================================================================================================
PhaseList::PhaseList(const PhaseList& right) :
    m_OrderedByDimensionality(true),
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
    m_totNumEdgeSpecies(0),
    SurPhaseList(0),
    SurPhaseXMLNodes(0),
    SurPhaseHasKinetics(0),
    PhaseList_(0),
    PhaseNames_(0),
    m_numElements(0),
    IOwnPhasePointers(true)
{
    /*
     * Call the assignment operator
     */
    PhaseList::operator=(right);
}
//==================================================================================================================================
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
        for (i = 0; i < PhaseList_.size(); i++) {
            PhaseList_[i] = 0;
        }
    }

    /*
     * Next copy over the size data
     */
    m_OrderedByDimensionality = right.m_OrderedByDimensionality;
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
    VolPhaseToGlobPhaseIndex_ = right.VolPhaseToGlobPhaseIndex_;

    m_NumSurPhases = right.m_NumSurPhases;
    m_NumEdgePhases = right.m_NumEdgePhases;
    m_totNumSurSpecies = right.m_totNumSurSpecies;
    m_totNumEdgeSpecies = right.m_totNumEdgeSpecies;

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
    SurPhaseToGlobPhaseIndex_ = right.SurPhaseToGlobPhaseIndex_;

    EdgePhaseList.resize(right.m_NumSurPhases, 0);
    if (IOwnPhasePointers) {
        for (i = 0; i < m_NumEdgePhases; i++) {
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
    EdgePhaseToGlobPhaseIndex_ = right.EdgePhaseToGlobPhaseIndex_;

    m_numElements = right.m_numElements;

    m_PhaseSpeciesStartIndex  = right.m_PhaseSpeciesStartIndex;

    IOwnPhasePointers = right.IOwnPhasePointers;

    m_GlobalElementObj = right.m_GlobalElementObj;
    m_atoms = right.m_atoms;
    m_globalToLocalEMap = right.m_globalToLocalEMap;
    m_localToGlobalEMap = right.m_localToGlobalEMap;

    PhaseList_.resize(m_NumTotPhases);
    PhaseNames_.resize(m_NumTotPhases);
    for (size_t i = 0; i < NumVolPhases_; i++) {
        PhaseList_[i] =  VolPhaseList[i];
        PhaseNames_[i] = VolPhaseList[i]->name();
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
    m_lastPhaseIndexAdded = right.m_lastPhaseIndexAdded;
    m_nSpLastPhaseAdded = right.m_nSpLastPhaseAdded;
    m_dimLastPhaseAdded = right.m_dimLastPhaseAdded;
    m_lastPhaseIndexDeleted = right.m_lastPhaseIndexDeleted;
    m_nSpLastPhaseDeleted = right.m_nSpLastPhaseDeleted;
    m_dimLastPhaseDeleted = right.m_dimLastPhaseDeleted;

    return *this;
}
//==================================================================================================================================
bool PhaseList::orderedByDims() const
{
    return m_OrderedByDimensionality;
}
//==================================================================================================================================
ThermoPhase* PhaseList::addVolPhase(const std::string& canteraFile, const std::string& phaseID, bool orderByDims)
{
    XML_Node* xroot = get_XML_File(canteraFile);
    if (!xroot) {
        return 0;
    }
    XML_Node* vPhase = findXMLPhase(xroot, phaseID);
    if (!vPhase) {
       return 0;
    }
    ThermoPhase *tp = newPhase(canteraFile, phaseID);
    if (tp) {
       addVolPhase(tp, vPhase, orderByDims);
    }
    return tp;
}
//===================================================================================================================================
ThermoPhase* PhaseList::addSurPhase(const std::string& canteraFile, const std::string& phaseID)
{
    XML_Node* xroot = get_XML_File(canteraFile);
    if (!xroot) {
        return 0;
    }
    XML_Node* vPhase = findXMLPhase(xroot, phaseID);
    if (!vPhase) {
       return 0;
    }
    ThermoPhase *tp = newPhase(canteraFile, phaseID);
    if (tp) {
        addSurPhase(tp, vPhase);
    }
    return tp;
}
//===================================================================================================================================
ThermoPhase* PhaseList::addEdgePhase(const std::string& canteraFile, const std::string& phaseID)
{
    XML_Node* xroot = get_XML_File(canteraFile);
    if (!xroot) {
        return 0;
    }
    XML_Node* vPhase = findXMLPhase(xroot, phaseID);
    if (!vPhase) {
       return 0;
    }
    ThermoPhase *tp = newPhase(canteraFile, phaseID);
    if (tp) {
        addEdgePhase(tp, vPhase);
    }
    return tp;
}
//==================================================================================================================================
ThermoPhase* PhaseList::addPhase(const std::string& canteraFile, const std::string& phaseID)
{
    XML_Node* xroot = get_XML_File(canteraFile);
    if (!xroot) {
        return 0;
    }
    XML_Node* vPhase = findXMLPhase(xroot, phaseID);
    if (!vPhase) {
       return 0;
    }
    ThermoPhase *tp = newPhase(canteraFile, phaseID);
    if (tp) {
        addPhase(tp, vPhase);
    }
    return tp;
}
//==================================================================================================================================
void PhaseList::addPhase(ThermoPhase* const vp, XML_Node* vPhase)
{
#ifdef DEBUG_MODE
    if (!vp) {
        throw CanteraError("PhaseList::addPhase()", "zero pointer for ThermoPhase");
    }
#endif
    if (vp->nDim() == 3) {
        addVolPhase(vp, vPhase);
    } else if (vp->nDim() == 2) {
        addSurPhase(vp, vPhase);
    } else if (vp->nDim() == 1) {
         addEdgePhase(vp, vPhase);
    } else {
          throw CanteraError("PhaseList::addPhase()", "unknown dimension");
    }
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
void PhaseList::addVolPhase(ThermoPhase* const vp, XML_Node* vPhase, bool orderedByDims)
{
    AssertThrow(vp!=0, "Volume Phase Pointer must be nonzero");

    // Check for incompatibilities
    if (m_NumTotPhases > 0) {
        std::string tname = vp->name();
        for (size_t k = 0; k < m_NumTotPhases; k++) {
            std::string pname = PhaseList_[k]->name();
            if (tname == pname) {
                throw CanteraError("PhaseList::addVolPhase()", "phase name, " + tname + " is a duplicated for different ThermoPhases\n");
            }
        }
    }
    if (!vPhase) {
      vPhase = &(vp->xml());
    }
    if (vp->nDim() != 3) {
        throw CanteraError("PhaseList::addVolPhase()", "number of dimensions isn't three");
    }

    m_lastPhaseIndexAdded = NumVolPhases_;
    m_dimLastPhaseAdded = 3; 
    VolPhaseToGlobPhaseIndex_.push_back(NumVolPhases_);
    NumVolPhases_++;
    VolPhaseList.resize(NumVolPhases_, 0);
    VolPhaseList[NumVolPhases_ - 1] = vp;

    size_t nSpecies = vp->nSpecies();
    m_NumTotPhases++;
    m_NumTotSpecies += nSpecies;
    m_totNumVolSpecies += nSpecies;
    m_nSpLastPhaseAdded = nSpecies;

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
        size_t indexP = NumVolPhases_ + m_NumSurPhases;
        for (size_t isp = 0; isp < m_NumSurPhases; isp++) {
            m_PhaseSpeciesStartIndex[indexP] = m_PhaseSpeciesStartIndex[indexP - 1] + nSpecies;
            indexP--;
        }
    }
    /*
     * Push back indecises of edge phases
     */
    if (m_NumEdgePhases > 0) {
        size_t indexP = NumVolPhases_ + m_NumSurPhases + m_NumEdgePhases;
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
    size_t numE = vp->nElements();
    for (size_t e = 0; e < numE; e++) {
        std::string symb1 = vp->elementName(e);
        doublevalue weight1 = vp->atomicWeight(e);
        int an = vp->atomicNumber(e);
        doublevalue entropy298 = vp->entropyElement298(e, true);
        int eType = vp->elementType(e);
        m_GlobalElementObj.addUniqueElement(symb1, weight1, an, entropy298, eType);
    }
    m_numElements = m_GlobalElementObj.nElements();

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

    calcElementMaps(true); 

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
void PhaseList::addSurPhase(ThermoPhase* const sp, XML_Node* sPhase)
{
    AssertThrow(sp != 0, "must be nonzero");
    // Check for incompatibilities
    if (m_NumTotPhases > 0) {
        std::string tname = sp->name();
        for (size_t k = 0; k < m_NumTotPhases; k++) {
            std::string pname = PhaseList_[k]->name();
            if (tname == pname) {
                throw CanteraError("PhaseList::addSurPhase()",
                                   "phase name, " + tname + " is a duplicated for different ThermoPhases\n");
            }
        }
    }
    // Get the storred phase XML tree
    if (!sPhase) {
      sPhase = &(sp->xml());
    }
    if (sp->nDim() >=  3) {
        throw CanteraError("PhaseList::addSurPhase()", "Number of dimensions is three or greater");
    }
    if (sp->nDim() <= 1) {
        throw CanteraError("PhaseList::addSurPhase()", "Number of dimensions is one or less");
    }
    m_lastPhaseIndexAdded = NumVolPhases_ + m_NumSurPhases;
    m_dimLastPhaseAdded = 2;
    m_NumSurPhases++;
    SurPhaseList.resize(m_NumSurPhases, 0);
    SurPhaseList[m_NumSurPhases - 1] = sp;
    SurPhaseToGlobPhaseIndex_.push_back(m_lastPhaseIndexAdded);

    size_t nSpecies = sp->nSpecies();
    m_NumTotPhases++;
    m_NumTotSpecies += nSpecies;
    m_totNumSurSpecies += nSpecies;
    m_nSpLastPhaseAdded = nSpecies;

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
        size_t indexP = NumVolPhases_ + m_NumSurPhases + m_NumEdgePhases;
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
    size_t numE = sp->nElements();
    for (size_t e = 0; e < numE; e++) {
        std::string symb1 = sp->elementName(e);
        doublevalue weight1 = sp->atomicWeight(e);
        int an = sp->atomicNumber(e);
        doublevalue entropy298 = sp->entropyElement298(e, true);
        int eType = sp->elementType(e);
        m_GlobalElementObj.addUniqueElement(symb1, weight1, an, entropy298, eType);
    }
    m_numElements = m_GlobalElementObj.nElements();

    /*
     * Check to see whether the phase has a kinetics object
     */
    SurPhaseHasKinetics.resize(m_NumSurPhases, 0);
    SurPhaseHasKinetics[m_NumSurPhases-1] = 0;
    if (sPhase->hasChild("kinetics")) {
        const XML_Node& kinNode = sPhase->child("kinetics");
        std::string smodel = kinNode["model"];
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

    calcElementMaps(true); 
}
//==================================================================================================================================
void PhaseList::addEdgePhase(ThermoPhase* const sp, XML_Node* ePhase)
{
    AssertThrow(sp != 0, "must be nonzero");
    // Check for incompatibilities
    if (m_NumTotPhases > 0) {
        std::string tname = sp->name();
        for (size_t k = 0; k < m_NumTotPhases; k++) {
            std::string pname = PhaseList_[k]->name();
            if (tname == pname) {
                throw CanteraError("PhaseList::addEdgePhase()",
                                   "phase name, " + tname + " is a duplicated for different ThermoPhases\n");
            }
        }
    }
    // Get the storred phase XML tree
    if (!ePhase) {
      ePhase = &(sp->xml());
    }
    if (sp->nDim() >=  2) {
        throw CanteraError("PhaseList::addEdgePhase()", "Number of dimensions is two or greater");
    }
  
    m_lastPhaseIndexAdded = NumVolPhases_ + m_NumSurPhases + m_NumEdgePhases;
    m_dimLastPhaseAdded = 1;
    m_NumEdgePhases++;
    EdgePhaseList.resize(m_NumEdgePhases, 0);
    EdgePhaseList[m_NumEdgePhases - 1] = sp;
    EdgePhaseToGlobPhaseIndex_.push_back(m_lastPhaseIndexAdded); 

    size_t nSpecies = sp->nSpecies();
    m_NumTotPhases++;
    m_NumTotSpecies += nSpecies;
    m_totNumEdgeSpecies += nSpecies;
    m_nSpLastPhaseAdded = nSpecies;

    /*
     * Store XML node
     */
    EdgePhaseXMLNodes.resize(m_NumEdgePhases, 0);
    AssertThrow(ePhase != 0,"must be nonzero");
    EdgePhaseXMLNodes[m_NumEdgePhases - 1] = ePhase;

    /*
     * Push back indecises of edge phases
     */
    m_PhaseSpeciesStartIndex.resize(m_NumTotPhases + 1, 0);

    size_t istart = NumVolPhases_ + m_NumSurPhases + m_NumEdgePhases;
    m_PhaseSpeciesStartIndex[istart] = m_PhaseSpeciesStartIndex[istart - 1] + nSpecies;

    /*
     * Check elements list -> enforce strict conformance.
     */
    size_t numE = sp->nElements();
    for (size_t e = 0; e < numE; e++) {
        std::string symb1 = sp->elementName(e);
        doublevalue weight1 = sp->atomicWeight(e);
        int an = sp->atomicNumber(e);
        doublevalue entropy298 = sp->entropyElement298(e, true);
        int eType = sp->elementType(e);
        m_GlobalElementObj.addUniqueElement(symb1, weight1, an, entropy298, eType);
    }
    m_numElements = m_GlobalElementObj.nElements();

    /*
     * Check to see whether the phase has a kinetics object
     */
    EdgePhaseHasKinetics.resize(m_NumEdgePhases, 0);
    EdgePhaseHasKinetics[m_NumEdgePhases-1] = 0;
    if (ePhase->hasChild("kinetics")) {
        const XML_Node& kinNode = ePhase->child("kinetics");
        std::string smodel = kinNode["model"];
        if (smodel != "" && smodel != "none" && smodel != "None") {
            EdgePhaseHasKinetics[m_NumEdgePhases-1] = 1;
        }
    }

    PhaseList_.resize(m_NumTotPhases);
    PhaseNames_.resize(m_NumTotPhases);
    istart = NumVolPhases_ + m_NumSurPhases;
    for (size_t i = 0; i < m_NumEdgePhases; i++) {
        PhaseList_[i + istart] = EdgePhaseList[i];
        PhaseNames_[i + istart] = EdgePhaseList[i]->name();
    }

    calcElementMaps(true); 
}
//==================================================================================================================================
size_t PhaseList::volPhaseIndex(const ThermoPhase* const vp) const
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
size_t PhaseList::surPhaseIndex(const ThermoPhase* const sp) const
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
size_t PhaseList::edgePhaseIndex(const ThermoPhase* const ep) const
{
    for (size_t i = 0; i < m_NumEdgePhases; i++) {
        const ThermoPhase* const temp = EdgePhaseList[i];
        if (temp == ep) {
            return i;
        }
    }
    return npos;
}
//==================================================================================================================================
std::string PhaseList::phase_name(size_t iphGlob) const
{
    return PhaseNames_[iphGlob];
}
//==================================================================================================================================
std::string PhaseList::phaseIDString(size_t iphGlob) const
{
    return thermo(iphGlob).phaseIDString();
}
//==================================================================================================================================
phaseID PhaseList::phaseIdentifier(size_t iphGlob) const
{
    return thermo(iphGlob).phaseIdentifier();
}
//==================================================================================================================================
std::string PhaseList::phase_id(size_t iphGlob) const
{
    return thermo(iphGlob).id();
}
//==================================================================================================================================
size_t PhaseList::globalPhaseIndex(const ThermoPhase* const tp) const
{
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        if (tp == PhaseList_[i]) {
            return i;
        }
    }
    return npos;
}
//==================================================================================================================================
size_t PhaseList::globalPhaseIndex(const std::string& pIDStr) const
{
#ifdef DEBUG_MODE
    phaseID pID(pIDStr);
    return globalPhaseIndex(pID);
#else
    return globalPhaseIndex(phaseID(pIDStr));
#endif
}
//==================================================================================================================================
size_t PhaseList::globalPhaseIndex(const phaseID& pID) const
{
    for (size_t ip = 0; ip < m_NumTotPhases; ++ip) {
        if (PhaseList_[ip]->canBeThisPhase(pID)) {
            return ip;
        }
    }
    return npos;
}
//==================================================================================================================================
size_t PhaseList::globalPhaseNameIndex(const std::string& phase_name) const
{
    for (size_t ip = 0; ip < m_NumTotPhases; ++ip) {
        if (PhaseNames_[ip] == phase_name) {
            return ip;
        }
    }
    return npos;
}
//==================================================================================================================================
size_t PhaseList::globalSpeciesIndex(const ThermoPhase* const ttp, size_t k) const
{
    size_t iphGlob;
#ifdef DEBUG_MODE
    iphGlob = globalPhaseIndex(ttp);
    if (iphGlob == npos) {
        return npos;
    }
#else
    if ((iphGlob = globalPhaseIndex(ttp)) == npos) return npos;
#endif
    return globalSpeciesIndex(iphGlob, k);
}
//==================================================================================================================================
size_t PhaseList::globalSpeciesIndex(const std::string& sIDStr, const std::string pname) const
{
    speciesID sID(sIDStr);
    size_t iphGlob;
    size_t k;
    if (pname == "") {
        return globalSpeciesIndex(sID);
    } else {
#ifdef DEBUG_MODE
        iphGlob = globalPhaseNameIndex(pname);
        if (iphGlob == npos) {
            return npos;
        }
        k = PhaseList_[iphGlob]->speciesIndex(sID);
        if (k == npos) {
            return npos;
        }
#else
        if ((iphGlob = globalPhaseNameIndex(pname))      == npos) return npos;
        if ((k = PhaseList_[iphGlob]->speciesIndex(sID)) == npos) return npos;
#endif
        return (m_PhaseSpeciesStartIndex[iphGlob] + k);
    }
    return npos;
}
//==================================================================================================================================
size_t PhaseList::globalSpeciesIndex(const speciesID& sID) const
{
    size_t k;
    for (size_t iphGlob = 0; iphGlob < m_NumTotPhases; iphGlob++) {
#ifdef DEBUG_MODE
        ThermoPhase* tp = PhaseList_[iphGlob];
        if (tp->canBeInThisPhase(sID)) {
            k = tp->speciesNameIndex(sID.m_species_name);
            if (k != npos) {
                return ( m_PhaseSpeciesStartIndex[iphGlob] + k);
            }
        }
#else
        if ((k = PhaseList_[iphGlob]->speciesIndex(sID)) != npos) return ( m_PhaseSpeciesStartIndex[iphGlob] + k);
#endif
    }
    return npos;
}
//==================================================================================================================================
size_t PhaseList::globalSpeciesIndex(size_t iphGlob, size_t k) const
{
#ifdef DEBUG_MODE
    AssertTrace(iphGlob < m_NumTotPhases);
    AssertTrace(k < PhaseList_[iphGlob]->nSpecies());
    size_t istart = m_PhaseSpeciesStartIndex[iphGlob];
    return (istart + k);
#else
    return (m_PhaseSpeciesStartIndex[iphGlob] + k);
#endif
}
//==================================================================================================================================
size_t PhaseList::globalSpeciesIndexVolPhaseIndex(size_t iphVol, size_t k) const
{
#ifdef DEBUG_MODE
    if ( ! m_OrderedByDimensionality) {
      throw CanteraError(" PhaseList::globalSpeciesIndexVolPhaseIndex", "m_OrderedByDimensionality not done yet");
    }
    AssertTrace(iphVol != npos);
    AssertTrace(k != npos);
    AssertTrace(iphVol < NumVolPhases_);
    const ThermoPhase* const tp = VolPhaseList[iphVol];
    AssertTrace(k < tp->nSpecies());
    size_t istart = m_PhaseSpeciesStartIndex[iphVol];
    return (istart + k);
#else
    return ( m_PhaseSpeciesStartIndex[iphVol] + k);
#endif
}
//==================================================================================================================================
size_t PhaseList::globalSpeciesIndexSurPhaseIndex(size_t iphSur, size_t k) const
{
#ifdef DEBUG_MODE
    if ( ! m_OrderedByDimensionality) {
      throw CanteraError(" PhaseList::globalSpeciesIndexSurPhaseIndex", "m_OrderedByDimensionality not done yet");
    }
    AssertTrace(iphSur != npos);
    AssertTrace(k != npos);
    AssertTrace(iphSur < m_NumSurPhases);
    const ThermoPhase* const tp = SurPhaseList[iphSur];
    AssertTrace(k < tp->nSpecies());
    size_t phaseIndex = NumVolPhases_ + iphSur;
    return ( m_PhaseSpeciesStartIndex[phaseIndex] + k);
#else
    return ( m_PhaseSpeciesStartIndex[NumVolPhases_ + iphSur] + k);
#endif
}
//==================================================================================================================================
size_t PhaseList::phaseIndexFromGlobalSpeciesIndex(size_t kGlob) const
{
#ifdef DEBUG_MODE
    AssertTrace(kGlob != npos);
    if (kGlob >= m_PhaseSpeciesStartIndex[m_NumTotPhases]) {
        throw CanteraError("PhaseList::phaseIndexFromGlobalSpeciesIndex ", "indexing error");
    }
#endif
    for (size_t iphGlob = 0; iphGlob < m_NumTotPhases; iphGlob++) {
        if (kGlob < m_PhaseSpeciesStartIndex[iphGlob + 1]) {
            return iphGlob;
        }
    }
    return npos;
}
//==================================================================================================================================
ThermoPhase* PhaseList::phasePtr(const char* const pName) const
{
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        //if phase names match, return this phase
        if (!PhaseNames_[i].compare(pName)) {
            return PhaseList_[i];
        }
    }
    // did not find matching phase
    return nullptr;
}
//==================================================================================================================================
ThermoPhase* PhaseList::phasePtr(const std::string& pname) const
{
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        //if phase names match, return this phase
        if (!PhaseNames_[i].compare(pname)) {
            return PhaseList_[i];
        }
    }
    return nullptr;
}
//====================================================================================================================================
size_t PhaseList::dimPhaseIndexFromGlobalPhaseIndex(size_t iphGlob, int *dimPtr)
{
    if (! m_OrderedByDimensionality) {
        throw CanteraError(" PhaseList:: dimPhaseIndexFromGlobalPhaseIndex", "not handled");
    }
    if (dimPtr) {
        *dimPtr = PhaseList_[iphGlob]->nDim();
    }
    if (iphGlob < NumVolPhases_) {
         return iphGlob; 
    }
    iphGlob -= NumVolPhases_;
    if (iphGlob < m_NumSurPhases) {
        return iphGlob;
    }
#ifdef DEBUG_MODE
    iphGlob -= m_NumSurPhases;
    if (iphGlob < m_NumEdgePhases) {
        return iphGlob;
    }
    return npos;
#else
    return iphGlob - m_NumSurPhases;
#endif
}
//====================================================================================================================================
void PhaseList::
getLocalIndecisesFromGlobalSpeciesIndex(size_t kGlob, size_t& iphGlob, size_t& k) const
{
#ifdef DEBUG_MODE
    AssertTrace(kGlob != npos);
    AssertTrace(kGlob < m_NumTotSpecies);
#endif
    iphGlob = phaseIndexFromGlobalSpeciesIndex(kGlob);
#ifdef DEBUG_MODE
    if (iphGlob == npos) {
         throw CanteraError("PhaseList::getLocalIndecisesFromGlobalSpeciesIndex()", "Could not fund phase index");
    }
#endif
    k = (kGlob - m_PhaseSpeciesStartIndex[iphGlob]);
}
//====================================================================================================================
//! Kernel function that checks consistency of the phase id() and the the name and number/name/order of species in a phase 
/*!
 *   It returns true if the phases are the same and have the same species in it in the same order.
 *   Note, the two phases may have a different phase name().
 *
 *   @param[in]      tpA                    Pointer to the first %ThermoPhase object
 *   @param[in]      tpB                    Pointer to the second ThermoPhase object
 *
 *   @return                                Returns true if the two ThermoPhase classes have the same eos and the same
 *                                          phase id() and the same species names. Returns false otherwise.
 */
static bool ThermoPhasesTheSameNames(const ThermoPhase* const tpA, const ThermoPhase* const tpB)
{
    // Check the id() attribute of the phase to see that it is the same
    std::string sna = tpA->id();
    std::string snb = tpB->id();
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
	std::string spA = tpA->speciesName(k);
	std::string spB = tpB->speciesName(k);
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
    //          volume phases are same number and are the same ThermoPhase id()
    //          surface phases are same number and are the same ThermoPhase id()
    int compGood = 1;
    if (m_NumTotPhases == plGuest->m_NumTotPhases) {
	if (NumVolPhases_ == plGuest->NumVolPhases_) {
	    for (size_t i = 0; i < NumVolPhases_; ++i) {
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
	    for (size_t i = 0; i < m_NumSurPhases; ++i) {
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
        if (m_NumEdgePhases == plGuest->m_NumEdgePhases) {
            for (size_t i = 0; i < m_NumEdgePhases; ++i) {
                const ThermoPhase* tpa = EdgePhaseList[i];
                const ThermoPhase* tpb = plGuest->EdgePhaseList[i];
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
    for (size_t i = 0; i < NumVolPhases_; ++i) {
	bool iFound = false;
	int numBestAgreement = 0;
	int numCurrentAgreement = 0;
	const ThermoPhase* tpa = VolPhaseList[i];
	std::string ida = tpa->id();
	std::string pna = tpa->name();
	int iEOSa = tpa->eosType();
	const std::vector<std::string>&  sNa_list = tpa->speciesNames();
	for (size_t iGuest =  0; i <  plGuest->NumVolPhases_; ++i) {
	    const ThermoPhase* tpb = plGuest->VolPhaseList[i];
	    numCurrentAgreement = 0;
	    if (ida == tpb->id()) {
		numCurrentAgreement++;
	    }
	    if (pna == tpb->name()) {
		numCurrentAgreement++;
	    }
	    for (size_t kGuest =  0; kGuest <  tpb->nSpecies(); ++kGuest) {
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
    for (size_t i = 0; i <  m_NumSurPhases; ++i) {
	bool iFound = false;
	int numBestAgreement = 0;
	int numCurrentAgreement = 0;
	const ThermoPhase* tpa = SurPhaseList[i];
	std::string ida = tpa->id();
	std::string pna = tpa->name();
	int iEOSa = tpa->eosType();
	const std::vector<std::string>&  sNa_list = tpa->speciesNames();
	for (size_t iGuest =  0; i <  plGuest->m_NumSurPhases; ++i) {
	    const ThermoPhase* tpb = plGuest->SurPhaseList[i];
	    numCurrentAgreement = 0;
	    if (ida == tpb->id()) {
		numCurrentAgreement++;
	    }
	    if (pna == tpb->name()) {
		numCurrentAgreement++;
	    }
	    for (size_t kGuest =  0; kGuest <  tpb->nSpecies(); ++kGuest) {
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
PhaseList::thermo(int iphGlob) const
{
#ifdef DEBUG_MODE
    if (iphGlob < 0) {
        throw CanteraError("PhaseList::thermo()", "index out of range");
    }
#endif
    return thermo((size_t) iphGlob);
}
//==================================================================================================================================
ThermoPhase&
PhaseList::thermo(size_t iphGlob) const
{
#ifdef DEBUG_MODE
    if (iphGlob == npos) {
        throw CanteraError("PhaseList::phase()", "index out of range");
    }
    if (iphGlob >= (NumVolPhases_ + m_NumSurPhases + m_NumEdgePhases)) {
        throw CanteraError("PhaseList::phase()", "index out of range");
    }
#endif
    return *(PhaseList_[iphGlob]);
}
//==================================================================================================================================
ThermoPhase&
PhaseList::thermo(const char* const pName) const
{
#ifdef DEBUG_MODE
    std::string pp(pName);
    return thermo(pp);
#else
    return thermo( std::string(pName) );
#endif
}
//==================================================================================================================================
ThermoPhase& PhaseList::thermo(const std::string& pname) const
{
    for (size_t ip = 0; ip <  m_NumTotPhases; ++ip) {
	if (pname == PhaseNames_[ip]) {
	    return *PhaseList_[ip];
	}
    }
    throw CanteraError("PhaseList::thermo()", "Can't find phase named, " + pname + ", in list of phases");
}
//==================================================================================================================================
const Elements* PhaseList::globalElements() const
{
    return &m_GlobalElementObj;
}
//==================================================================================================================================
std::string PhaseList::elementName(size_t mIndx) const
{
    return m_GlobalElementObj.elementName(mIndx);
}
//==================================================================================================================================
size_t PhaseList::elementIndex(const std::string& elemName) const
{
    return m_GlobalElementObj.elementIndex(elemName);
}
//==================================================================================================================================
size_t PhaseList::nElements() const
{
    return m_numElements;
}
//==================================================================================================================================
XML_Node* PhaseList::surPhaseXMLNode(size_t iphSur) const
{
    return SurPhaseXMLNodes[iphSur];
}
//==================================================================================================================================
XML_Node* PhaseList::edgePhaseXMLNode(size_t iphEdge) const
{
    return EdgePhaseXMLNodes[iphEdge];
}
//==================================================================================================================================
size_t PhaseList::localToGlobalElementIndex(size_t iphGlob, size_t mLocalIndex)
{
#ifdef DEBUG_MODE
    if (iphGlob >= PhaseList_.size()) {
        throw CanteraError("PhaseList::localToGlobalElementIndex()", "phase index is out of range");
    }
    if (mLocalIndex >= PhaseList_[iphGlob]->nElements()) {
        throw CanteraError("PhaseList::localToGlobalElementIndex()", "element index is out of range");
    }
#endif
    return m_localToGlobalEMap[iphGlob][mLocalIndex];
}
//==================================================================================================================================
size_t PhaseList::globalToLocalElementIndex(size_t iphGlob, size_t mGlobalIndex)
{
#ifdef DEBUG_MODE
    if (iphGlob >= PhaseList_.size()) {
        throw CanteraError("PhaseList::globalToLocalElementIndex()", "phase index is out of range");
    }
    if (mGlobalIndex >= m_GlobalElementObj.nElements()) {
        throw CanteraError("PhaseList::globalToLocalElementIndex()", "element index is out of range");
    }
#endif
    return m_globalToLocalEMap(mGlobalIndex, iphGlob);
}
//====================================================================================================================
doublevalue PhaseList::nAtoms(const size_t kGlob, const size_t eGlob) const
{
      return m_atoms(eGlob, kGlob);
}
//====================================================================================================================
XML_Node* PhaseList::volPhaseXMLNode(size_t iphVol) const
{
#ifdef DEBUG_MODE
    AssertTrace(iphVol < NumVolPhases_);
#endif
    return VolPhaseXMLNodes[iphVol];
}
//====================================================================================================================
ThermoPhase& PhaseList::volPhase(size_t iphVol)
{
#ifdef DEBUG_MODE
    AssertTrace(iphVol < NumVolPhases_);
#endif
    return *(VolPhaseList[iphVol]);
}
//====================================================================================================================
ThermoPhase& PhaseList::surPhase(size_t iphSur)
{
#ifdef DEBUG_MODE
    AssertTrace(iphSur < m_NumSurPhases);
#endif
    return *(SurPhaseList[iphSur]);
}
//====================================================================================================================
ThermoPhase& PhaseList::edgePhase(size_t iphEdge)
{
#ifdef DEBUG_MODE
    AssertTrace(iphEdge < m_NumEdgePhases);
#endif
    return *(EdgePhaseList[iphEdge]);
}
//====================================================================================================================
size_t PhaseList::nVolPhases() const
{
    return NumVolPhases_;
}
//====================================================================================================================
size_t PhaseList::nSurPhases() const
{
    return m_NumSurPhases;
}
//====================================================================================================================
size_t PhaseList::nEdgePhases() const
{
    return m_NumEdgePhases;
}
//====================================================================================================================
size_t PhaseList::nVolSpecies() const
{
    return m_totNumVolSpecies;
}
//====================================================================================================================
size_t PhaseList::nSurSpecies() const
{
    return m_totNumSurSpecies;
}
//====================================================================================================================
size_t PhaseList::nEdgeSpecies() const
{
    return m_totNumEdgeSpecies;
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
bool PhaseList::volPhaseHasKinetics(size_t iphVol) const
{
    return VolPhaseHasKinetics[iphVol];
}
//====================================================================================================================
bool PhaseList::surPhaseHasKinetics(size_t iphSur) const
{
    return SurPhaseHasKinetics[iphSur];
}
//====================================================================================================================
bool PhaseList::edgePhaseHasKinetics(size_t iphEdge) const
{
    return EdgePhaseHasKinetics[iphEdge];
}
//====================================================================================================================
bool PhaseList::phaseHasKinetics(size_t iphGlob) const
{
    if (iphGlob < NumVolPhases_) {
        return VolPhaseHasKinetics[iphGlob];
    }
    iphGlob -= NumVolPhases_;
    if (iphGlob < m_NumSurPhases) {
        return SurPhaseHasKinetics[iphGlob];
    }
    iphGlob -= m_NumSurPhases;
    return EdgePhaseHasKinetics[iphGlob];
}
//====================================================================================================================
std::string PhaseList::speciesName(size_t kGlob) const
{
#ifdef DEBUG_MODE
    AssertTrace(kGlob != npos);
    AssertTrace(kGlob < m_NumTotSpecies);
    size_t iphGlob = phaseIndexFromGlobalSpeciesIndex(kGlob);
    size_t kStart = m_PhaseSpeciesStartIndex[iphGlob];
    size_t kLocal = kGlob - kStart;
    ThermoPhase& tp = thermo(iphGlob);
    return tp.speciesName(kLocal);
#else
    size_t iphGlob = getPhaseIndexFromGlobalSpeciesIndex(kGlob);
    return thermo(iphGlob).speciesName( kGlob - m_PhaseSpeciesStartIndex[iphGlob] );
#endif
}
//======================================================================================================================
std::string PhaseList::speciesIDString(size_t kGlob) const
{
    size_t iphGlob = phaseIndexFromGlobalSpeciesIndex(kGlob);
    return thermo(iphGlob).speciesIDString( kGlob - m_PhaseSpeciesStartIndex[iphGlob] );
}
//======================================================================================================================
speciesID PhaseList::speciesIdentifier(size_t kGlob) const
{
    size_t iphGlob = phaseIndexFromGlobalSpeciesIndex(kGlob);
    return thermo(iphGlob).speciesIdentifier( kGlob - m_PhaseSpeciesStartIndex[iphGlob] );
}
//======================================================================================================================
void PhaseList::setState_TP(doublevalue temperature, doublevalue pressure)
{
    for (size_t i = 0; i <  NumVolPhases_; i++) {
        ThermoPhase *tp = VolPhaseList[i];
        tp->setState_TP(temperature, pressure);
    } 
    for (size_t i = 0; i < m_NumSurPhases; i++) {
        ThermoPhase *tp = SurPhaseList[i];
        tp->setState_TP(temperature, pressure);
    }
    for (size_t i = 0; i < m_NumEdgePhases; i++) {
        ThermoPhase *tp = EdgePhaseList[i];
        tp->setState_TP(temperature, pressure);
    }
}
//==================================================================================================================================
void PhaseList::calcElementMaps(bool forceAll)
{
    bool redoAll = forceAll;
    if (m_globalToLocalEMap.nRows() != m_numElements) {
        redoAll = true;
    }
    m_globalToLocalEMap.resize(m_numElements, m_NumTotPhases);
    m_localToGlobalEMap.resize(m_NumTotPhases);
    ThermoPhase* tp;
    std::string ename;
    size_t eGlob, eLoc;
    size_t eGtot = m_GlobalElementObj.nElements();

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        std::vector<size_t> localEmap;
        tp = PhaseList_[iph];
        for (eLoc = 0; eLoc < tp->nElements(); ++eLoc) {
            ename = tp->elementName(eLoc);
            eGlob = m_GlobalElementObj.elementIndex(ename); 
            localEmap.push_back(eGlob);
        }
        m_localToGlobalEMap[iph] = localEmap;
        for (size_t eGlob = 0; eGlob < eGtot; eGlob++) {
            ename = m_GlobalElementObj.elementName(eGlob);
            eLoc = tp->elementIndex(ename); 
            m_globalToLocalEMap(eGlob,iph) = eLoc;
        }
    }

    // Fill in the m_atoms matrix from the underlying ThermoPhases
    m_atoms.resize(m_numElements, m_NumTotSpecies, 0.0);
    for (eGlob = 0; eGlob < eGtot; eGlob++) {
        size_t kGlob = 0;
        // iterate over the phases
        for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
            tp =PhaseList_[iph];
            eLoc = m_globalToLocalEMap(eGlob, iph);
            for (size_t kp = 0; kp < tp->nSpecies(); kp++) {
                if (eLoc != npos) {
                    m_atoms(eGlob, kGlob) = tp->nAtoms(kp, eLoc);
                } else {
                    m_atoms(eGlob,kGlob) = 0.0;
                }
                kGlob++;
            }
        }
    }
}
//==================================================================================================================================
void PhaseList::resizeGlobalSpeciesVectorForAddedPhase(std::vector<doublevalue>& speciesProp)
{
   if (m_lastPhaseIndexAdded == npos) {
       return;
   } 
   size_t cursize = speciesProp.size();
   size_t newsize = cursize + m_nSpLastPhaseAdded;
   if (m_totNumSurSpecies != newsize) {
      throw CanteraError("PhaseList()::resizeGlobalSpeciesVectorForLastPhaseAdded()", " possible misapplied situation");
   }
   speciesProp.resize(newsize, 0.0);
   size_t nspck = m_PhaseSpeciesStartIndex[m_lastPhaseIndexAdded+1] - m_PhaseSpeciesStartIndex[m_lastPhaseIndexAdded];
   if (nspck !=  m_nSpLastPhaseAdded) {
       throw CanteraError("PhaseList()::resizeGlobalSpeciesVectorForLastPhaseAdded()", "species check failure");
   }
   size_t kNewStart = m_PhaseSpeciesStartIndex[m_NumTotPhases];
   size_t numSpeciesToMove = kNewStart - m_PhaseSpeciesStartIndex[m_lastPhaseIndexAdded+1];
   if (numSpeciesToMove > 0) {
       size_t kOldStart = kNewStart - m_nSpLastPhaseAdded;
       // have to do this by hand as there is overlapping memory 
       for (size_t k = 0; k < numSpeciesToMove; ++k) {
           speciesProp[kNewStart - k] = speciesProp[kOldStart - k];
       }
       // Zero the entries for the new phase, if we have to.
       if (m_lastPhaseIndexAdded != (m_NumTotPhases - 1)) {
           kNewStart = m_PhaseSpeciesStartIndex[m_lastPhaseIndexAdded];
           for (size_t k = 0; k < m_nSpLastPhaseAdded; ++k) { 
               speciesProp[kNewStart + k] = 0.0;
           }
       }
   }
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

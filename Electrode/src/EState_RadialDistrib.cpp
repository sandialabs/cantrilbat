/*
 * $Id: EState_RadialDiffusion.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */


#include "cantera/base/ctml.h"

#include "EState_RadialDistrib.h"

#include "Electrode_SimpleDiff.h"

#include <cstdio>
#include <cstdlib>
using namespace Cantera;
using namespace std;

namespace Cantera
{
//======================================================================================================================
/*
 * EState constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
EState_RadialDistrib::EState_RadialDistrib() :
    EState(),
    numRCells_(0),
    numKRSpecies_(0),
    numSPhases_(0),
    rnodePos_(0),
    cellBoundR_(0),
    rLatticeCBR_(0),
    concTot_SPhase_Cell_(0),
    concKRSpecies_Cell_(0),
    spMoles_KRsolid_Cell_(0),
    onRegionBoundary_(-1)
{
    EST_fileToBeWritten_ = EST_RADIALDISTRIB;
    electrodeTypeString_ = "SimpleDiff";
}
//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
EState_RadialDistrib::EState_RadialDistrib(const EState_RadialDistrib& right) :
    EState(),
    numRCells_(0),
    numKRSpecies_(0),
    numSPhases_(0),
    rnodePos_(0),
    cellBoundR_(0),
    rLatticeCBR_(0),
    concTot_SPhase_Cell_(0),
    concKRSpecies_Cell_(0),
    spMoles_KRsolid_Cell_(0),
    onRegionBoundary_(-1)
{
    EST_fileToBeWritten_ = EST_MULTIPLATEAU;
    electrodeTypeString_ = "SimpleDiff";

    /*
     * Call the assignment operator.
     */
    EState_RadialDistrib::operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
EState_RadialDistrib& EState_RadialDistrib::operator=(const EState_RadialDistrib& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    EState::operator=(right);

    numRCells_                             = right.numRCells_;
    numKRSpecies_                          = right.numKRSpecies_;
    numSPhases_                            = right.numSPhases_;
    rnodePos_                              = right.rnodePos_;
    cellBoundR_                            = right.cellBoundR_;
    rLatticeCBR_                           = right.rLatticeCBR_;
    concTot_SPhase_Cell_                   = right.concTot_SPhase_Cell_;
    concKRSpecies_Cell_                    = right.concKRSpecies_Cell_;
    spMoles_KRsolid_Cell_                  = right.spMoles_KRsolid_Cell_;
    onRegionBoundary_                      = right.onRegionBoundary_;

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
EState_RadialDistrib::~EState_RadialDistrib()
{
}
//======================================================================================================================
// Duplicator function for this class
/*
 *  @return Returns a duplication of the current state as a pointer to the base class
 */
EState* EState_RadialDistrib::duplMyselfAsEState() const
{
    EState_RadialDistrib* es = new EState_RadialDistrib(*this);
    return dynamic_cast<EState*>(es);
}
//======================================================================================================================
int EState_RadialDistrib::initialize(const Cantera::Electrode_SimpleDiff* const e)
{
    EState::initialize(e);

    numRCells_                   = e->numRCells_;
    numKRSpecies_                = e->numKRSpecies_;
    numSPhases_                  = e->numSPhases_;
    rnodePos_                    = e->rnodePos_final_;
    cellBoundR_                  = e->cellBoundR_final_;
    rLatticeCBR_                 = e->rLatticeCBR_final_;
    concTot_SPhase_Cell_         = e->concTot_SPhase_Cell_final_;
    concKRSpecies_Cell_          = e->concKRSpecies_Cell_final_;
    spMoles_KRsolid_Cell_        = e->spMoles_KRsolid_Cell_final_;

    return 1;
}
//======================================================================================================================
// Write the ElectrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node*   EState_RadialDistrib::writeIdentificationToXML() const
{
    XML_Node* x = new XML_Node("ElectrodeIdentification");

    ctml::addString(*x, "electrodeTypeString", electrodeTypeString_);
    ctml::addInteger(*x, "EState_Type",         EST_fileToBeWritten_);
    ctml::addInteger(*x, "electrodeModelType",  electrodeModelType_);
    ctml::addInteger(*x, "electrodeDomainNumber",  electrodeDomainNumber_);
    ctml::addInteger(*x, "electrodeCellNumber",  electrodeCellNumber_);

    return x;
}
//======================================================================================================================
// Write the ElectrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node* EState_RadialDistrib::writeStateToXML() const
{
    XML_Node* x = EState::writeStateToXML();

    ctml::addInteger(*x, "numRCells", numRCells_);
    ctml::addInteger(*x, "numKRSpecies", numKRSpecies_);
    ctml::addInteger(*x, "numSPhases", numSPhases_);

    ctml::addNamedFloatArray(*x, "rnodePos", numRCells_, DATA_PTR(rnodePos_), "m");
    ctml::addNamedFloatArray(*x, "cellBoundR", numRCells_, DATA_PTR(cellBoundR_), "m");
    ctml::addNamedFloatArray(*x, "rLatticeCBR", numRCells_, DATA_PTR(rLatticeCBR_), "m");
    ctml::addNamedFloatArray(*x, "concTot_SPhase_Cell", numRCells_ * numSPhases_, DATA_PTR(concTot_SPhase_Cell_), "kmol/m3");
    ctml::addNamedFloatArray(*x, "concKRSpecies_Cell", numRCells_ * numKRSpecies_, DATA_PTR(concKRSpecies_Cell_), "kmol/m3");
    ctml::addNamedFloatArray(*x, "spMoles_KRsolid_Cell", numRCells_ * numKRSpecies_, DATA_PTR(spMoles_KRsolid_Cell_), "kmol");

    return x;
}

//======================================================================================================================
//  Read the state from the XML_Node  given by the argument
/*
 *  @return pointer to the XML_Node tree
 */
void EState_RadialDistrib::readStateFromXML(const XML_Node& xmlEState)
{
    EState::readStateFromXML(xmlEState);

    numRCells_ = ctml::getInteger(xmlEState, "numRCells");
    numKRSpecies_ = ctml::getInteger(xmlEState, "numKRSpecies");
    numSPhases_ = ctml::getInteger(xmlEState, "numRCells");


    ctml::getFloatArray(xmlEState, rnodePos_, true, "", "rnodePos");
    ctml::getFloatArray(xmlEState, cellBoundR_, true, "", "cellBoundR"); 
    ctml::getFloatArray(xmlEState, rLatticeCBR_, true, "", "rLatticeCBR");

    ctml::getFloatArray(xmlEState, concTot_SPhase_Cell_, true, "", "concTot_SPhase_Cell");
    ctml::getFloatArray(xmlEState, concKRSpecies_Cell_, true, "", "concKRSpecies_Cell");
    ctml::getFloatArray(xmlEState, spMoles_KRsolid_Cell_, true, "", "spMoles_KRsolid_Cell");
}
//======================================================================================================================
// Set the State of this object from the state of the Electrode object
/*
 *  (virtual function)
 *
 *  virtual function, because the base class is called from general code, allowing the child classes to be invoked.
 *
 *  This function takes the electrode objects _final_ state and copies it into this object.
 *  This function must be carried out before the XML_Node tree is written.
 *
 *  @param e   Pointer to the Electrode object. Note, this class may use dynamic casting
 *             to choose a child object, and then may invoke an error if the match isn't
 *             correct.
 */
void EState_RadialDistrib::copyElectrode_intoState(const Cantera::Electrode* const e)
{
    EState::copyElectrode_intoState(e);

    const Cantera::Electrode_SimpleDiff* const emp = dynamic_cast<const Cantera::Electrode_SimpleDiff* const>(e);

    if (emp) {
	rnodePos_                    = emp->rnodePos_final_;
	cellBoundR_                  = emp->cellBoundR_final_;
	rLatticeCBR_                 = emp->rLatticeCBR_final_;
	concTot_SPhase_Cell_         = emp->concTot_SPhase_Cell_final_;
	concKRSpecies_Cell_          = emp->concKRSpecies_Cell_final_;
	spMoles_KRsolid_Cell_        = emp->spMoles_KRsolid_Cell_final_;
    } else {
        throw CanteraError("EState_RadialDistrib::copyElectrode_intoState","bad cast");
    }
}
//======================================================================================================================
//    Set the state of the Electrode from the state of this object
void EState_RadialDistrib::setStateElectrode_fromEState(Cantera::Electrode* const e) const
{
    EState::copyEState_toElectrode(e);

    Electrode_SimpleDiff* emp = dynamic_cast<Electrode_SimpleDiff*>(e);

    emp->rnodePos_final_             = 	rnodePos_;    
    emp->cellBoundR_final_           = 	cellBoundR_;
    emp->rLatticeCBR_final_          =  rLatticeCBR_;
    emp->concTot_SPhase_Cell_final_  = 	concTot_SPhase_Cell_;
    emp->concKRSpecies_Cell_final_   = 	concKRSpecies_Cell_;
    emp->spMoles_KRsolid_Cell_final_ = 	spMoles_KRsolid_Cell_;

    /*
     * Now we can do an update
     */
    emp->updateState();
    emp->stateToPhaseFlagsReconciliation(false);
    emp->setInitStateFromFinal(true);

}
//====================================================================================================================
} // End of namespace Cantera
//======================================================================================================================

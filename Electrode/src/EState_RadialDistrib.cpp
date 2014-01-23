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
    rnodePos_(0),
    cellBoundR_(0),
    rLatticeCBR_(0),
    concTot_SPhase_Cell_(0),
    concKRSpecies_Cell_(0),
    spMoles_KRsolid_Cell_(0),
    spMf_KRSpecies_Cell_(0),
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
   rnodePos_(0),
    cellBoundR_(0),
    rLatticeCBR_(0),
    concTot_SPhase_Cell_(0),
    concKRSpecies_Cell_(0),
    spMoles_KRsolid_Cell_(0),
    spMf_KRSpecies_Cell_(0),
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

    rnodePos_                              = right.rnodePos_;
    cellBoundR_                            = right.cellBoundR_;
    rLatticeCBR_                           = right.rLatticeCBR_;
    concTot_SPhase_Cell_                   = right.concTot_SPhase_Cell_;
    concKRSpecies_Cell_                    = right.concKRSpecies_Cell_;
    spMoles_KRsolid_Cell_                  = right.spMoles_KRsolid_Cell_;
    spMf_KRSpecies_Cell_                   = right.spMf_KRSpecies_Cell_;
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

    rnodePos_                    = e->rnodePos_final_;
    cellBoundR_                  = e->cellBoundR_final_;
    rLatticeCBR_                 = e->rLatticeCBR_final_;
    concTot_SPhase_Cell_         = e->concTot_SPhase_Cell_final_;
    concKRSpecies_Cell_          = e->concKRSpecies_Cell_final_;
    spMoles_KRsolid_Cell_        = e->spMoles_KRsolid_Cell_final_;
    spMf_KRSpecies_Cell_         = e->spMf_KRSpecies_Cell_final_;

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
XML_Node*   EState_RadialDistrib::writeStateToXML() const
{
    XML_Node* x = EState::writeStateToXML();

    // FILL IN WITH RADIAL STUFF

    return x;
}

//======================================================================================================================
//! Write the ElectrodeState to an XML_Node tree
/*!
 *  @return pointer to the XML_Node tree
 */
void EState_RadialDistrib::readStateFromXML(const XML_Node& xmlEState)
{
    EState::readStateFromXML(xmlEState);

    //inner_platNum_ = ctml::getInteger(xmlEState, "inner_platNum");

    // FILL IN WITH RADIAL STUFF
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
        // FILL IN WITH RADIAL stuff
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

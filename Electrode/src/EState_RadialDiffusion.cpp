/*
 * $Id: EState_RadialDiffusion.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>

#include "cantera/base/ctml.h"

#include "EState_RadialDiffusion.h"
#include "Electrode_RadialDiffusion_NoDiff.h"


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
EState_RadialDiffusion::EState_RadialDiffusion() :
    EState(),
{
    EST_fileToBeWritten_ = EST_RADIALDIFFUSION;
    electrodeTypeString_ = "RadialDiffusion_NoDiff";
}
//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
EState_RadialDiffusion::EState_RadialDiffusion(const EState_RadialDiffusion& right) :
    EState()
{
    EST_fileToBeWritten_ = EST_MULTIPLATEAU;
    electrodeTypeString_ = "RadialDiffusion_NoDiff";

    /*
     * Call the assignment operator.
     */
    EState_RadialDiffusion::operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
EState_RadialDiffusion& EState_RadialDiffusion::operator=(const EState_RadialDiffusion& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    EState::operator=(right);

    inner_platNum_   = right.inner_platNum_;
    second_platNum_  = right.second_platNum_;

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
EState_RadialDiffusion::~EState_RadialDiffusion()
{
}
//======================================================================================================================
// Duplicator function for this class
/*
 *  @return Returns a duplication of the current state as a pointer to the base class
 */
EState* EState_RadialDiffusion::duplMyselfAsEState() const
{
    EState_RadialDiffusion* es = new EState_RadialDiffusion(*this);
    return dynamic_cast<EState*>(es);
}
//======================================================================================================================
int EState_RadialDiffusion::initialize(const Cantera::Electrode_RadialDiffusion_NoDiff* const e)
{
    EState::initialize(e);

    return 1;
}
//======================================================================================================================
// Write the ElectrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node*   EState_RadialDiffusion::writeIdentificationToXML() const
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
XML_Node*   EState_RadialDiffusion::writeStateToXML() const
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
void EState_RadialDiffusion::readStateFromXML(const XML_Node& xmlEState)
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
void EState_RadialDiffusion::copyElectrode_intoState(const Cantera::Electrode* const e)
{
    EState::copyElectrode_intoState(e);

    const Cantera::Electrode_RadialDiffusion_NoDiff* const emp = dynamic_cast<const Cantera::Electrode_RadialDiffusion_NoDiff* const>(e);

    if (emp) {
        // FILL IN WITH RADIAL stuff
    } else {
        throw CanteraError("EState_RadialDiffusion::copyElectrode_intoState","bad cast");
    }
}
//======================================================================================================================
//    Set the state of the Electrode from the state of this object
void EState_RadialDiffusion::setStateElectrode_fromEState(Cantera::Electrode* const e) const
{

    EState::copyEState_toElectrode(e);

    Electrode_RadialDiffusion_NoDiff* emp = dynamic_cast<Electrode_RadialDiffusion_NoDiff*>(e);

    /*
     *  Set the plateaus according to what's storred in the EState file
     */
    emp->plateauSetup2(inner_platNum_, second_platNum_);

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
                                                                                           /

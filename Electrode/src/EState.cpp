/*
 * $Id: EState.cpp 584 2013-04-03 00:29:38Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>

#include "cantera/base/ctml.h"

#include "EState.h"
#include "Electrode.h"
#include "Electrode_MP_RxnExtent.h"
#include "Electrode_Factory.h"


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
EState::EState() :
    eRef_(0),
    electrodeTypeString_(""),
    EST_lastFileRead_(EST_CSTR),
    EST_Version_lastFileRead_(1),
    EST_fileToBeWritten_(EST_CSTR),
    EST_Version_fileToBeWritten_(1),
    spMoles_(0),
    phaseVoltages_(0),
    temperature_(298.15),
    pressure_(OneAtm),
    electrodeModelType_(0),
    electrodeDomainNumber_(-1),
    electrodeCellNumber_(-1),
    particleNumberToFollow_(0.0),
    electrodeSolidVolume_(0.0),
    grossVolume_(0.0),
    radiusExterior_(0.0),
    surfaceAreaRS_(0),
    electrodeMoles_(0.0),
    electrodeCapacityType_(CAPACITY_ANODE_ECT),
    capacityLeft_(0.0),
    capacityInitial_(0.0),
    depthOfDischarge_(0.0),
    depthOfDischargeStarting_(0.0),
    relativeElectronsDischargedPerMole_(0.0),
    relativeDepthOfDischarge_(0.0),
    capacityDischargedToDate_(0.0),
    deltaTsubcycle_init_next_(1.0E300)
{
}
//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
EState::EState(const EState& right) :
    eRef_(right.eRef_),
    electrodeTypeString_(""),
    EST_lastFileRead_(EST_CSTR),
    EST_Version_lastFileRead_(1),
    EST_fileToBeWritten_(EST_CSTR),
    EST_Version_fileToBeWritten_(1),
    spMoles_(0),
    phaseVoltages_(0),
    temperature_(298.15),
    pressure_(OneAtm),
    electrodeModelType_(0),
    electrodeDomainNumber_(-1),
    electrodeCellNumber_(-1),
    particleNumberToFollow_(0.0),
    electrodeSolidVolume_(0.0),
    grossVolume_(0.0),
    radiusExterior_(0.0),
    surfaceAreaRS_(0),
    electrodeMoles_(0.0),
    electrodeCapacityType_(CAPACITY_ANODE_ECT),
    capacityLeft_(0.0),
    capacityInitial_(0.0),
    depthOfDischarge_(0.0),
    depthOfDischargeStarting_(0.0),
    relativeElectronsDischargedPerMole_(0.0),
    relativeDepthOfDischarge_(0.0),
    capacityDischargedToDate_(0.0),
    deltaTsubcycle_init_next_(1.0E300)
{
    /*
     * Call the assignment operator.
     */
    EState::operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
EState& EState::operator=(const EState& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    eRef_                              = right.eRef_;
    electrodeTypeString_               = right.electrodeTypeString_;
    EST_lastFileRead_                  = right.EST_lastFileRead_;
    EST_Version_lastFileRead_          = right.EST_Version_lastFileRead_;
    spMoles_                           = right.spMoles_;
    phaseVoltages_                     = right.phaseVoltages_;
    temperature_                       = right.temperature_;
    pressure_                          = right.pressure_;
    electrodeModelType_                = right.electrodeModelType_;
    electrodeDomainNumber_             = right.electrodeDomainNumber_;
    electrodeCellNumber_               = right.electrodeCellNumber_;
    particleNumberToFollow_            = right.particleNumberToFollow_;
    electrodeSolidVolume_              = right.electrodeSolidVolume_;
    grossVolume_                       = right.grossVolume_;
    radiusExterior_                    = right.radiusExterior_;
    surfaceAreaRS_                     = right.surfaceAreaRS_;
    electrodeMoles_                    = right.electrodeMoles_;
    electrodeCapacityType_             = right.electrodeCapacityType_;
    capacityLeft_                      = right.capacityLeft_;
    capacityInitial_                   = right.capacityInitial_;
    depthOfDischarge_                  = right.depthOfDischarge_;
    depthOfDischargeStarting_          = right.depthOfDischargeStarting_;
    relativeElectronsDischargedPerMole_= right.relativeElectronsDischargedPerMole_;
    relativeDepthOfDischarge_          = right.relativeDepthOfDischarge_;
    capacityDischargedToDate_          = right.capacityDischargedToDate_;
    deltaTsubcycle_init_next_          = right.deltaTsubcycle_init_next_;

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
EState::~EState()
{
}
//======================================================================================================================
// Duplicator function for this class
/*
 *  @return Returns a duplication of the current state as a pointer to the base class
 */
EState* EState::duplMyselfAsEState() const
{
    EState* es = new EState(*this);
    return es;
}
//======================================================================================================================
int EState::initialize(const Cantera::Electrode* const e)
{
    eRef_ = e;
    electrodeTypeString_    = Electrode_Types_Enum_to_string(e->electrodeType());
    electrodeModelType_     = e->electrodeModelType_;
    electrodeDomainNumber_  = e->electrodeDomainNumber_;
    electrodeCellNumber_    = e->electrodeCellNumber_;

    copyElectrode_intoState(eRef_);
    return 1;
}
//======================================================================================================================
// Write the ElectrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node*   EState::writeIdentificationToXML() const
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
XML_Node*   EState::writeStateToXML() const
{
    XML_Node* x = new XML_Node("electrodeState");

    int nsp = spMoles_.size();
    ctml::addNamedFloatArray(*x, "spMoles", nsp, DATA_PTR(spMoles_), "kmol");

    int np = phaseVoltages_.size();
    ctml::addNamedFloatArray(*x, "phaseVoltages", np, DATA_PTR(phaseVoltages_), "volt");

    ctml::addFloat(*x, "temperature",  temperature_, "Kelvin");
    ctml::addFloat(*x, "pressure",  pressure_, "Pa");
    ctml::addInteger(*x, "electrodeModelType",  electrodeModelType_);
    ctml::addInteger(*x, "electrodeDomainNumber",  electrodeDomainNumber_);
    ctml::addInteger(*x, "electrodeCellNumber",  electrodeCellNumber_);
    ctml::addFloat(*x, "particleNumberToFollow",  particleNumberToFollow_, "");
    ctml::addFloat(*x, "electrodeSolidVolume", electrodeSolidVolume_, "m3");
    ctml::addFloat(*x, "grossVolume", grossVolume_, "m3");
    ctml::addFloat(*x, "radiusExterior", radiusExterior_, "m3");
    int ns = surfaceAreaRS_.size();
    ctml::addNamedFloatArray(*x, "surfaceAreaRS", ns, DATA_PTR(surfaceAreaRS_), "m2");
    ctml::addFloat(*x, "electrodeMoles", electrodeMoles_, "kmol");
    ctml::addInteger(*x, "electrodeCapacityType", (int) electrodeCapacityType_);
    ctml::addFloat(*x, "capacityLeft", capacityLeft_, "coulomb");
    ctml::addFloat(*x, "capacityInitial", capacityInitial_, "coulomb");
    ctml::addFloat(*x, "depthOfDischarge", depthOfDischarge_, "coulomb");
    ctml::addFloat(*x, "depthOfDischargeStarting", depthOfDischargeStarting_, "coulomb");
    ctml::addFloat(*x, "relativeElectronsDischargedPerMole", relativeElectronsDischargedPerMole_, "");
    ctml::addFloat(*x, "relativeDepthOfDischarge", relativeDepthOfDischarge_, "");
    ctml::addFloat(*x, "capacityDischargedToDate", capacityDischargedToDate_, "coulomb");
    ctml::addFloat(*x, "deltaTsubcycle_init_next", deltaTsubcycle_init_next_, "s");

    return x;
}
//======================================================================================================================
// Write the ElectrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
void EState::readStateFromXMLRoot(const XML_Node& xmlRoot)
{
    if (xmlRoot.hasChild("electrodeState")) {
        XML_Node& x = xmlRoot.child("electrodeState");
        readStateFromXML(x);
    }
}
//======================================================================================================================
void EState::readIdentificationFromXML(const XML_Node& xmlEState)
{
    std::string typeSS;
    ctml::getNamedStringValue(xmlEState,"electrodeTypeString", electrodeTypeString_ , typeSS);

    EST_lastFileRead_ = (EState_Type_Enum)  ctml::getInteger(xmlEState, "EState_Type");
    if (EST_lastFileRead_ !=  EST_CSTR) {
        throw CanteraError(" EState::readIdentificationFromXML()",
                           "read the wrong EState type - new situation");
    }

    EST_Version_lastFileRead_ = ctml::getInteger(xmlEState, "fileVersionNumber");
    if (EST_Version_lastFileRead_ != 1) {
        throw CanteraError(" EState::readIdentificationFromXML()",
                           "read the wrong EState version type - new situation");
    }

    electrodeCapacityType_ = (Cantera::Electrode_Capacity_Type_Enum) ctml::getInteger(xmlEState, "electrodeCapacityType");
    electrodeModelType_ = ctml::getInteger(xmlEState, "electrodeModelType");
    electrodeDomainNumber_ = ctml::getInteger(xmlEState, "electrodeDomainNumber");
    electrodeCellNumber_ = ctml::getInteger(xmlEState, "electrodeCellNumber");
}
//======================================================================================================================

//! Write the ElectrodeState to an XML_Node tree
/*!
 *  @return pointer to the XML_Node tree
 */
void EState::readStateFromXML(const XML_Node& xmlEState)
{
    ctml::getFloatArray(xmlEState, spMoles_, true, "", "spMoles");
    ctml::getFloatArray(xmlEState, phaseVoltages_, true, "volts", "phaseVoltages");
    temperature_ = ctml::getFloat(xmlEState, "temperature", "toSI");
    pressure_ = ctml::getFloat(xmlEState, "pressure", "toSI");
    electrodeModelType_ = ctml::getInteger(xmlEState, "electrodeModelType");
    electrodeDomainNumber_ = ctml::getInteger(xmlEState, "electrodeDomainNumber");
    electrodeCellNumber_ = ctml::getInteger(xmlEState, "electrodeCellNumber");
    particleNumberToFollow_ = ctml::getFloat(xmlEState, "particleNumberToFollow", "toSI");
    electrodeSolidVolume_ = ctml::getFloat(xmlEState, "electrodeSolidVolume", "toSI");
    grossVolume_ = ctml::getFloat(xmlEState, "grossVolume", "toSI");
    radiusExterior_ = ctml::getFloat(xmlEState, "radiusExterior", "toSI");
    ctml::getFloatArray(xmlEState, surfaceAreaRS_, true, "m2", "surfaceAreaRS");
    electrodeMoles_ = ctml::getFloat(xmlEState, "electrodeMoles", "toSI");
    electrodeCapacityType_ = (Cantera::Electrode_Capacity_Type_Enum) ctml::getInteger(xmlEState, "electrodeCapacityType");
    capacityLeft_ = ctml::getFloat(xmlEState, "capacityLeft", "toSI");
    capacityInitial_ = ctml::getFloat(xmlEState, "capacityInitial", "toSI");
    depthOfDischarge_ = ctml::getFloat(xmlEState, "depthOfDischarge", "toSI");
    depthOfDischargeStarting_ = ctml::getFloat(xmlEState, "depthOfDischargeStarting", "toSI");
    relativeElectronsDischargedPerMole_ = ctml::getFloat(xmlEState, "relativeElectronsDischargedPerMole", "toSI");
    relativeDepthOfDischarge_ = ctml::getFloat(xmlEState, "relativeDepthOfDischarge", "toSI");
    capacityDischargedToDate_ = ctml::getFloat(xmlEState, "capacityDischargedToDate", "toSI");
    deltaTsubcycle_init_next_ = ctml::getFloat(xmlEState, "deltaTsubcycle_init_next", "toSI");
}
//======================================================================================================================
// Set the State of this object from the state of the Electrode object
void EState::copyElectrode_intoState(const Cantera::Electrode* const e)
{
    eRef_                              = e;
    spMoles_                           = e->spMoles_final_;
    phaseVoltages_                     = e->phaseVoltages_;
    temperature_                       = e->temperature_;
    pressure_                          = e->pressure_;
    electrodeModelType_                = e->electrodeModelType_;
    electrodeDomainNumber_             = e->electrodeDomainNumber_;
    electrodeCellNumber_               = e->electrodeCellNumber_;
    particleNumberToFollow_            = e->particleNumberToFollow_;
    electrodeSolidVolume_              = e->ElectrodeSolidVolume_;
    double currentSolidVol = e->SolidVol();
    grossVolume_ = currentSolidVol / (1.0 - e->porosity_);
    radiusExterior_                    = e->Radius_exterior_final_;
    surfaceAreaRS_                     = e->surfaceAreaRS_final_;
    electrodeMoles_                    = e->SolidTotalMoles();
    electrodeCapacityType_             = e->capacityType();
    capacityLeft_                      = e->capacityLeft();
    capacityInitial_                   = e->capacityInitial();
    depthOfDischarge_                  = e->depthOfDischarge();
    depthOfDischargeStarting_          = e->depthOfDischargeStarting();

    relativeDepthOfDischarge_          = e->depthOfDischargeFraction();

    capacityDischargedToDate_          = e->capacityDischarged();

    deltaTsubcycle_init_next_          = e->deltaTsubcycle_init_next_;

    double capUsed = capacityInitial_ - capacityLeft_;

    double relCapUsed = capUsed / electrodeMoles_;

    const Electrode_MP_RxnExtent* emp = dynamic_cast<const Electrode_MP_RxnExtent*>(e);
    if (emp) {
        relativeElectronsDischargedPerMole_ = emp->RelativeExtentRxn_final_;
    } else {
        relativeElectronsDischargedPerMole_ = relCapUsed;
    }

}
//======================================================================================================================
//Set the state of the Electrode from the state of this object
void EState::setStateElectrode_fromEState(Cantera::Electrode* const e) const
{
    EState::copyEState_toElectrode(e);

    e->stateToPhaseFlagsReconciliation(false);
    e->updateState();
    e->setInitStateFromFinal(true);
}
//======================================================================================================================
// Set the state of the Electrode Class from the state of the EState object
/*
 *  This is not a virtual function
 */
void EState::copyEState_toElectrode(Cantera::Electrode* const e) const
{
    e->spMoles_final_                     = spMoles_;
    e->phaseVoltages_                     = phaseVoltages_;
    e->temperature_                       = temperature_;
    e->pressure_                          = pressure_;
    e->electrodeModelType_                = electrodeModelType_;
    e->electrodeDomainNumber_             = electrodeDomainNumber_;
    e->electrodeCellNumber_               = electrodeCellNumber_;
    e->particleNumberToFollow_            = particleNumberToFollow_;
    e->ElectrodeSolidVolume_              = electrodeSolidVolume_;
    e->Radius_exterior_final_             = radiusExterior_;
    e->surfaceAreaRS_final_               = surfaceAreaRS_;
    e->depthOfDischargeStarting_          = depthOfDischargeStarting_;
    e->electronKmolDischargedToDate_      = capacityDischargedToDate_ / Cantera::Faraday;

    e->setCapacityType(electrodeCapacityType_);

    for (int iph = 0; iph < e->m_NumTotPhases; iph++) {
        e->updatePhaseNumbers(iph);
    }

    e->deltaTsubcycle_init_next_          = deltaTsubcycle_init_next_;
    e->deltaTsubcycle_init_init_          = deltaTsubcycle_init_next_;


    Electrode_MP_RxnExtent* const emp = dynamic_cast<Electrode_MP_RxnExtent* const>(e);
    if (emp) {
        emp->RelativeExtentRxn_final_ = relativeElectronsDischargedPerMole_;
    }
}
//====================================================================================================================
} // End of namespace Cantera
//======================================================================================================================

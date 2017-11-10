/*
 * $Id: EState.cpp 584 2013-04-03 00:29:38Z hkmoffa $
 */

#include "EState.h"
#include "Electrode.h"
#include "Electrode_MP_RxnExtent.h"
#include "Electrode_CSTR.h"
#include "Electrode_Factory.h"
#include "EState_XML.h"

#include "mdp_allo.h"

#include <string>

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
EState_Identification::EState_Identification() :
    electrodeTypeString(""),
    EST_Type(EST_UNKNOWN_TYPE),
    EState_Type_String("Unknown"),
    EST_Version(1),
    electrodeChemistryModelType(0),
    electrodeDomainNumber(0),
    electrodeCellNumber(0),
    electrodeCapacityType(CAPACITY_ANODE_ECT)
{
}
//==================================================================================================================================
ZZCantera::XML_Node* EState_Identification::writeIdentificationToXML() const
{
    XML_Node* x = new XML_Node("ElectrodeIdentification");
    ZZctml::addNamedString(*x, "electrodeTypeString", electrodeTypeString);
    ZZctml::addInteger(*x, "EState_Type",         EST_Type);
    ZZctml::addNamedString(*x, "EState_Type_String", EState_Type_String);
    ZZctml::addInteger(*x, "fileVersionNumber",  EST_Version);
    ZZctml::addInteger(*x, "electrodeModelType",  electrodeChemistryModelType);
    ZZctml::addInteger(*x, "electrodeDomainNumber",  electrodeDomainNumber);
    ZZctml::addInteger(*x, "electrodeCellNumber",  electrodeCellNumber);
    if (electrodeCapacityType == CAPACITY_ANODE_ECT) {
	ZZctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Anode");
    } else if (electrodeCapacityType == CAPACITY_CATHODE_ECT) {
	ZZctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Cathode");
    } else {
	ZZctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Other");
    }
    return x;
}
//==================================================================================================================================
void EState_Identification::readIdentificationFromXML(const XML_Node& xmlEI)
{  
    std::string typeSS;
    std::string nn = xmlEI.name();
    const XML_Node* x = &xmlEI;
    if (nn != "ElectrodeIndentification") {
        x = xmlEI.findByName("ElectrodeIdentification");
        if (!x) {
            throw Electrode_Error("EState::readIdentificationFromXM",
                                  "Could not find the XML node named ElectrodeIdentification");
        }
    }

    ZZctml::getNamedStringValue(*x, "electrodeTypeString", electrodeTypeString , typeSS);

    EST_Type = (EState_Type_Enum)  ZZctml::getInteger(*x, "EState_Type");
    if (x->hasChild("EState_Type_String")) {
	ZZctml::getNamedStringValue(*x, "EState_Type_String", EState_Type_String , typeSS);
    } else {
        EState_Type_String = "EState_CSTR";
    }

    if (x->hasChild("fileVersionNumber")) {
	EST_Version = ZZctml::getInteger(*x, "fileVersionNumber");
    } else {
	EST_Version =  1;

    }

    if (x->hasChild("electrodeModelType")) {
        electrodeChemistryModelType = ZZctml::getInteger(*x, "electrodeModelType");
    } else {
        electrodeChemistryModelType = 0;
    }

    if (x->hasChild("electrodeDomainNumber")) {
        electrodeDomainNumber = ZZctml::getInteger(*x, "electrodeDomainNumber");
    } else {
        electrodeDomainNumber = 0;
    }

    if (x->hasChild("electrodeCellNumber")) {
     electrodeCellNumber = ZZctml::getInteger(*x, "electrodeCellNumber");
    } else {
     electrodeCellNumber = 0;
    }
 
    if (x->hasChild("electrodeCapacityType")) {
	ZZctml::getNamedStringValue(*x, "electrodeCapacityType", nn, typeSS);
	if (nn == "Capacity_Anode") {
	    electrodeCapacityType = CAPACITY_ANODE_ECT;
	} else if (nn == "Capacity_Cathode") {
	    electrodeCapacityType = CAPACITY_CATHODE_ECT;
	} else if (nn == "Capacity_Other") {
	    electrodeCapacityType = CAPACITY_OTHER_ECT;
	}
    } else {
	electrodeCapacityType = CAPACITY_ANODE_ECT;
    }
  
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
EState::EState(std::string EState_type_string) :
    eRef_(nullptr),
    es_id_(),
    electrodeTypeString_(""),
    EST_lastFileRead_(EST_CSTR),
    EST_Version_lastFileRead_(1),
    EST_fileToBeWritten_(EST_CSTR),
    EST_Version_fileToBeWritten_(1),
    phaseVoltages_(0),
    temperature_(298.15),
    pressure_(OneAtm),
    electrodeChemistryModelType_(0),
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
    electronKmolDischargedToDate_(0.0),
    deltaTsubcycle_init_next_(1.0E300)
{
    es_id_.EState_Type_String = EState_type_string;
}
//======================================================================================================================
EState::EState(const EState& right) :
    eRef_(right.eRef_),
    es_id_(right.es_id_),
    electrodeTypeString_(""),
    EST_lastFileRead_(EST_CSTR),
    EST_Version_lastFileRead_(1),
    EST_fileToBeWritten_(EST_CSTR),
    EST_Version_fileToBeWritten_(1),
    phaseVoltages_(0),
    temperature_(298.15),
    pressure_(OneAtm),
    electrodeChemistryModelType_(0),
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
    electronKmolDischargedToDate_(0.0),
    deltaTsubcycle_init_next_(1.0E300)
{
    /*
     * Call the assignment operator.
     */
    EState::operator=(right);
}
//==================================================================================================================================
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
    es_id_                             = right.es_id_;
    electrodeTypeString_               = right.electrodeTypeString_;
    EST_lastFileRead_                  = right.EST_lastFileRead_;
    EST_Version_lastFileRead_          = right.EST_Version_lastFileRead_;
    spMoles_                           = right.spMoles_;
    phaseVoltages_                     = right.phaseVoltages_;
    temperature_                       = right.temperature_;
    pressure_                          = right.pressure_;
    electrodeChemistryModelType_       = right.electrodeChemistryModelType_;
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
    electronKmolDischargedToDate_      = right.electronKmolDischargedToDate_;

    deltaTsubcycle_init_next_          = right.deltaTsubcycle_init_next_;
    solnDot_                           = right.solnDot_;

    return *this;
}
//==================================================================================================================================
EState::~EState()
{
}
//==================================================================================================================================
EState* EState::duplMyselfAsEState(Electrode* e) const
{
    EState* es = new EState(*this);
    if (e) {
        es->eRef_ = e;
    }
    return es;
}
//==================================================================================================================================
int EState::initialize(const ZZCantera::Electrode* const e)
{
    eRef_ = e;

    es_id_.electrodeTypeString = Electrode_Types_Enum_to_string(e->electrodeType());
    es_id_.electrodeChemistryModelType = e->electrodeChemistryModelType_;
    es_id_.electrodeDomainNumber =  e->electrodeDomainNumber_;
    es_id_.electrodeCellNumber = e->electrodeCellNumber_;

    electrodeTypeString_    = Electrode_Types_Enum_to_string(e->electrodeType());
    electrodeChemistryModelType_  = e->electrodeChemistryModelType_;
    electrodeDomainNumber_  = e->electrodeDomainNumber_;
    electrodeCellNumber_    = e->electrodeCellNumber_;

    copyElectrode_intoState(eRef_);
    return 0;
}
//==================================================================================================================================
const std::string& EState::electrodeType() const
{
    return electrodeTypeString_;
}
//==================================================================================================================================
// Write the Identification to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node* EState::writeIdentificationToXML() const
{
    XML_Node* x = new XML_Node("ElectrodeIdentification");

    ZZctml::addNamedString(*x, "electrodeTypeString", electrodeTypeString_);
    ZZctml::addInteger(*x, "EState_Type",         EST_fileToBeWritten_);
    ZZctml::addNamedString(*x, "EState_Type_String", esmodel::EState_Type_Enum_to_string(EST_fileToBeWritten_));
    ZZctml::addInteger(*x, "fileVersionNumber", EST_Version_lastFileRead_);
    ZZctml::addInteger(*x, "electrodeModelType",  electrodeChemistryModelType_);
    ZZctml::addInteger(*x, "electrodeDomainNumber",  electrodeDomainNumber_);
    ZZctml::addInteger(*x, "electrodeCellNumber",  electrodeCellNumber_);
    if (electrodeCapacityType_ == CAPACITY_ANODE_ECT) {
	ZZctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Anode");
    } else if (electrodeCapacityType_ == CAPACITY_CATHODE_ECT) {
	ZZctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Cathode");
    } else {
	ZZctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Other");
    }
    x->addChildToTree( write_PhaseListID_ToXML());

    return x;
}
//======================================================================================================================
// Write the PhaseList ID to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node* EState::write_PhaseListID_ToXML() const
{
    XML_Node* x = new XML_Node("PhaseList_Id");
    ZZctml::addInteger(*x, "PhaseList_DomainNumber",  electrodeDomainNumber_);
    ZZctml::addInteger(*x, "PhaseList_CellNumber",  electrodeCellNumber_);
    ZZctml::addNamedString(*x, "PhaseList_Name",  eRef_->electrodeName_);
    ZZctml::addInteger(*x, "nPhases",  eRef_->nPhases());
    std::vector<std::string> v;
    for (size_t i = 0; i < eRef_->nPhases(); ++i) {
        v.push_back(eRef_->phaseIDString(i));
    }
    ZZctml::addTokenArray(*x, "phaseIDStrings", v);

    std::vector<int> iv;
    for (size_t i = 0; i < eRef_->nPhases(); ++i) {
        iv.push_back(eRef_->thermo(i).nSpecies());
    }
    ZZctml::addNamedIntegerArray(*x, "phaseNumSpecies", eRef_->nPhases(), &iv[0]);

    v.clear();
    for (size_t i = 0; i < eRef_->nGlobalSpecies(); ++i) {
        v.push_back(eRef_->speciesName(i));
    }
    ZZctml::addTokenArray(*x, "speciesNames", v);

    return x;
}
//==================================================================================================================================
// Write the electrodeState to an XML_Node tree
XML_Node* EState::write_electrodeState_ToXML() const
{
    XML_Node* x = new XML_Node("electrodeState");

    size_t nsp = spMoles_.size();
    ZZctml::addNamedFloatArray(*x, "spMoles", nsp, spMoles_.data(), "kmol");

    size_t np = phaseVoltages_.size();
    ZZctml::addNamedFloatArray(*x, "phaseVoltages", np, phaseVoltages_.data(), "volt");

    ZZctml::addFloat(*x, "temperature",  temperature_, "Kelvin");
    ZZctml::addFloat(*x, "pressure",  pressure_, "Pa");
    ZZctml::addInteger(*x, "electrodeModelType",  electrodeChemistryModelType_);
    ZZctml::addInteger(*x, "electrodeDomainNumber",  electrodeDomainNumber_);
    ZZctml::addInteger(*x, "electrodeCellNumber",  electrodeCellNumber_);
    ZZctml::addFloat(*x, "particleNumberToFollow",  particleNumberToFollow_, "");
    ZZctml::addFloat(*x, "electrodeSolidVolume", electrodeSolidVolume_, "m3");
    ZZctml::addFloat(*x, "grossVolume", grossVolume_, "m3");
    ZZctml::addFloat(*x, "radiusExterior", radiusExterior_, "m3");
    int ns = surfaceAreaRS_.size();
    ZZctml::addNamedFloatArray(*x, "surfaceAreaRS", ns, DATA_PTR(surfaceAreaRS_), "m2");
    ZZctml::addFloat(*x, "electrodeMoles", electrodeMoles_, "kmol");
    ZZctml::addInteger(*x, "electrodeCapacityType", (int) electrodeCapacityType_);
    ZZctml::addFloat(*x, "capacityLeft", capacityLeft_, "coulomb");
    ZZctml::addFloat(*x, "capacityInitial", capacityInitial_, "coulomb");
    ZZctml::addFloat(*x, "depthOfDischarge", depthOfDischarge_, "coulomb");
    ZZctml::addFloat(*x, "depthOfDischargeStarting", depthOfDischargeStarting_, "coulomb");
    ZZctml::addFloat(*x, "relativeElectronsDischargedPerMole", relativeElectronsDischargedPerMole_, "");
    ZZctml::addFloat(*x, "relativeDepthOfDischarge", relativeDepthOfDischarge_, "");
    ZZctml::addFloat(*x, "capacityDischargedToDate", capacityDischargedToDate_, "coulomb");
    ZZctml::addFloat(*x, "deltaTsubcycle_init_next", deltaTsubcycle_init_next_, "s");
    
    size_t neq = solnDot_.size();
    if (neq > 0) {
        ZZctml::addNamedFloatArray(*x, "solnDot", neq, solnDot_.data(), "");
    } 
    return x;
}
//==================================================================================================================================
void EState::readStateFromXMLRoot(const XML_Node& xmlRoot)
{
    if (xmlRoot.hasChild("electrodeState")) {
        XML_Node& x = xmlRoot.child("electrodeState");
        readStateFromXML(x);
    }
}
//==================================================================================================================================
void EState::readIdentificationFromXML(const XML_Node& xmlEI)
{
    std::string typeSS;
    std::string nn = xmlEI.name();
    const XML_Node* x = &xmlEI;
    if (nn != "ElectrodeIndentification") {
	x = xmlEI.findByName("ElectrodeIdentification");
	if (!x) {
	    throw Electrode_Error("EState::readIdentificationFromXM", "Could not find the XML node named ElectrodeIdentification");
	}
    }
    ZZctml::getNamedStringValue(*x, "electrodeTypeString", electrodeTypeString_ , typeSS);
    EST_lastFileRead_ = (EState_Type_Enum) ZZctml::getInteger(*x, "EState_Type");

    std::string EState_Type_String;
    if (x->hasChild("EState_Type_String")) {
	ZZctml::getNamedStringValue(*x, "EState_Type_String", EState_Type_String, typeSS);
	ZZCantera::EState_Type_Enum echeck = esmodel::string_to_EState_Type_Enum(EState_Type_String);
	if (echeck != EST_lastFileRead_ ) {
	    throw Electrode_Error("EState::readIdentificationFromXM",
				  "Incompatibility between EState_Type and EState_Type_String ");
	}
    } 
   
    if (x->hasChild("fileVersionNumber")) {
	EST_Version_lastFileRead_ = ZZctml::getInteger(*x, "fileVersionNumber");
    }
    if (x->hasChild("electrodeCapacityType")) {
	electrodeCapacityType_ = (ZZCantera::Electrode_Capacity_Type_Enum) ZZctml::getInteger(*x, "electrodeCapacityType");
    }
    if (x->hasChild("electrodeModelType")) {
	electrodeChemistryModelType_ = ZZctml::getInteger(*x, "electrodeModelType");
    } 
    if (x->hasChild("electrodeDomainNumber")) {
	electrodeDomainNumber_ = ZZctml::getInteger(*x, "electrodeDomainNumber");
    }
    if (x->hasChild("electrodeDomainNumber")) {
     electrodeCellNumber_ = ZZctml::getInteger(*x, "electrodeDomainNumber");
    }

    if (x->hasChild("electrodeCapacityType")) {
	ZZctml::getNamedStringValue(*x, "electrodeCapacityType", nn, typeSS);
	if (nn == "Capacity_Anode") {
	    electrodeCapacityType_ = CAPACITY_ANODE_ECT;
	} else if (nn == "Capacity_Cathode") {
	    electrodeCapacityType_ = CAPACITY_CATHODE_ECT;
	} else if (nn == "Capacity_Other") {
	    electrodeCapacityType_ = CAPACITY_OTHER_ECT;
	}
    }

}
//==================================================================================================================================
// Read identification information from a struct EState_Identification object
void EState::readIdentificationFromStruct(const EState_ID_struct& es_ID)
{
    // @todo do checking on compatible values
    electrodeTypeString_ = es_ID.electrodeTypeString;
    
    EST_lastFileRead_ = es_ID.EST_Type;
    EST_Version_lastFileRead_ =  es_ID.EST_Version;
    electrodeCapacityType_ = es_ID.electrodeCapacityType;
    electrodeChemistryModelType_ = es_ID.electrodeChemistryModelType; 
    electrodeDomainNumber_ = es_ID.electrodeDomainNumber;
    electrodeCellNumber_ = es_ID.electrodeCellNumber;
}
//==================================================================================================================================
void EState::readStateFromXML(const XML_Node& xmlEState)
{
    std::string nodeName = xmlEState.name();
    if (nodeName != "electrodeState") {
	throw Electrode_Error(" EState::readStateFromXML",
			      " Name of the xml node should have been electrodeState. Instead it was " + nodeName);
    }
    ZZctml::getFloatArray(xmlEState, spMoles_, true, "", "spMoles");
    ZZctml::getFloatArray(xmlEState, phaseVoltages_, true, "volts", "phaseVoltages");
    temperature_ = ZZctml::getFloat(xmlEState, "temperature", "toSI");
    pressure_ = ZZctml::getFloat(xmlEState, "pressure", "toSI");
    electrodeChemistryModelType_ = ZZctml::getInteger(xmlEState, "electrodeModelType");
    electrodeDomainNumber_ = ZZctml::getInteger(xmlEState, "electrodeDomainNumber");
    electrodeCellNumber_ = ZZctml::getInteger(xmlEState, "electrodeCellNumber");
    particleNumberToFollow_ = ZZctml::getFloat(xmlEState, "particleNumberToFollow", "toSI");
    electrodeSolidVolume_ = ZZctml::getFloat(xmlEState, "electrodeSolidVolume", "toSI");
    grossVolume_ = ZZctml::getFloat(xmlEState, "grossVolume", "toSI");
    radiusExterior_ = ZZctml::getFloat(xmlEState, "radiusExterior", "toSI");
    ZZctml::getFloatArray(xmlEState, surfaceAreaRS_, true, "m2", "surfaceAreaRS");
    electrodeMoles_ = ZZctml::getFloat(xmlEState, "electrodeMoles", "toSI");
    electrodeCapacityType_ = (ZZCantera::Electrode_Capacity_Type_Enum) ZZctml::getInteger(xmlEState, "electrodeCapacityType");
    capacityLeft_ = ZZctml::getFloat(xmlEState, "capacityLeft", "toSI");
    capacityInitial_ = ZZctml::getFloat(xmlEState, "capacityInitial", "toSI");
    depthOfDischarge_ = ZZctml::getFloat(xmlEState, "depthOfDischarge", "toSI");
    depthOfDischargeStarting_ = ZZctml::getFloat(xmlEState, "depthOfDischargeStarting", "toSI");
    relativeElectronsDischargedPerMole_ = ZZctml::getFloat(xmlEState, "relativeElectronsDischargedPerMole", "toSI");
    relativeDepthOfDischarge_ = ZZctml::getFloat(xmlEState, "relativeDepthOfDischarge", "toSI");
    capacityDischargedToDate_ = ZZctml::getFloat(xmlEState, "capacityDischargedToDate", "toSI");
    electronKmolDischargedToDate_ =  capacityDischargedToDate_ / ZZCantera::Faraday;
    if (electrodeCapacityType_ == CAPACITY_CATHODE_ECT) {
	electronKmolDischargedToDate_ *= -1.0;
    }

    deltaTsubcycle_init_next_ = ZZctml::getFloat(xmlEState, "deltaTsubcycle_init_next", "toSI");

    if (xmlEState.hasChild("solnDot")) {
        ZZctml::getFloatArray(xmlEState, solnDot_, true, "", "solnDot");
    }
}
//==================================================================================================================================
// Set the State of this object from the state of the Electrode object
void EState::copyElectrode_intoState(const Electrode* const e, bool doFinal)
{
    eRef_                              = e;
    spMoles_                           = e->spMoles_final_;
    phaseVoltages_                     = e->phaseVoltages_;
    temperature_                       = e->temperature_;
    pressure_                          = e->pressure_;
    electrodeChemistryModelType_       = e->electrodeChemistryModelType_;
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

    electronKmolDischargedToDate_      = e->electronKmolDischargedToDate_;

    deltaTsubcycle_init_next_          = e->deltaTsubcycle_init_next_;

    double capUsed = capacityInitial_ - capacityLeft_;

    double relCapUsed = 1.0;
    if (electrodeMoles_ > 1.0E-200) {
        relCapUsed = capUsed / electrodeMoles_;
    } 

    const Electrode_MP_RxnExtent* emp = dynamic_cast<const Electrode_MP_RxnExtent*>(e);
    if (emp) {
        relativeElectronsDischargedPerMole_ = emp->RelativeExtentRxn_final_;
    } else {
        const Electrode_CSTR* ecstr = dynamic_cast<const Electrode_CSTR*>(e);
        if (ecstr) {
            relativeElectronsDischargedPerMole_ = ecstr->RelativeExtentRxn_final_;
        } else {
            relativeElectronsDischargedPerMole_ = relCapUsed;
        }
    }

    const Electrode_Integrator* ei = dynamic_cast<const Electrode_Integrator*>(e);
    if (ei) {
        if (doFinal) {
            solnDot_                  = ei->solnDot_final_;
        } else {
            solnDot_                  = ei->solnDot_init_;
        }
    }
}
//======================================================================================================================
//Set the state of the Electrode from the state of this object -> both init and final and init_init and final_final
void EState::setStateElectrode_fromEState(Electrode* const e) const
{
    EState::copyEState_toElectrode(e);

    e->stateToPhaseFlagsReconciliation(false);
    e->updateState();
    e->setInitStateFromFinal(true);
}
//======================================================================================================================
// Set the state of the Electrode Class from the state of the EState object
//   @deprecated
void EState::copyEState_toElectrode(Electrode* const e) const
{
    e->spMoles_final_                     = spMoles_;
    e->phaseVoltages_                     = phaseVoltages_;
    e->temperature_                       = temperature_;
    e->pressure_                          = pressure_;
    e->electrodeChemistryModelType_       = electrodeChemistryModelType_;
    e->electrodeDomainNumber_             = electrodeDomainNumber_;
    e->electrodeCellNumber_               = electrodeCellNumber_;
    e->particleNumberToFollow_            = particleNumberToFollow_;
    e->ElectrodeSolidVolume_              = electrodeSolidVolume_;
    // grossVolume_
    e->Radius_exterior_final_             = radiusExterior_;
    e->surfaceAreaRS_final_               = surfaceAreaRS_;
    // electrodeMoles_
    e->setCapacityType(electrodeCapacityType_);
    // capacityLeft_ -> ok no explicit storage of this quantity in Electrode object
    e->capacityInitialZeroDod_            = capacityInitial_;
    // depthOfDischarge_  -> ok no explicit storage of this quantity in Electrode object
    e->depthOfDischargeStarting_          = depthOfDischargeStarting_;
    // relativeElectronsDischargedPerMole_
    // relativeDeptOfDischarge_
    // capacityDischargedToDate_
    // e->electronKmolDischargedToDate_      = capacityDischargedToDate_ / ZZCantera::Faraday;
    e->electronKmolDischargedToDate_      = electronKmolDischargedToDate_;

    for (size_t iph = 0; iph < e->m_NumTotPhases; iph++) {
        e->updateState_Phase(iph);
    }

    e->deltaTsubcycle_init_next_          = deltaTsubcycle_init_next_;
    e->deltaTsubcycle_init_init_          = deltaTsubcycle_init_next_;

    Electrode_Integrator* const ei = dynamic_cast<Electrode_Integrator* const>(e);
    if (ei->solnDot_final_.size() == solnDot_.size()) {
        ei->solnDot_final_ = solnDot_;
        ei->solnDot_init_ = solnDot_;
        ei->solnDot_init_init_ = solnDot_;
    }

    Electrode_MP_RxnExtent* const emp = dynamic_cast<Electrode_MP_RxnExtent* const>(e);
    if (emp) {
        emp->RelativeExtentRxn_final_ = relativeElectronsDischargedPerMole_;
    } else {
        Electrode_CSTR* const ecstr = dynamic_cast<Electrode_CSTR* const>(e);
        if (ecstr) {
            ecstr->RelativeExtentRxn_final_ = relativeElectronsDischargedPerMole_;
        }
    }
}
//=================================================================================================================================
int EState::printHead(int printLvl) const
{
    static bool printHead = false;
    static int num = 0;
    if (printLvl == -1) {
        printHead = false;
        num = 0;
    } else {
    if (printLvl >= 2) {
        if (!printHead) {
            printf("\t\tEState:: FAILURE: differing quantities\n");
            printf("\t\t   Quantity           Index  Value            Guest_Value \n");
            printHead = true;
        }
    }
    num++; 
    }
    return num;
}
//=================================================================================================================================
// print a difference between two strings
void EState::printDiff(const std::string& vexp, bool significant,  const std::string& val, 
                       const std::string& gval, int printLvl) const
{
    std::string istr = "no ";
    if (significant) {
	istr = "yes";
    }
    int num = printHead(printLvl);
    if ((printLvl == 2 && num < 100) || (printLvl >= 3) ) {
	printf("\t\t   %-15.15s  %-3.3s   ", vexp.c_str(), istr.c_str());
	printf("%17s %17s\n", val.c_str(), gval.c_str());
    }
}
//==================================================================================================================================
void EState::printDiff(const std::string& vexp, int index, int val, int gval, int printLvl) const
{
    std::string istring = "   ";
    if (index >= 0) {
	istring = int2str(index, "%3d");
    }
    int num = printHead(printLvl);
    if (printLvl >= 2) {
        if ( (printLvl == 2 && num < 100) || (printLvl >= 3) ) {
	    printf("\t\t   %-15.15s  %-3.3s  ", vexp.c_str(), istring.c_str() );
	    if (val != MDP_INT_NOINIT) {
		printf("%15d ", val);
	    } else {
		printf("- NotAvail -    ");
	    }
	    if (gval != MDP_INT_NOINIT) {
		printf("%15d", gval);
	    } else {
		printf("-NotAvail-     ");
	    }
	    printf("\n");	
	}
    }
}
//==================================================================================================================================
void EState::printDiff(const std::string& vexp, int index, double val, double gval, int printLvl) const
{
    std::string istring = "   ";
    if (index >= 0) {
	istring = int2str(index, "%3d");
    }
    int num = printHead(printLvl);    
    if (printLvl >= 2) {
        if ( (printLvl == 2 && num < 15) || (printLvl >= 3) ) {
	    printf("\t\t   %-15.15s  %-3.3s    ",
		   vexp.c_str(), istring.c_str() );
	    if (val != MDP_DBL_NOINIT) {
		printf("%-15.7E ", val);
	    } else {
		printf("- NotAvail -    ");
	    }
	    if (gval != MDP_DBL_NOINIT) {
		printf("%15.7E", gval);
	    } else {
		printf("-NotAvail-     ");
	    }
	    printf("\n");	
	}
    }
}
//==================================================================================================================================
void EState::printVecDiff(const std::string& vexp, const std::vector<double>& val, const std::vector<double>& gval,
			  int printLvl) const
{
    std::string istring = "   ";
    size_t j1 = val.size();
    size_t j2 = gval.size();
    size_t jmax = std::max(j1, j2);
    printHead(printLvl);
    if (printLvl >= 2) {
	for (size_t j = 0; j < jmax; ++j) {
	    if (j > j2) {
		printDiff(vexp, j, val[j], MDP_DBL_NOINIT, printLvl);
	    } else if (j > j1) {
		printDiff(vexp, j, MDP_DBL_NOINIT, gval[j], printLvl);
	    } else {
		printDiff(vexp, j, val[j], gval[j], printLvl);
	    }
	
	}
    }
}
//==================================================================================================================================
/*
 *    printLvl settings
 *              All printing starts two tabs in
 *              0  Absolutely no printing is done
 *              1  One line is printed
 *              2  A heading and about 10 lines are printed
 *              3  A heading and unlimited lines are printed
 *
 */
bool EState::compareOtherState(const EState* const ESguest, double molarAtol, int nDigits, bool includeHist, int printLvl) const
{
    bool btotal = true;
    bool boolR = true;
    printHead(-1);
 
    if (printLvl > 1) {
	printf("\t\tEState::compareOtherState() Start comparison of types %s vs. %s\n",electrodeTypeString_.c_str(), 
	       ESguest->electrodeTypeString_.c_str());
    }

    // electrodeTypeString_
    /*
     *    Differences in electrodeTypeString_ don't cause errors. We can solve the problem with different models and 
     *    get the same answer.
     */
    boolR = (electrodeTypeString_ == ESguest->electrodeTypeString_);
    if (!boolR) {
	printDiff("electrodeTypeString_", false, electrodeTypeString_, ESguest->electrodeTypeString_, printLvl);
    }

    // EST_lastFileRead_  The file type read
    /*
     *  Presumably differences in this are fatal
     */
    boolR = (EST_lastFileRead_ == ESguest->EST_lastFileRead_);
    if (!boolR) {
	printDiff("EST_lastFileRead_", true, esmodel::EState_Type_Enum_to_string(EST_lastFileRead_),
		  esmodel::EState_Type_Enum_to_string(ESguest->EST_lastFileRead_), printLvl);
    }
    btotal = boolR && btotal;

    // EST_fileToBeWritten_  The file type that EState will write out
    /*
     *  Presumably differences in this not significant
     */
    boolR = (EST_fileToBeWritten_ == ESguest->EST_fileToBeWritten_);
    if (!boolR) {
	printDiff("EST_fileToBeWritten_", false, esmodel::EState_Type_Enum_to_string(EST_fileToBeWritten_),
		  esmodel::EState_Type_Enum_to_string(ESguest->EST_fileToBeWritten_), printLvl);
    }

    // EST_Version_fileToBeWritten_  The file type version that EState will write out
    /*
     *  Presumably differences in this not significant
     */
    boolR = (EST_Version_fileToBeWritten_ == ESguest->EST_Version_fileToBeWritten_);
    if (!boolR) {
	printDiff("EST_fileToBeWritten_", false, EST_Version_fileToBeWritten_ ,ESguest->EST_Version_fileToBeWritten_, printLvl);
    }

    // temperature
    boolR = esmodel::doubleEqual(temperature_, ESguest->temperature_, 0.0, nDigits);
    if (!boolR) {
	printDiff("Temperature", -1, temperature_, ESguest->temperature_,printLvl);
    }
    btotal = boolR && btotal;

    // Compare pressure
    boolR = esmodel::doubleEqual(pressure_, ESguest->pressure_, 0.0, nDigits);
    if (!boolR) {
	printDiff("Pressure", -1, pressure_, ESguest->pressure_, printLvl);
    }
    btotal = boolR && btotal;

    // Compare gross Volumes  
    // based on 55 kmol / m3 ( water)
    double volAtol =  molarAtol / 55.;
    boolR = esmodel::doubleEqual(grossVolume_, ESguest->grossVolume_, volAtol, nDigits);
    if (!boolR) {
	printDiff("grossVolume", -1,  grossVolume_, ESguest->grossVolume_, printLvl);
    }
    btotal = boolR && btotal;

    // Compare exterior radius
    double radiusAtol = pow(volAtol, 0.3333);
    boolR = esmodel::doubleEqual(radiusExterior_, ESguest->radiusExterior_, radiusAtol, nDigits);
    if (!boolR) {
	printDiff("radiusExterior_", -1, radiusExterior_, ESguest->radiusExterior_, printLvl);
    }
    btotal = boolR && btotal;

    // Compare surface area vector
    double surfaceAtol = 12 * radiusAtol * radiusAtol;
    boolR = esmodel::doubleVectorEqual(surfaceAreaRS_, ESguest->surfaceAreaRS_, surfaceAtol, nDigits);
    if (!boolR) {
	printVecDiff("surfaceAreaRS_", surfaceAreaRS_, ESguest->surfaceAreaRS_, printLvl);
    }
    btotal = boolR && btotal;
	
    
    // Compare spMoles vector
    boolR = esmodel::doubleVectorEqual(spMoles_, ESguest->spMoles_, molarAtol, nDigits);
    if (!boolR) {
	printVecDiff("speciesMoles", spMoles_, ESguest->spMoles_,  printLvl);
    }
    btotal = boolR && btotal;

 
    // Compare phase voltages
    boolR = esmodel::doubleVectorEqual(phaseVoltages_, ESguest->phaseVoltages_, molarAtol, nDigits);
    if (!boolR) {
	printVecDiff("phaseVoltages", phaseVoltages_, ESguest->phaseVoltages_,  printLvl);
    }
    btotal = boolR && btotal;

    // Electrode Model Chemistry Type
    /*
     *   Specific chemistry models
     */
    boolR = (electrodeChemistryModelType_ == ESguest->electrodeChemistryModelType_);
    if (!boolR) {
	printDiff("electrodeChemistryModelType_", -1, electrodeChemistryModelType_, ESguest->electrodeChemistryModelType_, 
		  printLvl);
    }
    btotal = boolR && btotal;

    // Compare particle number
    double numberAtol = volAtol / (4.0E-18);
    boolR = esmodel::doubleEqual(particleNumberToFollow_, ESguest->particleNumberToFollow_, numberAtol, nDigits);
    if (!boolR) {
	printDiff("particleNumberToFollow_", -1, particleNumberToFollow_, ESguest->particleNumberToFollow_,  printLvl);
    }
    btotal = boolR && btotal;

    // compare the electrode solid volume
    boolR = esmodel::doubleEqual(electrodeSolidVolume_, ESguest->electrodeSolidVolume_, volAtol, nDigits);
    if (!boolR) {
	printDiff("electrodeSolidVolume", -1, electrodeSolidVolume_, ESguest->electrodeSolidVolume_, printLvl);
    }
    btotal = boolR && btotal;


    // compare the total solid moles
    boolR = esmodel::doubleEqual(electrodeMoles_, ESguest->electrodeMoles_, molarAtol, nDigits);
    if (!boolR) {
	printDiff("electrodeMoles_", -1, electrodeMoles_, ESguest->electrodeMoles_, printLvl);
    }
    btotal = boolR && btotal;

    // Compare the electrode capacity type
    boolR = (electrodeCapacityType_ == ESguest->electrodeCapacityType_);
    if (!boolR) {
	printDiff("electrodeMoles_", -1, int(electrodeCapacityType_), int(ESguest->electrodeCapacityType_), printLvl);
    }
    btotal = boolR && btotal;

    // Compare the capacity Left
    double capAtol = molarAtol * Faraday;
    boolR = esmodel::doubleEqual(capacityLeft_, ESguest->capacityLeft_, capAtol, nDigits);
    if (!boolR) {
	printDiff("capacityLeft_", -1, capacityLeft_, ESguest->capacityLeft_, printLvl);
    }
    btotal = boolR && btotal;

    // Compare the capacity initial
    boolR = esmodel::doubleEqual(capacityInitial_, ESguest->capacityInitial_, capAtol, nDigits);
    if (!boolR) {
	printDiff("capacityInitial_", -1, capacityInitial_, ESguest->capacityInitial_, printLvl);
    }
    btotal = boolR && btotal;

    // Compare the depth of discharge
    boolR = esmodel::doubleEqual(depthOfDischarge_, ESguest->depthOfDischarge_, capAtol, nDigits);
    if (!boolR) {
	printDiff("depthOfDischarge_", -1, depthOfDischarge_, ESguest->depthOfDischarge_, printLvl);
    }
    btotal = boolR && btotal;

    // Compare the relative discharge
    double unitlessAtol = 1.0E-10;
    if (electrodeMoles_ > 1.0E-100) {
	unitlessAtol = molarAtol / electrodeMoles_;
    }
    boolR = esmodel::doubleEqual(relativeDepthOfDischarge_, ESguest->relativeDepthOfDischarge_, unitlessAtol, nDigits);
    if (!boolR) {
	printDiff("relativeDepthOfDischarge_", -1, relativeDepthOfDischarge_, ESguest->relativeDepthOfDischarge_, printLvl);
    }
    btotal = boolR && btotal;

     
    // Compare the capacity discharged to date -> this is a number that is dependent on the past time history of the simulation
    boolR = esmodel::doubleEqual(capacityDischargedToDate_, ESguest->capacityDischargedToDate_, capAtol, nDigits);
    if (!boolR) {
	printDiff("capacityDischargedToDate_", -1, capacityDischargedToDate_, ESguest->capacityDischargedToDate_, printLvl);
    }
    if (includeHist) {
	btotal = boolR && btotal;
    }

    // Compare the electron moles discharged to date -> this is a number that is dependent on the past time history of the simulation
    boolR = esmodel::doubleEqual(electronKmolDischargedToDate_, ESguest->electronKmolDischargedToDate_, molarAtol, nDigits);
    if (!boolR) {
	printDiff("electronKmolDischargedToDate_", -1, electronKmolDischargedToDate_, ESguest->electronKmolDischargedToDate_,
		  printLvl);
    }
    if (includeHist) {
	btotal = boolR && btotal;
    }

    // Compare the next deltaT -> this is a number that is dependent on the past time history of the simulation
    boolR = esmodel::doubleEqual(deltaTsubcycle_init_next_, ESguest->deltaTsubcycle_init_next_, unitlessAtol, nDigits);
    if (!boolR) {
	printDiff("deltaTsubcycle_init_next_", -1, deltaTsubcycle_init_next_, ESguest->deltaTsubcycle_init_next_,
		  printLvl);
    }
    if (includeHist) {
	btotal = boolR && btotal;
    }

    if (!btotal) {
	if (printLvl == 1) {
	    printf("\t\tEState:: ERROR:   States are not the same\n");
        }
    } else {
	if (printLvl >= 1) {
	    printf("\t\tEState:: SUCCESS: States are the same\n");
        }
    }
    return btotal;
}
//==================================================================================================================================
double EState::electrodeMoles() const
{
    return electrodeMoles_;
}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------


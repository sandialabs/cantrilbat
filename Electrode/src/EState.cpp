/*
 * $Id: EState.cpp 584 2013-04-03 00:29:38Z hkmoffa $
 */

#include "EState.h"
#include "Electrode.h"
#include "Electrode_MP_RxnExtent.h"
#include "Electrode_Factory.h"
#include "EState_XML.h"

#include "mdp_allo.h"

#include <string>

namespace Cantera
{
//==================================================================================================================================
EState_Identification::EState_Identification() :
    electrodeTypeString_(""),
    EST_Type_(EST_UNKNOWN_TYPE),
    EState_Type_String_("Unknown"),
    EST_Version_(1),
    electrodeChemistryModelType_(0),
    electrodeDomainNumber_(0),
    electrodeCellNumber_(0),
    electrodeCapacityType_(CAPACITY_ANODE_ECT)
{
}
//==================================================================================================================================
Cantera::XML_Node* EState_Identification::writeIdentificationToXML() const
{
    XML_Node* x = new XML_Node("ElectrodeIdentification");
    ctml::addNamedString(*x, "electrodeTypeString", electrodeTypeString_);
    ctml::addInteger(*x, "EState_Type",         EST_Type_);
    ctml::addNamedString(*x, "EState_Type_String", EState_Type_String_);
    ctml::addInteger(*x, "fileVersionNumber",  EST_Version_);
    ctml::addInteger(*x, "electrodeModelType",  electrodeChemistryModelType_);
    ctml::addInteger(*x, "electrodeDomainNumber",  electrodeDomainNumber_);
    ctml::addInteger(*x, "electrodeCellNumber",  electrodeCellNumber_);
    if (electrodeCapacityType_ == CAPACITY_ANODE_ECT) {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Anode");
    } else if (electrodeCapacityType_ == CAPACITY_CATHODE_ECT) {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Cathode");
    } else {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Other");
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

    ctml::getNamedStringValue(*x, "electrodeTypeString", electrodeTypeString_ , typeSS);
    EST_Type_ = (EState_Type_Enum)  ctml::getInteger(*x, "EState_Type");
    if (x->hasChild("EState_Type_String")) {
	ctml::getNamedStringValue(*x, "EState_Type_String", EState_Type_String_ , typeSS);
    } else {
        EState_Type_String_ = "EState_CSTR";
    }
    if (x->hasChild("fileVersionNumber")) {
	EST_Version_ = ctml::getInteger(*x, "fileVersionNumber");
    }
    if (x->hasChild("electrodeModelType")) {
        electrodeChemistryModelType_ = ctml::getInteger(*x, "electrodeModelType");
    }
    if (x->hasChild("electrodeDomainNumber")) {
        electrodeDomainNumber_ = ctml::getInteger(*x, "electrodeDomainNumber");
    }
    if (x->hasChild("electrodeDomainNumber")) {
     electrodeCellNumber_ = ctml::getInteger(*x, "electrodeDomainNumber");
    }
 
    if (x->hasChild("electrodeCapacityType")) {
	ctml::getNamedStringValue(*x, "electrodeCapacityType", nn, typeSS);
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
//===================================================================================================================================
int EState::initialize(const Cantera::Electrode* const e)
{
    eRef_ = e;
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
// Write the ElectrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node*   EState::writeIdentificationToXML() const
{
    XML_Node* x = new XML_Node("ElectrodeIdentification");

    ctml::addNamedString(*x, "electrodeTypeString", electrodeTypeString_);
    ctml::addInteger(*x, "EState_Type",         EST_fileToBeWritten_);
    ctml::addNamedString(*x, "EState_Type_String", esmodel::EState_Type_Enum_to_string(EST_fileToBeWritten_));
    ctml::addInteger(*x, "fileVersionNumber", EST_Version_lastFileRead_);
    ctml::addInteger(*x, "electrodeModelType",  electrodeChemistryModelType_);
    ctml::addInteger(*x, "electrodeDomainNumber",  electrodeDomainNumber_);
    ctml::addInteger(*x, "electrodeCellNumber",  electrodeCellNumber_);
    if (electrodeCapacityType_ == CAPACITY_ANODE_ECT) {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Anode");
    } else if (electrodeCapacityType_ == CAPACITY_CATHODE_ECT) {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Cathode");
    } else {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Other");
    }
    return x;
}
//======================================================================================================================
// Write the electrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node* EState::write_electrodeState_ToXML() const
{
    XML_Node* x = new XML_Node("electrodeState");

    int nsp = spMoles_.size();
    ctml::addNamedFloatArray(*x, "spMoles", nsp, DATA_PTR(spMoles_), "kmol");

    int np = phaseVoltages_.size();
    ctml::addNamedFloatArray(*x, "phaseVoltages", np, DATA_PTR(phaseVoltages_), "volt");

    ctml::addFloat(*x, "temperature",  temperature_, "Kelvin");
    ctml::addFloat(*x, "pressure",  pressure_, "Pa");
    ctml::addInteger(*x, "electrodeModelType",  electrodeChemistryModelType_);
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
    ctml::addFloat(*x, "electronKmolDischargedToDate", electronKmolDischargedToDate_, "kmol");
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
void EState::readIdentificationFromXML(const XML_Node& xmlEI)
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

    ctml::getNamedStringValue(*x, "electrodeTypeString", electrodeTypeString_ , typeSS);

    EST_lastFileRead_ = (EState_Type_Enum)  ctml::getInteger(*x, "EState_Type");

    std::string EState_Type_String;
    if (x->hasChild("EState_Type_String")) {

	ctml::getNamedStringValue(*x, "EState_Type_String", EState_Type_String, typeSS);
	Cantera::EState_Type_Enum echeck = esmodel::string_to_EState_Type_Enum(EState_Type_String);
	if (echeck !=  EST_lastFileRead_ ) {
	    throw Electrode_Error("EState::readIdentificationFromXM",
				  "Incompatibility between EState_Type and EState_Type_String ");
	}
    } 
   
    if (x->hasChild("fileVersionNumber")) {
	EST_Version_lastFileRead_ = ctml::getInteger(*x, "fileVersionNumber");
    }
   
    if (x->hasChild("electrodeCapacityType")) {
	electrodeCapacityType_ = (Cantera::Electrode_Capacity_Type_Enum) ctml::getInteger(*x, "electrodeCapacityType");
    }
    if (x->hasChild("electrodeModelType")) {
	electrodeChemistryModelType_ = ctml::getInteger(*x, "electrodeModelType");
    } 
    if (x->hasChild("electrodeDomainNumber")) {
	electrodeDomainNumber_ = ctml::getInteger(*x, "electrodeDomainNumber");
    }
    if (x->hasChild("electrodeDomainNumber")) {
     electrodeCellNumber_ = ctml::getInteger(*x, "electrodeDomainNumber");
    }

    if (x->hasChild("electrodeCapacityType")) {
	ctml::getNamedStringValue(*x, "electrodeCapacityType", nn, typeSS);
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
/*
 *
 */
void EState::readIdentificationFromStruct(const EState_ID_struct& es_ID)
{
    // @todo do checking on compatible values
    electrodeTypeString_ = es_ID.electrodeTypeString_;
    
    EST_lastFileRead_ = es_ID.EST_Type_;
    EST_Version_lastFileRead_ =  es_ID.EST_Version_;
    electrodeCapacityType_ = es_ID.electrodeCapacityType_;
    electrodeChemistryModelType_ = es_ID.electrodeChemistryModelType_; 
    electrodeDomainNumber_ = es_ID.electrodeDomainNumber_;
    electrodeCellNumber_ = es_ID.electrodeCellNumber_;
}
//======================================================================================================================
//  Read the  electrode state from  an XML_Node tree
/*
 *   
 */
void EState::readStateFromXML(const XML_Node& xmlEState)
{
    std::string nodeName = xmlEState.name();
    if (nodeName != "electrodeState") {
	throw Electrode_Error(" EState::readStateFromXML",
			      " Name of the xml node should have been electrodeState. Instead it was " + nodeName);
    }
    ctml::getFloatArray(xmlEState, spMoles_, true, "", "spMoles");
    ctml::getFloatArray(xmlEState, phaseVoltages_, true, "volts", "phaseVoltages");
    temperature_ = ctml::getFloat(xmlEState, "temperature", "toSI");
    pressure_ = ctml::getFloat(xmlEState, "pressure", "toSI");
    electrodeChemistryModelType_ = ctml::getInteger(xmlEState, "electrodeModelType");
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
    if (xmlEState.hasChild("electronKmolDischargedToDate")) {
	electronKmolDischargedToDate_ = ctml::getFloat(xmlEState, "electronKmolDischargedToDate", "toSI");
    } else {
	electronKmolDischargedToDate_ =  capacityDischargedToDate_ / Cantera::Faraday;
    }
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
    e->electrodeChemistryModelType_                = electrodeChemistryModelType_;
    e->electrodeDomainNumber_             = electrodeDomainNumber_;
    e->electrodeCellNumber_               = electrodeCellNumber_;
    e->particleNumberToFollow_            = particleNumberToFollow_;
    e->ElectrodeSolidVolume_              = electrodeSolidVolume_;
    e->Radius_exterior_final_             = radiusExterior_;
    e->surfaceAreaRS_final_               = surfaceAreaRS_;
    e->depthOfDischargeStarting_          = depthOfDischargeStarting_;
    e->electronKmolDischargedToDate_      = capacityDischargedToDate_ / Cantera::Faraday;
    e->electronKmolDischargedToDate_      = electronKmolDischargedToDate_;
    e->setCapacityType(electrodeCapacityType_);

    for (size_t iph = 0; iph < e->m_NumTotPhases; iph++) {
        e->updateState_Phase(iph);
    }

    e->deltaTsubcycle_init_next_          = deltaTsubcycle_init_next_;
    e->deltaTsubcycle_init_init_          = deltaTsubcycle_init_next_;


    Electrode_MP_RxnExtent* const emp = dynamic_cast<Electrode_MP_RxnExtent* const>(e);
    if (emp) {
        emp->RelativeExtentRxn_final_ = relativeElectronsDischargedPerMole_;
    }
}
//==================================================================================================================================
/*
static double setAtolEm40(double a1, double a2)
{
    double m1 = std::max(a1, a2);
    double m2 = std::max(m1, 1.0E-40);
    return m2;
}
*/
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
//! print a difference between two strings
void EState::printDiff(const std::string& vexp, bool significant,  const std::string& val, const std::string& gval, int printLvl) const
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
	    printf("\t\t   %-15.15s  %-3.3s  ",
		   vexp.c_str(), istring.c_str() );
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
// Create a new Electrode Object
/*
 * @param model  String to look up the model against
 *
 * @return    Returns a pointer to a new Electrode instance matching the definition within EState object. Returns NULL if
 *            something went wrong. Throws an exception if the electrodeType isn't covered.
 *
 *  - Can't do this because you need an underlying PhaseList object to have been formed.
 */
/*
Electrode* newElectrodeObject(const Cantera::EState& es, double currentTime, Cantera::Electrode_Factory* f)
{
    if (f == 0) {
        f = Electrode_Factory::factory();
    }

    const std::string smodel = es.electrodeType();

    Electrode* e = f->newElectrodeObject(smodel);

    es.copyEState_toElectrode(e);
    e->setTime(currentTime);

    return e;
}
*/
//==================================================================================================================================
} // End of namespace Cantera
//==================================================================================================================================

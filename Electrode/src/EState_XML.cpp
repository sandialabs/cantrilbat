/**
 *
 */

#include "EState_XML.h"

#include "cantera/base/xml.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include "Electrode_Exception.h"
#include "EState.h"
#include "EState_RadialDistrib.h"

#include "mdp_stringUtils.h"

#include <cstdio>
#include <fstream>

using namespace Cantera;
using namespace std;

//! esmodel stands for Electrode source model. It will expand to encompass everything in this directory
namespace esmodel
{
//==================================================================================================================================
Map_ESEnum_String gMap_ESEnum_String;

//==================================================================================================================================
static void create_string_maps()
{
    if (gMap_ESEnum_String.string_maps_created) {
        return;
    }
    std::map<Cantera::EState_Type_Enum, std::string>& estate_types_string = gMap_ESEnum_String.estate_types_string;
    std::map<std::string , EState_Type_Enum>& string_estate_types= gMap_ESEnum_String.string_estate_types;

    gMap_ESEnum_String.string_maps_created = true;
    estate_types_string[EST_CSTR]                       =  "EState_CSTR";
    estate_types_string[EST_MULTIPLATEAU]               =  "EState_MultiPlateau";
    estate_types_string[EST_RADIALDISTRIB]              =  "EState_RadialDistrib";

    // Invert the maps automatically.
    for (std::map<EState_Type_Enum, std::string>::iterator pos = estate_types_string.begin();
            pos != estate_types_string.end(); ++pos) {
        string_estate_types[pos->second] = pos->first;
        std::string lll =  Cantera::lowercase(pos->second);
        string_estate_types[lll] = pos->first;
    }

}
//=================================================================================================================================
Cantera::XML_Node* getElectrodeOutputFile(const std::string& fileName, int index)
{
    /*
     * Call the basic Cantera routine that reads the XML file. If the file isn't found a NULL is returned 
     * causing this program to also return NULL
     */
    Cantera::XML_Node* xSavedSoln = get_XML_File(fileName);
    if (!xSavedSoln) {
        return xSavedSoln;
    }
    /*
     *  Find the electrodeOutput XML element.
     *   
     */

    XML_Node* eRecord = xSavedSoln->findNameIDIndex("electrodeOutput", "", index);
    if (!eRecord) {
        ESModel_Warning("esmodel::getElectrodeOutputFile()",
	                "Could not find a node named electrodeOutput with index " + mdpUtil::int2str(index));
    }
    return eRecord;
}
//==================================================================================================================================
std::string EState_Type_Enum_to_string(const EState_Type_Enum& etype)
{
    if (!gMap_ESEnum_String.string_maps_created) {
        create_string_maps();
    }
    std::map<EState_Type_Enum, std::string>& estate_types_string = gMap_ESEnum_String.estate_types_string;

    std::map<EState_Type_Enum, std::string>::iterator pos = estate_types_string.find(etype);
    if (pos == estate_types_string.end()) {
        return "UnknownEStateType";
    }
    return pos->second;
}
//==================================================================================================================================
EState_Type_Enum string_to_EState_Type_Enum(const std::string& input_string)
{
    if (!gMap_ESEnum_String.string_maps_created) {
        create_string_maps();
    }
    std::map<std::string , EState_Type_Enum>& string_estate_types = gMap_ESEnum_String.string_estate_types;
    std::map<std::string, EState_Type_Enum>::iterator pos = string_estate_types.find(input_string);
    if (pos == string_estate_types.end())  {
        std::string iii = Cantera::lowercase(input_string);
        pos = string_estate_types.find(iii);
        if (pos == string_estate_types.end())  {
            return EST_UNKNOWN_TYPE;
        }
    }
    return pos->second;
}
//====================================================================================================================
/*
 *  Defining memory for a static member of the clase
 */
EState_Factory* EState_Factory::s_factory = 0;
#if defined(THREAD_SAFE_CANTERA)
boost::mutex EState_Factory::estate_mutex;
#endif

//==================================================================================================================================
// Private constructors prevents usage
EState_Factory::EState_Factory()
{
}
//==================================================================================================================================
EState_Factory::~EState_Factory()
{
}
//==================================================================================================================================
// Static function that creates a static instance of the factory.
EState_Factory* EState_Factory::factory()
{
#if defined(THREAD_SAFE_CANTERA)
    boost::mutex::scoped_lock lock(estate_mutex);
#endif
    if (!s_factory) {
        s_factory = new EState_Factory;
    }
    return s_factory;
}
//==================================================================================================================================
// delete the static instance of this factory
void  EState_Factory::deleteFactory()
{
#if defined(THREAD_SAFE_CANTERA)
    boost::mutex::scoped_lock lock(estate_mutex);
#endif
    if (s_factory) {
        delete s_factory;
        s_factory = 0;
    }
}
//==================================================================================================================================
ETimeState::ETimeState() :
    cellNumber_(0),
    domainNumber_(0),
    es_(0),
    stateType_("t_final"),
    timeIncrType_("global"),
    time_(0.0),
    iOwnES_(true)
{
}

//==================================================================================================================================
ETimeState::ETimeState(const ETimeState& r) :
    cellNumber_(r.cellNumber_),
    domainNumber_(r.domainNumber_),
    es_(0),
    stateType_(r.stateType_),
    timeIncrType_(r.timeIncrType_),
    time_(r.time_),
    iOwnES_(r.iOwnES_)
{
    if (iOwnES_) {
	es_ = (r.es_)->duplMyselfAsEState();
    } else {
	es_ = (r.es_)->duplMyselfAsEState();
	iOwnES_ = true;
    }
}
//==================================================================================================================================
ETimeState::ETimeState(const Cantera::XML_Node& xTimeState, const Cantera::EState_ID_struct& e_id) :
    cellNumber_(0),
    domainNumber_(0),
    es_(0),
    stateType_("t_final"),
    timeIncrType_("global"),
    time_(0.0),
    iOwnES_(true)
{
    cellNumber_ = e_id.electrodeCellNumber_;
    domainNumber_ = e_id.electrodeDomainNumber_;
    read_ETimeState_fromXML(xTimeState, e_id);
}
//==================================================================================================================================
ETimeState& ETimeState::operator=(const ETimeState& r)
{
    if (this == &r) return *this;
    if (iOwnES_) {
	if (es_) {
	    delete es_;
	}
    }
    cellNumber_ = r.cellNumber_;
    domainNumber_ = r.domainNumber_;
    stateType_ = r.stateType_;
    timeIncrType_ = r.timeIncrType_;
    time_ = r.time_;
    iOwnES_ = r.iOwnES_;
    if (iOwnES_) {
	es_ = (r.es_)->duplMyselfAsEState();
    } else {
	es_ = r.es_;
    }
    return *this;
}
//==================================================================================================================================
XML_Node* ETimeState::write_ETimeState_ToXML() const
 {
     const std::string fmt = "%22.14E";
     XML_Node* xes = es_->write_electrodeState_ToXML();
     if (!xes) {
	 return NULL;
     }
     XML_Node* xmi = new XML_Node("timeState");
     if (stateType_ == "t_init" || stateType_ == "t_final" ||  stateType_ == "t_intermediate") {
	 xmi->addAttribute("type", stateType_);
     } else {
	 throw Electrode_Error( "ETimeState::write_ETimeState_ToXML()",
			      " Unknown state " + stateType_);
     }
     xmi->addAttribute("domain", domainNumber_);
     xmi->addAttribute("cellNumber", cellNumber_);
     xmi->addChild("time", time_, fmt);
     xmi->addChild(*xes);
     return xmi;
 }
//==================================================================================================================================
void ETimeState::read_ETimeState_fromXML(const Cantera::XML_Node& xTimeState, const Cantera::EState_ID_struct& e_id)
{
    /*
     *   Check to see that we are in the right spot
     */
    if (xTimeState.name() != "timeState") {
	throw Electrode_Error("read_ETimeState_fromXML", "Error: expecting timeState node but got " + xTimeState.name());
    }

    string ss = xTimeState["cellNumber"];
    cellNumber_ = atoi(ss.c_str());
    if (cellNumber_ != e_id.electrodeCellNumber_) {
	throw Electrode_Error( " ETimeState::read_ETimeState_fromXML",  "different cellNumbers");
    }

    ss = xTimeState["domain"];
    domainNumber_ = atoi(ss.c_str());
    if (domainNumber_ !=  e_id.electrodeDomainNumber_) {
	throw Electrode_Error( " ETimeState::read_ETimeState_fromXML",  "different domainNumbers");
    }

    ss = xTimeState["type"];
    if (ss != "t_init" && ss != "t_final"  && ss != "t_intermediate") {
	throw Electrode_Error( " ETimeState::read_ETimeState_fromXML",  "unknown type: " + ss);
    } 
    stateType_ = ss;
    timeIncrType_ = "global";

    if (es_) {
	if (iOwnES_) {
	    delete es_;
	}
    }
    iOwnES_ = 1;
    
    string typeString;
    string timeValStr;
    ctml::getNamedStringValue(xTimeState, "time", timeValStr, typeString);
    time_ = fpValueCheck(timeValStr);

    const XML_Node* xEState =  xTimeState.findByName("electrodeName");

    es_ =  createEState_fromXML(*xEState, e_id);

}
//==================================================================================================================================
//  Compare the current state of this object against another guest state to see if they are the same
/*
 *    We compare the state of the solution up to a certain number of digits.
 *
 *     @param[in]       ESguest          Guest state object to be compared against
 *     @param[in]       molarAtol        Absolute tolerance of the molar numbers in the state.
 *                                       Note from this value, we can get all other absolute tolerance inputs.
 *     @param[in]       nDigits          Number of digits to compare against
 *     @param[in]       includeHist      Include capacityDischarged and nextDeltaT variables in final bool comparison
 *     @param[in]       printLvl         print level of the routine
 *
 *     @return                           Returns true
 */
bool ETimeState::compareOtherTimeState(const ETimeState* const ETSguest, double molarAtol, int nDigits,
				       bool includeHist, int printLvl) const
{
    EState* esGuest_ =  ETSguest->es_;
    
    if (printLvl > 1) {
	printf("   CompareotherTimeState:        state1            state2\n");
	printf("               Time:             %11.5E         %11.5E\n", time_, ETSguest->time_);
	printf("               type:             %15s               %15s \n",  stateType_.c_str(), ETSguest->stateType_.c_str());
	printf("               Cell#:            %3d                  %3d\n", cellNumber_, ETSguest->cellNumber_);
    }

    bool ok =  es_->compareOtherState(esGuest_, molarAtol, nDigits, includeHist, printLvl);
    return ok;
}
//==================================================================================================================================
ETimeState::~ETimeState()
{
    if (iOwnES_) {
	delete es_;
    }
}
//==================================================================================================================================
ETimeInterval::ETimeInterval() :
     intervalType_("global"),
     index_(-1),
     numIntegrationSubCycles_(1),
     etsList_(0),
     deltaTime_init_next_(1.0E-8)
{
}

//==================================================================================================================================
ETimeInterval::~ETimeInterval() 
{
     for (size_t k = 0; k < etsList_.size(); ++k) {
        ETimeState* ets = etsList_[k];
        if (ets) {
           delete ets;
        } 
     }
}
//==================================================================================================================================
ETimeInterval::ETimeInterval(const Cantera::XML_Node& xTimeInterval, const Cantera::EState_ID_struct& e_id) :
    intervalType_("global"),
     index_(-1),
     numIntegrationSubCycles_(1),
     etsList_(0),
     deltaTime_init_next_(1.0E-8)
{
    read_ETimeInterval_fromXML(xTimeInterval, e_id);
}
//==================================================================================================================================
ETimeInterval::ETimeInterval(const ETimeInterval& right) :
     intervalType_(right.intervalType_),
     index_(right.index_),
     numIntegrationSubCycles_(right.numIntegrationSubCycles_),
     etsList_(0),
     deltaTime_init_next_(right.deltaTime_init_next_)
{
     etsList_.resize(numIntegrationSubCycles_+1, 0);
     for (size_t k = 0; k < (size_t) (numIntegrationSubCycles_+1); ++k) {
         etsList_[k] = new ETimeState(*(right.etsList_[k]));
     }
}
//==================================================================================================================================
// Create/Malloc an XML Node containing the ETimeInterval data contained in this object
/*
 *   @return   Returns the malloced XML_Node with name globalTimeStep containing the information in this
 *             object. The calling program is responsible for freeing this.
 */
Cantera::XML_Node* ETimeInterval::write_ETimeInterval_ToXML(int index) const
{
     int ii = index_;
     if (index >= 0) {
        ii = index;
     } 
     XML_Node* xtg = new XML_Node("globalTimeStep");
     xtg->addAttribute("index", int2str(ii));
     ctml::addInteger(*xtg,"numIntegrationSubCycles", numIntegrationSubCycles_);
     ctml::addFloat(*xtg, "deltaTime_init_next", deltaTime_init_next_);
     XML_Node* xti = new XML_Node("timeIncrement");
     xti->addAttribute("type", "global");
     for (size_t k = 0; k < etsList_.size(); ++k) {
         XML_Node* xmi = etsList_[k]->write_ETimeState_ToXML();
         xti->addChild(*xmi); 
     }
     xtg->addChild(*xti);
     return xtg;
}
//==================================================================================================================================
void ETimeInterval::read_ETimeInterval_fromXML(const Cantera::XML_Node& xTimeInterval, const Cantera::EState_ID_struct& e_id)
{
    string nn = xTimeInterval.name();
    if (nn != "globalTimeStep") {
	throw Electrode_Error("ETimeInterval::read_ETimeInterval_fromXML()",
			      "Was expecting the name globalTimeStep, but got instead: " + nn);
    }
    nn = xTimeInterval["index"];
    index_ = atoi(nn.c_str());
    numIntegrationSubCycles_ = ctml::getInteger(xTimeInterval, "numberIntegrationSybCycles");
    deltaTime_init_next_ = ctml::getFloat(xTimeInterval, "deltaTime_init_next");
    const XML_Node* xTimeIncr = xTimeInterval.findByName("timeIncrement");
    std::vector<XML_Node*> xStatesList = xTimeIncr->getChildren("timeState");
    size_t num =  xStatesList.size();
    etsList_.resize(num, 0);
    for (size_t k = 0; k < num; ++k) {
	XML_Node* xState = xStatesList[k];
	etsList_[k] = new ETimeState(*xState, e_id);
    }
    intervalType_ = "global";
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
ElectrodeTimeEvolutionOutput::ElectrodeTimeEvolutionOutput() :
    index_(1),
    timeStamp_(""),
    e_ID_(),
    etiList_(0)
{
}
//==================================================================================================================================
ElectrodeTimeEvolutionOutput::~ElectrodeTimeEvolutionOutput()
{
    for (size_t k = 0; k < etiList_.size(); ++k) {
	delete etiList_[k];
    }
}
//==================================================================================================================================
ElectrodeTimeEvolutionOutput::ElectrodeTimeEvolutionOutput(const Cantera::XML_Node& xElectrodeOutput) :
    index_(1),
    timeStamp_(""),
    e_ID_(),
    etiList_(0)
{
     read_ElectrodeTimeEvolutionOutput_fromXML(xElectrodeOutput);
}
//==================================================================================================================================
ElectrodeTimeEvolutionOutput::ElectrodeTimeEvolutionOutput(const ElectrodeTimeEvolutionOutput& right) :
    index_(right.index_),
    timeStamp_(right.timeStamp_),
    e_ID_(right.e_ID_)
{
    etiList_.resize(right.etiList_.size(), 0);
    for (size_t k = 0; k < etiList_.size(); ++k) {
	etiList_[k] = new ETimeInterval(*(right.etiList_[k]));
    }
}
//==================================================================================================================================
Cantera::XML_Node* ElectrodeTimeEvolutionOutput::write_ElectrodeTimeEvolutionOutput_ToXML(int index) const
{
    XML_Node* xEO = 0; new XML_Node("electrodeOutput");
    (*xEO)["index"] = int2str(index);
    ctml::addString(*xEO, "timeStamp", timeStamp_);
    XML_Node* xID = e_ID_.writeIdentificationToXML();
    xEO->addChild(*xID);
    for (size_t k = 0; k < etiList_.size(); ++k) {
	ETimeInterval* eti = etiList_[k];
	int timeIndex =  eti->index_;
	XML_Node* xETI = eti->write_ETimeInterval_ToXML(timeIndex);
	xEO->addChild(*xETI);
    }
    return xEO;
}
//==================================================================================================================================
void ElectrodeTimeEvolutionOutput::read_ElectrodeTimeEvolutionOutput_fromXML(const Cantera::XML_Node& xElectrodeOutput)
{
    string valueString, typeString;
    string nn = xElectrodeOutput.name();
    if (nn != "electrodeOutput") {
	throw Electrode_Error("ElectrodeTimeEvolutionOutput::read_ElectrodeTimeEvolutionOutput_fromXML",
			      "Was expecting the name electrodeOutput, but got instead: " + nn);
    }
    nn = xElectrodeOutput["index"];
    index_ = atoi(nn.c_str());

    ctml::getNamedStringValue(xElectrodeOutput, "timeStamp", valueString, typeString);
    timeStamp_ = valueString;
    const XML_Node* xmlEI = xElectrodeOutput.findByName("ElectrodeIdentification");
    e_ID_.readIdentificationFromXML(*xmlEI);

    std::vector<XML_Node*> xGlobalTimeSteps = xElectrodeOutput.getChildren("globalTimeStep");
    numGlobalTimeIntervals_ =  xGlobalTimeSteps.size();
    etiList_.resize(numGlobalTimeIntervals_, 0);
    for (size_t k = 0; k <(size_t) numGlobalTimeIntervals_; ++k) {
	XML_Node* xState = xGlobalTimeSteps[k];
	etiList_[k] = new ETimeInterval(*xState, e_ID_);
    }
}
//==================================================================================================================================
// Create a new EState Object
/*
 * @param model  String to look up the model against
 *
 * @return    Returns a pointer to a new EState instance matching the  model string. Returns NULL if
 *            something went wrong. Throws an exception if the string wasn't matched.
 */
Cantera::EState* EState_Factory::newEStateObject(std::string model)
{
    /*
     *  Look up the string to find the enum
     */
    EState_Type_Enum ieos = string_to_EState_Type_Enum(model);
    EState* ee = 0;
    /*
     *  Do the object creation
     */
    switch (ieos) {
    case EST_CSTR:
        ee = new EState();
        break;
    case EST_MULTIPLATEAU:
	throw Electrode_Error("EState_Factory::newEStateObject()",
			      "Unknown EState model: " + model);
	break;
    case EST_RADIALDISTRIB:
	ee = new EState_RadialDistrib();
	break;
    default:
        throw Electrode_Error("EState_Factory::newEStateObject()",
			      "Unknown EState model: " + model);
    }
    return ee;
}
//==================================================================================================================================

Cantera::XML_Node* getElectrodeIndentification()
{

    return 0;
}

//==================================================================================================================================
Cantera::XML_Node* selectLastGlobalTimeStepInterval(Cantera::XML_Node* xSoln, int& globalTimeStepNum)
{
    /*
     *  Find the electrodeOutput XML element.
     *     -> Later we will generalize this to search amongst multiple electrodeOutput objects
     */
    XML_Node* eOutput = xSoln->findByName("electrodeOutput");
    if (!eOutput) {
        ESModel_Warning("Electrode::selectGlobalTimeStepIncrement()",
                        "could not find a node named electrodeOutput");
        globalTimeStepNum = 0;
        return eOutput;
    }
    /*
     *  Get pointers to all of the children representing global time steps.
     */
    std::vector<XML_Node*> xGlobalTimeSteps = eOutput->getChildren("globalTimeStep");

    size_t numG =  xGlobalTimeSteps.size();
    if (numG == 0) {
        ESModel_Warning("Electrode::selectGlobalTimeStepIncrement()",
                        "Could not find globalTimeStep children");
        globalTimeStepNum = 0;
        return NULL;
    }
    int maxG = 0;
    int currIndex = 0;
    size_t jMax = 0;
    for (size_t j = 0; j < numG; ++j) {
	XML_Node* nodeG = xGlobalTimeSteps[j];
	std::string currSIndex = nodeG->attrib("index");
	currIndex = intValue(currSIndex);
	if (currIndex >= maxG) {
	    jMax = j;
	    maxG = currIndex;
	}
    }
    
    globalTimeStepNum = maxG;
    return xGlobalTimeSteps[jMax];
}
//==================================================================================================================================
Cantera::XML_Node* locateTimeLast_GlobalTimeStepIntervalFromXML(const XML_Node& xmlGlobalTimeStep, double& timeVal, int printSteps)
{
    string typeString;
    string timeValStr;
   if (xmlGlobalTimeStep.name() != "globalTimeStep") {
        ESModel_Warning(" EState::readGlobalTimeStepIntervalFromXML",
                        " Name of the xml node should have been globalTimeStep. Instead it was " + xmlGlobalTimeStep.name());
        timeVal = 0.0;
        return NULL;
    }
    std::vector<XML_Node*>  xSteps = xmlGlobalTimeStep.getChildren("timeIncrement"); 
    XML_Node* lastTimeIncrement = xSteps.back();
    
    std::vector<XML_Node*> xTimeStates = lastTimeIncrement->getChildren("timeState");             
   
    if (printSteps) {
	for (size_t n = 0; n <  xTimeStates.size(); ++n) {
	    XML_Node* xTimeState = xTimeStates[n];
	    std::string tt = xTimeState->attrib("type");
	    ctml::getNamedStringValue(*xTimeState, "time", timeValStr, typeString);
	    timeVal = Cantera::fpValueCheck(timeValStr);
	    printf(" type = %s  time = %g\n", tt.c_str(), timeVal);   
	}
    }
    XML_Node* xTimeState = xTimeStates.back();
    if (xTimeState->attrib("type") != "t_final") {
      printf("confused");
      exit(-1);
    }

    // Read the time
   
    ctml::getNamedStringValue(*xTimeState, "time", timeValStr, typeString);
    timeVal = fpValueCheck(timeValStr);
    XML_Node* eStateX = xTimeState->findByName("electrodeState");

    return eStateX;
}
//====================================================================================================================================
bool get_Estate_Indentification(const Cantera::XML_Node& xSoln, EState_ID_struct & e_id)
{


    bool retn = true;
    std::string typeSS;

    /*
     *  Find the head XML node for the identification
     */
    const XML_Node* xID = xSoln.findByName("ElectrodeIdentification");
    if (!xID) {
        throw Electrode_Error("esmodel::get_Estate_Indentification",
                              "could not find a node named ElectrodeIndentification");
    }

   
    ctml::getNamedStringValue(*xID, "electrodeTypeString", e_id.electrodeTypeString_ , typeSS);

    e_id.EST_Type_ = (enum EState_Type_Enum)  ctml::getInteger(*xID, "EState_Type");
    e_id.EState_Type_String_ = EState_Type_Enum_to_string(e_id.EST_Type_ );

  
    ctml::getNamedStringValue(*xID,"electrodeTypeString", e_id.electrodeTypeString_ , typeSS);

    if (xID->hasChild("fileVersionNumber")) {
	e_id.EST_Version_ = ctml::getInteger(*xID, "fileVersionNumber");
    }
    if (e_id.EST_Version_ != 1) {
        throw Electrode_Error(" EState::readIdentificationFromXML()",
                           "read the wrong EState version type - new situation");
    }
    if (xID->hasChild("electrodeCapacityType")) {
	e_id.electrodeCapacityType_ = (Cantera::Electrode_Capacity_Type_Enum) ctml::getInteger(*xID, "electrodeCapacityType");
    }
    if (xID->hasChild("electrodeModelType")) {
	e_id.electrodeChemistryModelType_ = ctml::getInteger(*xID, "electrodeModelType");
    }
    if (xID->hasChild("electrodeDomainNumber")) {
	e_id.electrodeDomainNumber_ = ctml::getInteger(*xID, "electrodeDomainNumber");
    }
    if (xID->hasChild("electrodeCellNumber")) {
	e_id.electrodeCellNumber_ = ctml::getInteger(*xID, "electrodeCellNumber");
    }
    return retn;
}
//==================================================================================================================================
Cantera::EState* newEStateObject(std::string model, EState_Factory* f)
{
    if (f == 0) {
        f = EState_Factory::factory();
    }
    return f->newEStateObject(model);

}
//==================================================================================================================================
/*
 *  Read the last time step from an EState XML file and return the EState object
 */
Cantera::EState* readEStateFileLastStep(const std::string& XMLfileName, double& timeRead)
{
    timeRead = 0.0;
    /*
     *   Read the XML file name getting the first ElectrodeOutput XML tree.
     *    -> that's what the 1 argument is. May do more with this in the future.
     */
    Cantera::XML_Node* xEout = getElectrodeOutputFile(XMLfileName, 1);
    if (!xEout) {
	ESModel_Warning("getElectrodeOutputFile", "Error");
	return NULL;
    }
    /*
     *   Read the identification structure, putting that into e_id struct
     */
    EState_ID_struct e_id;
    get_Estate_Indentification(*xEout , e_id);
    /*
     *   Malloc the appropriate type of EState object
     */
    EState* es = newEStateObject(e_id.EState_Type_String_);
    if (!es) {
	return NULL;
    }
    /*
     *   Put the identification information back into the EState object
     */
    es->readIdentificationFromXML(*xEout);
    /*
     *   Select the last global time step increment
     */
    int globalTimeStepNum = 0;
    Cantera::XML_Node* x = selectLastGlobalTimeStepInterval(xEout, globalTimeStepNum);
    if (!x) {
	delete es;
	return NULL;
    }
    /*
     *  Select the t_final step from the last global time step
     */
    double timeVal;
    Cantera::XML_Node* xSt = locateTimeLast_GlobalTimeStepIntervalFromXML(*x, timeVal, 0);
    if (!xSt) {
	delete es;
	return NULL;
    }
    /*
     *  Read the t_final state into es and also return the time
     */
    timeRead = timeVal;
    es->readStateFromXML(*xSt);
    // printf("Read global time step num %d representing solution at time %g\n", globalTimeStepNum, timeVal);

    return es;
}
//==================================================================================================================================
Cantera::EState* createEState_fromXML(const Cantera::XML_Node& xEState, const Cantera::EState_ID_struct& e_id)
{
    /*
     *   Check to see that we are in the right spot
     */
    if (xEState.name() != "electrodeState") {
	ESModel_Warning("createEState_fromXML", "Error: expecting electrodeState node but got " + xEState.name());
	return NULL;
    }
    /*
     *   Malloc the appropriate type of EState object
     */
    EState* es = newEStateObject(e_id.EState_Type_String_);
    if (!es) {
	return NULL;
    }
    /*
     *   Put the identification information back into the EState object
     */
    es->readIdentificationFromStruct(e_id); 
    /*
     *   Read the state into es and return
     */
    es->readStateFromXML(xEState);
    return es;
}
//==================================================================================================================================
esmodel::ElectrodeTimeEvolutionOutput* readXMLElectrodeOutput(const Cantera::XML_Node& Xfile, int index)
{
    /*
     *  Find the electrodeOutput XML element.
     *   
     */
    XML_Node* xRecord = Xfile.findNameIDIndex("electrodeOutput", "", index);
    if (!xRecord) {
        ESModel_Warning("esmodel::getElectrodeOutputFile()",
	                "Could not find a node named electrodeOutput with index " + mdpUtil::int2str(index));
    }
    ElectrodeTimeEvolutionOutput* e_teo = new ElectrodeTimeEvolutionOutput(*xRecord);
    return e_teo;
}
//==================================================================================================================================
esmodel::ElectrodeTimeEvolutionOutput* readFileElectrodeOutput(const std::string& XMLfileName, int index)
{
    /*
     *   Read the XML file name getting the first ElectrodeOutput XML tree.
     *    -> that's what the 1 argument is. May do more with this in the future.
     */
    Cantera::XML_Node* xEout = getElectrodeOutputFile(XMLfileName, index);
    if (!xEout) {
	ESModel_Warning("getElectrodeOutputFile", "Error");
	return NULL;
    }
    ElectrodeTimeEvolutionOutput* e_teo = new ElectrodeTimeEvolutionOutput(*xEout);
    return e_teo;
}
//==================================================================================================================================
void writeElectrodeOutputFile(std::string fileName, const esmodel::ElectrodeTimeEvolutionOutput& e_teo)
{
     Cantera::XML_Node root("--");
     Cantera::XML_Node& ct = root.addChild("ctml");
     Cantera::XML_Node* x_eteo = e_teo.write_ElectrodeTimeEvolutionOutput_ToXML();
     ct.addChild(*x_eteo);
     std::fstream s(fileName.c_str(), fstream::in | fstream::out | fstream::trunc);
     root.write(s);
     s.close();
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------


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
//======================================================================================================================
Cantera::XML_Node* getElectrodeOutputFile(const std::string& fileName, int index)
{
    Cantera::XML_Node* xSavedSoln = get_XML_File(fileName);
    if (!xSavedSoln) {
	throw Electrode_Error(" getElectrodeOutputFile()",
			      " Could not open file, " + fileName);
    }
    /*
     *  Find the electrodeOutput XML element.
     *   
     */
    XML_Node* eOutput = xSavedSoln->findByName("electrodeOutput");
    if (!eOutput) {
        throw Electrode_Error("esmodel::getElectrodeOutputFile()",
			      "could not find a node named electrodeOutput");
    }

    XML_Node* eRecord = eOutput->findNameIDIndex("electrodeOutput", "", index);
    if (!eRecord) {
        throw Electrode_Error("esmodel::getElectrodeOutputFile()",
			      "could not find a node named electrodeOutput with index " + mdpUtil::int2str(index));
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

Cantera::XML_Node* selectLastGlobalTimeStepIncrement(Cantera::XML_Node* xSoln, int& globalTimeStepNum)
{
    /*
     *  Find the electrodeOutput XML element.
     *     -> Later we will generalize this to search amongst multiple electrodeOutput objects
     */
    XML_Node* eOutput = xSoln->findByName("electrodeOutput");
    if (!eOutput) {
        throw CanteraError("Electrode::selectGlobalTimeStepIncrement()",
                           "could not find a node named electrodeOutput");
    }
    std::vector<XML_Node*>  xGlobalTimeSteps = eOutput->getChildren("globalTimeStep");



    size_t numG =  xGlobalTimeSteps.size();
    if (numG == 0) {
	throw CanteraError("Electrode::selectGlobalTimeStepIncrement()",
                           "could not find a node named electrodeOutput");
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
Cantera::XML_Node* locateTimeLast_GlobalTimeStepFromXML(const XML_Node& xmlGlobalTimeStep, double& timeVal, int printSteps)
{
    string typeString;
    string timeValStr;
   if (xmlGlobalTimeStep.name() != "globalTimeStep") {
        throw Electrode_Error(" EState::readGlobalTimeStepFromXML",
                              " Name of the xml node should have been globalTimeStep. Instead it was " + xmlGlobalTimeStep.name());
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

    // locate the electrodeState

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
        throw CanteraError(" EState::readIdentificationFromXML()",
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

}
//----------------------------------------------------------------------------------------------------------------------------------







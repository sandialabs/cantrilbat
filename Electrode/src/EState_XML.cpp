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

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

using namespace std;
//----------------------------------------------------------------------------------------------------------------------------------
//! esmodel stands for Electrode source model. It will expand to encompass everything in this directory
namespace esmodel
{
//==================================================================================================================================
Map_ESEnum_String gMap_ESEnum_String;

//==================================================================================================================================
//! Create String maps for the Types of EState objects
/*!
 *  We automatically create the inversion of the maps
 */
static void create_string_maps()
{
    if (gMap_ESEnum_String.string_maps_created) {
        return;
    }
    std::map<ZZCantera::EState_Type_Enum, std::string>& estate_types_string = gMap_ESEnum_String.estate_types_string;
    std::map<std::string , EState_Type_Enum>& string_estate_types = gMap_ESEnum_String.string_estate_types;

    gMap_ESEnum_String.string_maps_created = true;
    estate_types_string[EST_CSTR]                       =  "EState_CSTR";
    estate_types_string[EST_CSTR_LiCoO2Cathode]         =  "EState_CSTR_LiCoO2Cathode";
    estate_types_string[EST_CSTR_MCMBAnode]             =  "EState_CSTR_MCMBAnode";
    estate_types_string[EST_MULTIPLATEAU]               =  "EState_MultiPlateau";
    estate_types_string[EST_RADIALDISTRIB]              =  "EState_RadialDistrib";

    // Invert the maps automatically.
    for (std::map<EState_Type_Enum, std::string>::iterator pos = estate_types_string.begin();
            pos != estate_types_string.end(); ++pos) {
        string_estate_types[pos->second] = pos->first;
        std::string lll =  ZZCantera::lowercase(pos->second);
        string_estate_types[lll] = pos->first;
    }
}
//=================================================================================================================================
ZZCantera::XML_Node* getElectrodeOutputFile(const std::string& fileName, int index)
{
    /*
     * Call the basic Cantera routine that reads the XML file. If the file isn't found a NULL is returned 
     * causing this program to also return NULL
     */
    ZZCantera::XML_Node* xSavedSoln = get_XML_File(fileName);
    if (!xSavedSoln) {
        return xSavedSoln;
    }
    XML_Node* xCTML = xSavedSoln->findByName("ctml");
    if (!xCTML) {
	ESModel_Warning("esmodel::getElectrodeOutputFile()", "Could not find a node named ctml");
	return xCTML;
    }
    /*
     *  Find the electrodeOutput XML element.
     *   
     */

    XML_Node* eRecord = xCTML->findNameIDIndex("electrodeOutput", "", index);
    if (!eRecord) {
        ESModel_Warning("esmodel::getElectrodeOutputFile()",
	                "Could not find a node named electrodeOutput with index " + mdpUtil::int2str(index));
    }
    return eRecord;
}
//==================================================================================================================================
std::string EState_Type_Enum_to_string(const ZZCantera::EState_Type_Enum& etype)
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
        std::string iii = ZZCantera::lowercase(input_string);
        pos = string_estate_types.find(iii);
        if (pos == string_estate_types.end())  {
            return EST_UNKNOWN_TYPE;
        }
    }
    return pos->second;
}
//==================================================================================================================================

EState_Factory* EState_Factory::s_factory = 0;

#if defined(THREAD_SAFE_CANTERA)
boost::mutex EState_Factory::state_mutex;
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
//==================================================================================================================================
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
ETimeState::ETimeState(const ZZCantera::XML_Node& xTimeState, const ZZCantera::EState_ID_struct& e_id) :
    cellNumber_(0),
    domainNumber_(0),
    es_(0),
    stateType_("t_final"),
    timeIncrType_("global"),
    time_(0.0),
    iOwnES_(true)
{
    cellNumber_ = e_id.electrodeCellNumber;
    domainNumber_ = e_id.electrodeDomainNumber;
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
     XML_Node* xmi = new XML_Node("timeState");
     if (stateType_ == "t_init" || stateType_ == "t_final" ||  stateType_ == "t_intermediate") {
	 xmi->addAttribute("type", stateType_);
     } else {
	 throw Electrode_Error( "ETimeState::write_ETimeState_ToXML()", " Unknown state " + stateType_);
     }
     xmi->addAttribute("domain", domainNumber_);
     xmi->addAttribute("cellNumber", cellNumber_);
     // write the time out with 15 digits of accuracy -> want a restart capability, but don't want roundoff diffs.
     const std::string fmt = "%22.14E";
     xmi->addChild("time", time_, fmt);

     // Create a new node with the solution written into it
     XML_Node* xes = es_->write_electrodeState_ToXML();
     if (!xes) {
	 return NULL;
     }
     // Then, add this node (without a copy) to the tree
     xmi->addChildToTree(xes);

     return xmi;
 }
//==================================================================================================================================
void ETimeState::read_ETimeState_fromXML(const ZZCantera::XML_Node& xTimeState, const ZZCantera::EState_ID_struct& e_id)
{
    /*
     *   Check to see that we are in the right spot
     */
    if (xTimeState.name() != "timeState") {
	throw Electrode_Error("read_ETimeState_fromXML", "Error: expecting timeState node but got " + xTimeState.name());
    }

    std::string ss = xTimeState["cellNumber"];
    cellNumber_ = atoi(ss.c_str());
    if (cellNumber_ != e_id.electrodeCellNumber) {
	throw Electrode_Error("ETimeState::read_ETimeState_fromXML",  "different cellNumbers");
    }

    ss = xTimeState["domain"];
    domainNumber_ = atoi(ss.c_str());
    if (domainNumber_ !=  e_id.electrodeDomainNumber) {
	throw Electrode_Error("ETimeState::read_ETimeState_fromXML",  "different domainNumbers");
    }

    ss = xTimeState["type"];
    if (ss != "t_init" && ss != "t_final"  && ss != "t_intermediate") {
	throw Electrode_Error("ETimeState::read_ETimeState_fromXML",  "unknown type: " + ss);
    } 
    stateType_ = ss;
    timeIncrType_ = "global";

    if (es_) {
	if (iOwnES_) {
	    delete es_;
	}
    }
    iOwnES_ = 1;
    
    std::string typeString;
    std::string timeValStr;
    ZZctml::getNamedStringValue(xTimeState, "time", timeValStr, typeString);
    time_ = fpValueCheck(timeValStr);

    const XML_Node* xEState = xTimeState.findByName("electrodeState");
    if (!xEState) {
	throw Electrode_Error("ETimeState::read_ETimeState_fromXML()", "Could not find the XML element electrodeState");
    }

    es_ = createEState_fromXML(*xEState, e_id);
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
 
    //  Note, we print out differences in time, but we don't throw an error if time_ isn't the same between the two structs
   
    if (printLvl > 1) {
	printf("   CompareOtherTimeState:        state1            state2\n");
	printf("                    Time:      %11.5E         %11.5E\n", time_, ETSguest->time_);
	printf("                    type:       %15s               %15s \n",  stateType_.c_str(), ETSguest->stateType_.c_str());
	printf("                   Cell#:        %3d                  %3d\n", cellNumber_, ETSguest->cellNumber_);
	printf("            timeIncrType:  %15s               %15s \n", timeIncrType_.c_str(), ETSguest->timeIncrType_.c_str());
    }
    bool ok =  es_->compareOtherState(esGuest_, molarAtol, nDigits, includeHist, printLvl);
    return ok;
}
//==================================================================================================================================
double ETimeState::electrodeMoles() const
{
    return es_->electrodeMoles();
}
//==================================================================================================================================
ETimeState::~ETimeState()
{
    if (iOwnES_) {
	delete es_;
    }
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
ETimeInterval::ETimeInterval() :
     intervalType_("global"),
     index_(-1),
     numIntegrationSubCycles_(1),
     etsList_(0),
     deltaTime_init_next_(1.0E-8),
     deltaTime_init_init_(1.0E-8)
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
ETimeInterval::ETimeInterval(const ZZCantera::XML_Node& xTimeInterval, const ZZCantera::EState_ID_struct& e_id) :
    ETimeInterval()
{
    read_ETimeInterval_fromXML(xTimeInterval, e_id);
}
//==================================================================================================================================
ETimeInterval::ETimeInterval(const ETimeInterval& right) :
     intervalType_(right.intervalType_),
     index_(right.index_),
     numIntegrationSubCycles_(right.numIntegrationSubCycles_),
     deltaTime_init_next_(right.deltaTime_init_next_),
     deltaTime_init_init_(right.deltaTime_init_init_)
{
     etsList_.resize(numIntegrationSubCycles_+1, 0);
     for (size_t k = 0; k < (size_t) (numIntegrationSubCycles_+1); ++k) {
         etsList_[k] = new ETimeState(*(right.etsList_[k]));
     }
}
//==================================================================================================================================
ETimeInterval& ETimeInterval::operator=(const ETimeInterval& right)
{
     if (this != &right) {
         intervalType_ = right.intervalType_;
         index_ = right.index_;
         numIntegrationSubCycles_ = right.numIntegrationSubCycles_;

         for (size_t k = 0; k < etsList_.size(); ++k) {
            if (etsList_[k]) {
                delete (etsList_[k]);
            }
         }
         etsList_.resize( right.etsList_.size(), 0);
         for (size_t k = 0; k < etsList_.size(); ++k) {
             if (right.etsList_[k]) {
                 etsList_[k] = new ETimeState( *(right.etsList_[k]) );
             } else {
                 etsList_[k] = nullptr;
             }
         }

         deltaTime_init_next_ = right.deltaTime_init_next_;
         deltaTime_init_init_ = right.deltaTime_init_init_;
     }
     return *this;
}
//==================================================================================================================================
// Create/Malloc an XML Node containing the ETimeInterval data contained in this object
/*
 *   @return   Returns the malloced XML_Node with name globalTimeStep containing the information in this
 *             object. The calling program is responsible for freeing this.
 */
ZZCantera::XML_Node* ETimeInterval::write_ETimeInterval_ToXML(int index, int windex ) const
{
     int ii = index_;
     if (index >= 0) {
        ii = index;
     } 

     std::string fmt = "%22.14E";

     XML_Node* xtg = new XML_Node("globalTimeStep");
     xtg->addAttribute("index", int2str(ii));
     xtg->addAttribute("windex", int2str(1));
     double t_init_init = startingTime();
     xtg->addAttribute("t_init_init", fp2str(t_init_init, fmt));
     double t_final_final = endingTime();
     xtg->addAttribute("t_final_final", fp2str(t_final_final, fmt));

     ZZctml::addFloat(*xtg, "deltaTime_init_next", deltaTime_init_next_);
     ZZctml::addFloat(*xtg, "deltaTime_init_init", deltaTime_init_init_);
     ZZctml::addInteger(*xtg,"numIntegrationSubCycles", numIntegrationSubCycles_);
     XML_Node* xti = new XML_Node("timeIncrement");
     xti->addAttribute("type", "global");
     for (size_t k = 0; k < etsList_.size(); ++k) {
         XML_Node* xmi = etsList_[k]->write_ETimeState_ToXML();
         xti->addChildToTree(xmi); 
     }
     xtg->addChildToTree(xti);
     return xtg;
}
//==================================================================================================================================
void ETimeInterval::read_ETimeInterval_fromXML(const ZZCantera::XML_Node& xTimeInterval, const ZZCantera::EState_ID_struct& e_id)
{
    std::string nn = xTimeInterval.name();
    if (nn != "globalTimeStep") {
	throw Electrode_Error("ETimeInterval::read_ETimeInterval_fromXML()",
			      "Was expecting the name \"globalTimeStep\", but got instead: " + nn);
    }
    nn = xTimeInterval["index"];
    index_ = atoi(nn.c_str());
    numIntegrationSubCycles_ = ZZctml::getInteger(xTimeInterval, "numIntegrationSubCycles");
    deltaTime_init_next_ = ZZctml::getFloat(xTimeInterval, "deltaTime_init_next");
    deltaTime_init_init_ = ZZctml::getFloat(xTimeInterval, "deltaTime_init_init");
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
//  Compare the current state of this object against another guest state to see if they are the same
/*
 *    We compare the state of the solution up to a certain number of digits.
 *
 *     @param[in]       ETIguest         Guest time interval to be compared against
 *     @param[in]       molarAtol        Absolute tolerance of the molar numbers in the state.
 *                                       Note from this value, we can get all other absolute tolerance inputs.
 *     @param[in]       nDigits          Number of digits to compare against
 *     @param[in]       includeHist      Include capacityDischarged and nextDeltaT variables in final bool comparison
 *     @param[in]       compareType      Comparison type:
 *                                           0 Intermediates and initial state has to be the same
 *     @param[in]       printLvl         print level of the routine
 *
 *     @return                           Returns true if the times are the same and the states are the same.
 */
bool ETimeInterval::compareOtherETimeInterval(const ETimeInterval* const ETIguest, double molarAtol, double unitlessAtol, int nDigits,
					      bool includeHist, int compareType, int printLvl) const
{
    bool total_ok = true, ok;
    ok = numIntegrationSubCycles_ == ETIguest->numIntegrationSubCycles_;
    total_ok = total_ok && ok;

    ok = etsList_.size() == ETIguest->etsList_.size();
    total_ok = total_ok && ok;

    int subPrint = printLvl - 1;
    if (subPrint < 0) subPrint = 0;

    if (ok) {
	if (subPrint == 1) {
	    printf("         CompareOtherETimeInterval:\n");
	    printf("               time1         time2            state1   state2         success\n");
	}
	for (size_t k = 0; k < etsList_.size(); ++k) {
	    const ETimeState* ets = etsList_[k];
	    const ETimeState* etsGuest =  ETIguest->etsList_[k];
	   
	    ok = doubleEqual(ets->time_, etsGuest->time_, unitlessAtol, nDigits);
	    if (printLvl) {
		if (printLvl == 1) {
		    printf("        %E15.8   %E15.8   %15s %15s   ", 
			   ets->time_, etsGuest->time_, ets->stateType_.c_str(), etsGuest->stateType_.c_str());	
		} else {
		    if (ok) {
			printf("CompareTimeIntervals: times are equal time = %g\n",  ets->time_);
		    } else {
			printf("CompareTimeIntervals: times different: %g %g\n",  ets->time_, etsGuest->time_);
		    }
		}	
	    }
	    total_ok = total_ok && ok;

	    ok = ets->compareOtherTimeState(etsGuest, molarAtol, nDigits, true, subPrint);
	    if (printLvl) {
		if (printLvl == 1) {
		    if (ok) {
			printf("    OK\n");
		    } else {
			printf("    Different\n");
		    }
		} else {
		    if (ok) {
			printf("CompareTimeIntervals: ETimeStates are equivalent\n");
		    } else {
			printf("CompareTimeIntervals: ETimeStates are different\n");
		    }
		}
	    }
	    total_ok = total_ok && ok;
	}
    
    }

    return total_ok;
}
//==================================================================================================================================
double ETimeInterval::startingTime() const
{
    ETimeState* ets = etsList_[0];
    return ets->time_;
}
//==================================================================================================================================
double ETimeInterval::endingTime() const
{
    ETimeState* ets = etsList_.back();
    return ets->time_;
}
//==================================================================================================================================
double ETimeInterval::electrodeInitialMoles() const
{
    ETimeState* ets = etsList_[0];
    return ets->electrodeMoles();
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
ElectrodeTimeEvolutionOutput::ElectrodeTimeEvolutionOutput(const ZZCantera::XML_Node& xElectrodeOutput) :
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
ZZCantera::XML_Node* ElectrodeTimeEvolutionOutput::write_ElectrodeTimeEvolutionOutput_ToXML(int index) const
{
    int windex = index_;
    if (index != -1) {
        windex = index;
    }
    XML_Node* xEO = new XML_Node("electrodeOutput");
    xEO->addAttribute("index",  int2str(windex));
    ZZctml::addString(*xEO, "timeStamp", timeStamp_);
    XML_Node* xID = e_ID_.writeIdentificationToXML();
    xEO->addChildToTree(xID);
    for (size_t k = 0; k < etiList_.size(); ++k) {
	ETimeInterval* eti = etiList_[k];
	int timeIndex =  eti->index_;
	XML_Node* xETI = eti->write_ETimeInterval_ToXML(timeIndex, windex);
	xEO->addChildToTree(xETI);
    }
    return xEO;
}
//==================================================================================================================================
void ElectrodeTimeEvolutionOutput::read_ElectrodeTimeEvolutionOutput_fromXML(const ZZCantera::XML_Node& xElectrodeOutput)
{
    std::string valueString, typeString;
    std::string nn = xElectrodeOutput.name();
    if (nn != "electrodeOutput") {
	throw Electrode_Error("ElectrodeTimeEvolutionOutput::read_ElectrodeTimeEvolutionOutput_fromXML",
			      "Was expecting the name electrodeOutput, but got instead: " + nn);
    }
    nn = xElectrodeOutput["index"];
    index_ = atoi(nn.c_str());

    ZZctml::getNamedStringValue(xElectrodeOutput, "timeStamp", valueString, typeString);
    timeStamp_ = valueString;
    const XML_Node* xmlEI = xElectrodeOutput.findByName("ElectrodeIdentification");
    if (! xmlEI) {
       throw Electrode_Error("ElectrodeTimeEvolutionOutput::read_ElectrodeTimeEvolutionOutput_fromXML",
                              "Was expecting to find the ElectrodeIdentification XML Element, but didn't. Something is seriously wrong");
    }
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
bool ElectrodeTimeEvolutionOutput::compareOtherTimeEvolution(const ElectrodeTimeEvolutionOutput* const ETOguest,
							     double molarAtol, double unitlessAtol, int nDigits,
							     bool includeHist, int compareType, int printLvl) const
{
    // We'll first code up a simple 1 - 1 comparison
    bool total_ok = true, ok;
    ok = numGlobalTimeIntervals_ == ETOguest->numGlobalTimeIntervals_;
    total_ok = total_ok && ok;

    ok = etiList_.size() == ETOguest->etiList_.size();
    total_ok = total_ok && ok;

    int subPrint = printLvl - 1;
    if (subPrint < 0) subPrint = 0;

    if (ok) {
        if (subPrint == 1) {
            printf("         CompareOtherETimeInterval:\n");
            printf("               time1 - time1Last        time2- time2Last              state1   state2         success\n");
        }
        for (size_t k = 0; k < etiList_.size(); ++k) {
            const ETimeInterval* ets = etiList_[k];
            const ETimeInterval* etsGuest = ETOguest->etiList_[k];

	    ETimeState* firstTimeState = ets->etsList_[0];
	    double firstTime =  firstTimeState->time_;

	    ETimeState* firstTimeStateG = etsGuest->etsList_[0];
	    double firstTimeG =  firstTimeStateG->time_;

	    size_t klast = ets->etsList_.size() - 1;
	    ETimeState* lastTimeState = ets->etsList_[klast];
	    double lastTime =  lastTimeState->time_;

	    size_t klastG = etsGuest->etsList_.size() - 1;
	    ETimeState* lastTimeStateG = etsGuest->etsList_[klastG];
	    double lastTimeG =  lastTimeStateG->time_;

            ok = doubleEqual(firstTime, firstTimeG, unitlessAtol, nDigits);
	    bool ok2 = doubleEqual(lastTime, lastTimeG, unitlessAtol, nDigits);
	    ok = ok && ok2;
            if (printLvl) {
                if (printLvl == 1) {
                    printf("        %E15.8   %E15.8  %E15.8   %E15.8   ",
                           firstTime, lastTime, firstTimeG, lastTimeG);
                } else {
                    if (ok) {
                        printf("CompareTimeIntervals: times are equal time = %g % g\n", firstTime, lastTime );
                    } else {
                        printf("CompareTimeIntervals: times different: %g - %g , %g - %g \n",
			       firstTime, lastTime, firstTimeG, lastTimeG);
                    }
                }
            }
            total_ok = total_ok && ok;
	    int compareType = 0;
            ok = ets->compareOtherETimeInterval(etsGuest, molarAtol, unitlessAtol, nDigits, true,  compareType, subPrint);
            if (printLvl) {
                if (printLvl == 1) {
                    if (ok) {
                        printf("    OK\n");
                    } else {
                        printf("    Different\n");
                    }
                } else {
                    if (ok) {
                        printf("CompareTimeIntervals: ETimeStates are equivalent\n");
                    } else {
                        printf("CompareTimeIntervals: ETimeStates are different\n");
                    }
                }
            }
            total_ok = total_ok && ok;
        }

    }
    return total_ok;
}
//==================================================================================================================================
bool ElectrodeTimeEvolutionOutput::compareOtherTimeEvolutionSub(const ElectrodeTimeEvolutionOutput* const ETOguest,
								int& numZonesNeededToPass,
								double molarAtol, double unitlessAtol, int nDigits,
								bool includeHist, int compareType, int printLvl) const
{
    // We'll first code up a simple 1 - 1 comparison
    bool total_ok = true, ok, ok2;
    int numChecked = 0;
    int numPassed = 0;
   
    std::vector<int> guestZoneMatch(numGlobalTimeIntervals_, -1);
    std::vector<int> hostZoneMatch(ETOguest->numGlobalTimeIntervals_, -1);


    int subPrint = printLvl - 1;
    if (subPrint < 0) subPrint = 0;

   
    if (subPrint == 1) {
	printf("         CompareOtherETimeInterval:\n");
	printf("            Zone1   time1 - time1Last    Zone2    time2- time2Last              state1   state2         success\n");
    }
    for (size_t k = 0; k < etiList_.size(); ++k) {
	const ETimeInterval* ets = etiList_[k];
	ETimeState* firstTimeState = ets->etsList_[0];
	double firstTime =  firstTimeState->time_;
	size_t klast = ets->etsList_.size() - 1;
	ETimeState* lastTimeState = ets->etsList_[klast];
	double lastTime =  lastTimeState->time_;

	for (size_t j = 0; j < (size_t) ETOguest->numGlobalTimeIntervals_; ++j) {
	    const ETimeInterval* etsGuest = ETOguest->etiList_[j];
	    ETimeState* firstTimeStateG = etsGuest->etsList_[0];
	    double firstTimeG =  firstTimeStateG->time_;
	    size_t klastG = etsGuest->etsList_.size() - 1;
	    ETimeState* lastTimeStateG = etsGuest->etsList_[klastG];
	    double lastTimeG =  lastTimeStateG->time_;

	    ok = doubleEqual(firstTime, firstTimeG, unitlessAtol, nDigits);
	    ok2 = doubleEqual(lastTime, lastTimeG, unitlessAtol, nDigits);
	    ok = ok && ok2;
	    /*
	     *  If beginning and ending time match, we'll call that a match
	     */
	    if (ok) {
		if (guestZoneMatch[k] != -1) {
		    throw Electrode_Error(" ElectrodeTimeEvolutionOutput::compareOtherTimeEvolutionSub",
					  " zones don't match up 1 to 1, there is a 1 to 2 match");
		}
		guestZoneMatch[k] = j;
		if (hostZoneMatch[j] != -1) {
		    throw Electrode_Error(" ElectrodeTimeEvolutionOutput::compareOtherTimeEvolutionSub",
					  " zones don't match up 1 to 1, there is a 2 to 1 match");
		}
		hostZoneMatch[j] = k;

		if (printLvl) {
		    if (printLvl == 1) {
			printf("    %lu    %E15.8   %E15.8 %lu %E15.8   %E15.8   ",
			       k, firstTime, lastTime,  j, firstTimeG, lastTimeG);
		    } else {
			printf("CompareTimeIntervals: times in zone1: %lu are equal to zone2: %lu  time = %g % g\n", 
			       k, j, firstTime, lastTime);
		    }
		}
		total_ok = total_ok && ok;
		int compareType = 0;
		numChecked++;
		ok = ets->compareOtherETimeInterval(etsGuest, molarAtol, unitlessAtol, nDigits, true,  compareType, subPrint);
                if (ok) {
                    ++numPassed;
                }
		if (printLvl) {
		    if (printLvl == 1) {
			if (ok) {
			    printf("    OK\n");
			} else {
			    printf("    Different\n");
			}
		    } else {
			if (ok) {
			    printf("CompareTimeIntervals: ETimeStates are equivalent\n");
			} else {
			    printf("CompareTimeIntervals: ETimeStates are different\n");
			}
		    }
		}
		total_ok = total_ok && ok;
	    }
	}
    }
    
    if (numChecked < numZonesNeededToPass) {
	if (printLvl) {
	    printf(" CompareTimeIntervals:  Not enough zones checked, %d, compared to requested, %d, to pass\n",
		   numChecked, numZonesNeededToPass);
	    if (numChecked > 0) {
		printf(" CompareTimeIntervals: The zones that were checked did pass!\n");
	    } else {
		printf(" CompareTimeIntervals:  No zones were checked at all\n");
	    }
	}
	total_ok = false;
    } else {
  	if (printLvl) {
	    printf(" CompareTimeIntervals:  Enough zones were checked, %d, compared to requested, %d, to pass\n",
		   numChecked, numZonesNeededToPass);
	    if (total_ok) {
		printf(" CompareTimeIntervals: PASSED!\n");
	    } else {
		printf(" CompareTimeIntervals: However, Some zones had different data in them, so FAILED\n");
	    }
	}
    }
    numZonesNeededToPass = numPassed;
    return total_ok;
}
//==================================================================================================================================
double ElectrodeTimeEvolutionOutput::electrodeInitialMoles() const
{
    ETimeInterval* eti = etiList_[0];
    return eti->electrodeInitialMoles();
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
// Create a new EState Object
/*
 * @param model  String to look up the model against
 *
 * @return    Returns a pointer to a new EState instance matching the  model string. Returns NULL if
 *            something went wrong. Throws an exception if the string wasn't matched.
 */
ZZCantera::EState* EState_Factory::newEStateObject(std::string model)
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
	ee = new EState(model);
	ee->electrodeTypeString_ = "CSTR";
	break;
    case EST_CSTR_LiCoO2Cathode:
	ee = new EState(model);
	ee->electrodeTypeString_ = "CSTR_LiCoO2Cathode";
	break;
    case EST_CSTR_MCMBAnode:
        ee = new EState(model);
	ee->electrodeTypeString_ = "CSTR_MCMBAnode";
	break;
        break;
    case EST_MULTIPLATEAU:
	throw Electrode_Error("EState_Factory::newEStateObject()",
			      "Unknown EState model: " + model);
	break;
    case EST_RADIALDISTRIB:
	ee = new EState_RadialDistrib(model);
	break;
    default:
        throw Electrode_Error("EState_Factory::newEStateObject()",
			      "Unknown EState model: " + model);
    }
    return ee;
}
//==================================================================================================================================
ZZCantera::XML_Node* selectLastGlobalTimeStepInterval(ZZCantera::XML_Node* xSoln, int& globalTimeStepNum)
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
ZZCantera::XML_Node* locateTimeLast_GlobalTimeStepIntervalFromXML(const ZZCantera::XML_Node& xmlGlobalTimeStep, double& timeVal, 
                                                                  int printSteps)
{
    std::string typeString;
    std::string timeValStr;
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
	    ZZctml::getNamedStringValue(*xTimeState, "time", timeValStr, typeString);
	    timeVal = ZZCantera::fpValueCheck(timeValStr);
	    printf(" type = %s  time = %g\n", tt.c_str(), timeVal);   
	}
    }
    XML_Node* xTimeState = xTimeStates.back();
    if (xTimeState->attrib("type") != "t_final") {
      printf("confused");
      exit(-1);
    }

    // Read the time
   
    ZZctml::getNamedStringValue(*xTimeState, "time", timeValStr, typeString);
    timeVal = fpValueCheck(timeValStr);
    XML_Node* eStateX = xTimeState->findByName("electrodeState");

    return eStateX;
}
//==================================================================================================================================
ZZCantera::EState* newEStateObject(std::string model, EState_Factory* f)
{
    if (f == 0) {
        f = EState_Factory::factory();
    }
    return f->newEStateObject(model);
}
//==================================================================================================================================
ZZCantera::EState* readEState_XMLFile_LastStep(const std::string& XMLfileName, double& timeRead)
{
    timeRead = 0.0;
    /*
     *   Read in the XML file using the XMLfile name as input. getting the first ElectrodeOutput XML tree.
     *    -> that's what the 1 argument is. May do more with this in the future.
     */
    ZZCantera::XML_Node* xEout = getElectrodeOutputFile(XMLfileName, 1);
    if (!xEout) {
	ESModel_Warning("getElectrodeOutputFile", "Error");
	return NULL;
    }
    /*
     *   Read the identification structure, putting that into e_id struct
     */
    EState_ID_struct e_id;
    e_id.readIdentificationFromXML(*xEout);

    /*
     *   Malloc the appropriate type of EState object
     */
    EState* es = newEStateObject(e_id.EState_Type_String);
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
    ZZCantera::XML_Node* x = selectLastGlobalTimeStepInterval(xEout, globalTimeStepNum);
    if (!x) {
	delete es;
	return NULL;
    }
    /*
     *  Select the t_final step from the last global time step
     */
    double timeVal;
    ZZCantera::XML_Node* xSt = locateTimeLast_GlobalTimeStepIntervalFromXML(*x, timeVal, 0);
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
ZZCantera::EState* createEState_fromXML(const ZZCantera::XML_Node& xEState, const ZZCantera::EState_ID_struct& e_id)
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
    EState* es = newEStateObject(e_id.EState_Type_String);
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
size_t reportXMLElectrodeOutput_NumRecords(const ZZCantera::XML_Node& Xfile)
{
    //
    //  Need to be at the ctml level for the next command to work
    //
    const XML_Node* xCTML = Xfile.findByName("ctml");
    if (!xCTML) {
        ESModel_Warning("esmodel::reportXMLElectrodeOutput_NumRecords()", "Could not find a node named \"ctml\"");
        return 0;
    }
    return xCTML->nChildren("electrodeOutput");
}
//==================================================================================================================================
esmodel::ElectrodeTimeEvolutionOutput* readXMLElectrodeOutput(const ZZCantera::XML_Node& Xfile, int index)
{
    //
    //  Need to be at the ctml level for the next command to work
    //
    const XML_Node* xCTML = Xfile.findByName("ctml");
    if (!xCTML) {
        ESModel_Warning("esmodel::readXMLElectrodeOutput()", "Could not find a node named \"ctml\"");
        return nullptr;
    }
    /*
     *  Find the correct electrodeOutput XML element.
     *   
     */
    const XML_Node* xRecord = xCTML->findNameIDIndex("electrodeOutput", "", index);
    if (!xRecord) {
        ESModel_Warning("esmodel::getElectrodeOutputFile()",
	                "Could not find a node named electrodeOutput with index " + mdpUtil::int2str(index));
	return nullptr;
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
    ZZCantera::XML_Node* xEout = getElectrodeOutputFile(XMLfileName, index);
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
     ZZCantera::XML_Node root("--");
     ZZCantera::XML_Node& ct = root.addChild("ctml");
     ZZCantera::XML_Node* x_eteo = e_teo.write_ElectrodeTimeEvolutionOutput_ToXML();
     ct.addChildToTree(x_eteo);
     std::fstream s(fileName.c_str(), fstream::in | fstream::out | fstream::trunc);
     root.write(s);
     s.close();
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------


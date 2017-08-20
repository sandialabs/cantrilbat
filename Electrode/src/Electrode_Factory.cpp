/*
 * @file Electrode_Factory.cpp
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "Electrode_Factory.h"

#include "Electrode_InfCapacity.h"
#include "Electrode_MP_RxnExtent.h"
#include "Electrode_MP_RxnExtent_FeS2.h"
#include "Electrode_SimpleDiff.h"
#include "Electrode_DiffTALE.h"
#include "Electrode_SimplePhaseChangeDiffusion.h"
#include "Electrode_CSTR.h"
#include "Electrode_CSTR_MCMBAnode.h"
#include "Electrode_CSTR_LiCoO2Cathode.h"
//#include "Electrode_SuccessiveSubstitution.h"
#include "Electrode_RadialDiffRegions.h"

#include "RSD_OCVmodel.h"


//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//====================================================================================================================
/*
 *  Defining memory for a static member of the clase
 */
Electrode_Factory* Electrode_Factory::s_factory = 0;
#if defined(THREAD_SAFE_CANTERA)
boost::mutex Electrode_Factory::electrode_mutex;
#endif
//====================================================================================================================

// This is the definition for this global struct. 
Map_ETEnum_String gMap_ETEnum_String;

//====================================================================================================================
//! Routine to create the maps used in the factory lookups
/*!
 *  If you are adding an Electrode type:
 *      Add your Name and 
 */
static void create_string_maps()
{
    if (gMap_ETEnum_String.string_maps_created) {
        return;
    }
    auto& electrode_types_string = gMap_ETEnum_String.electrode_types_string;
    auto& string_electrode_types= gMap_ETEnum_String.string_electrode_types;

    gMap_ETEnum_String.string_maps_created = true;
    electrode_types_string[BASE_TYPE_ET]                       =  "BaseType";
    electrode_types_string[INF_CAPACITY_ET]                    =  "InfCapacity";
    electrode_types_string[MP_RXNEXTENT_ET]                    =  "MP_RxnExtent";
    electrode_types_string[MULTIPLATEAU_NODIFF_ET]             =  "MultiPlateau_NoDiff";
    electrode_types_string[SIMPLE_DIFF_ET]                     =  "SimpleDiff";
    electrode_types_string[DIFF_TALE_ET]                       =  "DiffTALE";
    electrode_types_string[SIMPLE_PHASE_CHANGE_DIFFUSION_ET]   =  "SimplePhaseChangeDiffusion";
    electrode_types_string[CSTR_ET]                            =  "CSTR";
    electrode_types_string[CSTR_MCMB_ANODE_ET]                 =  "CSTR_MCMBAnode";
    electrode_types_string[CSTR_LICO2_CATHODE_ET]              =  "CSTR_LiCoO2Cathode";
    electrode_types_string[SUCCESSIVE_SUBSTITUTION_ET]         =  "SuccessiveSubstitution";
    electrode_types_string[MP_RXNEXTENT_FES2_ET]               =  "MP_RxnExtent_FeS2";
    electrode_types_string[MP_RXNEXTENT_LISI_ET]               =  "MP_RxnExtent_LiSi";
    electrode_types_string[RADIAL_DIFF_REGIONS_ET]             =  "Radial_Diff_Regions";

    // create the OCV model map
    createOCVmodel_map(gMap_ETEnum_String.OCVmodel_string);

    // Invert the maps automatically.
    for (auto pos = electrode_types_string.begin(); pos != electrode_types_string.end(); ++pos) {
        string_electrode_types[pos->second] = pos->first;
        std::string lll =  ZZCantera::lowercase(pos->second);
        string_electrode_types[lll] = pos->first;
    }

    // Invert the maps automatically.
    for (auto pos = gMap_ETEnum_String.OCVmodel_string.begin(); pos != gMap_ETEnum_String.OCVmodel_string.end(); ++pos) {
        gMap_ETEnum_String.string_OCVmodel[pos->second] = pos->first;
        std::string lll =  ZZCantera::lowercase(pos->second);
        gMap_ETEnum_String.string_OCVmodel[lll] = pos->first;
    }

}
//==================================================================================================================================
// Enum to String routine for the enum Electrode_Types_Enum
/*
 *  @param etype The model of the electrode
 *
 *  @return Returns the characteristic string for that Electrode Model
 */
std::string Electrode_Types_Enum_to_string(const Electrode_Types_Enum& etype)
{
    if (!gMap_ETEnum_String.string_maps_created) {
        create_string_maps();
    }
    std::map<Electrode_Types_Enum, std::string>& electrode_types_string = gMap_ETEnum_String.electrode_types_string;

    std::map<Electrode_Types_Enum, std::string>::iterator pos = electrode_types_string.find(etype);
    if (pos == electrode_types_string.end()) {
        return "UnknownElectrodeType";
    }
    return pos->second;
}
//====================================================================================================================
// String to Enum Routine for the enum Electrode_Types_Enum
/*
 *  Matches are first made using case. Then, they are made by ignoring case
 *
 *  @param       input_string
 *
 *  @return      Returns the Enum type for the string
 */
Electrode_Types_Enum string_to_Electrode_Types_Enum(const std::string& input_string)
{
    if (!gMap_ETEnum_String.string_maps_created) {
        create_string_maps();
    }
    std::map<std::string , Electrode_Types_Enum>& string_electrode_types = gMap_ETEnum_String.string_electrode_types;

    std::map<std::string, Electrode_Types_Enum>::iterator pos = string_electrode_types.find(input_string);
    if (pos == string_electrode_types.end())  {
        std::string iii = ZZCantera::lowercase(input_string);
        pos = string_electrode_types.find(iii);
        if (pos == string_electrode_types.end())  {
            return UNKNOWN_ET;
        }
    }
    return pos->second;
}
//===================================================================================================================
int stringName_RCD_OCVmodel_to_modelID(const std::string& input_string)
{
    if (!gMap_ETEnum_String.string_maps_created) {
        create_string_maps();
    }
    std::map<std::string, int>& string_OCVmodel = gMap_ETEnum_String.string_OCVmodel;
    std::map<std::string, int>::iterator pos    = string_OCVmodel.find(input_string);
    if (pos == string_OCVmodel.end())  {
        std::string iii = ZZCantera::lowercase(input_string);
        pos = string_OCVmodel.find(iii);
        if (pos == string_OCVmodel.end())  {
            return -1;
        }
    }
    return pos->second;
}
//===================================================================================================================
std::string modelID_to_stringName_RCD_OCVmodel(int modelID)
{
    if (!gMap_ETEnum_String.string_maps_created) {
        create_string_maps();
    }
    std::map<int, std::string>& OCVmodel_string = gMap_ETEnum_String.OCVmodel_string;
    std::map<int, std::string>::iterator pos = OCVmodel_string.find(modelID);
    if (pos == OCVmodel_string.end())  {
        return "UnknownModelType";
    }
    return pos->second;
}
//====================================================================================================================
// Private constructors prevents usage
Electrode_Factory::Electrode_Factory()
{
}
//====================================================================================================================
Electrode_Factory::~Electrode_Factory()
{
}
//====================================================================================================================
// Static function that creates a static instance of the factory.
Electrode_Factory* Electrode_Factory::factory()
{
#if defined(THREAD_SAFE_CANTERA)
    boost::mutex::scoped_lock lock(electrode_mutex);
#endif
    if (!s_factory) {
        s_factory = new Electrode_Factory;
    }
    return s_factory;
}
//==================================================================================================================================
void  Electrode_Factory::deleteFactory()
{
#if defined(THREAD_SAFE_CANTERA)
    boost::mutex::scoped_lock lock(electrode_mutex);
#endif
    if (s_factory) {
        delete s_factory;
        s_factory = 0;
    }
}
//==================================================================================================================================
Electrode* Electrode_Factory::newElectrodeObject(std::string model)
{
    /*
     *  Look up the string to find the enum
     */
    Electrode_Types_Enum ieos = string_to_Electrode_Types_Enum(model);
    Electrode* ee = 0;
    /*
     *  Do the object creation
     */
    switch (ieos) {
    case          INF_CAPACITY_ET:
        ee = new    Electrode_InfCapacity();
        break;
    case          MP_RXNEXTENT_ET:
        ee = new    Electrode_MP_RxnExtent();
        break;
    case          SIMPLE_DIFF_ET:
	ee = new Electrode_SimpleDiff();
        break;
    case          DIFF_TALE_ET:
	ee = new Electrode_DiffTALE();
        break;
    case          SIMPLE_PHASE_CHANGE_DIFFUSION_ET:
        ee = new    Electrode_SimplePhaseChangeDiffusion();
        break;
    case          CSTR_ET:
        ee = new    Electrode_CSTR();
        break;
    case          CSTR_MCMB_ANODE_ET:
        ee = new    Electrode_CSTR_MCMBAnode();
        break;
    case          CSTR_LICO2_CATHODE_ET:
        ee = new    Electrode_CSTR_LiCoO2Cathode();
        break;
    case          SUCCESSIVE_SUBSTITUTION_ET:
        //ee = new    Electrode_SuccessiveSubstitution();
        ee = 0 ;
        break;
    case          MP_RXNEXTENT_FES2_ET:
        ee = new    Electrode_MP_RxnExtent_FeS2();
        break;
    case          RADIAL_DIFF_REGIONS_ET:
        ee = new    Electrode_RadialDiffRegions();
        break;
    default:
        throw Electrode_Error("Electrode_Factory::newElectrodeObject()", "Unknown Electrode model: " + model);
    }
    return ee;
}
//==================================================================================================================================
//    Create a new ELECTRODE_KEY_INPUT Object given a model name
/*
 * @param model   String to look up the model
 * @param f       ThermoFactor instance to use in matching the string
 *
 * @return
 *   Returns a pointer to a new ELECTRODE_KEY_INPUT instance matching the
 *   model string for the Electrode. Returns NULL if something went wrong.
 *   Throws an exception  if the string
 *   wasn't matched.
 */
ELECTRODE_KEY_INPUT* Electrode_Factory::newElectrodeKeyInputObject(std::string model)
{
    /*
     *  Look up the string to find the enum
     */
    Electrode_Types_Enum ieos = string_to_Electrode_Types_Enum(model);
    ELECTRODE_KEY_INPUT* ei = 0;
    /*
     *  Do the object creation
     */
    switch (ieos) {
    case          MP_RXNEXTENT_ET:
        ei = new    ELECTRODE_MP_RxnExtent_KEY_INPUT();
        break;
    case          CSTR_ET:
    case          CSTR_MCMB_ANODE_ET:
    case          CSTR_LICO2_CATHODE_ET:
        ei = new    ELECTRODE_CSTR_KEY_INPUT();
        break;

    case          BASE_TYPE_ET:
    case          SUCCESSIVE_SUBSTITUTION_ET:
    case          INF_CAPACITY_ET:
    case          SIMPLE_PHASE_CHANGE_DIFFUSION_ET:
        ei = new    ELECTRODE_KEY_INPUT();
        break;

    case          SIMPLE_DIFF_ET:
    case          DIFF_TALE_ET:
    case          RADIAL_DIFF_REGIONS_ET:
        ei = new ELECTRODE_RadialDiffRegions_KEY_INPUT();
        break;

    default:
        throw Electrode_Error("Electrode_Factory::newElectrodeKeyInputObject()",
                           "Unknown Electrode model: " + model);
        break;
    }
    return ei;
}
//==================================================================================================================================
RSD_OCVmodel* Electrode_Factory::newRSD_OCVmodel(std::string smodel)
{
    int  ieos = stringName_RCD_OCVmodel_to_modelID(smodel);
    if (ieos == -1) {
	throw Electrode_Error("Electrode_Factory::newRSD_OCVmodel()",
                           "Unknown OCVoverride model: " + smodel);
    }
    RSD_OCVmodel* ei;
    switch (ieos) {
    case         OCVAnode_MCMB2528 :
        ei = new  RSD_OCVmodel(OCVAnode_MCMB2528);
        break;
    case         OCVAnode_MCMB2528_dualfoil :
        ei = new  RSD_OCVmodel(OCVAnode_MCMB2528_dualfoil);
        break;
    case         OCVCathode_CoO2_dualfoil :
        ei = new  RSD_OCVmodel(OCVCathode_CoO2_dualfoil);
        break;
    default:
	throw Electrode_Error("Electrode_Factory::newRSD_OCVmodel()",
                           "Unknown OCVoverride model: " + smodel);
        break;
    }
    return ei;
}
//==================================================================================================================================
//  Create a new  instance of an Electrode object
/*
 * @param model   String to look up the model against
 * @param f       ThermoFactor instance to use in matching the string, Defaults to 0
 *
 * @return
 *   Returns a pointer to a new Electrode instance matching the
 *   model string. Returns NULL if something went wrong.
 *   Throws an exception UnknownThermoPhaseModel if the string
 *   wasn't matched.
 */
Electrode* newElectrodeObject(std::string model, Electrode_Factory* f)
{
    if (f == 0) {
        f = Electrode_Factory::factory();
    }
    return f->newElectrodeObject(model);
}
//==================================================================================================================================
//  Create a new ELECTRODE_KEY_INPUT Object
/*
 * @param model   String to look up the model against
 * @param f       ThermoFactor instance to use in matching the string
 *
 * @return
 *   Returns a pointer to a new ELECTRODE_KEY_INPUT instance matching the
 *   model string for the Electrode. Returns NULL if something went wrong.
 *   Throws an exception  if the string
 *   wasn't matched.
 */
ELECTRODE_KEY_INPUT* newElectrodeKeyInputObject(std::string model, Electrode_Factory* f)
{
    if (f == 0) {
        f = Electrode_Factory::factory();
    }
    return f->newElectrodeKeyInputObject(model);
}
//==================================================================================================================================
//  Create a new RSD_OCVmodel Object
/*
 * @param model   String to look up the model against
 * @param f       Eledtrode Factory instance to use in matching the string
 *
 * @return
 *   Returns a pointer to a new RSD_OCVmodel  instance matching the
 *   model string for the RSD_OCVmodel. Returns NULL if something went wrong.
 *   Throws an exception  if the string wasn't matched.
 */
RSD_OCVmodel* newRSD_OCVmodel(std::string smodel, Electrode_Factory *f)
{
   if (f == 0) {
        f = Electrode_Factory::factory();
    }
   return f->newRSD_OCVmodel(smodel);
}

//==================================================================================================================================
} // End of ZZCantera Namespace
//----------------------------------------------------------------------------------------------------------------------------------

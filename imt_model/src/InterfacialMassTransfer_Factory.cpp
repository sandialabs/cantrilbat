
/*
 * $Id: InterfacialMassTransfer_Factory.cpp 507 2013-01-07 22:48:29Z hkmoffa $
 */



#include "InterfacialMassTransfer.h"
#include "InterfacialMassTransfer_Factory.h"

#include "imtPSS_NoSurf.h"
#include "imtPSS_NoSurf_DiffBL.h"
#include "InterfacialMassTransfer_1to1Distrib.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif


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
  InterfacialMassTransfer_Factory* InterfacialMassTransfer_Factory::s_factory = 0;
#if defined(THREAD_SAFE_CANTERA)
  boost::mutex InterfacialMassTransfer_Factory::imt_mutex;
#endif
 //====================================================================================================================

  static bool string_maps_created = false;
  static std::map<IMT_Types_Enum, std::string> IMT_types_string;
  static std::map<std::string , IMT_Types_Enum> string_IMT_types;

  //====================================================================================================================
  static void create_string_maps()
  {
    if (string_maps_created) {
      return;
    }
    string_maps_created = true;
    IMT_types_string[UNKNOWN_IMT]                        =  "UnknownType";
    IMT_types_string[BASE_TYPE_IMT]                      =  "BaseType";
    IMT_types_string[PSEUDOSS_NOSURF_IMT]                =  "PseudoSS_NoSurf";
    IMT_types_string[KIN_1TO1DISTRIB_IMT]                =  "Kin_1to1Distrib";
    IMT_types_string[PSEUDOSS_NOSURF_DIFFBL_IMT]         =  "PseudoSS_NoSurf_Diffusion_BL";

    // Invert the maps automatically. 
    for (std::map<IMT_Types_Enum, std::string>::iterator pos = IMT_types_string.begin();
         pos != IMT_types_string.end(); ++pos) {
      string_IMT_types[pos->second] = pos->first;
      std::string lll =  ZZCantera::lowercase(pos->second);
      string_IMT_types[lll] = pos->first;
    }
  }
  //====================================================================================================================
  // Enum to String routine for the enum IMT_Types_Enum 
  /*
   *  @param etype The model of the imt
   *
   *  @return Returns the characteristic string for that IMT Model
   */
  std::string to_string(const IMT_Types_Enum & etype)
  {
    if (!string_maps_created) {
      create_string_maps();
    }
    std::map<IMT_Types_Enum, std::string>::iterator pos = IMT_types_string.find(etype);
    if (pos == IMT_types_string.end()) {
      return  IMT_types_string[UNKNOWN_IMT];
    }
    return pos->second;
  }
  //====================================================================================================================
  // String to Enum Routine for the enum IMT_Types_Enum 
  /*
   *  Matches are first made using case. Then, they are made by ignoring case
   *
   *  @param       input_string
   *
   *  @return      Returns the Enum type for the string
   */
  IMT_Types_Enum string_to_IMT_Types_Enum(const std::string & input_string)
  {
    if (!string_maps_created) {
      create_string_maps();
    }
    std::map<std::string, IMT_Types_Enum>::iterator pos = string_IMT_types.find(input_string);
    if (pos == string_IMT_types.end())  {
      std::string iii = ZZCantera::lowercase(input_string);
      pos = string_IMT_types.find(iii);
      if (pos == string_IMT_types.end())  {
        return UNKNOWN_IMT;
      }
    }
    return pos->second;
  }
  //====================================================================================================================
  // Private constructors prevents usage
  InterfacialMassTransfer_Factory::InterfacialMassTransfer_Factory() :
    ZZCantera::FactoryBase()
  {
  }
  //====================================================================================================================
  InterfacialMassTransfer_Factory::~InterfacialMassTransfer_Factory()
  {
  } 
  //====================================================================================================================

  /*
   * This method returns a new instance of a subclass of ThermoPhase
   */
  InterfacialMassTransfer* InterfacialMassTransfer_Factory::newInterfacialMassTransferObject(std::string model) {
    
    /*
     *  Look up the string to find the enum
     */
    IMT_Types_Enum ieos = string_to_IMT_Types_Enum(model);
    InterfacialMassTransfer * ee = 0;


    switch (ieos) {
    case BASE_TYPE_IMT:
      ee = new InterfacialMassTransfer();
      break;
    case PSEUDOSS_NOSURF_IMT:
      ee = new imtPSS_NoSurf();
      break;

    case KIN_1TO1DISTRIB_IMT:
      ee = new InterfacialMassTransfer_1to1Distrib();
      break;

    case PSEUDOSS_NOSURF_DIFFBL_IMT:
      ee = new imtPSS_NoSurf_DiffBL();
 
    default: 
      throw CanteraError("IMT_Factory::newIMTObject()",
			 "Unknown IMT model: " + model);
    }
    return ee;
  }
  //====================================================================================================================
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
  IMT_KEY_INPUT * InterfacialMassTransfer_Factory::newIMTKeyInputObject(std::string model)
  {
    /*
     *  Look up the string to find the enum
     */
    IMT_Types_Enum ieos = string_to_IMT_Types_Enum(model);
    IMT_KEY_INPUT* ei = 0;
    /*
     *  Do the object creation
     */

    switch (ieos) {
    case BASE_TYPE_IMT:
      ei = new IMT_KEY_INPUT ();
      break;

    case PSEUDOSS_NOSURF_IMT:  
      ei = new IMT_KEY_INPUT ();
      break;
    
    case KIN_1TO1DISTRIB_IMT: 
      ei = new IMT_KEY_INPUT ();
      break;
 
    default: 
      throw CanteraError("IMT_Factory::newIMTObject()",
			 "Unknown IMT model: " + model);
    }

    return ei;
  }
  //====================================================================================================================
  //  Create a new IMT_KEY_INPUT Object
  /*
   * @param model   String to look up the model against
   * @param f        instance to use in matching the string
   *
   * @return 
   *   Returns a pointer to a new IMT_KEY_INPUT instance matching the
   *   model string for the IMT object. Returns NULL if something went wrong.
   *   Throws an exception  if the string
   *   wasn't matched.
   */ 
  IMT_KEY_INPUT * newIMTKeyInputObject(std::string model, InterfacialMassTransfer_Factory* f)
  {
    if (f == 0) {
      f = InterfacialMassTransfer_Factory::factory();
    }
    return f->newIMTKeyInputObject(model);
  }
  //====================================================================================================================
 

  
//====================================================================================================================
} // 
//======================================================================================================================

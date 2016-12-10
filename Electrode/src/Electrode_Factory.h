/**
 *  @file Electrode_Factory.h
 *     Headers for the factory class that can create known %Electrode objects
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_Factory Electrode_Factory\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef ELECTRODE_FACTORY_H
#define ELECTRODE_FACTORY_H

#include "Electrode.h"
#include "ApplBase_print.h"
#include "cantera/base/FactoryBase.h"

#if defined(THREAD_SAFE_CANTERA)
#include <boost/thread/mutex.hpp>
#endif

#include <string>
#include <map>
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

class RSD_OCVmodel;

//==================================================================================================================================
//! Structure to hold the mappings between the Electrode_Types_Enum and its string representation
/*!
 *  We also keep the OCVModel mappings in this struct as well.
 *  (Tidy way to keep the material together)
 */
struct Map_ETEnum_String {

    //! default constructors 
    Map_ETEnum_String() :
        string_maps_created(false)
    {
    }
    //! true if the maps have been created
    bool string_maps_created;

    //! Map with key, Electrode_Types_Enum, and the value being the string name for the key
    std::map<Electrode_Types_Enum, std::string> electrode_types_string;

    //! This is the reverse map of the above. Maps the string as key to the Electrode_Types_Enum
    std::map<std::string, Electrode_Types_Enum> string_electrode_types;

    //! Map with key, OCVmodel number, and the value being the string representation
    std::map<int, std::string> OCVmodel_string;

    //! Map with value being string representation of OCVmodel and the value is the int representation
    std::map<std::string, int> string_OCVmodel;
};
//==================================================================================================================================
//! We define the global gMap_ETEnum_String as the struct that contains these maps.
/*!
 *  This is the declaration. The definition is in the cpp file
 */
extern Map_ETEnum_String gMap_ETEnum_String;

//==================================================================================================================================
//! Enum to String routine for the enum Electrode_Types_Enum
/*!
 *  @param[in]               etype               The model of the electrode given as an Electrode_Types_Enum
 *
 *  @return                                      Returns the characteristic string for that Electrode Model
 */
std::string Electrode_Types_Enum_to_string(const Electrode_Types_Enum& etype);

//==================================================================================================================================
//! String to Enum Routine for the enum Electrode_Types_Enum
/*!
 *  Matches are first made using case. Then, they are made by ignoring case
 *
 *  @param[in]               input_string        String representing the %Electrode type
 *
 *  @return                                      Returns the Enum type for the string
 */
Electrode_Types_Enum string_to_Electrode_Types_Enum(const std::string& input_string);

//==================================================================================================================================
//! String to int routine for the OCV override model types
/*!
 *  
 *
 *   @param[in]              input_string        The string representation for the  OCV override model.
 *
 *   @return                                     Returns the int type for the string.  Unknown models return a value of -1.
 */
int stringName_RCD_OCVmodel_to_modelID(const std::string& input_string);

//==================================================================================================================================
//!  int type to string routine for OCV model types
/*!
 * 
 *  @param[in]      modelID                      int type for the string OCV model
 *
 *  @return                                      Returns the string representation for the string OCV model.
 */
std::string modelID_to_stringName_RCD_OCVmodel(int modelID);

//==================================================================================================================================
//! Factory class for Electrode simulations managers
/*!
 *  This class keeps a list of the known Electrode simulation classes, and is used to create new instances of these classes.
 *  In order to create a new class, the user can either add to this class or inherit from this class to create their
 *  own proprietary classes.
 */
class Electrode_Factory : public ZZCantera::FactoryBase
{

public:

    //! Static function that creates a static instance of the factory.
    /*!
     *  @return                                  Returns an instance of the factory type
     */
    static Electrode_Factory* factory();

    //! Delete the static instance of this factory
    virtual void deleteFactory();

    //! Destructor doesn't do anything.
    /*!
     *  We do not delete statically created single instance of this class here, because it would create an infinite
     *  loop if the destructor is called for that single instance.
     */
    virtual ~Electrode_Factory();

    //! Create a new Electrode Object
    /*!
     * @param model  String to look up the model against
     *
     * @return    Returns a pointer to a new Electrode instance matching the  model string. Returns NULL if
     *            something went wrong. Throws an exception if the string wasn't matched.
     */
    virtual Electrode* newElectrodeObject(std::string model);

    //! Create a new ELECTRODE_KEY_INPUT Object
    /*!
     *  @param[in]           model               String to look up the model against
     *
     *  @return                                  Returns a pointer to a new ELECTRODE_KEY_INPUT instance matching the
     *                                           model string for the Electrode. Returns NULL if something went wrong.
     *                                           Throws an exception  if the string wasn't matched.
     */
    virtual ELECTRODE_KEY_INPUT* newElectrodeKeyInputObject(std::string model);

    //! Create a new RSD_OCVmodel object
    /*!
     *  @param[in]           model               String to look up the model against
     *
     *  @return                                  Returns a pointer to a new ELECTRODE_KEY_INPUT instance matching the
     *                                           model string for the Electrode. Returns NULL if something went wrong.
     *                                           Throws an exception  if the string  wasn't matched.
     */
    virtual RSD_OCVmodel* newRSD_OCVmodel(std::string model);

private:
    //! Static member of a single instance
    static Electrode_Factory* s_factory;

protected:
    //! Private constructors prevents usage
    Electrode_Factory();

private:
#if defined(THREAD_SAFE_CANTERA)
    //! Declaration for locking mutex for electrode factory singelton
    static boost::mutex electrode_mutex;
#endif
};
//==================================================================================================================================
//! Create a new Electrode object
/*!
 *  @param[in]               model               String to look up the model against
 *  @param[in]               f                   Electrode Factory instance to use in matching the string
 *
 *  @return                                      Returns a pointer to a new ThermoPhase instance matching the model string. 
 *                                               Returns NULL if something went wrong. 
 *                                               Throws an exception if the string wasn't matched.
 */
Electrode* newElectrodeObject(std::string model, Electrode_Factory* f = 0);

//==================================================================================================================================
//! Create a new ELECTRODE_KEY_INPUT object
/*!
 *  @param[in]               model               String to look up the model against
 *  @param[in]               f                   Electrode Factory instance to use in matching the string
 *
 *  @return                                      Returns a pointer to a new ELECTRODE_KEY_INPUT instance matching the
 *                                               model string for the Electrode. Returns NULL if something went wrong.
 *                                               Throws an exception if the string wasn't matched.
 */
ELECTRODE_KEY_INPUT* newElectrodeKeyInputObject(std::string model, Electrode_Factory* f = 0);

//==================================================================================================================================
//! Create a new RSD_OCVmodel object
/*!
 *  @param[in]               model               String to look up the model against
 *  @param[in]               f                   Electrode Factory instance to use in matching the string
 *
 *  @return                                      Returns a pointer to a new RSD_OCVmodel  instance matching the
 *                                               model string for the RSD_OCVmodel. Returns NULL if something went wrong.
 *                                               Throws an exception if the string wasn't matched.
 */
RSD_OCVmodel* newRSD_OCVmodel(std::string model, Electrode_Factory *f = 0);

//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

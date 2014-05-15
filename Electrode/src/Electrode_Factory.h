/**
 *  @file ThermoFactory.h
 *     Headers for the factory class that can create known %ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 *
 */

/*
 * $Author: hkmoffa $
 * $Revision: 584 $
 * $Date: 2013-04-02 18:29:38 -0600 (Tue, 02 Apr 2013) $
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

namespace Cantera
{

struct Map_ETEnum_String {
    Map_ETEnum_String() :
        string_maps_created(false)
    {
    }
    bool string_maps_created;
    std::map<Electrode_Types_Enum, std::string> electrode_types_string;
    std::map<std::string , Electrode_Types_Enum> string_electrode_types;
};
extern Map_ETEnum_String gMap_ETEnum_String;

//! Enum to String routine for the enum Electrode_Types_Enum
/*!
 *  @param etype The model of the electrode
 *
 *  @return Returns the characteristic string for that Electrode Model
 */
std::string Electrode_Types_Enum_to_string(const Electrode_Types_Enum& etype);

//! String to Enum Routine for the enum Electrode_Types_Enum
/*!
 *  Matches are first made using case. Then, they are made by ignoring case
 *
 *  @param       input_string
 *
 *  @return      Returns the Enum type for the string
 */
Electrode_Types_Enum string_to_Electrode_Types_Enum(const std::string& input_string);

//! Factory class for thermodynamic property managers.
/*!
 * This class keeps a list of the known ThermoPhase classes, and is
 * used to create new instances of these classes.
 */
class Electrode_Factory : public Cantera::FactoryBase
{

public:

    //! Static function that creates a static instance of the factory.
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

    //!  Create a new ELECTRODE_KEY_INPUT Object
    /*!
     * @param model   String to look up the model against
     * @param f       ThermoFactor instance to use in matching the string
     *
     * @return
     *   Returns a pointer to a new ELECTRODE_KEY_INPUT instance matching the
     *   model string for the Electrode. Returns NULL if something went wrong.
     *   Throws an exception  if the string
     *   wasn't matched.
     */
    virtual ELECTRODE_KEY_INPUT* newElectrodeKeyInputObject(std::string model);

private:
    //! static member of a single instance
    static Electrode_Factory* s_factory;

protected:
    //! Private constructors prevents usage
    Electrode_Factory();

private:
#if defined(THREAD_SAFE_CANTERA)
    //! Decl for locking mutex for electrode factory singelton
    static boost::mutex electrode_mutex;
#endif

};


//!  Create a new Electrode Object
/*!
 * @param model   String to look up the model against
 * @param f       ThermoFactor instance to use in matching the string
 *
 * @return
 *   Returns a pointer to a new ThermoPhase instance matching the
 *   model string. Returns NULL if something went wrong.
 *   Throws an exception UnknownThermoPhaseModel if the string
 *   wasn't matched.
 */
Electrode* newElectrodeObject(std::string model, Electrode_Factory* f = 0);

//!  Create a new ELECTRODE_KEY_INPUT Object
/*!
 * @param model   String to look up the model against
 * @param f       ThermoFactor instance to use in matching the string
 *
 * @return
 *   Returns a pointer to a new ELECTRODE_KEY_INPUT instance matching the
 *   model string for the Electrode. Returns NULL if something went wrong.
 *   Throws an exception  if the string
 *   wasn't matched.
 */
ELECTRODE_KEY_INPUT* newElectrodeKeyInputObject(std::string model, Electrode_Factory* f = 0);

}
#endif

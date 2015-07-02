/**
 *  @file EState_XML.h
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#ifndef _ESTATE_XML_H
#define _ESTATE_XML_H

#include "Electrode_defs.h"
#include "EState.h"

#include "cantera/base/xml.h"
#include "cantera/base/FactoryBase.h"
#include <string>
#include <vector>
//---------------------------------------------------------------------------------------------------------------------------------
namespace Cantera {
class Electrode;

}
//----------------------------------------------------------------------------------------------------------------------------------
namespace esmodel
{

//! This structure will hold globally defined variables associated with maps for EState types
/*!
 *  There will only be one of these
 */
struct Map_ESEnum_String {
    //! constructor
    Map_ESEnum_String() :
        string_maps_created(false)
    {
    }
    //! boolean indicating that we have created the variable
    bool string_maps_created;

    //! map between EState and string
    std::map<Cantera::EState_Type_Enum, std::string> estate_types_string;

    //! Map between string and EState_Type_Enum
    std::map<std::string , Cantera::EState_Type_Enum> string_estate_types;
};





//! Singleton variable for this structure.
/*!
 *   We use the prefix gMap to indicate that it is global
 */
extern Map_ESEnum_String gMap_ESEnum_String;


//! Enum to String routine for the enum EState_Type_Enum
/*!
 *  @param[in]  estype The model of the electrode state
 *
 *  @return Returns the characteristic string for that EState Model
 */
std::string EState_Type_Enum_to_string(const Cantera::EState_Type_Enum& estype);

//! String to Enum Routine for the enum EState_Type_Enum
/*!
 *  Matches are first made using case. Then, they are made by ignoring case
 *
 *  @param       input_string
 *
 *  @return      Returns the Enum type for the string
 */
Cantera::EState_Type_Enum string_to_EState_Type_Enum(const std::string& input_string);

//==================================================================================================================================
//! Factory class for EState object managers.
/*!
 * This class keeps a list of the known ThermoPhase classes, and is
 * used to create new instances of these classes.
 */
class EState_Factory : public Cantera::FactoryBase
{

public:

    //! Static function that creates a static instance of the factory.
    static EState_Factory* factory();

    //! Delete the static instance of this factory
    virtual void deleteFactory();

    //! Destructor doesn't do anything.
    /*!
     *  We do not delete statically created single instance of this class here, because it would create an infinite
     *  loop if the destructor is called for that single instance.
     */
    virtual ~EState_Factory();

    //! Create a new EState Object
    /*!
     * @param model  String to look up the model against
     *
     * @return    Returns a pointer to a new EState instance matching the  model string. Returns NULL if
     *            something went wrong. Throws an exception if the string wasn't matched.
     */
    virtual Cantera::EState* newEStateObject(std::string model);

private:
    //! static member of a single instance
    static EState_Factory* s_factory;

protected:
    //! Protected default constructor
    EState_Factory();

private:
#if defined(THREAD_SAFE_CANTERA)
    //! Declaration for locking mutex for electrode factory singelton
    static boost::mutex state_mutex;
#endif
};
//==================================================================================================================================
//!  Create a new EState Object
/*!
 * @param model   String to look up the model against
 * @param f       EState Factory instance to use in matching the string
 *
 * @return             Returns a pointer to a new EState instance matching the   model string. Returns NULL if something went wrong.
 *                     Throws an exception if the string wasn't matched.
 */
Cantera::EState* newEStateObject(std::string model, EState_Factory* f = 0);

//====================================================================================================================================


Cantera::XML_Node* getElectrodeOutputFile(const std::string& fileName, int index);


//! Given an Electrode solution file, select a particular global time step number given by the index
/*!
 *
 *  
 */
Cantera::XML_Node* selectLastGlobalTimeStepIncrement(Cantera::XML_Node* xSoln, int& globalTimeStepNum);



Cantera::XML_Node* locateTimeLast_GlobalTimeStepFromXML(const Cantera::XML_Node& xmlGlobalTimeStep, double& timeVal,
							int printSteps = 0);

bool get_Estate_Indentification(const Cantera::XML_Node& xSoln, Cantera::EState_ID_struct & e_id);

}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

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

class ETimeState {

public:

    ETimeState();

    ETimeState(const ETimeState& r);

    //! Destructor
    ~ETimeState();

    ETimeState& operator=(const ETimeState& r);

    //!  Compare the current state of this object against another guest state to see if they are the same
    /*!
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
    bool compareOtherTimeState(const ETimeState* const ETSguest, double molarAtol, int nDigits,
			       bool includeHist = false, int printLvl = 0) const;


    // Cell number of the electrode object
    int cellNumber_;

    //! domain number of the electrode object
    int domainNumber_;

    Cantera::EState* es_;

    std::string stateType_;

    //! type of the time increment
    /*!
     *  Two possible values: "global" or "local"
     */
    std::string timeIncrType_;

    //! Time within the simulation corresponding to the state
    double time_;
    
    //! whether I own the EState object
    bool iOwnES_;
    
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
//!  Read an XML Electrode output file and create an XML tree structure
/*!
 *
 *    @param[in]         fileName            File name of 
 *
 *
 *    @return          Returns the pointer to the XML tree. If the file can't be found or is the wrong type, a NULL pointer is
 *                     returned.
 */
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

//!   Read an EState XMLfile returning the last time step and its corresponding time
/*!
 *    @param[in]    XMLfileName     File name of the XML electrode solution object
 * 
 *    @param[out]   timeRead        corresponding time of the solution, read from the file
 *
 *    @return                       Return a malloced EState object of the appropriate form
 *                                  with the state of the elctrode at t_final in that Estate object.
 */
Cantera::EState* readEStateFileLastStep(const std::string& XMLfileName, double& timeRead);

}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

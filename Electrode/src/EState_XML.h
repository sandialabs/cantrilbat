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

//!  Class representing an EState time state object
/*!
 *
 */
class ETimeState {

public:

    //! Constructor
    ETimeState();

    //! Copy constructor
    /*!
     *   @param[in]      r   Object to be copied
     */
    ETimeState(const ETimeState& r);

    //! Destructor
    ~ETimeState();

    //! Assignment operator
    /*!
     *   @param[in]      r   Object to be copied
     *  
     *   @return             Returns a reference to the current object
     */
    ETimeState& operator=(const ETimeState& r);

    //! Create/Malloc an XML Node containing the timeState data contained in this object
    /*!
     *   @return   Returns the malloced XML_Node with name timeState containing the information in this
     *             object. The calling program is responsible for freeing this.
     */
    Cantera::XML_Node* write_ETimeState_ToXML() const;

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


    //! Cell number of the electrode object
    int cellNumber_;

    //! domain number of the electrode object
    int domainNumber_;

    //! Base class pointer for the solution at each time
    /*!
     *          Note, this may be overwritten by child objects
     */
    Cantera::EState* es_;

    //! Type of the state
    /*!
     *      Possible types
     *                t_init
     *                t_intermediate
     *                t_final
     */
    std::string stateType_;

    //! Type of the time increment
    /*!
     *  Two possible values: "global" or "local"
     *
     *   Right now, we are practically only using global
     */
    std::string timeIncrType_;

    //! Time within the simulation corresponding to the state
    double time_;
    
    //! Boolean indicating whether I own the EState object
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

//!  Read an XML Electrode output file and create an XML tree structure
/*!
 *    File doesn't throw on failure. Instead it returns a NULL pointer.
 *
 *    @param[in]         fileName            File name of XML file
 *    @param[in]         index               Index of the electrodeOutput node that is requested. The default is "1", which
 *                                           is the first index written.
 *
 *    @return          Returns the pointer to the XML tree. If the file can't be found or is the wrong type, a NULL pointer is
 *                     returned. The top XML_Node is set at the electrodeOutput node corresponding to the specified input index
 */
Cantera::XML_Node* getElectrodeOutputFile(const std::string& fileName, int index);

//! Given an Electrode solution file, select the last global time step number returning its XML element
/*!
 *   @param[in]      xElectrodeOutput   Input XML tree containing the electrodeOutput XML_Node.
 *   @param[out]     globalTimeStepNum  Returns the global time step number index selected. This is the last global time step
 *                                      number in the input file
 *
 *   @return                            Returns the XML_Node corresponding to the last global time step in the solution file.
 *                                      The node will have a named called, globalTimeStep, and will have the largest index
 *                                      in the electrodeOutput XML element.                  
 */
Cantera::XML_Node* selectLastGlobalTimeStepInterval(Cantera::XML_Node* xElectrodeOutput, int& globalTimeStepNum);

//! Given a global time step interval XML tree, this routine will locate the t_final Electrode Time State
/*!
 *   @param[in]   xmlGlobalTimeStep     XML node corresponding to the global time step interval
 *   @param[out]  timeVal               Value of the time read
 *   @param[in]   printSteps            If nonzero it prints the type of steps and the time value. THe default is zero.
 *    
 *   @return                            Returns a pointer to the  t_final Electrode Time State xml tree. If there is problem,
 *                                      it returns NULL.
 */
Cantera::XML_Node* locateTimeLast_GlobalTimeStepIntervalFromXML(const Cantera::XML_Node& xmlGlobalTimeStep, double& timeVal,
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

//!   Create a new EState object from an XML_Node and EState id structure
Cantera::EState* newEStatefromXML(const Cantera::XML_Node& XeState, const Cantera::EState_ID_struct& e_id);


}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

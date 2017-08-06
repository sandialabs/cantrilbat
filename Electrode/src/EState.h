/**
 *  @file EState.h
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#ifndef _ESTATE_H
#define _ESTATE_H

#include "Electrode_defs.h"


#include <string>
#include <vector>
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

class Electrode;
class XML_Node;
class Electrode_Factory;

//! Enum indicating the type of Electrode object output file that is written.
enum EState_Type_Enum {
    EST_UNKNOWN_TYPE = -1,
    //! Electrode State file containing CSTR information for the explicit State of the electrode.
    /*!
     * The chief state variable is the number of moles of each species in the electrode phase.
     * Other properties include states of charge variables and volumetric variables.
     * The idea is that they all represent a single state that are compatible with another
     * and the idea is that particular models that implement a CSTR-type model can pick and
     * choose which ones to impliment as state variables.
     */
    EST_CSTR = 0,
    EST_CSTR_LiCoO2Cathode,
    EST_CSTR_MCMBAnode,
    EST_MULTIPLATEAU,
    EST_RADIALDISTRIB
};
//==================================================================================================================================
//!  Structure that holds the identification of the Electrode object
/*!
 *   This structure contains the information about the type of electrode object and the file type for the EState object.
 *   It also contains the information about the domain number and cell number, which uniquely identifies the electrode
 *   object in a multi-cell simulation.
 */
struct EState_Identification 
{
    //! Constructor
    EState_Identification();

    //! Write the XML_Node element that contains this structure's information
    /*!
     *   @return   Returns a pointer to the malloced XML_Node containing the structure's information
     */
    ZZCantera::XML_Node* writeIdentificationToXML() const;

    //! Read this structure's information from the XML_Node
    /*!
     *   @param[in]     xmlEI                    Reference to the XML Node which stores the Identification information
     */
    void readIdentificationFromXML(const XML_Node& xmlEI);

    /*  ------------------------------------------- Data ----------------------------------------- */

    //! Electrode Type String
    /*!
     *  This string is the one used in the Electrode Factory routines. It identifies the model
     *  used for the Electrode Object uniquely.
     */
    std::string electrodeTypeString_;

    //! enum type for the EState. This is used in the factory routine
    enum ZZCantera::EState_Type_Enum  EST_Type_;

    //! String used to identify the EState type. This is used in the factory routine for the EState object.
    /*!
     *  Note this can be different than the electrodeTypeString_
     */
    std::string EState_Type_String_;

    //! Version number of the file
    int EST_Version_;

    //! Integer describing a version model number for a chemistry model
    int electrodeChemistryModelType_;
    
    //! Domain number of the electrode (starts from 0 with the usual convention that 0 refers to anode and 2 refers to cathode)
    int electrodeDomainNumber_;

    //! Cell number of the electrode (unique within each domain)
    int electrodeCellNumber_;

    //! Capacity type of the electrode. 
    /*!
     *    There is a basic sign issue with capacity that differs between anodes and cathodes.
     *    Capacity left goes does when an anodic electrode gives off an electron but increases for a cathode electrode
     */
    ZZCantera::Electrode_Capacity_Type_Enum  electrodeCapacityType_;
};

//! Typedef for the EState_ID structure
typedef struct EState_Identification EState_ID_struct;

//==================================================================================================================================
//! Base Class for the Electrode State class concept. We define the state of the electrode here
//! which can be used to set the Electrode object classes and can be used to write out to an XML file.
//! This is the main class involved with saves and restarts of the Electrode object.
/*!
 *     This is a glue class that allows one to write out states of Electrode objects.
 *     There is a direct correspondence with members of this class with an XML_Node object that can be written out to a file.
 *
 *     The class also keeps track of the identification information for the electrode object.
 *     The identification information can be used as header information for the electrode object.
 *     There is an instantiation of the header XML_Node block within this class as well.
 *
 *     Unless otherwise stated the units of all terms are in MKS.
 *
 *     How this object interacts with the electrode object
 *    ----------------------------------------------------
 *
 *    The function EState::initialize(Electrode *e) will initialize this object.
 *
 *    This object is a friend to the Electrode objects that it writes into. Therefore,
 *    it can use and call protected functions from these Electrode objects.
 *
 *    This function will set the _final_ state of the Electrode object.
 *    This function will also set all of the internal degrees of freedom such that
 *    the next time step taken by the object will result in the same values and time step
 *    that would have occurred if there hadn't been a restart. In particular this requirement
 *    means that the state information contain information required to choose the size of the
 *    next time step.
 *
 *    It also contains information about the initial amount discharged, even though that
 *    actually doesn't influence the state or the step calculations.
 *
 *     Several concepts of version control
 *     -----------------------------------------
 *
 *     The file format for this class is given in the pair of EST_fileToBeWritten_
 *     and  const integer EST_Version_fileToBeWritten_.
 *     This identifies the fields and XML format of the written information.
 *
 *     XML files that are just read have the variables  EST_lastFileRead_ and
 *     EST_Version_lastFileRead_ assigned from the files. These values are reconciled with the
 *     EST_fileToBeWritten_ and  EST_Version_fileToBeWritten_ pair to determine compatibility.
 *     Right now, if they are not the same, an error is generated.
 *
 *     The string variable electrodeTypeString_ contains the Factory method string for the
 *     instanteation of the Electrode object used to create the state information. This is
 *     different than the EST information because the file format services more than one
 *     Electrode object class.
 *
 */
class EState
{

public:

    //! Default constructor for the base object
    EState();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    EState(const EState& right);

    //! Destructor is virtual
    virtual ~EState();

    //! Assignment operator
    /*!
     *  @param right  object to be duplicated
     *
     *  @return                                  Returns a reference to the current object
     */
    EState& operator=(const EState& right);

    //! Duplicator function for this class
    /*!
     *  @return Returns a duplication of the current state as a pointer to the base class
     */
    virtual EState* duplMyselfAsEState() const;

    //! Initialize the EState object based on an electrode Base class
    /*!
     *   This call will initialize all of the arrays within this class.
     *   All of the species and phase identification information is created and the class is
     *   readied for use as a state maintainer.
     *
     *  @param[in]           e                   Pointer to the electrode base class
     *
     *  @return                                  Returns 0
     */
    virtual int initialize(const ZZCantera::Electrode* const e);

    //! Returns a string representing the electrode type 
    /*!
     *  Strings are listed in Electrode_Factory.cpp and other child classes of that.
     *
     *  @return                                  Returns a string representing the Electrode type
     */
    virtual const std::string& electrodeType() const;

    //! Create an indentification XML_Node element for this Electrode EState object
    /*!
     *  @return                                  Returns a malloced XML_Node tree containing the identification information
     */
    virtual XML_Node* writeIdentificationToXML() const;

    //! Create an indentification XML_Node element for the PhaseList Object
    /*!
     *  @return                                  Returns a  malloced XML_Node tree containing the identification information
     */
    virtual XML_Node* write_PhaseListID_ToXML() const;

    //! Write the electrodeState contained within the EState object to a new malloced XML_Node tree
    /*!
     *  (virtual function from EState)
     *
     *  This function creates/mallocs an XML_Node tree containing the contents of the electrodeState.
     *  This is a virtual function, because the base class is called from general code, allowing
     *  the child classes to be invoked.
     *
     *  @return       Pointer to the XML_Node tree that is malloced. The calling routine is responsible for
     *                freeing the XML_Node tree. 
     */
    virtual XML_Node* write_electrodeState_ToXML() const;

    //! Read the state from the XML_Node tree given by the argument
    /*!
     *  @param[in]           xmlRoot             Root of the xml tree to get the information from
     */
    void readStateFromXMLRoot(const XML_Node& xmlRoot);

    //! Read the state from the XML_Node given by the argument
    /*!
     *  (Virtual function from EState) -> main way to get the process going and the child functions called.
     *
     *  @param[in]           xmlEState           electrodeState XML element to be read from.
     *                                           The electrodeState XML element contains the state of the electrode
     *                                           at a particular time, though the time is not part of this record.
     */
    virtual void readStateFromXML(const XML_Node& xmlEState);

    //! Read identification information from the XML header record
    /*!
     *  @param[in]           xmlEState           reference to the XML node.
     */
    void readIdentificationFromXML(const XML_Node& xmlEState);

    //! Read identification information from a struct EState_Identification object
    /*!
     *  @param[in]           es_ID               Reference to the EState ID structure 
     */
    void readIdentificationFromStruct(const EState_ID_struct& es_ID);

    //! Set the State of this object from the state of the Electrode object
    /*!
     *  (virtual function)
     *
     *  virtual function, because the base class is called from general code, allowing
     *      the child classes to be invoked.
     *
     *  This function takes the electrode objects _final_ state and copies it into this object.
     *  This function must be carried out before the XML_Node tree is written.
     *
     *  @param e   Pointer to the Electrode object. Note, this class may use dynamic casting
     *             to choose a child object, and then may invoke an error if the match isn't
     *             correct.
     */
    virtual void copyElectrode_intoState(const ZZCantera::Electrode* const e);

    //! Set the state of the Electrode from the state of this object
    /*!
     *  (virtual function)
     *
     *   Virtual function -> main way to get the process going and the child functions called.
     *
     *   This function is called from the Electrode object or a battery simulator to
     *   initialize the Electrode object from an input file.
     *
     *   This function will set the _final_ state of the Electrode object.
     *   This function will also set all of the internal degrees of freedom such that
     *   the next time step taken by the object will result in the same values and time step
     *   that would have occurred if there hadn't been a restart.
     *
     *  @param  e  Changeable pointer to the base class Electrode object. This function may
     *             do dynamic casting to get the correct child Electrode object.
     */
    virtual void setStateElectrode_fromEState(ZZCantera::Electrode* const e) const;

    //! Set the state of the Electrode Class from the state of the EState object
    /*!
     *  (virtual function)
     *
     *  virtual function, because the base class is called from general code, allowing
     *      the child classes to be invoked.
     *
     *  This is not a virtual function.  It copies information from one fixed class to another
     *  using the "friend" paradigm.
     *
     *  @param e   Pointer to the Electrode object. Note, this class may use dynamic casting
     *             to choose a child object, and then may invoke an error if the match isn't
     *             correct.
     */
    void copyEState_toElectrode(ZZCantera::Electrode* const e) const;

    //! Print a heading
    /*!
     *   @param[in]   printLvl            Level of printing
     *
     *    @return                         Returns the number of lines printed.
     */
    int printHead(int printLvl) const;

    //! Print a difference between two strings
    /*!
     *   @param[in]        vexp                String containing the expression
     *   @param[in]        significant         Causes a comparison failure
     *   @param[in]        val                 String containing the val
     *   @param[in]        gval                String containing the guest val
     *   @param[in]        printLvl            print lvl
     */
    void printDiff(const std::string& vexp, bool significant,  const std::string& val, const std::string& gval, int printLvl) const;

    //! Print a difference between two int values, with the string representation of the variable
    /*!
     *   @param[in]        vexp                String containing the expression
     *   @param[in]        index               Int containing the index pertaining to the expression (used in printing)
     *   @param[in]        val                 int containing the val
     *   @param[in]        gval                int containing the guest val
     *   @param[in]        printLvl            print lvl
     */
    void printDiff(const std::string& vexp, int index, int val, int gval, int printLvl) const;

    //! Print a difference between two double values, with the string representation of the variable
    /*!
     *   @param[in]        vexp                String containing the expression
     *   @param[in]        index               Int containing the index pertaining to the expression (used in printing)
     *   @param[in]        val                 double containing the val
     *   @param[in]        gval                double containing the guest val
     *   @param[in]        printLvl            print lvl
     */
    void printDiff(const std::string& vexp, int index, double val, double gval, int printLvl) const;

    //! Print a difference between two double vector values, with a string representation of the variable
    /*!
     *   @param[in]        vexp                String containing the expression
     *   @param[in]        val                 vector of double containing the val
     *   @param[in]        gval                vector of double containing the guest val
     *   @param[in]        printLvl            print lvl
     */
    void printVecDiff(const std::string& vexp, const std::vector<double>& val, const std::vector<double>& gval, int printLvl) const;

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
    virtual bool compareOtherState(const EState* const ESguest, double molarAtol, int nDigits, 
				   bool includeHist = false, int printLvl = 0) const;

    //! Get the number of moles in the electrode in kmol
    /*!
     *  @return                                  Returns the number of moles
     */
    virtual double electrodeMoles() const;

    /* ------------------------------------------------ D A T A ----------------------------------------------------------------- */

protected:

    //! Constant reference to the Electrode object that this object refers to.
    const Electrode* eRef_;

public:
    //! Electrode Model
    /*!
     *    String identification of the model used to instantiate the Electrode object.
     *    This is the string used in the factory method for creating the Electrode object.
     */
    std::string electrodeTypeString_;

protected:
    //! Value of the Electrode State type file that was previously read
    /*!
     *  In almost all cases this is equal to the value to be written.
     */
    enum EState_Type_Enum  EST_lastFileRead_;

    //! Version number of the file that was last read
    int EST_Version_lastFileRead_;

    //! Electrode State type file to be written
    /*!
     *   Currently this is EST_CSTR
     */
    enum EState_Type_Enum EST_fileToBeWritten_;

    //! Version number of the file related to this Object.
    /*!
     *  Current this is 1
     */
    const int EST_Version_fileToBeWritten_;

    //  STATE INFORMATION

    //! Number of moles of each species in each phase at the end of each
    //! subcycle of the integration step
    /*!
     *  INDEPENDENT STATE VARIABLE FOR THE ELECTRODE PROBLEM
     *   Number of moles of each species in each phase at the end of each
     *   subcycle of the integration step.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMoles_;

    //! Vector of phase voltages
    /*!
     *  length = Number of total phases in PhaseList
     *  indexing of PhaseList
     */
    std::vector<double> phaseVoltages_;

    //! Current value of the temperature (Kelvin)
    double temperature_;

    //! Pressure (Pa)
    double pressure_;

    //! Electrode model type
    /*!
     *   Electrode model type is used to set the capacity coefficients
     *
     *      0  Undetermined type
     *      1  LiSi anode with an initial composition of pure Li13Si4(S) and a final state of pure Si(s).
     *      2  FeS2 anode with an initial composition of pure FeS2(S) and a final state of Li2S(s) and Fe(s).
     *      3  LiSi anode with one plateau and interstitial diffusion model
     *      4  FeS2 cathode with a simplified thermo
     */
    int electrodeChemistryModelType_;

    //! Domain number of the electrode object
    /*!
     *  This number refers to the domain number within 1DElectrode.
     */
    int electrodeDomainNumber_;

    //! Cell number within the domain
    int electrodeCellNumber_;

    //! Number of Particles to follow
    /*!
     *   All the extrinsic properties of the object are multiplied by this value.
     *   This is the number of particles that are in the electrode object
     */
    double particleNumberToFollow_;

    //! Total volume of the electrode solid phase (compared to the electrolyte)
    /*!
     * units of m**3
     */
    double electrodeSolidVolume_;

    //! Gross volume of the electrode (m3)
    /*!
     * Total volume of the electrode.
     */
    double grossVolume_;

    //!  Radius of the exterior of the particle
    /*!
     *  This is at the start of the global time step
     */
    double radiusExterior_;

    //! Vector of the surface area for each Reacting Surface in the electrode.
    /*!
     *  Each surface is assumed to have a surface area. The surface area is a
     *  constitutive function  of the composition (i.e., State of Charge) of the electrode.
     *  In most constitutive models only one external reacting surface will be present at
     *  any one time.  This vector is over internal and external surfaces.
     *
     *  length = number of external surfaces that may be present
     *  units m**2
     *
     *  This is the final value at each time step
     */
    std::vector<double> surfaceAreaRS_;

    //! Total number of moles in the solid phases of the electrode.
    double electrodeMoles_;

    //! Capacity type of the electrode.
    Electrode_Capacity_Type_Enum  electrodeCapacityType_;

    //! Total capacity left in the battery (coulombs)
    double capacityLeft_;

    //! Initial capacity (coulombs)
    double capacityInitial_;

    //! Depth of discharge (coulombs)
    double depthOfDischarge_;

    //! Starting depth of Discharge;
    double depthOfDischargeStarting_;

    //! Initial conditions for the electrode in terms of the initial discharged capacity
    /*!
     *  The default value for this is -1. If default, then this value is used.
     *  If not default, then this is the relative amount of electrons that are discharged
     *  as a function of the initial mole number of active species in the electrode.
     *  Some objects don't support setting the initial conditions by this method.
     *
     *  This translates directly into the value, e RelativeExtentRxn_ used in the MP_RxnExtent object; it's the same thing.
     *
     *  There are no units for this, as
     */
    double relativeElectronsDischargedPerMole_;

    //! Relative depth of discharge (unitless)
    double relativeDepthOfDischarge_;

    //! Capacity discharged to date -> this is a number that is dependent on the past time history of the simulation
    double capacityDischargedToDate_;

    //! kmol of electrons that are discharged to date 
    /*!
     *   This is a number that is dependent on the past time history of the simulation. So past simulations must provide
     *   initial values of this to the current simulation
     */
    double electronKmolDischargedToDate_;

    //! Initial value of the next subcycle deltaT
    double deltaTsubcycle_init_next_;

    //! Statement that the Electrode class can access any information in this class
    /*!
     *  NOTE, I'm not sure that this direction of access is needed ATM.
     */
    friend class ZZCantera::Electrode;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

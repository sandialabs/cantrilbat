/**
 *  @file Electrode_input.h 
 *     Headers for the declarations of the base Electrode_iput class, used to accumulate input information read in
 *     from an Electrode input file before transfering it to the Electrode object
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_input Electrode_input\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_INPUT_H
#define _ELECTRODE_INPUT_H

#include "cantera/base/config.h"
#include "Electrode_defs.h"

#include "cantera/multiphase/PhaseList.h"
#include "BE_BlockEntry.h"

#include "cantera/base/ctml.h"
#include <string>
#include <vector>
/*
 *-----------------------------------------------------------------------------
 *
 * Include file containing constant declarations for inputs to
 * mpequil
 *
 *-----------------------------------------------------------------------------
 */
#define MPEQUIL_MAX_NAME_LEN_P1 81
#define MPEQUIL_MAX_NAME_LEN    80

//-----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

class Electrode;
class EGRInput;

//===================================================================================================================================
//! This class contain information about the electrode input from the input file
/*!
 *   We use this as a temporary holding area for variable inputs from the input file. The ELECTRODE_KEY_INPUT class
 *   then shorts it out.
 */
class ElectrodeBath
{
public:

    //! Pointer to the PhaseList
    PhaseList* m_pl;

    //! Species mole fractions in one phase
    double* XmolPLSpecVec;

    //! Species molalities in one phase
    double* MolalitiesPLSpecVec;

    //! Capacity Left coefficients in each phase of the phase list
    double** CapLeftCoeffPhases;

    //! Capacity Left coefficients in one phase
    double* CapLeftCoeffSpecVec;

    //! Capacity  coefficients in each phase of the phase list
    double** CapZeroDoDCoeffPhases;

    //! Capacity Left coefficients in one phase
    double* CapZeroDoDCoeffSpecVec;

    //! total moles of given phase
    std::vector<double> PhaseMoles;

    //! total mass of given phase
    std::vector<double> PhaseMass;

    //! Default constructor
    /*!
     *  @param[in]           pl                  Reference to the PhaseList object
     */
    ElectrodeBath(PhaseList* pl_ptr = nullptr);

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    ElectrodeBath(const ElectrodeBath &right);

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *  @return                                  Returns a reference to the current object
     */
    ElectrodeBath& operator=(const ElectrodeBath& right);

    //! Destructor
    ~ElectrodeBath();
};
//===================================================================================================================================
//! Structure for storring the input for OCVoverride models
/*!
 *
 */
struct OCV_Override_input {

    //! Default constructor
    OCV_Override_input();

    //! Copy constructor
    /*!
     *  @param[in]           right               Object to be copied
     */
    OCV_Override_input(const OCV_Override_input& right);

    //! Assignment operator
    /*!
     *  @param[in]           right               Object to be copied
     *  @return                                  Returns a reference to the current object
     */
    OCV_Override_input& operator=(const OCV_Override_input& right);

    //! Destructor
    ~OCV_Override_input();

    //! Number of times
    int numTimes;

    //! PhaseID of the surface, where the reactions are taking place
    int surfacePhaseID;

    //! String name of the OCV model
    std::string OCVModel;

    //! String name of the bulk species whose thermo is to be replaced
    std::string replacedSpeciesName;
   
    //! The global species id for the species whose thermo will be replaced
    int replacedGlobalSpeciesID;

    //! Local value of the species ID in the ThermoPhase object whose thermo is being modified
    int replacedLocalSpeciesID;

    //! PhaseID within the ReactingSurface object of the phase containing the replaced species
    int replacedSpeciesPhaseID;

    //! OCV_Format defaults to 0
    int OCV_Format_;

    //! Name of the species which will be used as a surrogate for the depth of discharge.
    std::string DoDSurrogateSpeciesName;

    //! Identity of the species whose mole fraction will be identified as being equivalent to the Depth of discharge variable
    size_t MF_DoD_LocalSpeciesID;

    //! ReactionID of the reaction whose Delta Thermo has been fit to experiment
    int rxnID;

    //! Reaction whose deltaS value is modeled by the dOCVdT calculation. This is usually the
    //! open circuit voltage for the Full-cell reaction or the half cell consisting of this reaction
    //! and the half-cell for the reference electrode reaction.
    int rxnID_deltaS;

    //! Temperature derivative type
    int temperatureDerivType;

    //! Value of the temperature at which this formula is based
    double temperatureBase;

    //! Value for the dOCVdT. Sometimes this is just set to a constant value
    double temperatureDerivValue;

    //! String name for the temperature derivative model for the OCV
    std::string OCVTempDerivModel;
};
//===================================================================================================================================
//! Storage for Command file input
/*!
 *  This is the current command file specification of the problem statement.
 *  We use the class to accumulate input information read in from an Electrode input file before transfering it to the 
 *  Electrode object.
 *
 *  All member data are public, so that we can use this data source easily
 *
 *  The  basic idea is to read the input file three times. 
 *    The first time we read in the cantera files.
 *    The second time we initialize the ThermoPhases
 *    The third time, we set up the input deck so that we can read in all of the information and do error checking
 *    on the phases and names that are input.
 */
class ELECTRODE_KEY_INPUT
{
public:

    //! Constructor
    /*!
     *  @param[in]           printLvl            print level. This defaults to 0
     */
    ELECTRODE_KEY_INPUT(int printLvl = 0);

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    ELECTRODE_KEY_INPUT(const ELECTRODE_KEY_INPUT &right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     *  @return                                  Returns a reference to the current object
     */
    ELECTRODE_KEY_INPUT& operator=(const ELECTRODE_KEY_INPUT& right);

    //! Virtual destructor
    virtual ~ELECTRODE_KEY_INPUT();

    //! Initialize the fields within the Electrode for reading
    /*!
     *
     *  In this routine, we fill in the fields in the ELECTRODE_KEY_INPUT
     *  structure that will be needed to parse the keyword input.
     *  We also size the vectors in the structure appropriately, now
     *  we know the extent of the problem.
     *  The information to do this has already been gathered into the PhaseList
     *  structure
     *
     *  These are:
     *         nTotPhases
     *         nTotSpecies
     *         nTotElements
     *         SpeciesNames
     *         PhaseNames
     *         ElementNames
     *         CanteraFNSurface
     *
     *  The sized fields include:
     *         PhaseInclude
     *         MoleNumber
     *         MoleNumberIG
     *         ElementAbundances
     *
     *  @param pl   Pointer to the phase list
     */
    virtual void InitForInput(const ZZCantera::PhaseList* const pl);

    //! setup for pass1
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be set up by this routine
     */
    virtual void setup_input_pass1(BEInput::BlockEntry* cf);

    //! setup for pass2
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be added onto by this routine
     */
    virtual void setup_input_pass2(BEInput::BlockEntry* cf);

    //! setup for pass3
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be added onto by this routine
     */
    virtual void setup_input_pass3(BEInput::BlockEntry* cf);

    //!   Initialize some of the fields in the ELECTRODE_KEY_INPUT structure by
    //!   post processing some of the data entry
    /*!
     *     Fields that are filled in and changed by this routine
     *         MoleNumber[]
     *         MoleFraction[]
     *
     *  @param cf  Block entry used to parse the input file.
     *
     *   @return   Returns 0
     */
    int post_input_pass3(const BEInput::BlockEntry* cf);

    //! setup for child1
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be added onto by this routine
     */
    virtual void setup_input_child1(BEInput::BlockEntry* cf);

    //! setup for child2
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be added onto by this routine
     */
    virtual void setup_input_child2(BEInput::BlockEntry* cf);

    //! setup for child3
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be added onto by this routine
     */
    virtual void setup_input_child3(BEInput::BlockEntry* cf);


    //! Post processing done after the input
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be added onto by this routine
     */
    virtual void post_input_child1(BEInput::BlockEntry* cf);

    //! Post processing done after the input
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be added onto by this routine
     */
    virtual void post_input_child2(BEInput::BlockEntry* cf);

    //! Post processing done after the input
    /*!
     *  Virtual function that may get added onto
     *
     *  @param cf Pointer to the BlockEntry to be added onto by this routine
     */
    virtual void post_input_child3(BEInput::BlockEntry* cf);

    //! Parse the input file and fill up the current structure with the result
    /*!
     *   Calls the following functions in order
     *      setup_input_pass1
     *      process_electrode_input
     *      setup_input_pass2
     *      process_electrode_input
     *      InitForInput(pl);
     *      setup_input_pass3(cf);
     *      process_electrode_input()
     *      electrode_model_init(cf);
     *
     *   This fills in all of the Base Class information
     *
     *  @param commandFile  Current command file
     *  @param cf           Block entry input
     *
     *  @return                                Returns 0 for success
     *                                         returns -1 for failure 
     */
    int electrode_input(std::string commandFile, BEInput::BlockEntry* cf);

    //! Parse the input file and fill up the current structure with the result
    /*!
     *    Calls the following functions in order
     *      setup_input_pass1
     *      process_electrode_input
     *      setup_input_pass2
     *      process_electrode_input
     *      InitForInput(pl);
     *      setup_input_pass3(cf);
     *      process_electrode_input()
     *      electrode_model_init(cf);
     *      setup_input_child1(cf);
     *      process_electrode_input();
     *      post_input_child1();
     *
     *  @param commandFile  Current command file
     *  @param cf           Block entry input
     *
     *  @return                                Returns 0 for success
     *                                         returns -1 for failure 
     */
    virtual int electrode_input_child(std::string commandFile, BEInput::BlockEntry* cf);

    // -------------------------------------------------- D A T A -----------------------------------------------------

    //! PhaseList object
    /*!
     *   This includes all of the phases, "period". In particular this includes the surface phases
     */
    ZZCantera::PhaseList* m_pl;

    //! Print level
    /*!
     *     0 no printout
     *     1 error printout
     *     2 processed lines
     *     3 last pass printout
     *     4 processed lines and last pass printout
     *     5 all pass printout
     */
    int printLvl_;

    //! Command file that is used to initialize the object
    std::string commandFile_;

    //! Shallow pointer not owned
    BEInput::BlockEntry* lastBlockEntryPtr_;

    //! This is assigned to the first surface phase found in the PhaseList
    //std::string CanteraFNSurface;

    //! Number of Cantera files to be read in
    int NumberCanteraFiles;

    //! Pointer vector of Cantera file names
    char** CanteraFileNames;

    //! Temperature of the electrode (kelvin)
    double Temperature;

    //! Pressure of the electrode (Pascal)
    double Pressure;

    //! Electrode Bath Structure
    ElectrodeBath m_BG;

    //! Vector of mole numbers of the species
    std::vector<double> MoleNumber;

    //! Vector of mole fractions of the species
    std::vector<double> MoleFraction;

    //! Vector of the electric potential of the PhaseList Phases
    std::vector<double> PotentialPLPhases;

    //! Problem type
    int ProblemType;

    //! Vector of species names
    /*!
     *  c string format
     *  length: nTotSpecies
     */
    char** SpeciesNames;

    //! Names of the phases
    char** PhaseNames;

    //!  Names of the elements
    /*!
     *  List of the names of the elements in C-style, nullterminated strings
     */
    char** ElementNames;

    //! Storage for the OCV override information
    /*
     *   This is stored as a series of 
     */
    std::vector<ZZCantera::OCV_Override_input *> OCVoverride_ptrList;

    //! level of the xml State information created
    /*!
     *       0 - none (default)
     *       1 - minimal -> writes the end of the run and can initialize
     *       2 - Write global time step information
     *       3 - Can write intermediate and global time step information
     */
    int xmlStateInfoLevel;

    //! Electrode Model Name
    /*!
     * This is used in the factory routine to instanteate the obejct
     */
    std::string electrodeModelName;

    //! Electrode type
    /*!
     *  0 is an anode and 1 is a cathode
     */
    int electrodeCapacityType;

    //!  Name of the electrode when used in printouts
    /*!
     * This is the name of the electrode to be used in printouts to identify
     * the electrode model and chemistry
     */
    std::string electrodeName;

    //! Characteristic particle diameter
    /*!
     *  Characteristic particle diameter of particle in the electrode.
     *  This will yield the surface area.
     *  units = m
     */
    double particleDiameter;

    //! Integral number of particles to follow that make up this electrode
    /*!
     *  This is actually an integer.
     */
    double particleNumberToFollow;

    //! Electrode Gross Area
    /*!
     *   Gross area of the electrode to consider. This will be combined with the width
     *   and the porosity information to calculate a total gross volume that electrode occupies
     *
     * units = m**2
     */
    double electrodeGrossArea;

    //! Electrode Diameter
    /*!
     *   Diameter of the electrode to consider. This will be combined with the width
     *   and the porosity information to calculate a net amount of electrode material
     *
     * units = m**2
     */
    double electrodeGrossDiameter;

    //! Electrode Gross thickness
    /*!
     *   total macroscopic width of the electrode. This will be combined with the gross area
     *   and the porosity information to calculate a total gross volume that the electrode
     *   domain occupies.
     *
     * units = m
     */
    double electrodeGrossThickness;

    //! Porosity
    /*!
     *  Fraction of the available volume filled with electrolyte. The rest is assumed
     *  to be filled with electrode material. We don't have the capability yet to have
     *  multiple types of electrode materials.
     */
    double porosity;

    //! Maximum reaction extent at the top 
    /*!
     *  Maximum amount of the RxnExt variable that the electrode can be charged to
     *  -> defaults to -1.0, which indicates that there is no limit on this variable
     */
    double RxnExtTopLimit;

    //! Maximum reaction extent at the bottom
    /*!
     *  Maximum amount of the RxnExt variable that the electrode can be discharged to
     *  -> defaults to -1.0, which indicates that there is no limit on this variable
     */
    double RxnExtBotLimit;

    //! Number of extra global reactions
    int numExtraGlobalRxns;

    //! Vector of ExtraGlobalReactions that will be filled in by the MultiBlockVec input structure
    /*!
     *   The pointers will be allocated and filled in by the BE_MultiBlockVec routine
     */
    std::vector<EGRInput*> m_EGRList;

    //! Number of total phases in the phase list
    size_t nTotPhases;

    //! Total number of species in the Phase List
    size_t nTotSpecies;

    //! Total number of elements
    size_t nTotElements;

    //! Initial conditions for the electrode in terms of the initial discharged capacity
    /*!
     *  The default value for this is -1. If default, then this value is used.
     *  If not default, then this is the relative amount of electrons that are discharged
     *  as a function of the initial mole number of active species in the electrode.
     *  Some objects don't support setting the initial conditions by this method.
     */
    double RelativeCapacityDischargedPerMole;

    //! Maximum number of subGlobal time step iterations
    int maxNumberSubGlobalTimeSteps;

    //! relative minimum time step ratio
    double relativeLocalToGlobalTimeStepMinimum;
};
//===================================================================================================================================
//! Utility class to store the reaction index and reaction multiplier for an entry in formulating a global reaction
/*!
 *  Global reactions are made up of a linear combination of local reactions. Here we specify the contribution from
 *  one elementary reaction.
 *  
 *  Destructor, Copy constructor and assigment operator all use the default compiler ones.
 */
struct ERSSpec
{
public:
    //! Index value of the reaction index that is contributing to the global reaction
    int m_reactionIndex;

    //! Reaction multiplier for the reaction index
    double m_reactionMultiplier;

    //! Default constructor
    ERSSpec() :
        m_reactionIndex(-1),
        m_reactionMultiplier(0.0) {
    }
};
//===================================================================================================================================
//! Small all-public struct that constructs a global reaction out of a linear combination of single-step reactions
struct EGRInput
{
public:

    //! Index of the Reacting Surface domain where the kinetics object is located
    size_t m_RSD_index;
   
    //! Kinetic species that indicates the species to follow through the reaction
    /*!
     * 
     */
    int m_SS_KinSpeciesKindex;

    //! Number of single-step reactions in the vector
    int m_numElemReactions;

    //! List of ERSSpec structs which specifies which elementary reaction contributes to the
    //! global vector and what the reaction multiplier is for each reaction
    /*!
     *   Length:  m_numElemReactions
     */
    std::vector<ERSSpec*> m_ERSList;

    //! Constructor
    EGRInput();

    //! Copy constuctor
    /*!
     *  @param[in]           right               Object to be copied
     */
    EGRInput(const EGRInput& right);

    //! assignment operator
    /*!
     *  @param[in]           right               Object to be copied
     *  @return                                  Returns a reference to the object
     */
    EGRInput& operator=(const EGRInput& right);

    //! destructor
    ~EGRInput();
};
//===================================================================================================================================
//! Set up the default conditions for a single ThermoPhases within the Electrode object
/*!
 *  @param[in]               tp                  Reference to the ThermoPhase 
 *  @param[in]               ei                  reference to the ELECTRODE_KEY_INPUT object
 *  @param[in]               BG                  Reference to the ElectrodeBath object which contains default initial conditions
 *  @param[in]               iph                 Phase index within the PhaseList
 *  @param[in]               printLvl            Print level   
 */
void setElectrodeBathSpeciesConditions(ZZCantera::thermo_t_double& tp, ELECTRODE_KEY_INPUT& ei, ElectrodeBath& BG, 
                                       size_t iph, int printLvl);
//===================================================================================================================================
//! Process the electrode input deck using a file
/*!
 *  We process the input deck from a file. We catch Block input errors. Therefore, the return bool should be
 *  checked for error status. 
 *
 *  @param[in]               cf                  Pointer to the block entry structure used to parse the input file
 *  @param[in]               fileName            file name
 *  @param[in]               printFlag           print flag
 *  @param[in]               pass                Current value of the number of passes through the input file
 *
 *  @return                                      Returns true if the processing didn't produce an error
 *                                               Returns false if there was an error.
 */
bool process_electrode_input(BEInput::BlockEntry* cf, std::string fileName, int printFlag = 0, int pass = 0);

//===================================================================================================================================
} // end of namespace
//-----------------------------------------------------------------------------------------------------------------------------------
#endif


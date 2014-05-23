/*
 * $Id: Electrode_input.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_INPUT_H
#define _ELECTRODE_INPUT_H

#include "PhaseList.h"
//  adding this file here because most applications will use the object
#include "BE_BlockEntry.h"

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

#define MPEQUIL_SUCCESS 0

namespace Cantera
{
class Electrode;
}


//! This class contain information about the electrode input from the input file
/*!
*   TODO: There is no reason for this class. It should be transfered into the ELECTRODE_KEY_INPUT class
*/
class ElectrodeBath
{
public:
    //! species mole fractions in one phase
    double* XmolPLSpecVec;
    //! species mole fractions in each phase of phase list
    double** XmolPLPhases;
    //! species molalities in one phase
    double* MolalitiesPLSpecVec;
    //! species molalities in each phase of phase list
    double** MolalitiesPLPhases;

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

    ElectrodeBath();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    ElectrodeBath(const ElectrodeBath &right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    ElectrodeBath& operator=(const ElectrodeBath& right);

    ~ElectrodeBath();
};

namespace Cantera {
//! Structure for storring the input for OCVoverride models
struct OCV_Override_input {
    OCV_Override_input();
    OCV_Override_input(const OCV_Override_input& right);
    OCV_Override_input& operator=(const OCV_Override_input& right);
    ~OCV_Override_input();

    int numTimes;
    int surfacePhaseID;
    std::string surfacePhaseName;
    std::string OCVModel;
    std::string replacedSpeciesName;
    //! the global species id for the species whose thermo will be replaced
    int replacedGlobalSpeciesID;
    int replacedLocalSpeciesID;
    int replacedSpeciesPhaseID;
    std::string DoDSurrogateSpeciesName;
    int MF_DoD_LocalSpeciesID;
    int rxnID;
    int temperatureDerivType;
    double temperatureBase;
    std::string OCVTempDerivModel;
};
}

class EGRInput;


//! storage for Command file input
/*!
 * This is the current command file specification of the problem statement.
 */
class ELECTRODE_KEY_INPUT
{
public:

    //! Constructor
    ELECTRODE_KEY_INPUT(int printLvl = 0);

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    ELECTRODE_KEY_INPUT(const ELECTRODE_KEY_INPUT &right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
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
    virtual void InitForInput(const Cantera::PhaseList*   const pl);

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
     */
    virtual int electrode_input_child(std::string commandFile, BEInput::BlockEntry* cf);

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
    std::string CanteraFNSurface;

    //! Number of Cantera files to be read in
    int NumberCanteraFiles;

    //! Pointer vector of Cantera file names
    char** CanteraFileNames;

    //! Temperature of the electrode (kelvin)
    double Temperature;

    //! Pressure of the electrode (Pascal)
    double Pressure;

    double Vol;

    //! Electrode Bath Gas Structure
    ElectrodeBath* m_BG;


    double* MoleNumber;
    double* MoleFraction;
    double* PotentialPLPhases;
    int*    PhaseInclude;
    int    ProblemType;

    //! Vector of species names
    /*!
     *  c string format
     *  length nTotSpecies
     */
    char** SpeciesNames;

    char** PhaseNames;

    //!  Names of the elements
    /*!
     *  List of the names of the elements in C-style, nullterminated strings
     */
    char** ElementNames;
    double* ElementAbundances;
    bool  specifiedElementAbundances;
    int   specifiedBlockKmolSpecies;

    //! Storage for the OCV override information
    /*
     *   This is stored as a series of 
     */
    std::vector<Cantera::OCV_Override_input *> OCVoverride_ptrList;

    //! level of the xml State information created
    /*!
     *       0 - none (default)
     *       1 - minimal -> writes the end of the run and can initialize
     *       2 - Write global time step information
     *       3 - Can write intermediate and global time step information
     */
    int   xmlStateInfoLevel;

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


    double RxnExtTopLimit;
    double RxnExtBotLimit;

    int numExtraGlobalRxns;
    EGRInput** m_EGRList;

    //! Number of total phases in the phase list
    int nTotPhases;

    //! Total number of species in the Phase List
    int nTotSpecies;

    //! Total number of elements
    int nTotElements;

    //! Initial conditions for the electrode in terms of the initial discharged capacity
    /*!
     *  The default value for this is -1. If default, then this value is used.
     *  If not default, then this is the relative amount of electrons that are discharged
     *  as a function of the initial mole number of active species in the electrode.
     *  Some objects don't support setting the initial conditions by this method.
     */
    double RelativeCapacityDischargedPerMole;

    //! PhaseList object
    /*!
     * this includes all of the phases, "period".
     *  In particular this includes the surface phases
     */
    Cantera::PhaseList* m_pl;


};


class ERSSpec
{
public:
    int m_reactionIndex;
    double m_reactionMultiplier;
    ERSSpec() :
        m_reactionIndex(-1),
        m_reactionMultiplier(0.0) {
    }

    ~ERSSpec() {
    }

    ERSSpec(const ERSSpec& right) {
        m_reactionIndex      = right.m_reactionIndex;
        m_reactionMultiplier = right.m_reactionMultiplier;
    }
    ERSSpec& operator=(const ERSSpec& right) {
        if (&right == this) {
            return *this;
        }
        m_reactionIndex      = right.m_reactionIndex;
        m_reactionMultiplier = right.m_reactionMultiplier;
        return *this;
    }
};

class EGRInput
{
public:
    int m_SS_KinSpeciesKindex;
    int m_numElemReactions;
    ERSSpec** m_ERSList;

    EGRInput();

    EGRInput(const EGRInput& right);

    EGRInput& operator=(const EGRInput& right);
    ~EGRInput();
};

void setElectrodeBathSpeciesConditions(Cantera::ThermoPhase& g,
                                       ELECTRODE_KEY_INPUT& EI, ElectrodeBath& BG, int iph, int printLvl);



bool process_electrode_input(BEInput::BlockEntry* cf, std::string fileName, int printFlag = 0,
                             int pass = 0);

#endif
/*****************************************************************************/

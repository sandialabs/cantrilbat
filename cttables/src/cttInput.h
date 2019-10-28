/**
 *  @file cttInput.h
 *
 */


/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CTTINPUT_H
#define CTTINPUT_H

#include "zuzax/kinetics.h"
#include "mdp_allo.h"
#include "zuzax/multiphase/PhaseList.h"

namespace BEInput
{
class BlockEntry;
}

/*
 * Use Zuzax for namespace identification
 */

#define UNITS_KCAL_CGS 0
#define UNITS_KJOULE   1
#define UNITS_CGS      2

void setup_input_pass1(BEInput::BlockEntry* cf,
                       Zuzax::Kinetics* g_kin_ptr,
                       Zuzax::ThermoPhase* g_ptr);

void setup_input_pass2(BEInput::BlockEntry* cf,
                       Zuzax::Kinetics* g_kin_ptr,
                       Zuzax::ThermoPhase* g_ptr);

void setup_input_pass3(BEInput::BlockEntry* cf,
                       Zuzax::Kinetics* g_kin_ptr,
                       Zuzax::ThermoPhase* g_ptr,
                       Zuzax::PhaseList* pl);

int process_input(BEInput::BlockEntry* cf, std::string commandFile,
                  Zuzax::Kinetics* g_kin_ptr,
                  Zuzax::ThermoPhase* g_ptr ,
                  Zuzax::PhaseList* pl);


struct ERSSpec {
    int m_reactionIndex;
    double m_reactionMultiplier;
    ERSSpec() :
        m_reactionIndex(-1),
        m_reactionMultiplier(0.0)
    {
    }
};

struct EGRInput {
    int m_SS_KinSpeciesKindex;
    int m_numElemReactions;
    struct ERSSpec** m_ERSList;
    EGRInput() :
        m_SS_KinSpeciesKindex(0),
        m_numElemReactions(0),
        m_ERSList(0)
    {
        m_ERSList = (struct ERSSpec**) mdpUtil::mdp_alloc_ptr_1(2);
        m_ERSList[0] = new ERSSpec();
    }
    ~EGRInput()
    {
        struct ERSSpec** ptr;
        for (ptr = m_ERSList; *ptr != 0; ptr++) {
            delete *ptr;
        }
        mdpUtil::mdp_safe_free((void**) &m_ERSList);
    }
};

class IOoptions
{
public:
    /*
     * Constructor():
     */
    IOoptions();
    /*
     * Destructor():
     */
    ~IOoptions();

    void reprep();

    /*
     * This routine sets up the dimensions for arrays within the
     * structure. Initializations are done as well.
     */
    void Initialize(int nSpecies)
    {
        PrintThermoTable = mdpUtil::mdp_alloc_int_1(nSpecies, ProcessAll);
    }

    /**
     * -------------- DATA -------------------------
     */

    int NumberCanteraFiles;

    char** CanteraFileNames;

    /*
     * This is list of options and parameters that are specified
     * through the input deck
     *
     * ProcessAll: set the global processing default.
     *             If true, the default is to process all options.
     *             If false, the default is to do nothing.
     */
    int ProcessAll;
    /*
     *  Boolean to print a thermo table for each species
     *   (length = number of species)
     */
    int* PrintThermoTable;
    /*
     * OutputUnits
     */
    int OutputUnits;

    /*
     * ChemicalPotColumn
     *  If true, thermo tables have an extra column containing the
     *  absolute value of the chemical potential
     */
    bool ChemPotColumn;

    /*
     * IntEngColumn
     * If true, thermo tables have an extra column containing the
     * relative value of the internal energy
     */
    bool IntEnergyColumn;

    //! Skip the transport calculations if true
    /*!
     * The default is false
     */
    bool SkipTransport;

    /*
     * Number of points in the temperature table, not including
     * the explicitly added points
     */
    int    m_TTnpts;
    /*
     * DeltaT for the points in the temperature table
     */
    double m_TTDeltaT;
    bool   TTinc298;
    double m_TTTlow;
    int    NumAddedTemps;
    double* AddedTemperatures;

    /*
     * Number of points in the voltage table, not including
     * the explicitly added points
     */
    int    m_VVnpts;
    /*
     * DeltaV for the points in the voltage table
     */
    double m_VVDeltaV;
    bool   VVincZero;
    bool   VVincEzero;
    bool   VVincEeq;
    double m_VVVlow;
    int    NumAddedVoltages;
    double* AddedVoltages;

    /*
     * Boolean to indicate whether reference pressure or bath pressure
     * is used in the calculation of the Thermodynamics Tables.
     */
    bool UseRefPressureInThermoTables;

    int* PhaseInclude;
    double* MoleNumber;
    char** PhaseNames;
    char** ElementNames;
    double* ElementAbundances;

    //! List of Kinetics Species Lists.
    /*!
     * length is total number of species, nTotSpecies
     */
    char** SpeciesNames;

    int nTotPhases;

    //! Total number of species in the kinetic species list
    int nTotSpecies;
    int nTotElements;

    void InitForInput(Zuzax::PhaseList* pl);

    int numExtraGlobalRxns;
    struct EGRInput** m_EGRList;
};

extern IOoptions IOO;


#endif

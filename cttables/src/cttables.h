/**
 *  @file cttables.h
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CTTABLES_H
#define CTTABLES_H

#include "cantera/base/config.h"
#include "cantera/transport.h"
#include "cttInput.h"

#include "cantera/thermo.h"

#include "TemperatureTable.h"
#include "VoltageTable.h"

#include <cstdio>


#ifdef useZuzaxNamepace
#ifndef ZZCantera
#define ZZCantera Zuzax
#endif
#else
#ifndef ZZCantera
#define ZZCantera Cantera
#endif
#endif

const double R_kcalmol = ZZCantera::GasConstant / 4.184E6;
const double R_Jgmol = ZZCantera::GasConstant * 1.0E-3;
const double R_kJgmol = ZZCantera::GasConstant * 1.0E-6;

extern ZZCantera::Transport* GTran;

/*
 * Turn on debug printing from the command line.
 */
extern int DebugPrinting;

struct Bath {
    double Temperature;
    double Pressure;
    double* Xmol;
    double* XmolPLSpecVec;
    double** XmolPLPhases;
    double* Molalities;
    double* MolalitiesPLSpecVec;
    double** MolalitiesPLPhases;
    double* PotentialPLPhases;
    std::vector<double> PhaseMoles;
    int BathSpeciesID;
    int* BathSpeciesIDVec;
    std::string BathSpeciesName;
    std::vector<std::string> BathSpeciesNameVec;
    int MajorGasID;
    std::string MajorGasName;
    double Density;
    Bath() :
        Temperature(298.15),
        Pressure(1.01325E5),
        Xmol(0),
        XmolPLSpecVec(0),
        XmolPLPhases(0),
        Molalities(0),
        MolalitiesPLSpecVec(0),
        MolalitiesPLPhases(0),
        PotentialPLPhases(0),
        BathSpeciesID(0),
        BathSpeciesIDVec(0),
        MajorGasID(1),
        Density(0.0)
    {

    }
    ~Bath()
    {
        mdpUtil::mdp_safe_free((void**) &Xmol);
        mdpUtil::mdp_safe_free((void**) &XmolPLSpecVec);
        mdpUtil::mdp_safe_free((void**) &XmolPLPhases);
        mdpUtil::mdp_safe_free((void**) &Molalities);
        mdpUtil::mdp_safe_free((void**) &MolalitiesPLSpecVec);
        mdpUtil::mdp_safe_free((void**) &MolalitiesPLPhases);
        mdpUtil::mdp_safe_free((void**) &PotentialPLPhases);
        mdpUtil::mdp_safe_free((void**) &BathSpeciesIDVec);
    }
} ;
extern Bath BG;
extern TemperatureTable* TT_ptr;
extern VoltageTable* VV_ptr;

extern bool TopIsKineticsObject;

struct UnitsIO {
    std::string sGibbs;
    double mGibbs;
    std::string sEntropy;
    double mEntropy;
    int unitDef;

    void setup(int unit)
    {
        unitDef = unit;
        if (unit == UNITS_KCAL_CGS) {
            sGibbs = "kcal/gmol";
            mGibbs = R_kcalmol;

            sEntropy = "cal/mol*K";
            mEntropy = R_kcalmol * 1.0E3;

        } else if (unit == UNITS_KJOULE) {

            sGibbs   = " kJ/gmol ";
            mGibbs   = R_kJgmol;

            sEntropy = "J/gmolK";
            mEntropy = R_kJgmol * 1.0E3;

        } else {
            printf("ERROR unknown units");
            exit(-1);
        }
    }
};
extern UnitsIO UIO;

//==================================================================================================================================

//! Reference state Entropy at 298 for the current ThermoPhase
/*!
 *  Length: Number of species in ThermoPhase
 */
extern std::vector<double> g_S298_refThermo;

//==================================================================================================================================


/**
 *  This routine will print out a table of information about a species in an ideal gas thermo phse. It explicitly
 *  assumes that a multitransport object has been created for  the phase, and it presumes a NASA polynomial form for the
 *  species thermodynamics.
 */
extern void printIdealGasSpeciesTable(Zuzax::ThermoPhase& g, int k, TemperatureTable& TT, Zuzax::DenseMatrix& Cp_Table,
                               Zuzax::DenseMatrix& Hrel_Table, Zuzax::DenseMatrix& Grel_Table, Zuzax::DenseMatrix& S_Table,
                               bool haveSpeciesTransportProps, Zuzax::DenseMatrix& Visc_Table, Zuzax::DenseMatrix& Cond_Table,
                               Zuzax::DenseMatrix& Diff_Table, std::vector<double> H298, Zuzax::DenseMatrix& Urel_Table,
                               std::vector<double> U298);


/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

extern void print_char(const char c, const int nTimes);

extern void pr_sf(const std::string& s, const int w);
extern void pr_df(const double d, const int w, const int p);
extern void pr_dfp(const double d, const int p);
extern void pr_de(const double d, const int w, const int p);
extern void pr_dg(const double d, const int w, const int p);
extern void dnt(const int i);

extern void print_map(const std::map<std::string,double>& m,
                      const std::string& prefix);

extern void setAllBathSpeciesConditions(ZZCantera::PhaseList* pl);
extern void printAllBathSpeciesConditions(ZZCantera::PhaseList* pl);
extern void setBathSpeciesConditions(ZZCantera::ThermoPhase& g,
                                     ZZCantera::PhaseList* pl,
                                     int printLvl);
extern void printBathSpeciesConditions(ZZCantera::ThermoPhase& g,
                                       ZZCantera::PhaseList* pl,
                                       int printLvl);

// Functions in cttables_thermo.cpp

//! Gather the entropy of the elements of a species at 298 K.
/*!
 *  This is useful for going back and forth from the gibbs free energy of formation and the absolute gibbs free energy 
 *  in NIST format. 
 *  This routine calls entropyElement298(m, true) to get the absolute entropy of the element in its standard state at 298 K 
 *  and 1 bar
 *
 *  @param[in]               g_ptr               ThermoPhase pointer
 *  @param[in]               k                   Species index
 *
 *  @return                                      Total entropy of the elements that make up the ThermoPhase
 *                                               Returns ENTROPY298_UNKNOWN if any of the elements that make up the species 
 *                                               has an unknown entropy.
 */
double entropyElem298(ZZCantera::ThermoPhase* g_ptr, size_t k);

//==================================================================================================================================
//! Print out volume information about the Standard State
/*!
 *  @param[in]               tp                  Pointer to the ThermoPhase
 *  @param[in]               k                   species Index
 */
void printVolSpecies(ZZCantera::ThermoPhase* tp, size_t k);
//==================================================================================================================================

extern void printThermoCoeffSpecies(ZZCantera::ThermoPhase*, int);
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#endif

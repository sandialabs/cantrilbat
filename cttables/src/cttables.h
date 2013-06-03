/**
 *  @file cttables.h
 *
 * $Id: cttables.h 497 2013-01-07 21:17:04Z hkmoffa $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CTTABLES_H
#define CTTABLES_H

#include "IdealReactingGas.h"
#include "cantera/transport.h"
#include "cttInput.h"
#include "mdp_allo.h"

#include "TemperatureTable.h"
#include "VoltageTable.h"

const double R_kcalmol = 8.314472E7 /  4.184E7 * 1.0E-3;
const double R_Jgmol = 8.314472;
const double R_kJgmol = 8.314472 * 1.0E-3;


extern Cantera::Transport *GTran;


#include <vector>

/*
 * Turn on debug printing from the command line.
 */
extern int DebugPrinting;

struct Bath {
  double Temperature;
  double Pressure;
  double * Xmol;
  double * XmolPLSpecVec;
  double **XmolPLPhases;
  double * Molalities;
  double * MolalitiesPLSpecVec;
  double **MolalitiesPLPhases;
  double *PotentialPLPhases;
  std::vector<double> PhaseMoles;
  int BathSpeciesID;
  int * BathSpeciesIDVec;
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
  ~Bath() {
    mdpUtil::mdp_safe_free((void **) &Xmol);
    mdpUtil::mdp_safe_free((void **) &XmolPLSpecVec);
    mdpUtil::mdp_safe_free((void **) &XmolPLPhases);
    mdpUtil::mdp_safe_free((void **) &Molalities);
    mdpUtil::mdp_safe_free((void **) &MolalitiesPLSpecVec);
    mdpUtil::mdp_safe_free((void **) &MolalitiesPLPhases);
    mdpUtil::mdp_safe_free((void **) &PotentialPLPhases);
    mdpUtil::mdp_safe_free((void **) &BathSpeciesIDVec);
  }
} ;
extern Bath BG;
extern TemperatureTable *TT_ptr;
extern VoltageTable *VV_ptr;

extern bool TopIsKineticsObject;

struct UnitsIO {
    std::string sGibbs;
    double mGibbs;
    std::string sEntropy;
    double mEntropy;
    int unitDef;

    void setup(int unit) {
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

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

extern void print_char(const char c, const int nTimes);

extern void print_bool(const bool b);

extern void pr_if(const int i, const int w);

extern void pr_sf(const std::string s, const int w);
extern void pr_sf_lj(const std::string s, const int w, const int crop = 0);
extern void pr_df(const double d, const int w, const int p);
extern void pr_dfp(const double d, const int p);
extern void pr_de(const double d, const int w, const int p);
extern void pr_dg(const double d, const int w, const int p);
extern void dnt(const int i);

extern void print_map(const std::map<std::string,double>& m, 
		      const std::string& prefix);

extern void setAllBathSpeciesConditions(Cantera::PhaseList *pl);
extern void printAllBathSpeciesConditions(Cantera::PhaseList *pl);
extern void setBathSpeciesConditions(Cantera::ThermoPhase& g, 
				     Cantera::PhaseList *pl,
				     int printLvl);
extern void printBathSpeciesConditions(Cantera::ThermoPhase& g, 
		  		       Cantera::PhaseList *pl,
				       int printLvl);
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


extern void printVolSpecies(Cantera::ThermoPhase *, int);

extern void printThermoCoeffSpecies(Cantera::ThermoPhase *, int);
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#endif

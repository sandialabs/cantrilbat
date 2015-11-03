/**
 *  @file Cu2S_models.h
 */
/*
 *  $Id: Cu2S_models.h 506 2013-01-07 22:43:59Z hkmoffa $
 */
#ifndef CU2S_MODELS_H
#define CU2S_MODELS_H

#include "cantera/transport.h"


#include "ReactionRate.h"

#include "md_timer.h"
#include "mdp_allo.h"

#include "m1d_ProblemResidEval.h"

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdio>
#include <fstream>

#include <limits.h>
#include <math.h>

#define SAFE_DELETE(a) if (a) { delete (a); a = 0; }
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))

using namespace std;
//using namespace BEInput;
//using namespace TKInput;
using namespace mdpUtil;
/*
 *************** Global objects ********************************
 *
 * Ideal gas mixture
 */

/*
 * Setup options for the program. -> These get initialized
 * via interaction with the user through the input file.
 */
struct ProgramOptions
{

  /*
   *  Model Parameters Block
   */
  int Dim;
  int NumNodes;
  int IncludeGasSpecies;
  /*
   * BCTop:
   *        0 -> nondimensional linear kinetics growth
   *        1 -> Cu2S growth kinetics.
   */
  int BCTop;
  /*
   * BCBot:
   *        0 -> Dirichlet condition -> set C(x=xbot0) = C0
   *        1 -> Cu2S growth kinetics.
   */
  int BCBot;
  /*
   *  ProblemType -
   *        0 -> nondimensional linear kinetics growth
   *        1 -> Cu2S growth kinetics.
   *        3
   *        4 -> default 
   */
  int ProblemType;

  /*
   * Time Step Parameters Block
   */
  int Debugging;
  double DeltaTInitial;
  double TFinal;
  double OutputDeltaTime;
  int PrintFlag;
  int MaxNumTimeSteps;
  double RTol;
  double ATol;
  int numInitialConstantDeltaTSteps;
  int printSolnSteps;
  int printSolnInterval;
  int printSolnFirstSteps;

  /*
   * Physical Parameters Block
   */
  double Temperature;
  // Pressure of the gas -> currently not used
  double Pressure;
  double Xbot0; // Initial radius
  double Xtop0; // Initial radius of the top of the domain
  double Thickness0; // Initial Thickness. 
  double Cv_init;
  double k1;
  double K1_equil;
  double Conc_Ce;
  double Diff_Coeff;
  double Conc_Solid;
  double CV0_xbot;

  double P1_A1; // Prob 1 preexponential factor for rxn 1
  double P1_Aneg1; // Prob 1 preexponential factor for rxn -1
  double P1_A2; // Prob 1 preexponential factor for rxn 2
  double P1_XH2S; // Prob 1 H2S mole fraction at the surface

  double P3_AuThick; // Prob 3 thickness of gold plating
  double P3_AuPoreRadius; // Prob 3 - Pore radius 
  double MPS_LPoreSize; // MPS -> Low Pore size
  double MPS_HPoreSize; // MPS -> High Pore size

  int MPS_NumPoreCalc; // MPS -> Number of pore calculations.
  int OutputGRModel2;
  string SolnFile;
  double KirkBeta; // proportionality for probability of
  // Kirkendall Voiding.
  /*
   * Default values for the options
   */
  ProgramOptions();

  ~ProgramOptions();
};

struct TableHeader
{
  double *PoreRadii;
  int NumPoreSizes;
  TableHeader(int num = 1) :
    PoreRadii(0), NumPoreSizes(num)
  {
    PoreRadii = mdp_alloc_dbl_1(num, 0.0);
  }
  ~TableHeader()
  {
    mdp_safe_free((void**) &PoreRadii);
  }
};

struct TableItems
{
  double Time;
  double *HeightBlooms;
  int NumPoreSizes;
  TableItems(int num = 1) :
    Time(0.0), HeightBlooms(0), NumPoreSizes(num)
  {
    HeightBlooms = mdp_alloc_dbl_1(num, 0.0);
  }
  ~TableItems()
  {
    mdp_safe_free((void **) &HeightBlooms);
  }
};
//=====================================================================================================================
namespace CanteraLite
{
//=====================================================================================================================
/*******************************************************************/
/*  Some Specific Reaction Rate Models                             */
/*******************************************************************/
/*
 * Copper Sulfide Engineering Reaction Rate -> quadratic form.
 */
class Cu2SReactionRate : public ReactionRate
{
public:
  Cu2SReactionRate(double x_h2s = 1.41E-7,
                   double a1 = 2.71E3,
                   double e1 = 6.30,
                   double aneg1 = 1.88E3,
                   double eneg1 = 6.3,
                   double pres_atm = 1.0,
                   double temperature = 300.);

  ~Cu2SReactionRate();

  virtual double
  calculateRate(const double * const y);

public:
  double m_X_H2S_g;
  double m_C_H2S_g;
  double m_A1;
  double m_E1; // units = kcal/mole
  double m_Aneg1;
  double m_Eneg1;
};
//=====================================================================================================================
/*
 * Copper Sulfide - Copper Reaction Rate -> Quadratic irreversible
 *                  form
 */
class Cu2SCuRR : public ReactionRate
{
public:
  Cu2SCuRR(double a2 = 9.97, double e2 = 0.0, double temperature = 300.);

  ~Cu2SCuRR();

  virtual double
  calculateRate(const double * const y);

  double m_A2;
  double m_E2;
};

/*
 * Copper Sulfide - Copper Reaction Rate -> Quadratic irreversible
 *                  form With Diffusion through Pore
 */
class Cu2SCuRRPoreDiff : public ReactionRate
{
public:

  Cu2SCuRRPoreDiff(double a2 = 9.97,
                   double e2 = 0.0,
                   double temperature = 300.,
                   double Diff_Coeff = 2.81E-11,
                   double AuPoreThickness = 1.0E-4,
                   double AuPoreRadius = 5.0E-5);

  ~Cu2SCuRRPoreDiff();
  virtual double
  calculateRate(const double * const y);

  double m_A2;
  double m_E2;
  double m_Diff_Coeff;
  double m_AuPoreThickness;
  double m_AuPoreRadius;
  double m_DL;
};
//=====================================================================================================================
/*
 * Simple constant reaction rate
 */
class ConstantReactionRate : public ReactionRate
{
public:

  //! Constructor
  /*!
   *
   * @param sdd   Contains the surface domain description.
   */
  ConstantReactionRate(double k1, double temperature = 300.);

  //ConstantReactionRate(double k1, double temperature = 300.);

  ~ConstantReactionRate();

  virtual double
  calculateRate(const double * const y);

  double m_k1;
};

/*
 * Simple linear reaction rate, with a constant reaction
 * rate constant
 */
class ConstantLinearReactionRate : public ReactionRate
{
public:
  ConstantLinearReactionRate(double k1, double temperature = 300.);

  ~ConstantLinearReactionRate();

  virtual double
  calculateRate(const double * const y);

  double m_k1;
};

/*******************************************************************/
/*******************************************************************/
//=====================================================================================================================
}  // End of namespace CanteraLite
//=====================================================================================================================


//=====================================================================================================================
#endif
//=====================================================================================================================

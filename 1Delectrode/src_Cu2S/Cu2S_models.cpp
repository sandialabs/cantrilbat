/**
 *  @file cu2s_1d_mps.cpp
 */
/*
 *  $Id: Cu2S_models.cpp 506 2013-01-07 22:43:59Z hkmoffa $
 */

#include "cantera/transport.h"



#include "LE_PickList.h"
#include "BE_MoleComp.h"
#include "BE_UnitConversionPressure.h"
#include "LE_OneDblUnits.h"
#include "LE_OneStr.h"
#include "LE_OneBoolInt.h"
#include "LE_OneDbl.h"
#include "LE_OneInt.h"
#include "md_timer.h"
#include "mdp_allo.h"

#include "Cu2S_models.h"

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
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))

using namespace CanteraLite;
using namespace std;
using namespace BEInput;
using namespace TKInput;
using namespace mdpUtil;

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
/*
 *************** Global objects ********************************
 *
 * Ideal gas mixture
 */

ProgramOptions PO;

/*
 * Default values for the options
 */
ProgramOptions::ProgramOptions() :
  Dim(3), NumNodes(10), IncludeGasSpecies(0),
      BCTop(0),
      BCBot(0),
      ProblemType(4), // ProblemType always set to 4 for this calculation
      Debugging(0),
      DeltaTInitial(1.0E-5),
      TFinal(5000.0), // Time in seconds
      OutputDeltaTime(0.0), // Intermediate delta T for output of info
      PrintFlag(1), MaxNumTimeSteps(INT_MAX), RTol(1.0E-3), ATol(1.0E-8),
      numInitialConstantDeltaTSteps(0), printSolnSteps(1),
      printSolnInterval(0), printSolnFirstSteps(0), Temperature(298.15), // Temperature in Kelvin
      Pressure(101325.), // Pressure in pascals

      //! Initial value of the bottom of the computational domain
      Xbot0(1.0),
      //! Initial Value of the top of the computational domain
      //! could be a radius or a height depending on the coordinate system
      Xtop0(1.1),
      //! Initial thickness of the computational domain = Xtop0 - Xbot0
      Thickness0(0.1),

      Cv_init(0.1), k1(1.0), K1_equil(1.0), Conc_Ce(1.0),

      //! diffusion coefficient for vacancy transport (cm2 s-1)
      Diff_Coeff(1.0), Conc_Solid(1.0), CV0_xbot(1.0),

      P1_A1(2.71E5), // Prob 1 preexponential factor for rxn 1
      // units are cm  s-1
      P1_Aneg1(1.01E5), // Prob 1 preexponential factor for rxn -1
      // Units are cm4/mol/s
      P1_A2(9.97), // Prob 1 preexponential factor for rxn 2
      // Units are cm4/mol/s
      P1_XH2S(1.41E-7), P3_AuThick(0.0),
      //! Current value of the pore radius (cm)
      //! (calculated) -> it varies with the current pore size.
      P3_AuPoreRadius(-1.0),
      //! Minimum size of the pore radius (cm) (input)
      MPS_LPoreSize(1.0E-5),
      //! Maximum size of the pore radius (cm) (input)
      MPS_HPoreSize(5.0E-4),
      //! Number of pore calculations in the MPS (multiple pore solutions)
      //! runs. (input)
      MPS_NumPoreCalc(10), OutputGRModel2(0), KirkBeta(0.0) // Default is to not include any cutoff
// due to Kirkendall voiding.
{

}

ProgramOptions::~ProgramOptions()
{
}

bool
fp_compare(double a1, double a2, double tol)
{
  static double atol = 1.0E-10;
  if (fabs(a1) < atol) {
    if (fabs(a2) < tol * atol)
      return true;
    return false;
  }
  if (fabs(a2) < atol) {
    if (fabs(a1) < tol * atol)
      return true;
    return false;
  }
  double a = (fabs(a1) + fabs(a2)) * 0.5;
  if ((fabs(a1 - a2) / a) < tol)
    return true;
  return false;
}

struct TableHeader *TH_ptr = 0;

vector<struct TableItems *> TimeTable;
double TimeTableNext = 0.0;
/*
 * Global file pointer for model 2 pointer.
 */
FILE *ifp2 = 0;

/*
 * ievent = 1 : normal output of a time step. 
 */
void
addGREntry(int ievent,
           double t,
           double thickness,
           double nominalHeight,
           int numProblem)
{
  static double t_old = 0.0;
  static double L_old = -100000.;
  static double Psi = 0.0;
  static double KirkFactor = 1.0;
  static double KirkFactor_old = 1.0;
  static bool wStart = false;
  static double thickness_start = 0.0;

  if (ifp2 == 0) {
    ifp2 = fopen("mps_model2.txt", "w");
    if (!ifp2) {
      printf("Can't open mps_model2.txt file\n");
      exit(-1);
    }
    t_old = t;
    L_old = nominalHeight;
    Psi = 0.0;
    fprintf(ifp2, "Number of radii = %d\n", PO.MPS_NumPoreCalc);
    fprintf(ifp2, "Final Time = %g\n", PO.TFinal);
    fprintf(ifp2, "Intermediate Output Delta Time = %g\n", PO.OutputDeltaTime);
  }
  /*
   * Print out the header for a new Radius calculation
   */
  if (ievent == 0) {
    if (numProblem == 0) {
      wStart = true;
      Psi = 0.0;
      t_old = t;
      thickness_start = thickness;
      L_old = nominalHeight;
      KirkFactor = 1.0;
      KirkFactor_old = 1.0;
      fprintf(ifp2, "Radius = %13.6e\n", PO.P3_AuPoreRadius);

      fprintf(ifp2, "    time   \t    thickness \t    L  "
        "    \t    GR   \t    Psi     \t     ProbNum  \t    K\n");

    }
  }
  if (ievent == 2) {
    if (numProblem == 1) {
      //fprintf(ifp2," END OF INTEGRATION\n");
    }
    printf(" END OF INTEGRATION\n");
  }

  /*
   * Print out normal a row in the table
   */
  if (ievent == 1) {

    double L = nominalHeight;
    bool ldiff = !(fp_compare(L, L_old, 5.0E-4));
    bool tdiff = !(fp_compare(t, t_old, 5.0E-4));
    bool atEnd = false;
    if (fabs((t - PO.TFinal) / (t + PO.TFinal)) < 0.001) {
      atEnd = true;
    }
    if (atEnd || (tdiff && ldiff)) {
      double GR = (L - L_old) / (t - t_old);
      if (wStart) {

        fprintf(ifp2, "%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%8d\t%13.6e\n",
                t_old, thickness_start, L_old, GR, Psi, 0, 1.0);
        wStart = false;
      }
      Psi = Psi + (L - L_old) / GR;
      if (numProblem == 0) {
        double fac = (L - L_old) * ((L + L_old) / 2.0 + PO.P3_AuThick);
        double flog = (PO.KirkBeta * Pi * PO.P3_AuPoreRadius
            * PO.P3_AuPoreRadius * (fac) / (PO.P3_AuPoreRadius * GR));
        if (flog < 0.0) {
          printf("We found an error\n");
          exit(-1);
        }
        KirkFactor = KirkFactor_old * exp(-flog);
        if (KirkFactor > KirkFactor_old) {
          printf("We found an error\n");
          exit(-1);
        }
      } else {
        double hplug = PO.P3_AuThick + PO.P3_AuPoreRadius;
        double volplug = (Pi * PO.P3_AuPoreRadius * PO.P3_AuPoreRadius * hplug);
        double flog = (PO.KirkBeta * volplug * (L - L_old)
            / (PO.P3_AuPoreRadius * GR));
        double rb4n = L * L * L * L;
        double rb4nm1 = L_old * L_old * L_old * L_old;
        double rpore3 = (PO.P3_AuPoreRadius * PO.P3_AuPoreRadius
            * PO.P3_AuPoreRadius);
        double fac = rb4n - rb4nm1 - 4 * rpore3 * (L - L_old);
        flog += (PO.KirkBeta * Pi / 6.0 * fac / (PO.P3_AuPoreRadius * GR));
        if (flog < 0.0) {
          printf("We found an error\n");
        }
        KirkFactor = KirkFactor_old * exp(-flog);
      }
      fprintf(ifp2, "%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%8d\t%13.6e\n", t,
              thickness, L, GR, Psi, numProblem, KirkFactor);
      fflush(ifp2);
      L_old = L;
      KirkFactor_old = KirkFactor;
      t_old = t;
    }
  }

}

void
AddTableEntry(double t, int index, double value)
{
  int num = TH_ptr->NumPoreSizes;
  int iii = (int) TimeTable.size();
  struct TableItems *iti = 0;
  int i, iFound = -1;
  if (iii == 0) {
    iti = new struct TableItems(num);
    TimeTable.push_back(iti);
    iti->Time = 0.0;
    iii = 1;
  }

  bool foundit = false;
  for (i = 0; i < iii; i++) {
    iti = TimeTable[i];
    if (t > iti->Time) {
      continue;
    } else if (t == iti->Time) {
      foundit = true;
      iti->HeightBlooms[index] = value;
      iFound = i;
      TimeTableNext = iti->Time + PO.OutputDeltaTime;
      break;
    } else {
      printf("AddTableEntry: Possible problem t = %g, itiTime = %g\n", t,
             iti->Time);
      for (int j = 0; j < iii; j++) {
        struct TableItems *jti = TimeTable[j];
        printf(" %d time = %g\n", j, jti->Time);
      }
      exit(-1);
    }

  }
  if (!foundit) {
    iti = new struct TableItems(num);
    TimeTable.push_back(iti);
    iti->Time = t;
    iti->HeightBlooms[index] = value;
    iFound = TimeTable.size() - 1;

    TimeTableNext = iti->Time + PO.OutputDeltaTime;
    //printf("AddTableEntry -> adding time = %g\n", t);
  }

  iti = TimeTable[0];
  double v1 = iti->HeightBlooms[index];
  for (i = 1; i <= iFound; i++) {
    iti = TimeTable[i];
    double v2 = iti->HeightBlooms[index];
    if (v2 < v1) {
      printf("We have a problem\n");
    }
  }

  iii = TimeTable.size();
  //for (int j = 0; j < iii; j++) {
  //struct TableItems *jti =  TimeTable[j];
  //printf(" %d time = %g target time = %g\n", j, jti->Time,
  //       TimeTableNext );
  //}
}

//=====================================================================================================================
namespace CanteraLite
{
//=====================================================================================================================
ConstantLinearReactionRate::ConstantLinearReactionRate(double k1,
                                                       double temperature) :
  ReactionRate(temperature), m_k1(k1)
{
  /*
   * Cv is reactant
   * Csolid is product.
   */
  m_NumSpecies = 2;
  m_stoichCoeff.resize(m_NumSpecies, 0.0);
  m_stoichCoeff[0] = -1.0;
  m_stoichCoeff[1] = 1.0;
}
//=====================================================================================================================
ConstantLinearReactionRate::~ConstantLinearReactionRate()
{
}
//=====================================================================================================================
double
ConstantLinearReactionRate::calculateRate(const double * const y)
{
  return (m_k1 * y[0]);
}
//=====================================================================================================================
/*
 * Copper Sulfide Engineering Reaction Rate -> quadratic form.
 */

Cu2SReactionRate::Cu2SReactionRate(double x_h2s,
                                   double a1,
                                   double e1,
                                   double aneg1,
                                   double eneg1,
                                   double pres_atm,
                                   double temperature) :
  ReactionRate(temperature), m_X_H2S_g(x_h2s), m_C_H2S_g(0), m_A1(a1),
      m_E1(e1), m_Aneg1(aneg1), m_Eneg1(eneg1)
{
  /*
   * Cv is reactant
   * Csolid is product.
   * Both are products of the reaction.
   */
  m_NumSpecies = 2;
  m_stoichCoeff.resize(m_NumSpecies, 0.0);
  m_stoichCoeff[0] = 2.0;
  m_stoichCoeff[1] = 1.0;
  /*
   * Calculate the concentration of C_H2S_g
   */
  double ctot = pres_atm / (82.05 * temperature);
  m_C_H2S_g = m_X_H2S_g * ctot;
}
//=====================================================================================================================
Cu2SReactionRate::~Cu2SReactionRate()
{
}
//=====================================================================================================================
double
Cu2SReactionRate::calculateRate(const double * const y)
{
  const double Rkcal = 1.987E-3;
  double fratec = m_A1 * exp(-m_E1 / (Rkcal * m_Temperature));
  double rratec = m_Aneg1 * exp(-m_Eneg1 / (Rkcal * m_Temperature));
  double c = MAX(y[0], 0.0);
  double rate = fratec * m_C_H2S_g - rratec * c * c;
  if (rate < 0.0) {
    printf("Warning rate is less than zero\n");
    exit(-1);
  }
  return rate;
}
//=====================================================================================================================
Cu2SCuRR::Cu2SCuRR(double a2, double e2, double temperature) :
  ReactionRate(temperature), m_A2(a2), m_E2(e2)
{
  /*
   * Cv is reactant
   * Csolid is product.
   * Cv is a quadratic reactant in the reaction.
   * Cs doesn't participate.
   */
  m_NumSpecies = 2;
  m_stoichCoeff.resize(m_NumSpecies, 0.0);
  m_stoichCoeff[0] = -1.0;
  //m_stoichCoeff[0] = -2.0;
  m_stoichCoeff[1] = 0.0;
}
//=====================================================================================================================
Cu2SCuRR::~Cu2SCuRR()
{
}
//=====================================================================================================================
double
Cu2SCuRR::calculateRate(const double * const y)
{
  const double Rkcal = 1.987E-3;
  double fratec = m_A2 * exp(-m_E2 / (Rkcal * m_Temperature));
  double c = MAX(y[0], 0.0);
  double rate = fratec * c * c;
  return rate;
}
//=====================================================================================================================
Cu2SCuRRPoreDiff::Cu2SCuRRPoreDiff(double a2,
                                   double e2,
                                   double temperature,
                                   double Diff_Coeff,
                                   double AuPoreThickness,
                                   double AuPoreRadius) :
  ReactionRate(temperature), m_A2(a2), m_E2(e2), m_Diff_Coeff(Diff_Coeff),
      m_AuPoreThickness(AuPoreThickness), m_AuPoreRadius(AuPoreRadius),
      m_DL(0.0)
{
  /*
   * Cv is reactant
   * Csolid is product.
   * Cv is a quadratic reactant in the reaction.
   * Cs doesn't participate.
   */
  m_NumSpecies = 2;
  m_stoichCoeff.resize(m_NumSpecies, 0.0);
  m_stoichCoeff[0] = -1.0;
  //m_stoichCoeff[0] = -2.0;
  m_stoichCoeff[1] = 0.0;
  if (AuPoreThickness > 0.0) {
    m_DL = m_Diff_Coeff / m_AuPoreThickness;
  }
}
//=====================================================================================================================
Cu2SCuRRPoreDiff::~Cu2SCuRRPoreDiff()
{
}
//=====================================================================================================================
double
Cu2SCuRRPoreDiff::calculateRate(const double * const y)
{
  double rate = 0.0;
  const double Rkcal = 1.987E-3;
  double fratec = m_A2 * exp(-m_E2 / (Rkcal * m_Temperature));
  double c0 = MAX(y[0], 0.0);
  if (c0 <= 0.0)
    return 0.0;
  if (fratec <= 0.0)
    return 0.0;
  if (m_Diff_Coeff <= 0.0)
    return 0.0;
  if (m_AuPoreRadius <= 0.0)
    return 0.0;
  if (m_AuPoreThickness <= 0.0) {
    rate = fratec * c0 * c0;
    return rate;
  }
  double tmp = m_DL * (m_DL + 4.0 * c0 * fratec);
  double c1 = (-m_DL + sqrt(tmp)) / (2.0 * fratec);
  rate = m_DL * (c0 - c1);
  // The 0.5 is due to the difference between the area of a
  // circle and the area of the sphere.
  return (0.5 * rate);
}
//=====================================================================================================================
/*
 * Simple constant reaction rate
 */
ConstantReactionRate::ConstantReactionRate(double k1, double temperature) :
  ReactionRate(temperature), m_k1(k1)
{
  /*
   * Cv is reactant
   * Csolid is product.
   */
  m_NumSpecies = 2;
  m_stoichCoeff.resize(m_NumSpecies);
  m_stoichCoeff[0] = 1.0;
  m_stoichCoeff[1] = 1.0;
}
//=====================================================================================================================
ConstantReactionRate::~ConstantReactionRate()
{
}
//=====================================================================================================================
double
ConstantReactionRate::calculateRate(const double * const y)
{
  return (m_k1);
}
//=====================================================================================================================
}
//=====================================================================================================================

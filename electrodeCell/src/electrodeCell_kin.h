/**
 *  @file electrodeCell_kin.h
 *
 */
/*
 * $Id: electrodeCell_kin.h 504 2013-01-07 22:32:48Z hkmoffa $
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef ELECTRODECELL_KIN_H
#define ELECTRODECELL_KIN_H

#include "cantera/equilibrium.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/numerics/DenseMatrix.h"

#include "TemperatureTable.h"
#include "Electrode.h"
#include "ExtraGlobalRxn.h"

#include <vector>

using namespace Cantera;

//#include "PhaseList.h"
#include "cantera/kinetics/RxnMolChange.h"
 
namespace Cantera {
  class ExtraGlobalRxn;
}


class RxnTempTableStuff {
public:
  RxnTempTableStuff(int irxn, int irxnGE = -1);
  ~RxnTempTableStuff();
  std::vector<double> kfwd_Table;
  std::vector<double> krev_Table;
  std::vector<double> deltaGss_Table;
  std::vector<double> deltaHss_Table;
  std::vector<double> deltaSss_Table;
  std::vector<double> deltaG_Table;
  std::vector<double> deltaH_Table;
  std::vector<double> deltaS_Table;
  std::vector<double> Afwd_Table;
  std::vector<double> EafwddivR_Table;
  std::vector<double> Arev_Table;
  std::vector<double> EarevdivR_Table;
  std::vector<double> kfwdPrime_Table;
  std::vector<double> krevPrime_Table;

  std::vector<double> NetROP_Table;
  std::vector<double> FwdROP_Table;
  std::vector<double> RevROP_Table;
  std::vector<double> Anet_Table;
  std::vector<double> EanetdivR_Table;
  int m_irxn;
  int m_irxnGlobalExtra;
};


class TemperatureTable;

namespace Cantera {
  class PhaseList;
  class Kinetics;
}

double findV(Cantera::Electrode *electrode, doublereal Itarget,
	     double Elow, double Ehigh, int printLvl, doublereal err, int maxsteps);

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void doKineticsTablesHomog(Cantera::Electrode *electrode,
			   Cantera::Kinetics *gKinetics,
			   TemperatureTable &TT);

void doKineticsTablesHetero(Cantera::Electrode *electrode,
			    Cantera::InterfaceKinetics *gKinetics,   
			    TemperatureTable &TT);

void processCurrentVsPotTable(Cantera::RxnMolChange *rmc,
			      Cantera::Electrode *electrode, int irxn,
			      TemperatureTable &TT,
			      Cantera::Kinetics &kin, 
			      Cantera::DenseMatrix& kfwd_Table, 
			      Cantera::DenseMatrix& krev_Table,
			      Cantera::DenseMatrix& deltaG_Table,
			      Cantera::DenseMatrix& deltaH_Table,
			      Cantera::DenseMatrix& deltaS_Table,
			      Cantera::DenseMatrix& Afwd_Table,
			      Cantera::DenseMatrix& EafwddivR_Table,
			      Cantera::DenseMatrix& Arev_Table,
			      Cantera::DenseMatrix& EarevdivR_Table,
			      Cantera::DenseMatrix& kfwdPrime_Table, 
			      Cantera::DenseMatrix& krevPrime_Table);

void
getGERKineticsTables(TemperatureTable& TT, Cantera::Electrode *electrode,
		     Cantera::Kinetics &kin,
		     Cantera::ExtraGlobalRxn &egr,
		     RxnTempTableStuff &rts);

void printGERKineticsTable(Cantera::Electrode *electrode, int j,
			   TemperatureTable& TT,
			   Cantera::Kinetics &kin,
			   Cantera::ExtraGlobalRxn &egr,
			   Cantera::RxnMolChange *rmc,
			   RxnTempTableStuff &rts);

double processGERCurrent(Cantera::RxnMolChange *rmc,
			 Cantera::Electrode *electrode, int iGERrxn,
			 Kinetics &kin,
			 Cantera::ExtraGlobalRxn &egr);

void processGERCurrentVsPotTable(Cantera::RxnMolChange *rmc, Cantera::Electrode *electrode,
				 int irxn,
				 TemperatureTable &TT,
				 Cantera::Kinetics &kin,  
				 Cantera::ExtraGlobalRxn &egr,
				 RxnTempTableStuff &rts);
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#endif

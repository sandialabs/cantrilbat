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

#include "zuzax/equilibrium.h"
#include "zuzax/kinetics.h"
#include "zuzax/kinetics/InterfaceKinetics.h"
#include "zuzax/thermo/SurfPhase.h"
#include "zuzax/numerics/DenseMatrix.h"
#include "zuzax/kinetics/ExtraGlobalRxn.h"

#include "TemperatureTable.h"
#include "Electrode.h"

#include <vector>

//#include "PhaseList.h"
#include "zuzax/kinetics/RxnMolChange.h"
 


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

namespace Zuzax {
  class PhaseList;
  class Kinetics;
}

double findV(Zuzax::Electrode *electrode, double Itarget, double Elow, double Ehigh, int printLvl, double err, int maxsteps);

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void doKineticsTablesHomog(Zuzax::Electrode *electrode,
			   Zuzax::Kinetics *gKinetics,
			   TemperatureTable &TT);

void doKineticsTablesHetero(Zuzax::Electrode *electrode,
			    Zuzax::InterfaceKinetics *gKinetics,   
			    TemperatureTable &TT);

void processCurrentVsPotTable(Zuzax::RxnMolChange *rmc,
			      Zuzax::Electrode *electrode, int irxn,
			      TemperatureTable &TT,
			      Zuzax::Kinetics &kin, 
			      Zuzax::DenseMatrix& kfwd_Table, 
			      Zuzax::DenseMatrix& krev_Table,
			      Zuzax::DenseMatrix& deltaG_Table,
			      Zuzax::DenseMatrix& deltaH_Table,
			      Zuzax::DenseMatrix& deltaS_Table,
			      Zuzax::DenseMatrix& Afwd_Table,
			      Zuzax::DenseMatrix& EafwddivR_Table,
			      Zuzax::DenseMatrix& Arev_Table,
			      Zuzax::DenseMatrix& EarevdivR_Table,
			      Zuzax::DenseMatrix& kfwdPrime_Table, 
			      Zuzax::DenseMatrix& krevPrime_Table);

void
getGERKineticsTables(TemperatureTable& TT, Zuzax::Electrode *electrode,
		     Zuzax::Kinetics &kin,
		     Zuzax::ExtraGlobalRxn &egr,
		     RxnTempTableStuff &rts);

void printGERKineticsTable(Zuzax::Electrode *electrode, int j,
			   TemperatureTable& TT,
			   Zuzax::Kinetics &kin,
			   Zuzax::ExtraGlobalRxn &egr,
			   Zuzax::RxnMolChange *rmc,
			   RxnTempTableStuff &rts);

double processGERCurrent(Zuzax::RxnMolChange *rmc,
			 Zuzax::Electrode *electrode, int iGERrxn,
			 Zuzax::Kinetics &kin,
			 Zuzax::ExtraGlobalRxn &egr);

void processGERCurrentVsPotTable(Zuzax::RxnMolChange *rmc, Zuzax::Electrode *electrode,
				 int irxn,
				 TemperatureTable &TT,
				 Zuzax::Kinetics &kin,  
				 Zuzax::ExtraGlobalRxn &egr,
				 RxnTempTableStuff &rts);
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#endif

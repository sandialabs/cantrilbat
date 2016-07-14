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
#include "cantera/kinetics/ExtraGlobalRxn.h"

#include "TemperatureTable.h"
#include "Electrode.h"

#include <vector>

//#include "PhaseList.h"
#include "cantera/kinetics/RxnMolChange.h"
 


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

double findV(ZZCantera::Electrode *electrode, doublereal Itarget,
	     double Elow, double Ehigh, int printLvl, doublereal err, int maxsteps);

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void doKineticsTablesHomog(ZZCantera::Electrode *electrode,
			   ZZCantera::Kinetics *gKinetics,
			   TemperatureTable &TT);

void doKineticsTablesHetero(ZZCantera::Electrode *electrode,
			    ZZCantera::InterfaceKinetics *gKinetics,   
			    TemperatureTable &TT);

void processCurrentVsPotTable(ZZCantera::RxnMolChange *rmc,
			      ZZCantera::Electrode *electrode, int irxn,
			      TemperatureTable &TT,
			      ZZCantera::Kinetics &kin, 
			      ZZCantera::DenseMatrix& kfwd_Table, 
			      ZZCantera::DenseMatrix& krev_Table,
			      ZZCantera::DenseMatrix& deltaG_Table,
			      ZZCantera::DenseMatrix& deltaH_Table,
			      ZZCantera::DenseMatrix& deltaS_Table,
			      ZZCantera::DenseMatrix& Afwd_Table,
			      ZZCantera::DenseMatrix& EafwddivR_Table,
			      ZZCantera::DenseMatrix& Arev_Table,
			      ZZCantera::DenseMatrix& EarevdivR_Table,
			      ZZCantera::DenseMatrix& kfwdPrime_Table, 
			      ZZCantera::DenseMatrix& krevPrime_Table);

void
getGERKineticsTables(TemperatureTable& TT, ZZCantera::Electrode *electrode,
		     ZZCantera::Kinetics &kin,
		     ZZCantera::ExtraGlobalRxn &egr,
		     RxnTempTableStuff &rts);

void printGERKineticsTable(ZZCantera::Electrode *electrode, int j,
			   TemperatureTable& TT,
			   ZZCantera::Kinetics &kin,
			   ZZCantera::ExtraGlobalRxn &egr,
			   ZZCantera::RxnMolChange *rmc,
			   RxnTempTableStuff &rts);

double processGERCurrent(ZZCantera::RxnMolChange *rmc,
			 ZZCantera::Electrode *electrode, int iGERrxn,
			 ZZCantera::Kinetics &kin,
			 ZZCantera::ExtraGlobalRxn &egr);

void processGERCurrentVsPotTable(ZZCantera::RxnMolChange *rmc, ZZCantera::Electrode *electrode,
				 int irxn,
				 TemperatureTable &TT,
				 ZZCantera::Kinetics &kin,  
				 ZZCantera::ExtraGlobalRxn &egr,
				 RxnTempTableStuff &rts);
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#endif

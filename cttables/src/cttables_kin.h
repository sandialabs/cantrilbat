/**
 *  @file cttables_kin.h
 *
 */
/*
 * $Id: cttables_kin.h 497 2013-01-07 21:17:04Z hkmoffa $
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CTTABLES_KIN_H
#define CTTABLES_KIN_H

#include "cttables.h"
#include "TemperatureTable.h"
#include "PhaseList.h"
namespace Cantera {
  class ExtraGlobalRxn;
}
class RxnMolChange {
public:
  RxnMolChange(Cantera::Kinetics *kinPtr, int irxn);
  virtual ~RxnMolChange();
  RxnMolChange(Cantera::Kinetics *kinPtr, Cantera::ExtraGlobalRxn *egr);

  std::vector<double> m_phaseMoleChange;
  std::vector<double> m_phaseMassChange;
  std::vector<double> m_phaseChargeChange;
  std::vector<double> m_phasePotentials;
  std::vector<int>    m_phaseTypes;
  std::vector<int>    m_phaseDims;

  int m_nPhases;
  Cantera::Kinetics *m_kinBase;
  int m_iRxn;
  double m_ChargeTransferInRxn;
  double m_beta;
  Cantera::ExtraGlobalRxn *m_egr;
};

class RxnTempTableStuff {
public:
  RxnTempTableStuff(int irxn, int irxnGE = -1);
  ~RxnTempTableStuff();
  vector<double> kfwd_Table;
  vector<double> krev_Table;
  vector<double> deltaGss_Table;
  vector<double> deltaHss_Table;
  vector<double> deltaSss_Table;
  vector<double> deltaG_Table;
  vector<double> deltaH_Table;
  vector<double> deltaS_Table;
  vector<double> Afwd_Table;
  vector<double> EafwddivR_Table;
  vector<double> Arev_Table;
  vector<double> EarevdivR_Table;
  vector<double> kfwdPrime_Table;
  vector<double> krevPrime_Table;

  vector<double> NetROP_Table;
  vector<double> FwdROP_Table;
  vector<double> RevROP_Table;
  vector<double> Anet_Table;
  vector<double> EanetdivR_Table;
  int m_irxn;
  int m_irxnGlobalExtra;
};


class TemperatureTable;

//class DenseMatrix;
#include "cantera/kinetics.h"
#include "cantera/numerics/DenseMatrix.h"
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void doKineticsTablesHomog(Cantera::PhaseList *pl,
			   Cantera::Kinetics *gKinetics,
			   TemperatureTable &TT);

void doKineticsTablesHetero(Cantera::PhaseList *pl, 
			    Cantera::InterfaceKinetics *gKinetics,   
			    TemperatureTable &TT);

void processCurrentVsPotTable(RxnMolChange *rmc,
			      Cantera::PhaseList *pl, int irxn,
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
getGERKineticsTables(TemperatureTable& TT, Cantera::PhaseList *pl,
		     Cantera::Kinetics &kin,
		     Cantera::ExtraGlobalRxn &egr,
		     RxnTempTableStuff &rts);

void printGERKineticsTable(Cantera::PhaseList *pl, int j,
			   TemperatureTable& TT,
			   Cantera::Kinetics &kin,
			   Cantera::ExtraGlobalRxn &egr,
			   RxnMolChange *rmc,
			   RxnTempTableStuff &rts);

void processGERCurrentVsPotTable(RxnMolChange *rmc,
				 Cantera::PhaseList *pl, int irxn,
				 TemperatureTable &TT,
				 Cantera::Kinetics &kin,  
				 Cantera::ExtraGlobalRxn &egr,
				 RxnTempTableStuff &rts);
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#endif

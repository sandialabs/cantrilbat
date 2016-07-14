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
#include "cantera/kinetics/RxnMolChange.h"
#include "cantera/kinetics/ExtraGlobalRxn.h"

//--------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{
  class ExtraGlobalRxn;
}
//--------------------------------------------------------------------------------------------------------------------------------

class RxnMolChangeLocal 
{
public:
  RxnMolChangeLocal(ZZCantera::Kinetics *kinPtr, int irxn);
  virtual ~RxnMolChangeLocal();
  RxnMolChangeLocal(ZZCantera::Kinetics *kinPtr, ZZCantera::ExtraGlobalRxn *egr);

  std::vector<double> m_phaseMoleChange;
  std::vector<double> m_phaseMassChange;
  std::vector<double> m_phaseChargeChange;
  std::vector<double> m_phasePotentials;
  std::vector<int>    m_phaseTypes;
  std::vector<int>    m_phaseDims;

  int m_nPhases;
  ZZCantera::Kinetics *m_kinBase;
  int m_iRxn;
  double m_ChargeTransferInRxn;
  double m_beta;
  ZZCantera::ExtraGlobalRxn *m_egr;
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

void doKineticsTablesHomog(ZZCantera::PhaseList *pl,
			   ZZCantera::Kinetics *gKinetics,
			   TemperatureTable &TT);

void doKineticsTablesHetero(ZZCantera::PhaseList *pl, 
			    ZZCantera::InterfaceKinetics *gKinetics,   
			    TemperatureTable &TT);

void processCurrentVsPotTable(ZZCantera::RxnMolChange *rmc,
			      ZZCantera::PhaseList *pl, int irxn,
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
getGERKineticsTables(TemperatureTable& TT, ZZCantera::PhaseList *pl,
		     ZZCantera::Kinetics &kin,
		     ZZCantera::ExtraGlobalRxn &egr,
		     RxnTempTableStuff &rts);

void printGERKineticsTable(ZZCantera::PhaseList *pl, int j,
			   TemperatureTable& TT,
			   ZZCantera::Kinetics &kin,
			   ZZCantera::ExtraGlobalRxn &egr,
			   ZZCantera::RxnMolChange *rmc,
			   RxnTempTableStuff &rts);

void processGERCurrentVsPotTable(ZZCantera::RxnMolChange *rmc,
				 ZZCantera::PhaseList *pl, int irxn,
				 TemperatureTable &TT,
				 ZZCantera::Kinetics &kin,  
				 ZZCantera::ExtraGlobalRxn &egr,
				 RxnTempTableStuff &rts);

void printAffinityHeader(ZZCantera::RxnMolChange *rmc,
                              ZZCantera::PhaseList *pl, int irxn,
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
                              ZZCantera::DenseMatrix& krevPrime_Table, double *unitsfwd, double* unitsRev);

void processAffinityTable(ZZCantera::RxnMolChange *rmc,
                              ZZCantera::PhaseList *pl, int irxn,
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


/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#endif

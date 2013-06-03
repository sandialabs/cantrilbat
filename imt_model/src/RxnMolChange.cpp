/**
 *  @file example2.cpp
 *
 * $Id: RxnMolChange.cpp 507 2013-01-07 22:48:29Z hkmoffa $
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "RxnMolChange.h"

// Kinetics includes
#include "cantera/kinetics.h"

#include "ExtraGlobalRxn.h"

#include <iostream>
#include <new>

using namespace Cantera;
using namespace std;

RxnMolChange::RxnMolChange(Kinetics *kinPtr, int irxn) :
  m_nPhases(0),
  m_kinBase(kinPtr),
  m_iRxn(irxn),
  m_ChargeTransferInRxn(0.0),
  m_beta(0.0),
  m_egr(0)
{
  int nReac = kinPtr->nReactions();
  int iph;
  AssertTrace(irxn >= 0);
  AssertTrace(irxn < nReac);

  m_nPhases = kinPtr->nPhases();

  m_phaseMoleChange.resize(m_nPhases, 0.0);
  m_phaseMassChange.resize(m_nPhases, 0.0);
  m_phaseChargeChange.resize(m_nPhases, 0.0);
  m_phasePotentials.resize(m_nPhases, 0.0);
  m_phaseTypes.resize(m_nPhases, 0);
  m_phaseDims.resize(m_nPhases, 0);

  int m_kk = kinPtr->nTotalSpecies();

  for (int kKin = 0; kKin < m_kk; kKin++) {
   iph =  m_kinBase->speciesPhaseIndex(kKin);
   ThermoPhase &tpRef = m_kinBase->thermo(iph);
   int kLoc = kKin - m_kinBase->kineticsSpeciesIndex(0, iph);
   double rsc = m_kinBase->reactantStoichCoeff(kKin, irxn);
   double psc = m_kinBase->productStoichCoeff(kKin, irxn);
   double nsc = psc - rsc;
   m_phaseMoleChange[iph] += (nsc);
   double mw = tpRef.molecularWeight(kLoc);
   m_phaseMassChange[iph] += (nsc) * mw;
   double chg = tpRef.charge(kLoc);
   m_phaseChargeChange[iph] += nsc * chg;
  }

  for (iph = 0; iph < m_nPhases; iph++) {
    ThermoPhase &tpRef = m_kinBase->thermo(iph);
    m_phasePotentials[iph] =  tpRef.electricPotential();
    m_phaseDims[iph] = tpRef.nDim();
    m_phaseTypes[iph] = tpRef.eosType();
    if (m_phaseChargeChange[iph] != 0.0) {
      double tmp = fabs(m_phaseChargeChange[iph]);
      if (tmp >  m_ChargeTransferInRxn) {
	m_ChargeTransferInRxn = tmp;
      }
    }
  }

  if (m_ChargeTransferInRxn) {
    InterfaceKinetics *iK = dynamic_cast<InterfaceKinetics *>(kinPtr);
    if (iK) {
      m_beta = iK->electrochem_beta(irxn);
    } else {
      throw CanteraError("RxnMolChange", "unknown condition on charge");
    }
  }

}

RxnMolChange::RxnMolChange(Kinetics *kinPtr, ExtraGlobalRxn *egr) :
  m_nPhases(0),
  m_kinBase(kinPtr),
  m_iRxn(-1),
  m_ChargeTransferInRxn(0.0),
  m_beta(0.0),
  m_egr(egr)
{
  int iph;
  AssertTrace(egr != 0);

  m_nPhases = kinPtr->nPhases();

  m_phaseMoleChange.resize(m_nPhases, 0.0);
  m_phaseMassChange.resize(m_nPhases, 0.0);
  m_phaseChargeChange.resize(m_nPhases, 0.0);
  m_phasePotentials.resize(m_nPhases, 0.0);
  m_phaseTypes.resize(m_nPhases, 0);
  m_phaseDims.resize(m_nPhases, 0);

  int m_kk = kinPtr->nTotalSpecies();

  for (int kKin = 0; kKin < m_kk; kKin++) {
   iph =  m_kinBase->speciesPhaseIndex(kKin);
   ThermoPhase &tpRef = m_kinBase->thermo(iph);
   int kLoc = kKin - m_kinBase->kineticsSpeciesIndex(0, iph);
   double rsc = egr->reactantStoichCoeff(kKin);
   double psc = egr->productStoichCoeff(kKin);
   double nsc = psc - rsc;
   m_phaseMoleChange[iph] += (nsc);
   double mw = tpRef.molecularWeight(kLoc);
   m_phaseMassChange[iph] += (nsc) * mw;
   double chg = tpRef.charge(kLoc);
   m_phaseChargeChange[iph] += nsc * chg;
  }

  for (iph = 0; iph < m_nPhases; iph++) {
    ThermoPhase &tpRef = m_kinBase->thermo(iph);
    m_phasePotentials[iph] =  tpRef.electricPotential();
    m_phaseDims[iph] = tpRef.nDim();
    m_phaseTypes[iph] = tpRef.eosType();
    if ( m_phaseChargeChange[iph] != 0.0) {
      double tmp = fabs(m_phaseChargeChange[iph]);
      if (tmp >  m_ChargeTransferInRxn) {
	m_ChargeTransferInRxn = tmp;
      }
    }
  }

  if (m_ChargeTransferInRxn) {
    InterfaceKinetics *iK = dynamic_cast<InterfaceKinetics *>(kinPtr);
    if (iK) {
      m_beta = 0.0;
    } else {
      throw CanteraError("RxnMolChange", "unknown condition on charge");
    }
  }

}


RxnMolChange::~RxnMolChange() {

}

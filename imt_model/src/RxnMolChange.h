/**
 *  @file cttables_kin.h
 *
 */
/*
 * $Id: RxnMolChange.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef RXNMOLCHANGE_H
#define RXNMOLCHANGE_H


#include <vector>

namespace Cantera {
  class ExtraGlobalRxn;
  class Kinetics;
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

#endif

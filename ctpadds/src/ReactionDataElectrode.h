/**
 *  @file ReactionDataElectrode.h
 *
 */
/*
 * $Author: hkmoffa $
 * $Revision: 508 $
 * $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_REACTION_DATA_ELECTRODE_H
#define CT_REACTION_DATA_ELECTRODE_H

#include "cantera/kinetics/ReactionData.h"
#include "cantera/base/xml.h"


namespace Cantera {

  class XML_Node;
  class Kinetics;
  class ThermoPhase;

  class ReactionDataElectrode : public ReactionData {
  public:
    ReactionDataElectrode();
    virtual ~ReactionDataElectrode(); 

    void getCoverageDependence(const XML_Node& node,
			       ThermoPhase& surfphase);
    void getStick(const XML_Node& node, Kinetics& kin,
		  doublereal& A, doublereal& b, doublereal& E);

    void getRateCoefficient(const XML_Node &kf, Kinetics & kin, 
			    int negA);
  
  };
}

#endif

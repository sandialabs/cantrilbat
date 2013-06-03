/**
 * @file ReactingVolDomain.h
 *
 */
/*
 * $Author: hkmoffa $
 * $Revision: 497 $
 * $Date: 2013-01-07 14:17:04 -0700 (Mon, 07 Jan 2013) $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef REACTINGVOLDOMAIN_H
#define REACTINGVOLDOMAIN_H
 
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/transport.h"
 
#include <string>
#include <iostream>

 
namespace Cantera {
 
  class Transport;
  class PhaseList;

  class ReactingVolDomain
  {
  public:
    ReactingVolDomain();
    /*ReactingVolDomain(XML_Node& root, std::string id);*/
    virtual ~ReactingVolDomain();

    /**
     * Import the ReactingVolDomain from a single 
     * XML source.
     */
    bool importFromXML(XML_Node& root);

    bool importFromPL(Cantera::PhaseList *pl, int ivkin, int iskin);

    bool ready() { return m_ok; }
    friend std::ostream& operator<<(std::ostream& s,
				    ReactingVolDomain &vd);

    ThermoPhase& thermo(int n=0) {
      return (m_InterfaceKinetics->thermo(n));
    }

    
    //! This kinetics operator is associated with just one
    //! homogeneous phase, associated with tpList[0] phase
    /*!
     * This object owns the Kinetics object
     */
    Kinetics *m_kinetics;
   
    //! This kinetics operator is associated with multiple
    //! homogeneous and surface phases.
    /*!
     * This object owns the Kinetics object
     */
    InterfaceKinetics *m_InterfaceKinetics;
    Transport* m_transport;
    int numPhases;
    std::vector<ThermoPhase *>tpList;
    std::vector<XML_Node *> xmlList;
      
    std::vector<int> kinOrder;

    int iphaseKin;

    std::vector<int> tplRead;

    std::string transportModel;

    bool m_DoSurfKinetics;
    bool m_DoHomogKinetics;

  protected:
    bool m_ok;
    XML_Node* m_XMLPhaseTree;
  private:
    bool m_XMLTree_owned;
    bool m_I_Own_Thermo;
  };
}
#endif

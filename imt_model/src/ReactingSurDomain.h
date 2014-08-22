/**
 * @file ReactingVolDomain.h
 *
 */
/*
 * $Author: hkmoffa $
 * $Revision: 507 $
 * $Date: 2013-01-07 15:48:29 -0700 (Mon, 07 Jan 2013) $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef REACTINGSURDOMAIN_H
#define REACTINGSURDOMAIN_H
 
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/RxnMolChange.h"
#include "cantera/transport.h"
 
#include <string>
#include <iostream>
 
namespace Cantera {
 
  class Transport;
  class PhaseList;

  //! ReactingSurDomain is a class that is a wrapper around InterfaceKinetics that
  //! knows about the PhaseList object
  /*!
   *    There is one surface associated with the object and one interfacial kinetics object.
   *
   */
  class ReactingSurDomain : public InterfaceKinetics
  {
  public:
    //! Default constructor
    ReactingSurDomain();

    //! Copy Constructor for the %Kinetics object.
    /*!
     * Currently, this is not fully implemented. If called it will
     * throw an exception.
     */
    ReactingSurDomain(const ReactingSurDomain &right);

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %Kinetics object to be copied into the
     *                 current one.
     */
    ReactingSurDomain & operator=(const ReactingSurDomain &right);

    //! Default destructor
    virtual ~ReactingSurDomain();

    //! Duplication routine for objects which inherit from Kinetics
    /*!
     *  This virtual routine can be used to duplicate %Kinetics objects
     *  inherited from %Kinetics even if the application only has
     *  a pointer to %Kinetics to work with.
     *
     *  These routines are basically wrappers around the derived copy  constructor.
     *
     * @param  tpVector Vector of shallow pointers to ThermoPhase objects. this is the
     *                  m_thermo vector within this object
     */
    virtual Kinetics *duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

    //!   Import all the phases from a PhaseList and initialize the
    //!   object
    bool importFromPL(Cantera::PhaseList *pl, int ivkin, int iskin);

    //! Returns a reference to the calculated production rates of species per surface area
    /*!
     *   This routine calls thet getNetProductionRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species net 
     *         production rates per surface area (kmol m-2 s-1)
     */
    const std::vector<double> & calcNetProductionRates();

    //! Returns a reference to the calculated creation rates of species
    /*!
     *   This routine calls thet getCreationRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species creation rates
     */
    const std::vector<double> & calcCreationRates();

    //! Returns a reference to the calculated destruction rates of species
    /*!
     *   This routine calls thet getDestructionRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species destruction rates
     */
    const std::vector<double> & calcDestructionRates();

    //! Declare a printing routine as a friend to this class
    friend std::ostream& operator<<(std::ostream& s, ReactingSurDomain &vd);
 
    /*
     *
     */
 

    //! Pointer to the transport operator
    Transport* m_transport;

    //! Vector of additional information about each reaction
    std::vector<Cantera::RxnMolChange *> rmcVector;
   
    //! Number of phases within the class
    int numPhases;

  public:
    //!  Vector of pointers to xml trees
    std::vector<XML_Node *> xmlList;

    //! Mapping between the phase order in the interface object
    //! and the overall order in the PhaseList object 
    /*!
     *  Note in the phase list object, surface phases are listed last.
     *  Length is the number of phases in the interface kinetics object
     *  value is the id of the phase in the PhaseList object.
     *
     *   kinOrder[kph] = iph;
     *        kph = phase index in the interface kinetics object
     *        iph = phase index in the PhaseList object
     */
    std::vector<int> kinOrder;

    //! Vector of the indexes of each phase in the ReactionSurfaceDomain object
    //! given the index withint the PhaseList object
    /*!
     *       jph = PLtoKinPhaseIndex_[iph];
     *
     *          iph refers to the index of the phase in the Electrode_Model object
     *          jph refers to the index of the phase in the heterogeneous kinetics object
     *
     *  Length = number of phases in the PhaseList
     *
     *  A value of -1 in this slot means that the phase doesn't participate in the 
     *  current ReactingSurDomain object
     */
    std::vector<int> PLtoKinPhaseIndex_;

    std::vector<int> PLtoKinSpeciesIndex_;

    //! ID of the phase in the PhaseList object that has the kinetics
    //! object
    int iphaseKin;

    std::vector<int> tplRead;

    bool m_DoSurfKinetics;
 
    //! Vector that will expose the species production rates for this kinetics object
    std::vector<double> speciesProductionRates_;

    //! Vector that will expose the species creation rates for this kinetics object
    std::vector<double> speciesCreationRates_;

    //! Vector that will expose the species destruction rates for this kinetics object
    std::vector<double> speciesDestructionRates_;


    //! Pointer to the phaselist object that contains the ThermoPhase objects
    Cantera::PhaseList *m_pl;

  protected:
    bool m_ok;
    XML_Node* m_XMLPhaseTree;
 
  };
}
#endif

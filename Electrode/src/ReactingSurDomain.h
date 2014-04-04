/**
 * @file ReactingVolDomain.h
 *
 */
/*
 * $Author: hkmoffa $
 * $Revision: 571 $
 * $Date: 2013-03-26 10:44:21 -0600 (Tue, 26 Mar 2013) $
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
#include "cantera/transport.h"

#include <string>
#include <iostream>
class RxnMolChange;

namespace Cantera
{

class Transport;
class PhaseList;


class ReactingSurDomain : public Cantera::InterfaceKinetics
{
public:
    //! Default constructor
    ReactingSurDomain();

    //! Copy Constructor for the %Kinetics object.
    /*!
     * Currently, this is not fully implemented. If called it will
     * throw an exception.
     */
    ReactingSurDomain(const ReactingSurDomain& right);

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %Kinetics object to be copied into the
     *                 current one.
     */
    ReactingSurDomain& operator=(const ReactingSurDomain& right);

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
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*>& tpVector) const;

    //!   Import all the phases from a PhaseList and initialize the
    //!   object
    bool importFromPL(Cantera::PhaseList* pl, int ivkin, int iskin);

    //! Returns a reference to the calculated production rates of species
    /*!
     *   This routine calls thet getNetProductionRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species net
     *         production rates
     */
    const std::vector<double>& calcNetProductionRates();

    //! Returns a reference to the calculated creation rates of species
    /*!
     *   This routine calls thet getCreationRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species creation rates
     */
    const std::vector<double>& calcCreationRates();

    //! Returns a reference to the calculated destruction rates of species
    /*!
     *   This routine calls thet getDestructionRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species destruction rates
     */
    const std::vector<double>& calcDestructionRates();


    //! Get the exchange current density formulation for the current reaction rate
    /*!
     *  Returns the parameterization for a single reaction that defines the exchange
     *  current density formulation for that reaction. Note the exchange current density depends
     *  on the activity concentrations of the reactants and the products.
     *
     *         i  = io [ exp( (Betaf  Nu  nSt  F)/ RT) - exp(- ((1 - Betaf) Nu nSt F)/ RT) ]
     *
     *              where Nu = E - OCV is the overpotential
     *
     *              and E = Phi(Metal) - Phi(Soln)
     *
     *   where the current is defined as
     *
     *         i = - z_e- * F * d[e-]/dt
     *
     *    and
     *
     *         d[e-]/dt = (n) ROP_net
     *
     *   A note about signs
     *
     *    i is defined here as being positive if there is a net current going from the metal into solution.
     *    This in turn means that if there is a positive generation of electrons, then the current
     *    is positive.
     *
     *    In this formulation the stoichiometric number of electrons can be positive or negative
     *    It will be negative if the electrons appears on the reactant side of the equation.
     *
     *    In this formulation, io can be positive or negative. It's sign will be determined by the
     *    sign of the nSt, the stoichiometric electrons.
     *
     *   @param   irxn     Reaction id
     *   @param   nSt      Number of stoichiometric electrons transferred (can be negative)
     *   @param   OCV      Open circuit voltage (volts)
     *   @param   io       Exchange Current density value (can be negative)
     *                        units coulombs / sec / m^2
     *   @param   nu       Overpotential for the reaction (can be positive or negative)
     *   @param   beta     Symmetry factor
     *
     *  @return  returns the current for the reaction
     */
    double getExchangeCurrentFormulation(int irxn,  doublereal* nStoich, doublereal* OCV,
                                       doublereal* io, doublereal* nu, doublereal *beta);

    double calcCurrent(double, double, double, double, double) const;

    //!  Identify the metal phase and the electrons species
    /*!
     *   We fill in the internal variables, metalPhaseRS_ and kElectronRS_ here
     */
    void identifyMetalPhase();

public:
    //! Declare a printing routine as a friend to this class
    friend std::ostream& operator<<(std::ostream& s, ReactingSurDomain& vd);

    /*
     *
     */


    //! Pointer to the transport operator
    Transport* m_transport;

    //! Vector of additional information about each reaction
    std::vector<RxnMolChange*> rmcVector;

    //! Number of phases within the class
    int numPhases;


    //!  Vector of pointers to xml trees
    std::vector<XML_Node*> xmlList;

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

    //! Temp vector that may be eliminated in the future.
    std::vector<int> tplRead;

    bool m_DoSurfKinetics;

    //! Vector that will expose the species production rates for this kinetics object
    std::vector<double> speciesProductionRates_;

    //! Vector that will expose the species creation rates for this kinetics object
    std::vector<double> speciesCreationRates_;

    //! Vector that will expose the species destruction rates for this kinetics object
    std::vector<double> speciesDestructionRates_;

    //! Pointer to the phaselist object that contains the ThermoPhase objects.
    /*!
     *  This object doesn't own this
     */
    Cantera::PhaseList* m_pl;

    //! index of the metal phase in the list of phases for this surface
    int metalPhaseRS_;

    //! Index of the electrons species in the list of species for this surface, if none set it to -1
    int kElectronRS_;

    //! Index of the solution phase in the list of phases for this surface
    int solnPhaseRS_;

protected:


};
}
#endif

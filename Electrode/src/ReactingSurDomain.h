/**
 * @file ReactingVolDomain.h
 *
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef REACTINGSURDOMAIN_H
#define REACTINGSURDOMAIN_H

#include "cantera/kinetics/InterfaceKinetics.h"

#include "RSD_OCVmodel.h"

#include <string>
#include <vector>

class RxnMolChange;

namespace Cantera
{
class PhaseList;
class OCV_Override_input;

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

    //!   Import all the phases from a PhaseList and initialize the kinetics for this object
    /*!
     *
     */
    bool importFromPL(Cantera::PhaseList* pl, int iskin);

    //! Returns a reference to the calculated production rates of species
    /*!
     *   This routine calls thet getNetProductionRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species net
     *         production rates (kmol s-1 m-2)
     */
    const std::vector<double>& calcNetSurfaceProductionRateDensities();

    const std::vector<double>& calcNetSurfaceROP();

    //! Returns a reference to the calculated creation rates of species
    /*!
     *   This routine calls thet getCreationRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species creation rates
     */
    const std::vector<double>& calcSurfaceCreationRateDensities();

    //! Returns a reference to the calculated destruction rates of species
    /*!
     *   This routine calls thet getDestructionRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species destruction rates
     */
    const std::vector<double>& calcSurfaceDestructionRateDensities();


    //! Get the net current for the set of reactions on this surface
    /*!
     *       
     */
    double getCurrentDensityRxn(double *currentDensityRxn = 0);

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
     *   where the current density (amps m-2) is defined as
     *
     *         i = - z_e- * F * d[e-]/dt
     *
     *    and
     *
     *         d[e-]/dt = (n) ROP_net
     *
     *   A note about signs
     *
     *    i is defined here as being positive if there is a net current density going from the metal into solution.
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
     *  @return  returns the current density for the reaction (amps m-2)
     */
    double getExchangeCurrentDensityFormulation(int irxn,  doublereal* nStoich, doublereal* OCV,
                                                doublereal* io, doublereal* nu, doublereal *beta);

    //! Utility routine to calculate the current density given the parameters for
    //! an exchange current density formulation of the reaction rate
    /*!
     *
     *  @param returns the current density (amps m-2)
     */
    double calcCurrentDensity(double nu, double nStoich, double io, double beta, 
                              double temp) const;

    //!  Identify the metal phase and the electrons species
    /*!
     *   We fill in the internal variables, metalPhaseRS_ and kElectronRS_ here
     */
    void identifyMetalPhase();

    //! Add an open circuit voltage override feature to the current reacting surface
    /*!
     *     The open circuit voltage will replace the deltaG calculation for this reaction.
     *     It will do this by replacing the thermo functions for one of the species in the mechanism.
     *
     *      @param ocv_ptr    The open circuit voltage feature is described by this pointer.
     *
     */
    void addOCVoverride(OCV_Override_input *ocv_ptr);

    void deriveEffectiveChemPot();

public:
    //! Declare a printing routine as a friend to this class
    friend std::ostream& operator<<(std::ostream& s, ReactingSurDomain& vd);

    //! Vector of additional information about each reaction
    std::vector<RxnMolChange*> rmcVector;

    //! Number of phases within the class
    int numPhases_;

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
    int iphaseKin_;

    //! Temp vector that may be eliminated in the future.
    std::vector<int> tplRead;

    bool m_DoSurfKinetics;

    //! Vector that will expose the species production rates for this kinetics object
    /*!
     *  Length is the number of species in the kinetics vector
     */
    std::vector<double> speciesProductionRates_;

    //! Vector that will expose the species creation rates for this kinetics object
    /*!
     *  Length is the number of species in the kinetics vector
     */
    std::vector<double> speciesCreationRates_;

    //! Vector that will expose the species destruction rates for this kinetics object
    /*!
     *  Length is the number of species in the kinetics vector
     */
    std::vector<double> speciesDestructionRates_;

    std::vector<double> deltaGRxn_;

    //! Pointer to the phaselist object that contains the ThermoPhase objects.
    /*!
     *  This object doesn't own this. However, it uses this heavily.
     */
    Cantera::PhaseList* m_pl;

    //! index of the metal phase in the list of phases for this surface
    int metalPhaseRS_;

    //! Index of the electrons species in the list of species for this surface, if none set it to -1
    int kElectronRS_;

    //! Index of the solution phase in the list of phases for this surface
    int solnPhaseRS_;

    OCV_Override_input *ocv_ptr_;

    RSD_OCVmodel* OCVmodel_;

protected:

};

}
#endif

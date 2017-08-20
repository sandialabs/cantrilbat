/**
 *  @file  Electrode_FuncCurrent.h Declarations for utility residual functions used in calculating
 *         constant current conditions and in calculating open circuit potentials for nontrivial Electrochemical
 *         surface mechanisms.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_FUNCCURRENT_H
#define _ELECTRODE_FUNCCURRENT_H

#include "Electrode_Integrator.h"
#include "cantera/numerics/ResidEval.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
//! Rootfinder class that calculates the current given the voltage
/*!
 *  This is inherited from ResidEval in order to use the rootfinder there.
 *
 *  This class is used to integrate the Electrode object across the global time interval.
 *  Its main function is to calculate the current. This value is then used within the root finder to find the voltage that
 *  produces a particular current.
 *
 *  It takes as its primary object, a pointer to an Electrode object
 */
class Electrode_ECurr : public ZZCantera::ResidEval
{
public:

    //! Constructor
    /*!
     *  @param[in]           ee                  Pointer to the Electrode object
     *  @param[in]           deltaT              Value of the time step to be used
     */
    Electrode_ECurr(Electrode* ee, double deltaT) :
        ResidEval(),
        m_ee(ee),
        m_deltaT(deltaT),
        printLvl_(0),
        numStepsLast_(0)
    {
    }
    
    //! Returns the number of equations in the system
    /*!
     *  (virtual from ResidEval)
     *  @return                                  Returns the number of equations in the system
     */ 
    virtual int nEquations() const override {
        return 1;
    }

    //! Returns the number of integration steps used in the last calculation
    /*!
     *  @return                                  Returns the number of integration steps used in the last calculation
     */
    int nIntegrationSteps() const {
        return numStepsLast_;    
    }
   
    //! Evaluate the functional in the current interval, which in this case means integrate the Electrode object across
    //! the global time step
    /*!
     *  Note, there are some tolerance inputs for the integration, that should be open to be changed!
     *
     *  @param[in]           t                   Current time. Note, this parameter is unused.
     *  @param[in]           x                   Vector of state inputs. In this case, the first position is the voltage
     *                                           across the interface. This is the only input to the routine.
     *  @param[out]          r                   Vector of residual outputs. r[0] is defined as the amps produced by the 
     *                                           Electrode object during the time interval
     *
     *  @return                                  Returns 0 if the evaluation was successful
     */
    virtual int evalResidSS(const double t, const double* const x, double* const r) override {
        // set the metal and the electrolyte voltage
        m_ee->setVoltages(x[0], 0.0);
        if (printLvl_ > 0) {
            printf("Electrode_ECurr: call integrate with voltage = %20.13g\n", x[0]);
        }
        /*
         * Integrate the electrode for a certain time, m_deltaT,
         *   The tolerance is set at 1.0E-4.
         */
        numStepsLast_ = m_ee->integrate(m_deltaT, 1.0E-4);
        /*
         *  Get the amps produced by the integration.
         */
        double amps = m_ee->getIntegratedProductionRatesCurrent(srcNet);
        r[0] = amps;
#ifdef DEBUG_ELECTRODE_MODE
        int nSubs = m_ee->integrate(m_deltaT);
        if (nSubs > 1) {
            FILE* fp = fopen("iv.txt", "w");
            for (int n = 0; n < 200; n++) {
                double volts = x[0] - 0.001 + 0.00001 * n;
                m_ee->setVoltages(volts, 0.0);
                nSubs = m_ee->integrate(m_deltaT, 1.0E-5);
                amps = m_ee->getIntegratedProductionRatesCurrent(srcNet);
                fprintf(fp," %20.13g  %d %20.13g \n", volts, nSubs, amps);
            }
            fclose(fp);
            exit(-1);
        }
#endif
        if (printLvl_ > 0) {
            writelogf("Electrode_ECurr: Curr(voltage = %20.13g) = %20.13g   nsteps = %d\n", x[0], amps, numStepsLast_);
        }
	Electrode_Integrator* eeI = dynamic_cast<Electrode_Integrator*>(m_ee);
	if (eeI) {
	    if (printLvl_ > 1) {
		SubIntegrationHistory &sih = eeI->timeHistory();
		sih.print(3);
	    }
	}
        return 1;
    }

    //! Set the delta time for the interval
    /*!
     *  @param[in]           deltaT              Value of the delta time for the time interval
     */
    void set_deltaT(double deltaT) {
        m_deltaT = deltaT;
    }
private:
    //! Pointer to the underlying Electrode object
    Electrode* m_ee;

    //! Value of delta T
    double m_deltaT;

    //! Vector Value of the srcNet used to hold source terms
    double srcNet[50];
public:
    //! Print level for the object only.
    int printLvl_;
private:
    //! Number of steps taken on last integration
    int numStepsLast_;
};

//==================================================================================================================================
//! Rootfinder residual evaluation  class calculates the total electron production from a reacting surface domain
//! as a functional, given the electric potential of the metalPhase as the single unknown.
/*!
 *  This calcualtion is used in finding the open circuit potential of a reacting
 *  surface with multiple electron reactions occurring. The main purpose it to calculate the net electron production
 *  rate as a function of the applied voltage. This functional is then used within the root finder to locate
 *  the open circuit voltage.
 *
 */
class RSD_ElectronProduction : public ZZCantera::ResidEval
{
public:

    //! Constructor
    /*!
     *  @param[in]           rsd                 Pointer to the reacting surface domain
     *  @param[in]           metalPhase          Phase index of the metal phase within the PhaseList
     *  @param[in]           deltaT              Value of the time interval
     */
    RSD_ElectronProduction(ReactingSurDomain* rsd, size_t metalPhase, double deltaT) :
        ResidEval(),
        m_rsd(rsd),
        metalPhase_(metalPhase),
        m_deltaT(deltaT),
        printLvl_(0),
        electronSpeciesIndex_(npos) 
    {
        ikMetalPhase_ = rsd->PLtoKinPhaseIndex_[metalPhase];
        int nsp = rsd->nTotalSpecies();
        spNet.resize(nsp, 0.0);
        electronSpeciesIndex_ = rsd->kineticsSpeciesIndex(0, ikMetalPhase_);
    }

    //! Returns the number of equations in the system
    /*!
     *  (virtual from ResidEval)
     *  @return                                  Returns the number of equations in the system
     */
    virtual int nEquations() const override {
        return 1;
    }

    //! Evalulate the functional for the residual used in the root finder.
    /*!
     *
     *  @param[in]           t                   Current time. Note, this parameter is unused.
     *  @param[in]           x                   Vector of state inputs. In this case, the first position is the voltage
     *                                           across the interface. This is the only input to the routine.
     *  @param[out]          r                   Vector of residual outputs. r[0] is defined as the  rate of production
     *                                           of electrons from the surface.
     *
     *  @return                                  Returns 1 if the evaluation was successful
     */
    virtual int evalResidSS(const double t, const double* const x, double* const r) override {
        // set the metal and the electrolyte voltage

        ThermoPhase& mtp = m_rsd->thermo(ikMetalPhase_);
        mtp.setElectricPotential(x[0]);

        m_rsd->getNetProductionRates(&(spNet[0]));
        double eProd = spNet[electronSpeciesIndex_];
        r[0] = eProd;
        if (printLvl_) {
            writelogf("Electrode_ECurr: electronProd(voltage = %20.13g) = %20.13g\n", x[0], eProd);
        }
        return 1;
    }

    //! Set the delta time for the interval
    /*!
     *  @param[in]           deltaT              Value of the delta time for the time interval
     */
    void set_deltaT(double deltaT) {
        m_deltaT = deltaT;
    }
private:
    //! Pointer to the reacting surface domain
    ReactingSurDomain* m_rsd;
 
    //! Metal phase where the electrons reside
    /*!
     *  This is the index number within the PhaseList
     */
    size_t metalPhase_;

    //! Value of the delta time
    double m_deltaT;
public:
    //! Print level
    int printLvl_;
private:
    //! Phase index within the ReactingSurfaceDomain for the metal phase
    /*!
     *  This is the phase index within the InterfaceKinetics object for the electrons
     */
    size_t ikMetalPhase_;

    //! Vector of production rates for the species on the surface
    /*!
     *  Length: Number of species in the InterfaceKinetics object
     *  units:  kmol/m2/s
     */
    std::vector<double> spNet;
 
    //! Kinetics species Index number for the electrons species within the InterfaceKinetics object
    size_t  electronSpeciesIndex_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

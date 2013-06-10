
/*
 * $Id: FuncElectrodeCurrent.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#ifndef _ELECTRODE_FUNCCURRENT_H
#define _ELECTRODE_FUNCCURRENT_H

#include "cantera/numerics/ResidEval.h"
#include "cantera/numerics/RootFind.h"

#include "Electrode.h"

#include "ReactingSurDomain.h"
#include "Electrode_Integrator.h"

using namespace std;
using namespace Cantera;
using namespace VCSnonideal;
using namespace mdpUtil;

namespace Cantera
{

//!  Rootfinder class that calculates the current given the voltage
class Electrode_ECurr : public Cantera::ResidEval
{
public:

    //! constructor
    Electrode_ECurr(Electrode* ee, double deltaT) :
        ResidEval(),
        m_ee(ee),
        m_deltaT(deltaT),
        printLvl_(0),
        numStepsLast_(0)
    {
    }

    int nEquations() const {
        return 1;
    }
    int nIntegrationSteps() const {
        return numStepsLast_;    
    }

    virtual int evalSS(const doublereal t, const doublereal* const x,
                       doublereal* const r) {
        // set the metal and the electrolyte voltage
        m_ee->setVoltages(x[0], 0.0);
        if (printLvl_ > 2) {
            printf("Electrode_ECurr: call integrate with voltage = %20.13g\n", x[0]);
        }
        /*
         * Integrate the electrode for a certain time, m_deltaT,
         *   The tolerance is set at 1.0E-4.
         */
        numStepsLast_ =  m_ee->integrate(m_deltaT, 1.0E-4);
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
	Electrode_Integrator* eeI = dynamic_cast<Electrode_Integrator*>(m_ee);

        if (printLvl_ > 2) {
            printf("Electrode_ECurr: Curr(voltage = %20.13g) = %20.13g   nsteps = %d\n", x[0], amps, numStepsLast_);
	  
        }
	if (eeI) {
	    if (printLvl_ > 3) {
		SubIntegrationHistory &sih = eeI->timeHistory();
		sih.print(3);
	    }
	}
        return 0;
    }

    void set_deltaT(double deltaT) {
        m_deltaT = deltaT;
    }

    Electrode* m_ee;
    double m_deltaT;
    double srcNet[50];
    int printLvl_;
    int numStepsLast_;
};



//! This class calculates the total electron production from a reacting surface domain
//! as a functional, given the electric potential of the metalPhase
/*!
 *  This calcualtion is used in finding the open circuit potential of a reacting
 *  surface with multiple electron reactions occurring.
 */
class RSD_ElectronProduction : public Cantera::ResidEval
{
public:

    RSD_ElectronProduction(ReactingSurDomain* rsd, int metalPhase, double deltaT) :
        ResidEval(),
        m_rsd(rsd),
        metalPhase_(metalPhase),
        m_deltaT(deltaT),
        printLvl_(0),
        electronSpeciesIndex_(-1) {
        ikMetalPhase_ =  rsd->PLtoKinPhaseIndex_[metalPhase];
        int nsp = rsd->nTotalSpecies();
        spNet.resize(nsp, 0.0);
        electronSpeciesIndex_ = rsd->kineticsSpeciesIndex(0, ikMetalPhase_);
    }

    int nEquations() const {
        return 1;
    }

    virtual int evalSS(const doublereal t, const doublereal* const x,
                       doublereal* const r) {
        // set the metal and the electrolyte voltage

        ThermoPhase& mtp = m_rsd->thermo(ikMetalPhase_);
        mtp.setElectricPotential(x[0]);

        m_rsd->getNetProductionRates(&(spNet[0]));
        double eProd = spNet[electronSpeciesIndex_];
        r[0] = eProd;
        if (printLvl_) {
            printf("Electrode_ECurr: electronProd(voltage = %20.13g) = %20.13g\n", x[0], eProd);
        }
        return 0;
    }

    void set_deltaT(double deltaT) {
        m_deltaT = deltaT;
    }

    ReactingSurDomain* m_rsd;
    int metalPhase_;
    double m_deltaT;
    int printLvl_;
    int ikMetalPhase_;
    std::vector<double> spNet;
    int  electronSpeciesIndex_;
};

}

#endif

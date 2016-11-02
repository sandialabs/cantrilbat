/**
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_OCVMOD_REDLICHKISTERVPSSTP_H
#define CT_OCVMOD_REDLICHKISTERVPSSTP_H

#include "cantera/thermo/RedlichKisterVPSSTP.h"

#include "ReactingSurDomain.h"
#include "RSD_OCVmodel.h"

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

/**
 */

//!  RedlichKisterVPSSTP is a derived class of GibbsExcessVPSSTP that employs
//!  an overriden deltaG relationship
/*!
 *   HKM - DEAD END CURRENTLY
 *   This will turn into a templated class
 */
class OCVMOD_RedlichKisterVPSSTP : public RedlichKisterVPSSTP
{
public:
    //! Constructor
    /*!
     * This doesn't do much more than initialize constants with
     * default values.
     */
    OCVMOD_RedlichKisterVPSSTP();

    //! Construct and initialize a RedlichKisterVPSSTP ThermoPhase object
    //! directly from an xml input file
    /*!
     *
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    OCVMOD_RedlichKisterVPSSTP(const std::string& inputFile, const std::string& id = "");

    //! Construct and initialize a RedlichKisterVPSSTP ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    OCVMOD_RedlichKisterVPSSTP(XML_Node& phaseRef, const std::string& id = "");


    //! Copy constructor
    /*!
     * @param b class to be copied
     */
    OCVMOD_RedlichKisterVPSSTP(const OCVMOD_RedlichKisterVPSSTP& b);

    //! Assignment operator
    /*!
     * @param b class to be copied.
     */
    OCVMOD_RedlichKisterVPSSTP& operator=(const OCVMOD_RedlichKisterVPSSTP& b);
    //
    //  ----------------------------------------------------------------------------------------------------
    // 
    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(double* mu) const;


    //
    //  ----------------------------------------------------------------------------------------------------
    //                      EXTRA FUNCTIONS
    //  ----------------------------------------------------------------------------------------------------
    //


    void setup_Start(std::vector<double>& stoichVector, RSD_OCVmodel* OCVmodel);


    void setup_AddThermoPhase(std::vector<double>& stoichVector, ZZCantera::ThermoPhase *tp);
    
    //!
    /*!
     *    Has to be const because base class calls are from const functions
     */
    void deriveGibbsCorrection() const;

    //
    //  ----------------------------------------------------------------------------------------------------
    //                      NEW MEMBER DATA
    //  ----------------------------------------------------------------------------------------------------
    // 

    ReactingSurDomain*  rsd_;


    RSD_OCVmodel* OCVmodel_;

    std::vector<double> stoichRxn_;

    std::vector<ZZCantera::ThermoPhase*> extraTPList_;
    std::vector< std::vector<double> > extraStoichRxn_;

    size_t replacedLocalSpeciesID;

    //! Vector of gibbs values
    mutable std::vector<double> m_mu;

    mutable std::vector< std::vector<double> > m_muExtraList_;

};

}

#endif

/**
 * @file m1d_BDD_porousFlow.h
 */


#ifndef M1D_BDD_POROUSFLOW_H_
#define M1D_BDD_POROUSFLOW_H_

#include "m1d_BulkDomainDescription.h"

namespace Cantera {
  class Electrode;
  class Transport;
  class ThermoPhase;
}

namespace m1d
{
class ExtraPhase;

//! This class consists of multiple species diffusing in a time
//! dependent manner.  There is a net flow and a net electric current.
/*!
 * 
 */
class BDD_porousFlow : public BulkDomainDescription
{
public:

    //! Constructor
    /*!
     * This constructor constructs the bulk domain description from a DomainLayout object.
     *
     * In the constructor, we have typically been laying out what the unknowns are
     * and what the equations are, that are solved within the domain.
     *
     * @param dl_ptr   Pointer to the domain layout object
     */
    BDD_porousFlow(DomainLayout *dl_ptr, std::string domainFunctionName = "", std::string domainName = "");

    //! Destructor
    virtual ~BDD_porousFlow(); 

    //! Copy Constructor
    /*!
     * @param r Object to be copied
     */
    BDD_porousFlow(const BDD_porousFlow &r);

    //! Assignment operator
    /*!
     * @param r    Object to be copied
     * @return     Returns a changeable reference to the current object
     */
    BDD_porousFlow &operator=(const BDD_porousFlow &r);

    //! Read in the possible models for each domain
    /*!
     *  This procedure is done before the Equations anv variable list are set up.
     *  Needed information about what is possible is input here.
     *  We read the Cantera ThermoPhase and transport object into DomainDescriptions here.
     *
     *   We loop over volume and then surface domains.
     */
    virtual void ReadModelDescriptions();

    //! Determine the list of Equations and Variables
    /*!
     *  This routine is responsible for setting the variables:
     *    - VariableNameList
     *    - EquationNameList
     */
    virtual void SetEquationsVariablesList();


    //! Malloc and Return the object that will calculate the residual efficiently
    /*!
     * @return  Returns a pointer to the object that will calculate the residual efficiently
     */
    virtual BulkDomain1D *mallocDomain1D();

    //! This is done after the equations are set up
    /*!
     *  We loop over volume and then surface domains here.
     */
    virtual void DetermineConstitutiveModels();

    // --------------------------------------------------------------------------------------------
    //            DATA
    // --------------------------------------------------------------------------------------------

    //! Pointer to the thermo object for the electrolyte
    /*!
     *   We own this object
     */
    Cantera::ThermoPhase *ionicLiquid_;

    //! Pointer to the transport object for the electrolyte
    /*!
     * we own this object
     */
    Cantera::Transport* trans_;

    //! Pointer to the thermo object for the porous solid comprising the solid skeleton phase
    /*!
     * We own this object
     */
    Cantera::ThermoPhase* solidSkeleton_;

    //!  Total number of species in the electrolyte
    size_t nSpeciesElectrolyte_;

    //!  Identity of the species to be used for Mole Fraction summation equation
    size_t iMFS_index_;

    //!  Identity of the species index to be used for charge neutrality summation equation
    size_t iCN_index_;

    //! Vector of extra phases
    std::vector<ExtraPhase*> ExtraPhaseList_;

    //! Porosity Equation Problem Type
    /*! 
     *   List of available options : 
     *    0 ->  Constant                   
     *    1 ->  Calculated Out Of Equation System                    
     *    2 ->  Calculated in equation system    
     *    3 ->  Calculated in equation System as  part Of Mechanics system
     *          Stress equation causes change in strains causing change of geometry
     */
    int Porosity_prob_type_;

};
//=====================================================================================================================
}
//=====================================================================================================================
#endif 

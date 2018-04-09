/**
 * @file m1d_DomainDescription.h
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_DOMAINDESCRIPTION_H
#define M1D_DOMAINDESCRIPTION_H

#include "m1d_defs.h"
#include "m1d_EqnVarTypes.h"

#include <vector>
#include <string>
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

//forward declarations
class SurfDomainDescription;
class SurDomain1D;
class BulkDomain1D;
class GlobalIndices;
class DomainLayout;

//==================================================================================================================================
//! This is a light weight base class that describes a surface or bulk domain using global indexing, which will be the same for
//! all processors.
/*!
 *  The base class for the heavyweight objects that calculate the actual residuals  is called Domain1D.
 */
class DomainDescription
{

public:
    //! Constructor
    /*!
     *  Functional names of domain may be used to define objects that may be mapped to input file parameters
     *  For example, functional names for bulk domains are "anode", "separator", and "cathode".
     *
     *  @param[in]           dl_ptr              Domain layout pointer
     *  @param[in]           domainFunctionName  Functional name of domain
     *  @param[in]           domainName          name of domain
     */
    DomainDescription(DomainLayout* dl_ptr, std::string domainFunctionName = "", std::string domainName = "");

    //! Copy Constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    DomainDescription(const DomainDescription& r);


    //! Virtual destructor
    virtual ~DomainDescription();

    //! Assignment operator
    /*!
     *  @param[in]           r                   object to be copied.
     *  @return                                  Returns a reference to the current object
     */
    DomainDescription& operator=(const DomainDescription& r);

    //! Read in the possible models for each domain
    /*!
     *  (virtual from DomainDescription)
     *
     *  This procedure is done before the Equations anv variable list are set up.
     *  Needed information about what is possible is input here.
     *  We read the Cantera ThermoPhase and transport object into DomainDescriptions here.
     *
     *   We loop over volume and then surface domains.
     */
    virtual void ReadModelDescriptions();

    //! Determine the list of Equations and Variables
    /*!
     *  (virtual from DomainDescription)
     *
     *  This routine is responsible for setting the variables:
     *    - VariableNameList
     *    - EquationNameList
     */
    virtual void SetEquationsVariablesList();

    //! Set the equation descriptions. Determine connectivity with surroundings.
    /*!
     *  (virtual from DomainDescription)
     *
     *  This routine is responsible for setting the variables:
     *    - NumEquationsPerNode
     *    - EquationIndexStart_EqName
     *    - VariableIndexStart_VarName
     */
    virtual void SetEquationDescription();

    //! Finish setting up the Constitutive Models specific to this class
    /*!
     *  (virtual from DomainDescription)
     *
     *  For example, we may instantiate ThermoPhase or Transport classes here
     */
    virtual void DetermineConstitutiveModels();

    //! Reorder the variables and equations on this domain
    /*!
     *  (virtual from DomainDescription)
     *
     *  This is a necessary step for the nodal vars step.
     */
    virtual void ReorderVariablesEquations();

    // -------------------------------------------------------------------------------
    //               DATA
    // -------------------------------------------------------------------------------

    //! Number of equations
    /*!
     *   This is the number of equations per grid point
     *   within the domain. For surface domains, this is the number of
     *   unique additional equations at the node, over and above those
     *   represented by bulk domains.
     */
    int NumEquationsPerNode;

    //! Vector containing the variable names as they appear in the unknown solution for this domain
    /*!
     * This vector defines the ordering of the equations on the domain. However, this vector
     * does not define the ordering within the solution vector.
     *
     * Length = number of equations defined on the domain
     */
    std::vector<VarType> VariableNameList;

    //! Listing of the equations types as a vector
    /*!
     * This vector defines the ordering of the equations on the domain. However, this vector
     * does not define the ordering within the solution vector.
     *
     * Length = number of equations defined on the domain
     */
    std::vector<EqnType> EquationNameList;

    //! Listing of the index of the first equation
    //! of a particular type within the domain
    /*!
     *  Length is the length of the EQ_Name_Enum structure
     */
    std::vector<int> EquationIndexStart_EqName;

    //! Listing of the index of the first equation
    //! of a particular type within the domain
    /*!
     *  Length is the length of the EQ_Name_Enum structure
     *
     *  NOT_RECCOMMENDED FOR Bulk domains
     *               This indexing is not sufficiently complex for situations
     *               at the boundaries of the domains. The variables within a
     *               domain are not contiguous within the solution vector !!!!
     */
    std::vector<int> VariableIndexStart_VarName;

    //! Vector containing the info on whether these variables are algebraic constraints
    //! or whether they have a time derivative
    /*!
     * Length = number of equations defined on the domain
     */
    std::vector<int> IsAlgebraic_NE;

    //! Vector containing the info on whether these variables are arithmetically scaled variables
    /*!
     * Length = number of equations defined on the domain
     */
    std::vector<int> IsArithmeticScaled_NE;

    //! Name of the domain.
    std::string DomainName;

    //! Functional name of the domain
    std::string DomainFunctionName_;

    //! Shallow pointer to the domain layout for this Domain Description
    DomainLayout* DL_ptr_;

    //! Print level that is set through the input file.
    /*!
     *   0 -> Don't print anything
     *   1 -> Print only about significant issues going on
     *   2 -> Print status information at regular intervals.
     *   3 -> Print ShowSolution at regular intervals
     *   4 -> Print ShowSolution at all successful time steps
     *   5 -> Print additional information at each time step
     *   6 -> Print some information about each electrode object at each time step
     *   7 -> Print a lot of information about each electrode object at each time step
     */
    int SolutionBehavior_printLvl_;

    //! Level of residual information printing done to stdout
    /*!
     *   0 -> Don't print anything
     *   1 -> Print only about significant issues going on
     *   2 -> Print status information at regular intervals.
     *   3 -> Print ShowResidual at regular intervals
     *   4 -> Print ShowResidual at all successful time steps
     *   5 -> Print additional information when ShowSolution is called.
     *   6 -> Print additional information when any residual is called.
     *   7 -> Print a lot of information about each when ShowSolution is called.
     *   8 -> Print a lot of information when any base or show residual is called
     *   9 -> Print a lot of information when any residual is called
     */
    int Residual_printLvl_;

    //! Porosity equation type
    /*!
     *  This turns on the calculation of the porosity volume fraction in the equation system.
     *  This also
     *
     *  None     =                   0x00,
     *  Constant =                   0x01,
     *  CalculatedOutOfEqnSystem =   0x02,
     *  CalculatedInEqnSystem  =     0x04,
     *  PartOfMechanics =            0x08,
     *  AddedPhasesInEqnSystem =     0x16
     */
    int porosityEquationProbType_;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

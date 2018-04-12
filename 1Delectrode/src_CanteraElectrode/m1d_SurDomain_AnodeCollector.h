/**
 * @file m1d_SurDomain_AnodeCollector.h
 *  Basic object to calculate the surface residuals for surface domains.
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#ifndef M1D_SURDOMAIN_ANODECOLLECTOR_H_
#define M1D_SURDOMAIN_ANODECOLLECTOR_H_
//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_SurDomain1D.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
// Forward declarations
class NodalVars;
class LocalNodeIndices;
class BulkDomainDescription;

//==================================================================================================================================
//! Specification of a set of boundary conditions on the top of the Cu2S surface
/*!
 *
 */
class SurDomain_AnodeCollector : public SurBC_Dirichlet
{
public:
    //! Constructor
    /*!
     *
     * @param sdd   Contains the surface domain description.
     */
    SurDomain_AnodeCollector(m1d::SurfDomainDescription& sdd, int problemType);

    //! Copy Constructor
    /*!
     *
     * @param r  Item to be copied
     */
    SurDomain_AnodeCollector(const SurDomain_AnodeCollector& r);

    //! Destructor
    virtual
    ~SurDomain_AnodeCollector();

    //! Assignment operator
    /*!
     * @param r      Object to be copied into the current object
     * @return       Returns a changeable reference to the current object
     */
    SurDomain_AnodeCollector&
    operator=(const SurDomain_AnodeCollector& r);

    //! Prepare all of the indices for fast calculation of the residual
    /*!
     *  Here we collect all of the information necessary to
     *  speedily implement SpecFlag_NE and Value_NE within the
     *  residual calculation.
     *  We transfer the information from SDT_Dirichlet structure to
     * this structure for quick processing.
     */
    virtual void
    domain_prep(LocalNodeIndices* li_ptr);

    //! Basic function to calculate the residual for the domain.
    /*!
     *  We calculate the additions and/or replacement of the
     *  residual here for the equations that this dirichlet condition
     *  is responsible for.
     *
     * @param res           Output vector containing the residual
     * @param doTimeDependentResid  boolean indicating whether the time
     *                         dependent residual is requested
     * @param soln_ptr     solution vector at which the residual should be
     *                     evaluated
     * @param solnDot_ptr  solution dot vector at which the residual should
     *                     be evaluated.
     *  @param t           time
     *  @param rdelta_t    inverse of delta_t
     *  @param residType   Residual evaluation type
     */
    virtual void
    residEval(Epetra_Vector& res,
              const bool doTimeDependentResid,
              const Epetra_Vector* soln_ptr,
              const Epetra_Vector* solnDot_ptr,
              const Epetra_Vector* solnOld_ptr,
              const double t,
              const double rdelta_t,
              const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
              const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve);

    //!  Get the voltages for the metal and solution
    /*!
     *  @param solnElectrolyte  solution at the current node
     */
    void
    getVoltages(const double* const solnElectrolyte);

    //! Base class for writing the solution on the domain to a logfile.
    /*!
     *
     * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
     * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
     * @param soln_ptr             Pointer to the solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the solution vector
     * @param solnOld_ptr          Pointer to the solution vector at the old time step
     * @param residInternal _ptr   Pointer to the current value of the residual just calculated
     *                             by a special call to the residEval()
     * @param t                    time
     * @param rdelta_t             The inverse of the value of delta_t
     * @param indentSpaces         Indentation that all output should have as a starter
     * @param duplicateOnAllProcs  If this is true, all processors will include
     *                             the same log information as proc 0. If
     *                             false, the loginfo will only exist on proc 0.
     */
    virtual void
    showSolution(const Epetra_Vector* soln_GlAll_ptr,
                 const Epetra_Vector* solnDot_GlAll_ptr,
                 const Epetra_Vector* soln_ptr,
                 const Epetra_Vector* solnDot_ptr,
                 const Epetra_Vector* solnOld_ptr,
                 const Epetra_Vector_Owned* residInternal_ptr,
                 const double t,
                 const double rdelta_t,
                 int indentSpaces,
                 bool duplicateOnAllProcs = false);

    virtual void 
    eval_HeatBalance(const int ifunc,
		     const double t,
		     const double deltaT,
		     const Epetra_Vector *soln_ptr,
		     const Epetra_Vector *solnDot_ptr,
		     const Epetra_Vector *solnOld_ptr,
		     struct globalHeatBalVals& dVals);

    //! Generate the initial conditions
    /*!
     *   For surface dirichlet conditions, we impose the t = 0- condition.
     *
     * @param doTimeDependentResid    Boolean indicating whether we should
     *                                formulate the time dependent residual
     * @param soln                    Solution vector. This is the input to
     *                                the residual calculation.
     * @param solnDot                 Solution vector. This is the input to
     *                                the residual calculation.
     * @param t                       Time
     * @param delta_t                 delta_t for the initial time step
     */
    virtual void
    initialConditions(const bool doTimeDependentResid,
                      Epetra_Vector* soln,
                      Epetra_Vector* solnDot,
                      const double t,
                      const double delta_t);

    //! Class for saving the solution on the domain in an xml node.
    /*!
     *
     * @param oNode                Reference to the XML_Node
     * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
     * @param t                    time
     *
     * @param duplicateOnAllProcs  If this is true, all processors will include
     *                             the same XML_Node information as proc 0. If
     *                             false, the xml_node info will only exist on proc 0.
     */
    virtual void
    saveDomain(ZZCantera::XML_Node& oNode,
               const Epetra_Vector* soln_GlAll_ptr,
               const Epetra_Vector* solnDot_GlAll_ptr,
               const double t,
               bool duplicateOnAllProcs = false);


    // ******************************************************************************
    //  Member Data for this boundary condition
    // ******************************************************************************
protected:

    //! Pointer to the bulk domain description object
    //! for the porous anode-electrolyte bulk region
    BulkDomainDescription* bedd_;

public:
    //! voltage electrolyte
    double phiElectrolyte_;

    //! voltage of the electrode-collector interface
    double phiAnode_;

    //! Voltage at the edge of the anode current collector
    double phiAnodeCC_;

    //! Current at the collector
    double icurrCollector_;

    //! Thickness of the copper anode current collector
    double CCThickness_;

    //! Extra anode resistance put onto the whole battery
    /*!
     *  This is not on a per area basis.
     *
     *  Units: ohms
     */
    double extraAnodeResistance_;

    //! Temperature of the collector (Kelvin)
    double TempCollector_;
};
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
#endif

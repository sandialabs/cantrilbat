/*
 * m1d_infPorousLiKCl_LiSiAnode_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#ifndef M1D_INFPOROUSLIKCL_LISIANODE_DOM1D_H_
#define M1D_INFPOROUSLIKCL_LISIANODE_DOM1D_H_

//! This is derived class  provides the function
//! evaluation for a porous electrolyte bulk domain.
/*!
 * The porous electrolyte domain is characterized by a 
 * current conservation equation and several species 
 * conservation equations describing the electrolyte.
 * A porosity/tortuosity is also associated with the domain.
 */

#include <cantera/transport.h>      // transport properties
#include <cantera/thermo.h>      // transport properties
#include <cantera/thermo/IonsFromNeutralVPSSTP.h>  // ion properties
#include "m1d_DomainDescription.h"
#include "m1d_BulkDomain1D.h"
#include "m1d_porousElectrode_dom1D.h"
#include "Electrode.h"
#include "m1d_BDT_infPorAnode_LiKCl.h"

//======================================================================================================================
namespace Cantera
{
class Electrode;
}
//======================================================================================================================
namespace m1d
{
class LocalNodeIndices;
//======================================================================================================================
//! Class for solving residuals for bulk domains
/*!
 * There is a 1 to 1 mapping between the local control volume indexing and
 * the Local Consecutive Ordering indexing
 *
 * There is a 1 to 1 mapping between the global control volume indexing
 * and the Global node number indexing that is given by a single offset.
 */
class infPorousLiKCl_LiSiAnode_dom1D : public porousElectrode_dom1D
{

public:

    //! Constructor
    /*!
     * @param bdd   Contains the bulk domain description.
     */
    infPorousLiKCl_LiSiAnode_dom1D(m1d::BDT_infPorAnode_LiKCl* bdd_anode_ptr);

  //! Copy constructor
  /*!
   * @param r      Object to be copied into the current object
   */
  infPorousLiKCl_LiSiAnode_dom1D(const infPorousLiKCl_LiSiAnode_dom1D &r);

  //! Destructor
  virtual
  ~infPorousLiKCl_LiSiAnode_dom1D();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  infPorousLiKCl_LiSiAnode_dom1D &
  operator=(const infPorousLiKCl_LiSiAnode_dom1D &r);

  //! Prepare all of the indices for fast calculation of the residual
  /*!
   *  Ok, at this point, we will have figured out the number of equations
   *  to be calculated at each node point. The object NodalVars will have
   *  been fully formed.
   *
   *  We use this to figure out what local node numbers/ cell numbers are
   *  needed and to set up indices for their efficient calling.
   *
   *  Child objects of this one will normally call this routine in a
   *  recursive fashion.
   */
  virtual void
  domain_prep(LocalNodeIndices *li_ptr);

   virtual void instantiateElectrodeCells();

  virtual void advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                                   const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr,
                                   const double t, const double t_old);


  //! Basic function to calculate the residual for the domain.
  /*!
   *  All residual terms are written with the following sign convention
   *  based on keeping the time derivative term positive.
   *
   *       res = dcdt - dc2 /dx2 - src = 0
   *
   * @param res  Output vector containing the residual
   * @param doTimeDependentResid  boolean indicating whether the time
   *                         dependent residual is requested
   * @param soln_ptr     solution vector at which the residual should be
   *                     evaluated
   * @param solnDot_ptr  solution dot vector at which the residual should
   *                     be evaluated.
   *  @param t           time
   *  @param rdelta_t    inverse of delta_t
   */
  virtual void
  residEval(Epetra_Vector &res,
            const bool doTimeDependentResid,
            const Epetra_Vector *soln_ptr,
            const Epetra_Vector *solnDot_ptr,
            const Epetra_Vector *solnOld_ptr,
            const double t,
            const double rdelta_t,
            const ResidEval_Type_Enum residType = Base_ResidEval,
	    const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

  //!  Calculate the electrode reaction rates and store it in internal variables
  void
  calcElectrode();

  //!  Setup shop at a particular point in the domain, calculating intermediate quantites
  //!  and updating Cantera's objects
  /*!
   *  All member data with the suffix, _Curr_, are updated by this function.
   *
   * @param solnElectrolyte_Curr  Current value of the solution vector
   * @param type                  Type of call
   *                              0 - at the current cell center
   */
  void
  SetupThermoShop1Old(const double* const solnElectrolyte_Curr, int type);
  void
  SetupThermoShop1(const NodalVars* const nv, const double* const soln_Curr);


  //!  Setup shop at a particular point in the domain, calculating intermediate quantites
  //!  and updating Cantera's objects
  /*!
   *  This routine will set up shop at a point intermediate to a left and a right point
   *  All member data with the suffix, _Curr_, are updated by this function.
   *
   * @param solnElectrolyte_CurrL  Current value of the solution vector at the left side
   * @param solnElectrolyte_CurrR  Current value of the solution vector at the right side
   * @param type                  Type of call
   *                              0 - at the left cell boundary
   *                              1 - at the right cell boundary
   */
  void
  SetupThermoShop2Old(const double* const solnElectrolyte_CurrL, const double* const solnElectrolyte_CurrR, int type);

  //! Calculate gradients and fluxes at the current point
  /*!
   * @param xdel            size of the cell
   * @param type            Type of call
   *                              0 - at the left cell boundary
   *                              1 - at the right cell boundary
   */
  void
  SetupTranShop(const double xdel, const int type);

  //! Function updates the ThermoPhase object for the electrolyte
  //! given the solution vector
  /*!
   *
   * @param solnElectrolyte
   */
  void
  updateElectrolyteOld(const double* const solnElectrolyte);
  void
  updateElectrolyte(const NodalVars* const nv, const double* const solnElectrolyte);


  //! Functions updates the Electrode object from the current values that are storred within the object
  void
  updateElectrode();

  //! Retrieves the voltages from the solution vector and puts them into local storage
  /*!
   * @param solnElectrolyte start of the solution vector at the current node
   */
  void
  getVoltagesOld(const double * const solnElectrolyte);

  void
  getMFElectrolyte_solnOld(const double* const solnElectrolyte);

  //! Base class for saving the solution on the domain in an xml node.
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
             const Epetra_Vector *soln_GlAll_ptr,
             const Epetra_Vector *solnDot_GlAll_ptr,
             const double t,
             bool duplicateOnAllProcs = false);

  // Method for writing the header for the surface domain to a tecplot file.
  /*
   * Only proc0 will write tecplot files.
   */
  void writeSolutionTecplotHeader();

  // Method for writing the solution on the surface domain to a tecplot file.
  /*
   * Only proc0 will write tecplot files.
   *
   * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
   * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
   * @param t                    time
   *
   */
  void writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr,
			    const Epetra_Vector *solnDot_GlAll_ptr,
			    const double t );

  //! Class for writing the solution on the domain to a logfile.
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
  showSolution(const Epetra_Vector *soln_GlAll_ptr,
               const Epetra_Vector *solnDot_GlAll_ptr,
               const Epetra_Vector *soln_ptr,
               const Epetra_Vector *solnDot_ptr,
               const Epetra_Vector *solnOld_ptr,
               const Epetra_Vector_Owned *residInternal_ptr,
               const double t,
               const double rdelta_t,
               int indentSpaces,
               bool duplicateOnAllProcs = false);

  //! Generate the initial conditions
  /*!
   *   The basic algorithm is to loop over the volume domains.
   *   Then, we loop over the surface domains
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
                    Epetra_Vector *soln,
                    Epetra_Vector *solnDot,
                    const double t,
                    const double delta_t);

 //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
  /*!
   * @param atolDefault             Default atol value
   * @param soln                    Solution vector. This is a constant
   *                                the residual calculation.
   * @param atolVector              (OUTPUT) Reference for the atol vector to fill up
   */
  virtual void setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
			     Epetra_Vector_Ghosted & atolVector,
			     const Epetra_Vector_Ghosted * const atolV = 0);

  //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
  /*!
   * @param atolDefault             Default atol value
   * @param soln                    Solution vector. This is a constant
   *                                the residual calculation.
   * @param atolVector              (OUTPUT) Reference for the atol vector to fill up
   */
  virtual void setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
				     const Epetra_Vector_Ghosted & solnDot,
				     Epetra_Vector_Ghosted & atolVector_DAEInit,
				     const Epetra_Vector_Ghosted * const atolV = 0);

  //! Evaluates the atol vector used in the delta damping process.
  /*!
   *   @param relcoeff     Relative constant to multiply all terms by
   *   @param soln         current solution vector.
   *   @param atolDeltaDamping      If non-zero, this copies the vector into the object as input
   *                      The default is zero.
   */
  virtual void
  setAtolDeltaDamping(double atolDefault, double relcoeff, 
		      const Epetra_Vector_Ghosted & soln, 
		      Epetra_Vector_Ghosted & atolDeltaDamping,
		      const Epetra_Vector_Ghosted * const atolV = 0);
		    

 //! Evaluates the atol vector used in the delta damping process for the DAE problem
  /*!
   *   @param relcoeff     Relative constant to multiply all terms by
   *   @param soln         current solution vector.
   *   @param solnDot      Current solutionDot vector.
   *   @param atolDeltaDamping       If non-zero, this copies the vector into the object as input
   *                       The default is zero.
   */
  virtual void
  setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
			      const Epetra_Vector_Ghosted & soln,
			      const Epetra_Vector_Ghosted & solnDot,
			      Epetra_Vector_Ghosted & atolDeltaDamping,
			      const Epetra_Vector_Ghosted * const atolV = 0);


  /**
   * Method to check for precipitation of the salts.  
   * Returns index of offending cation or -1 if no precipitation
   */
  int checkPrecipitation(  );

protected:

   m1d::BDT_infPorAnode_LiKCl* BDT_anode_ptr_;

  //! Pointer to the thermo object for the molten salt
  /*!
   *   We do not own this object
   */
  //ZZCantera::IonsFromNeutralVPSSTP *ionicLiquid_;

  //! Pointer to the transport object for the molten salt
  /*!
   * We do not own this object
   */
  //ZZCantera::Transport* trans_;

  //! Pointer to the electrode object
  /*!
   *   We do not own the electrode object
   */
  ZZCantera::Electrode *Electrode_;

  //! Number of phases solved
  int nph_;

  //! number of species solved
  int nsp_;
  //total concentration

  double concTot_cent_;
  double concTot_cent_old_;

  //Need to define a void fraction variable here.
  //The void fraction does not change with time
  //for the electrolyte, but it will for the electrode.
  //Since I'm not sure how I want to define this porosity
  //at this instant, I will set it to a constant for now.

  // ------------------------------------------------------------------------

  //! Volume Fraction of the electrolyte within the cell
  /*!
   *  This is now a constant.
   *  Length is number of cells on the processor.
   */
  // std::vector<double> porosity_Cell_;

  //! Volume Fraction of the electrolyte within the cell at the previous time step
  /*!
   *  This is now a constant.
   *  Length is number of cells on the processor.
   */
  //std::vector<double> porosity_Cell_old_;

  //! Surface Area of the electrolyte - electrode interface within the cell
  /*!
   *  Length is number of cells on the processor.
   *  units = m2 / m2. (surface area per cross-sectional area
   */
  std::vector<double> surfaceAreaDensity_Cell_;

  //! Electrode Current per surface area of the electrode
  /*!
   *  The surface area is specifically defined as the external surface of the electrode
   *  Length is number of cells on the processor
   *  units = amps / m2
   */
  std::vector<double> icurrInterfacePerSurfaceArea_Cell_;


  //! Total concentration of the electrolyte at cell centers
  /*!
   *  Length is number of cells on the processor.
   */
  std::vector<double> concTot_Cell_;

  //!Total concentration of the electrolyte at cell centers
  /*!
   *  Length is number of cells on the processor.
   */
  std::vector<double> concTot_Cell_old_;


  std::vector<double> capacityDischarged_Cell_;

  std::vector<double> depthOfDischarge_Cell_;

  std::vector<double> capacityLeft_Cell_;

  std::vector<double> capacityZeroDoD_Cell_;


  // ------------------------------------------------------------------------
  //!  Cell storage -> storage of cell related quantities

  //! Cell index number
  // int cIndex_cc_;

  //! Axial velocity - left cell boundary
  double Fleft_cc_;
  //! Axial Velocity - right cell boundary
  double Fright_cc_;

  //! Electrostatic potential - Left cell
  double Vleft_cc_;
  //! Electrostatic potential - center cell
  double Vcent_cc_;
  //! Electrostatic potential - right cell
  double Vright_cc_;

  //! Electrostatic potential in the electrode - Left cell
  double VElectrodeLeft_cc_;
  //! Electrostatic potential in the electrode  - center cell
  double VElectrodeCent_cc_;
  //! Electrostatic potential  in the electrode - right cell
  double VElectrodeRight_cc_;

  std::vector<double> Xleft_cc_;

  std::vector<double> Xcent_cc_;

  std::vector<double> Xright_cc_;

  std::vector<double> spCharge_;

  // -----------------------------------------------------------------------
  //!  Current Thermo value of quantities at the current point

  //! Local value of the temperature
  double temp_Curr_;

  //! Local value of the pressure
  double pres_Curr_;

  //!  Current value of the voltage
  double phiElectrolyte_Curr_;

  //!  Current value of the voltage in the anode
  double phiElectrode_Curr_;

  //!  Total concentration
  double concTot_Curr_;

  //! Current porosity
  double porosity_Curr_;

  //! Electrical conductivity of the electrode
  /*!
   *   units are S m-1
   *   The default is 1 E6
   */
  double conductivityElectrode_;

  // --------------------------------------------------------------------------
  //!  Current transport values of quantities at the current point

  //! Gradient of the temperature
  double gradT_trCurr_;

  // Gradient of the Electric Potential
  double gradV_trCurr_;

  //! Gradient of the potential in the electrode phase
  double gradVElectrode_trCurr_;

  // Gradient of the Mole fraction
  std::vector<double> gradX_trCurr_;

  std::vector<double> Vdiff_trCurr_;

  //! diffusive flux
  /*!
   *
   */
  std::vector<double> jFlux_trCurr_;

  //! current flow in the electrode due to conduction
  double icurrElectrode_trCurr_;

  //! Electrode production rates for all species in the Electrode object
  /*!
   *   Length is the number of species in the electrode object
   */
  std::vector<double> electrodeSpeciesProdRates_;

  //! Current value of the interface current from the electrode in the cell
  /*!
   *   Units are amps / m2
   */
  double icurrInterface_Curr_;

  //! Phase mole fluxes from the electrode reactions
  /*!
   *  Vector of mole fluxes for each phase in the electrode
   *  units = kmol /m2 sec
   */
  std::vector<double> phaseMoleFlux_;

  //! soln Phase mole fluxes from the electrode reactions
  /*!
   *  units = kmol /m2 sec
   */
  double solnMoleFluxInterface_Curr_;

  //! Electrode Current at the cell boundaries
  std::vector<double> icurrElectrode_CBL_;
  std::vector<double> icurrElectrode_CBR_;

  //! Electrolyte Current at the cell boundaries - left
  std::vector<double> icurrElectrolyte_CBL_;

  //! Electrolyte Current at the cell boundaries - right
  std::vector<double> icurrElectrolyte_CBR_;

  //! Electrolyte Current at the cell boundaries - right
  std::vector<double> deltaV_Cell_;
  std::vector<double> Ess_Cell_;
  std::vector<double> overpotential_Cell_;
  std::vector<double> icurrRxn_Cell_;
  std::vector<double> LiFlux_Cell_;

  // --------------------------------------------------------------------------

  std::vector<double> solnTemp;

private:
  void
  err(const char *msg);

public:

};
//======================================================================================================================
}
//======================================================================================================================
#endif // M1D_POROUSLIKCL_LISIANODE_DOM1D_H_

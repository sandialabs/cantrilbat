/**
 * @file m1d_porousFlow_dom1D.h
 */

/*
 *   $Id: m1d_porousFlow_dom1D.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_POROUSFLOW_DOM1D_H_
#define M1D_POROUSFLOW_DOM1D_H_

#include "m1d_BulkDomain1D.h"

#include <cantera/transport.h>    



namespace m1d
{
class LocalNodeIndices;
class cellTmps;

//======================================================================================================================
//! This is derived class  provides the function
//! evaluation for a porous electrode.
/*!
 * The porous electrolyte domain is characterized by a 
 * current conservation equation and several species 
 * conservation equations describing the electrolyte.
 * A porosity/tortuosity is also associated with the domain.
 *
 * There is a 1 to 1 mapping between the local control volume indexing and
 * the Local Consecutive Ordering indexing
 *
 * There is a 1 to 1 mapping between the global control volume indexing
 * and the Global node number indexing that is given by a single offset.
 */
class porousFlow_dom1D : public BulkDomain1D
{

public:

  //! Constructor
  /*!
   * @param bdd   Contains the bulk domain description.
   */
  porousFlow_dom1D(m1d::BulkDomainDescription &bdd);

  //! Copy constructor
  /*!
   * @param r      Object to be copied into the current object
   */
  porousFlow_dom1D(const porousFlow_dom1D &r);

  //! Destructor
  virtual  ~porousFlow_dom1D();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  porousFlow_dom1D&
  operator=(const porousFlow_dom1D&r);

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

  //! Function that gets called at end the start of every time step
  /*!
   *  This function provides a hook for a residual that gets called whenever a
   *  time step has been accepted and we are about to move on to the next time step.
   *  The call is made with the current time as the time
   *  that is accepted. The old time may be obtained from t and rdelta_t_accepted.
   *
   *  After this call interrogation of the previous time step's results will not be
   *  valid.
   *
   *   @param  doTimeDependentResid  This is true if we are solving a time dependent
   *                                 problem.
   *   @param  soln_ptr              Solution value at the current time
   *   @param  solnDot_ptr           derivative of the solution at the current time.
   *   @param  solnOld_ptr           Solution value at the old time step, n-1
   *   @param  t                     current time to be accepted, n
   *   @param  t_old                 previous time step value
   */
  virtual void
  advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                      const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr,
                      const double t, const double t_old);


  virtual double heatSourceLastStep() const;
  
  virtual double heatSourceAccumulated() const;

  virtual void heatSourceZeroAccumulated() const;

  //! Set up tmps for quick calculation of residuals
  void residSetupTmps();

  //! Calculates and returns an estimate of the effective areal resistance of the layer
  /*!
   *   (virtual in porous_flow_dom1D)
   *
   *   resistance = ((potCathodic - potAnodic) - voltOCV) / current
   *
   *  @param potAnodic potential in the anodic direction. If the anode, this returns the potential of the
   *                   solid in the anode next to the anode current collector.
   *  @param potCathodic potential in the cathode direction. If the anode, this returns the potential of the
   *                   electrolyte in the anode next to the separator.
   *  @param voltOCV  OCV calculated in a quick manner. 
   *  @param current   Returns the total current going through the domain ( amps m-2)
   *  
   *  @return returns the effective resistance of the layer (ohm m2)
   */
  virtual double effResistanceLayer(double &potAnodic, double  &potCathodic, double &voltOCV, double &current);


  // -------------------------------------------------------------------------------------------------------------------
  // -----------------------------------------   DATA   ----------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------

  // ------------------  Options for Processing ------------------------------------------------------------------------

  //! Boolean indicating whether we are solving the enthalpy conservation equation on the domain
  int doEnthalpyEquation_;

  // ------------------- Thermodynamics quantities on the domain -------------------------------------------------------

  //!  Partial molar Heat Capacity  of the electrolyte species located in all of the cells
  /*!
   *   Vector of partial molar heat capacity const press (KRSpecies, iCell)
   *   Units of Joules/(kmol K)
   */
  std::vector<doublereal> CpPM_lyte_Cell_;

  //!  Partial molar Enthalpy  of the electrolyte species located in all of the cells
  /*!
   *   Vector of partial molar enthalpy  (KRSpecies, iCell)
   *   Units of Joules/(kmol)
   */
  std::vector<doublereal> EnthPM_lyte_Cell_;

  //
  // ------------------- Porosity of the Domain -----------------------------------------------------------------------

protected:
  //! Volume Fraction of the electrolyte within each control volume
  /*!
   * (change to CV)
   *  Length is number of cells on the processor.
   */
  std::vector<double> porosity_Cell_;

  //! Volume Fraction of the electrolyte within the cell at the previous time step
  /*!
   *  Length is number of cells on the processor.
   */
  std::vector<double> porosity_Cell_old_;

  //
  // ------------------- Locally derived quantities that are valid at the point of current interest --------------------
  //                         ( these are intermediate values and all have the suffix _Curr_ )
  //                         ( these all refer to the new time value)
  //                         ( these are all calculated by the routine SetupShopThermo1()

  //! Temperature at the current point (Kelvin)
  double temp_Curr_;

  //! Local value of the pressure (Pascal)
  double pres_Curr_;

  //!  Total concentration of the electrolyte at the current position (kmol m-3)
  double concTot_Curr_;

  //!  Current value of the electrolyte voltage (volts) 
  double phiElectrolyte_Curr_;

  //! Current porosity
  double porosity_Curr_;

  //! Vector of temporary indexing quantities for each cell
  /*!
   * These are calculated once at the start of the program
   */
  std::vector<cellTmps> cellTmpsVect_Cell_;

  //
  // ---------------------------
  //

  //! Source of heat during the current time step
  /*!
   *  qSource is the heat generated in each cell per time for the current time step.
   *  Therefore, the source term is integrated in the axial direction and it is
   *  integrated wrt to the time interval.
   *
   *  Length = number of local cells
   *  units = Joules / m2
   */
  mutable std::vector<double> qSource_Cell_curr_;

  //! Source of heat during the global interval
  /*!
   *  qSource is the heat generated in each cell per time and accumulated over 
   *  the time steps that make up the current interval.
   *  Therefore, the source term is integrated in the axial direction and it is
   *  integrated wrt to the time interval.
   *
   *  Length = number of local cells
   *  units = Joules / m2
   */
  mutable std::vector<double> qSource_Cell_accumul_;

  //! Source of joule heating in the electrolyte during the current time interval for the current cell
  /*!
   *  This is the heat generated due to joule heating in the electrolyte in each cell per time for the current time step
   *  Therefore, the source term is integrated in the axial direction and it is
   *  integrated wrt to the time interval.
   *
   *  Length = number of local cells
   *  units = Joules / m2
   */
  std::vector<double> jouleHeat_lyte_Cell_curr_;

  //! Source of joule heating in the solid during the current time interval for the current cell
  /*!
   *  This is the heat generated due to joule heating in the electrolyte in each cell per time for the current time step
   *  Therefore, the source term is integrated in the axial direction and it is
   *  integrated wrt to the time interval.
   *
   *  Length = number of local cells
   *  units = Joules / m2
   */
  std::vector<double> jouleHeat_solid_Cell_curr_;

  //! Source of heat release in the electrode during the current time interval for the current cell
  /*!
   *  This is the heat generated due to diffusion and reaction in the electrode
   *  in each cell per time for the current time step
   *  Therefore, the source term is integrated in the axial direction and it is
   *  integrated wrt to the time interval.
   *
   *  Length = number of local cells
   *  units = Joules / m2
   */
  std::vector<double> electrodeHeat_Cell_curr_; 

  //! Source of overpotential heat release in the electrode during the current time interval for the current cell
  /*!
   *  This is the heat generated due to irreversible reaction and diffusion processes in the electrode
   *  in each cell per time for the current time step
   *  Therefore, the source term is integrated in the axial direction and it is
   *  integrated wrt to the time interval.
   *
   *  Length = number of local cells
   *  units = Joules / m2
   */
  std::vector<double> overPotentialHeat_Cell_curr_;

  //! Source of reversible entropic heat release in the electrode during the current time interval for the current cell
  /*!
   *  This is the heat generated due to reversible reaction and diffusion processes in the electrode
   *  in each cell per time for the current time step
   *  Therefore, the source term is integrated in the axial direction and it is
   *  integrated wrt to the time interval.
   *
   *  Length = number of local cells
   *  units = Joules / m2
   */
  std::vector<double> deltaSHeat_Cell_curr_;

  //!  Anodic electric potential of the domain
  /*!
   *  This is the potential that is closest to the anode side of the domain.
   *  If we are in the anode, this is the potential of the elctrode near the anode
   *  current collector.
   *
   *  Quantity is evaluated in post processor
   */
  double potentialAnodic_;

  //!  Cathodic electric potential of the domain
  /*!
   *  This is the potential that is closest to the cathode side of the domain.
   *
   *  Quantity is evaluated in post processor
   */
  double potentialCathodic_;

  //!  Total enthalpy within each cell at the current conditions (Joules/m2)
  /*!
   *   This only depends on the current conditions of temperatue, pressure, mole numbers of species, and volume
   *   fractions. It is an extensive quantity.
   *
   *   This quantity is malloced always. But it is only calculated when doEnthalpyEquation_ is true.
   *
   *   Length = total number of cells.
   */
  std::vector<double> nEnthalpy_New_Cell_;

  //!  Total enthalpy within each cell at the previous time step (Joules/m2)
  /*!
   *   This only depends on the previous time step conditions of temperatue, pressure, mole numbers of species, and volume
   *   fractions. It is an extensive quantity.
   *
   *   This quantity is malloced always. But it is only calculated when doEnthalpyEquation_ is true.
   *
   *   Length = total number of cells.
   */
  std::vector<double> nEnthalpy_Old_Cell_;

  //! Velocity basis of the transport equations
  Cantera::VelocityBasis ivb_;

};
//======================================================================================================================
}
//======================================================================================================================
#endif 

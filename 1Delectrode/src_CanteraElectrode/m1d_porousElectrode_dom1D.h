/**
 * @file m1d_porousElectrode_dom1D.h
 */

#ifndef M1D_POROUSELECTRODE_DOM1D_H_
#define M1D_POROUSELECTRODE_DOM1D_H_

#include "m1d_porousFlow_dom1D.h"
#include "m1d_BDD_porousElectrode.h"
#include "cantera/base/Array.h"

//#include <cantera/transport.h>    

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
  class Electrode;
}

namespace m1d
{
class LocalNodeIndices;
class BDD_porousElectrode;

//=================================================================================================================================
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
class porousElectrode_dom1D : public porousFlow_dom1D
{

public:

    //! Constructor
    /*!
     * @param bdd   Contains the bulk domain description.
     */
    porousElectrode_dom1D(BDD_porousElectrode* bdd_pe_ptr_);
    
  //! Copy constructor
  /*!
   * @param r      Object to be copied into the current object
   */
  porousElectrode_dom1D(const porousElectrode_dom1D &r);

  //! Destructor
  virtual  ~porousElectrode_dom1D();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  porousElectrode_dom1D& operator=(const porousElectrode_dom1D&r);

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
  virtual void domain_prep(LocalNodeIndices *li_ptr);


  //!  An electrode object must be created and initialized for every cell in the domain
  /*!
   *      Create electrode objects for every cell by calling the electrode factory. Other related
   *      quantities are also initialized such as the porosity within the cell, the solid composition
   *      which is storred in the electrode object, and the cell size. The electrode cross sectional
   *      area of the cell is calculated.
   *
   *      Correct the volume and number of moles of  active material within each of these electrode 
   *      objects to correspond to the discretized volume.
   *
   *      This object is built and initialized according to the cathode.inp or anode.inp input files.
   *      They will be later changed to different initial conditions when we set the initial conditions.
   *
   *      Values filled in:
   *
   *              porosity_Cell_[iCell]   Initial porosity of the cell
   *              Electrode_Cell_[iCell]  Pointer to the electrode object which is initialized
   *                                      within this routine.
   *              xdelCell_Cell_[iCell]   Thickness of the cell
   */
  virtual void instantiateElectrodeCells();

  //! Return the  Maximum number of normal electrode subgrid integration steps taken in the last base residual
  /*!
   *   (birth and deaths of phases aren't counted)
   */
  virtual int getMaxSubGridTimeSteps() const;

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


    //! Returns the total capacity of the electrode in Amp seconds per cross-sectional area
    /*!
     *  Returns the capacity of the electrode in Amps seconds m-2.
     *  The PA stands for "per cross-sectional area".
     *  This is the same as the number of coulombs that can be delivered at any voltage.
     *  Note, this number differs from the capacity of electrodes that is usually quoted for
     *  a battery. That number depends on the rate of discharge and also depends on the
     *  specification of a cutoff voltage. Here, we dispense with both of these specifications.
     *  So, it should be considered a theoretical capacity at zero current and minimal cutoff voltage
     *  considering the current state of the battery. The initial theoretical capacity given
     *  ideal conditions is given by capacityInitial().
     *
     *  It will also include all plateaus that are defined by the electrode object.
     *
     *  This capacity may change as degradation mechanisms cause the electrode to lose capability.
     *  Therefore, the capacity will be a function of time.
     *  At all times the following relation holds:
     *
     *  capacity() = capacityDischarged() + capacityLeft().
     *
     *  The algorithm that is used is to sum up the individual cell electrode capacity() calculations.
     *  Then, divide by the cross sectional area.
     *
     *  @param platNum   Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity and
     *                   the capacity is stated for that plateau.
     *
     *  @return returns the theoretical capacity of the electrode in Amp seconds m-2 = coulombs m-2
     */
    virtual double capacityPA(int platNum = -1) const;

    //! Amount of charge that the electrode has discharged up to this point (coulombs) per cross-sectional area
    /*!
     *   We report the number in terms of Amp seconds = coulombs.
     *   Note the capacity discharged for cathodes will be defined as the negative of the electron
     *   source term, as this refers to the forward discharge of a cathode.
     *   This definition is necessary for the following to be true.
     *
     *         capacity() = capacityDischarged() + capacityLeft()
     *
     *   Note, the current is defined as the amount
     *   of positive charge that goes from the solid into the electrolyte.
     *
     *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity.
     */
    virtual double capacityDischargedPA(int platNum = -1) const;
 
    //! Amount of charge that the electrode that has available to be discharged per cross-sectional area
    /*!
     *  We report the number in terms of Amp seconds = coulombs. This accounts for loss mechanisms.
     *
     *        At all times the following relation holds:
     *
     *  capacity() = capacityDischarged() + capacityLeft() + depthOfDischargeStarting().
     *
     *    If there is capacity lost, this loss is reflected both in the capacityLeft() and depthOfDischargeStarting()
     *    quantities so that the above relation holds.
     *
     *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity.
     */
    virtual double capacityLeftPA(int platNum = -1, double voltsMax = 50.0, double voltsMin = -50.0) const;

    //! Report the current depth of discharge in Amp seconds per cross-sectional area
    /*!
     *  Report the current depth of discharge. This is roughly equal to the total
     *  number of electrons that has been theoretically discharged from a fully charged state.
     *  For multiple cycles, this becomes the true electron counter for the electrode.
     *
     *  Usually this is reported as a function of the discharge rate and there is a
     *  cutoff voltage at which the electron counting is turned off. Neither of these
     *  concepts is employed here.
     *
     *  The depth of discharge may be modified when there is capacity lost.
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *
     *  @return  returns the depth of discharge in Amp seconds m-2
     */
    virtual double depthOfDischargePA(int platNum = -1) const;

    //! Initial starting depth of discharge in coulombs per cross sectional area
    /*!
     *   When there is capacity lost, this number may be modified.
     *
     *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity.
     */
    virtual double depthOfDischargeStartingPA(int platNum = -1) const;

    //! Reset the counters that keep track of the amount of discharge to date
    virtual void resetCapacityDischargedToDate();

    //! Return a value for the open circuit potential without doing a formally correct calculation
    /*!
     *  Currently this is defined as the open circuit potential on the outside electrode.
     *
     *   @return return the open circuit potential 
     */ 
    virtual double openCircuitPotentialQuick() const;

    //! Calculate the porosity
    virtual double calcPorosity(size_t iCell);

    //! Add terms to the polarization analysis
    /*!
     *  @param[in]         phiCurrentCollector   Value of the Electric potential in the solid current collector
     *  @param[in]         region                Region of the analysis:
     *                                                0 = anode
     *                                                2 = cathode
     */
    virtual void doPolarizationAdditions(double phiCurrentCollector, int region);

    // -------------------------------------------------------------------------------------------------------------------
    // -----------------------------------------   DATA   ----------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------------------
protected:

    //! Pointer to the Most Derived Domain Description
    BDD_porousElectrode* BDD_PE_ptr_;
    
    //! Electrode Cell data for the anode/cathode cells
    /*!
     *  Length is the number of owned Cells on the processor
     *   -> not sure what owned here means. This will have to be figured out in the future.
     *   Note, in contrast to the other quantities in the calculation, the electrode object is extrinsic
     *   it has a volume associated with it.  The volume is equal to the control volume thickness
     *   multiplied by the cross-sectional area.
     */
    std::vector<ZZCantera::Electrode*> Electrode_Cell_;

    //! Maximum number of normal electrode subgrid integration steps taken in the last base residual
    //!  (birth and deaths of phases aren't counted)
    int maxElectrodeSubIntegrationSteps_;

    //! Surface area of the electrolyte - electrode interface within the cell
    /*!
     *  This is only used for printing output from showSolution.
     *  Length is number of cells on the processor.
     *  units = m2 / m2. surface area per cross sectional area in each cell
     */
    std::vector<double> surfaceArea_Cell_;

    //!  Total enthalpy within each electrode in each cell at the current conditions (Joules/m2)
    /*!
     *   This only depends on the current conditions of temperatue, pressure, mole numbers of species, and volume
     *   fractions. It is an extensive quantity.
     *
     *   Length = total number of cells on proc
     */
    std::vector<double> nEnthalpy_Electrode_New_Cell_;

    //!  Total enthalpy within each electrode in each cell at the previous time step (Joules/m2)
    /*!
     *   This only depends on the previous time step conditions of temperatue, pressure, mole numbers of species, and volume
     *   fractions. It is an extensive quantity.
     *
     *   Length: total number of cells on proc
     */
    std::vector<double> nEnthalpy_Electrode_Old_Cell_;

    //!  Total volume of the electrode in each cell at the no stress condition (m3)
    /*!
     *   This only depends on the current conditions of temperatue, pressure, mole numbers of species.
     *   It is an extensive quantity. It is calculated from electrode->SolidVol()
     *
     *   Length: total number of cells on proc
     */
    std::vector<double> nVol_zeroStress_Electrode_Cell_;

    //
    // ------------------- Locally derived quantities that are valid at the point of current interest --------------------
    //
    //                         ( these are intermediate values and all have the suffix _Curr_ )
    //                         ( these all refer to the new time value)
    //                         ( these are all calculated by the routine SetupShopThermo1() or SetupShopThermo2()
    //                         ( these are either node quantities or boundary quantities )


    //! Heat flux of the Enhanced Enthalpy through the solid at the current point
    /*!
     *  Units = Watts / m2
     */
    double jFlux_EnthalpyPhi_metal_trCurr_;

    //!  Partial molar Enthalpy of the metal electron species at the Curr point with Phi enhancement
    /*!
     *   Vector of partial molar enthalpy  (KRSpecies)
     *   Units of Joules/(kmol)
     */
    std::vector<double> EnthalpyPhiPM_metal_Curr_;

    //! Counter to calculate the number of electrode subcycles
    /*!
     *  Vector is over the number of electrode objects
     */
    std::vector<int> numElectrodeSubCycles_Cell_;

    //!  Counter that stores the moles of elements in each of the electrode objects
    /*!
     *
     */
    ZZCantera::Array2D elem_Solid_Old_Cell_;

    //! Pointer to the metal phase that does electrical conduction within the solid
    /*!
     *  We do not own this object
     */
    ZZCantera::ThermoPhase* metalPhase_;

    //!  Total volume of the electrode in each cell at the no stress condition from old time step (m3)
    /*!
     *   This only depends on the current conditions of temperatue, pressure, mole numbers of species.
     *   It is an extensive quantity. It is calculated from electrode->SolidVol()
     *
     *   Length = total number of cells.
     */
    std::vector<double> nVol_zeroStress_Electrode_Old_Cell_;
    // Debugging -> last part keeps getting written into
    size_t wBufff_[5];

};
//======================================================================================================================
}
//======================================================================================================================
#endif 

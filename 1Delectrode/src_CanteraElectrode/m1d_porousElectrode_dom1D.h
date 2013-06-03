/**
 * @file m1d_porousElectrode_dom1D.h
 */

/*
 *   $Id: m1d_porousElectrode_dom1D.h 564 2013-03-08 23:35:51Z hkmoffa $
 */

#ifndef M1D_POROUSELECTRODE_DOM1D_H_
#define M1D_POROUSELECTRODE_DOM1D_H_

#include "m1d_porousFlow_dom1D.h"

#include <cantera/transport.h>    

namespace Cantera
{
  class Electrode;
}

namespace m1d
{
class LocalNodeIndices;

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
class porousElectrode_dom1D : public porousFlow_dom1D
{

public:

  //! Constructor
  /*!
   * @param bdd   Contains the bulk domain description.
   */
  porousElectrode_dom1D(m1d::BulkDomainDescription &bdd);

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
  porousElectrode_dom1D&
  operator=(const porousElectrode_dom1D&r);

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
   *              electrodeCrossSectionalArea_
   *              xdelCell_Cell_[iCell]   Thickness of the cell
   */
  virtual void instantiateElectrodeCells();

    //! Return the  Maximum number of normal electrode subgrid integration steps taken in the last base residual
    /*!
     *   (birth and deaths of phases aren't counted)
     */
    virtual int getMaxSubGridTimeSteps() const;

protected:

    //! Maximum number of normal electrode subgrid integration steps taken in the last base residual
    //!  (birth and deaths of phases aren't counted)
    int maxElectrodeSubIntegrationSteps_;


};
//======================================================================================================================
}
//======================================================================================================================
#endif 

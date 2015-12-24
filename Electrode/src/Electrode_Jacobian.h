/*
 * $Id: Electrode_Jacobian.h
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_JACOBIAN_H
#define _ELECTRODE_JACOBIAN_H

#include "Electrode.h"

namespace Cantera {

//! Base class for computing Electrode object sensitivities
/*!
 *  This is a base class defining an interface that is used to get sensitivity
 *  values of source terms provided by an Electrode object with respect to the
 *  various possible dofs.
 *  The different available source terms are specified using the SOURCES enum,
 *  similarly the different dofs are specified using the DOFS enum.
 *  The basic usage is as follows:
 *    1) Create an Electrode_Jacobian object for the electrode that you need
 *       sensitivities to.
 *    2) Specify what sensitivities will be needed using add_entries_to_compute,
 *       add_entry_to_compute, and remove_entry_to_compute.
 *    3) As needed call compute_jacobian to update the sensitivity values at
 *       a given set of dof values.
 *    4) Access the calculated jacobian values using get_jacobian_value
 */
class Electrode_Jacobian {

public:

  //! This pair definition marries an independent variable specififed by a DOF enum with a source term specified by the SOURCES enum.
  //!  The two of them together signifies a Jacobian term (i.e., an entry in a 2D matrix). 
  typedef std::pair<DOFS, SOURCES> DOF_SOURCE_PAIR;

  //! Constructor
  /*!
   *  @param[in]   elect               Takes a pointer to the underlying electrode object
   */
  Electrode_Jacobian(Electrode* elect);

  //! Copy Constructor
  Electrode_Jacobian(const Electrode_Jacobian& right);

  //! Assignment Operator
  Electrode_Jacobian& operator=(const Electrode_Jacobian& right);

  //! Destructor
  virtual ~Electrode_Jacobian();

  std::string dofsString(enum DOFS dd) const;

  //! Compute the Jacobian at the point specified by centerpoint
  /*!
   *  The array centerpoint should contain the value of each dof where
   *  centerpoint[DOFS] = dof_value, where DOFS refers to a value from the DOFS enum.
   *  centerpoint.size() == (MAX_DOF + n_species - 1) (i.e. a value for all possible dofs must be specified even if
   *  Jacobian entries are not being computed for that dof)
   *  Species mole fractions are specified starting at centerpoint[SPECIES] and should be in the order expected by
   *  electrode->setElectrolyteMoleNumbers
   */
  virtual void compute_jacobian(const std::vector<double> & centerpoint, const double dt) = 0;

  //! Print the jacobian out as a table to stdout
  /*!
   *  @param[in]  indentSpaces      Number of indent spaces
   */
  virtual void print_jacobian(int indentSpaces = 0) const = 0;

  /*!
   * Return the partial derivative of the requested source term with respect to the requested dof,
   * requested source and dof are specifed by dof_source_pair
   *
   *  @param[in]   dof_source_pair                  DOF_SOURCE_PAIR pair representing the row and variable
   */
  double get_jacobian_value(const DOF_SOURCE_PAIR &dof_source_pair) const
  {
      std::map< DOF_SOURCE_PAIR, double >::const_iterator it = jacobian.find(dof_source_pair);
      if (it == jacobian.end()) {
	  throw CanteraError("Electrode_Jacobian::get_jacobian_value", "Jacobian Entry not computed");
      }
      return it->second;
  }

  // These 3 functions enable the user to specify which Jacobian entries need to be
  // calculated.
  virtual void add_entries_to_compute(const std::vector<DOF_SOURCE_PAIR> &entries);
  virtual void add_entry_to_compute(DOF_SOURCE_PAIR entry);
  virtual void remove_entry_to_compute(DOF_SOURCE_PAIR entry);

protected:

  //! Pointer to the electrode object
  Electrode*  electrode;

  //! Index within the electrode object for the start of the electrolyte species
  int electrolytePhaseSpeciesStart;

  // Store the desired Jacobian contributions as a map from [dof, source] -> result
  std::map< DOF_SOURCE_PAIR, double > jacobian;


  //! Jacobian centerpoint storage
  std::vector<double> jac_centerpoint;

  //! Jacobian delta_t
  double jac_dt;

  //! Jacobian t_init_init;
  double jac_t_init_init;

  ThermoPhase* tp_solnPhase;

  int jac_numSubs_Max;
  int jac_numSubs_Min;

  double jac_energySource;
  double jac_electrolytePhaseSource;
  double  jac_electronSource;
  std::vector<double> jac_lyteSpeciesSource;


private:
 
};

}
#endif
/*****************************************************************************/


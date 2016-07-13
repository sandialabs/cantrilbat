/*
 * Electrode_FD_Jacobian.h
 *
 *  Created on: Jun 10, 2013
 *      Author: vebruni
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef ELECTRODE_FD_JACOBIAN_H_
#define ELECTRODE_FD_JACOBIAN_H_

#include "Electrode_Jacobian.h"

#include <list>
#include <map>

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

//! Implementation of the Electrode_Jacobian interface that uses numerical finite differencing to calculate source term 
//! sensitivity values.
class Electrode_FD_Jacobian : public Electrode_Jacobian 
{
public:

    //! Constructor
    /*!
     *   @param[in]    elect                Pointer to the electrode object
     *   @param[in]    baseRelDelta         double for the calculation of the base delta - relative delta
     */
    Electrode_FD_Jacobian(Electrode* elect, double baseRelDelta);

    Electrode_FD_Jacobian(const Electrode_FD_Jacobian& right);

    Electrode_FD_Jacobian& operator=(const Electrode_FD_Jacobian& right);

    //! Destructor
    virtual ~Electrode_FD_Jacobian();

    //! Create a default setup of the Jacobian
    /*!
     *  @param[out]   centerpoint            Value of the dofs
     */
    virtual void default_setup(std::vector<double>& centerpoint);
 
    //!  Create a vector of DOFS that is used in the jacobian calculation
    /*!
     *  @param[out]      centerpoint             On exit, this contains the external variables at t = t_final that
     *                                           are DOFS.
     */
    virtual void default_dofs_fill(std::vector<double>& centerpoint);
  
    //! Calculate a two sided jacobian and store it in the object
    /*!
     *    @param[in]       centerpoint           External varialbes to be used for the calculation
     *    @param[in]       dt                    Delta T
     *    @param[in,out]   dof_Deltas            Input deltas for the dofs, or output dofs.
     *    @param[in]       useDefaultDeltas      Boolean indicating whether deltas are computed or input
     */
    virtual void compute_jacobian(const std::vector<double>& centerpoint, const double dt, double* dof_Deltas = 0,
				  bool useDefaultDeltas = true);

    //! Calculate a one sided jacobian and store it in the object
    /*!
     *    @param[in]       centerpoint           External varialbes to be used for the calculation
     *    @param[in]       dt                    Delta T
     */
    virtual void compute_oneSided_jacobian(const std::vector<double>& centerpoint, const double dt, double* dof_Deltas = 0,
					   bool useDefaultDeltas = true, bool baseAlreadyCalculated = false);

    //! Store the base calculation if the Jacobian object is being used as a storage vehicle
    //! for an externally evaluated jacobian calculation
    /*!
     *
     */
    virtual void store_base_calculation(const std::vector<double>& centerpoint, const double dt,
					const std::vector<double>& speciesSources, double const enthalpySrc);

    //! Add entry to compute
    virtual void add_entry_to_compute(DOF_SOURCE_PAIR entry);

    //! remove entry to compute
    virtual void remove_entry_to_compute(DOF_SOURCE_PAIR entry);

    //! Set the atol used in the perturbation calculation to a number
    /*!
     *   The default atols used in the perturbation calculation.
     *
     *  @param[in]         dof_Atol              Vector of atols length = number of DOF.
     */
    void set_Atol(const std::vector<double>& dof_Atol);

    //! Calculate default Atol values
    /*!
     *   @param[in]      centerpoint                   Default DOF values
     *   @param[out]     dof_Atol                      calculated default values returned to user.
     *                                                 The default is zero, in which case the values are not returned to the user
     */
    void calc_dof_Atol(const std::vector<double>& centerpoint, double* const dof_Atol = 0);

    //! Calculate the perturbations that should be used in the calculation of the jacobian
    virtual void calc_Perturbations(const std::vector<double>& centerpoint, std::vector<double>& dof_Deltas, 
                                    double base_RelDelta = 1.0E-4, double* dof_Atol = 0);

    //! @copydoc Electrode_Jacobian::print_jacobian
    virtual void print_jacobian(int indentSpaces = 0) const;

private:

    // This helper function handles running a single electrode->integrate call at
    // a given set of dof values.
    /*!
     *  @return              Returns the number of subintegration steps
     */
    int run_electrode_integration(const std::vector<double> & dof_values, double dt, bool base = true);

    // This helper function sets a specific jacobian entry using a centered difference formula
    // based on the specified source values and delta in the dof.
    void set_jacobian_entry(const DOF_SOURCE_PAIR & entry, double source_values[2], double delta);

protected:

    //! These variables are used to determine which dofs we need to finite difference with respect to.
    std::list<DOFS> dofs_to_fd;

    //! Map between the DOFS identification and the number of jacobians needed
    std::map<DOFS, int> num_sources_using_dof;

    //! Relative delta to use to create the jacobian
    double base_RelDelta;

    //! Deltas used on the last jacobian
    std::vector<double> jac_Delta;

    //! Atols used on the last external jacobian delta calculation
    std::vector<double> jac_dof_Atol;
};
//==================================================================================================================================
} // end of namespace
//----------------------------------------------------------------------------------------------------------------------------------
#endif /* ELECTRODE_FD_JACOBIAN_H_ */

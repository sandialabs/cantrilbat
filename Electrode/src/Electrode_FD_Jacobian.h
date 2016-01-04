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

namespace Cantera {

//! Implementation of the Electrode_Jacobian interface that uses numerical finite differencing to calculate source term 
//! sensitivity values.
class Electrode_FD_Jacobian : public Electrode_Jacobian 
{
public:

    //! Constructor
    /*!
     *   @param[in]    elect                Pointer to the electrode object
     *   @param[in]    baseRelDelta         Double for the 
     */
    Electrode_FD_Jacobian(Electrode* elect, double baseRelDelta);

    Electrode_FD_Jacobian(const Electrode_FD_Jacobian& right);

    Electrode_FD_Jacobian& operator=(const Electrode_FD_Jacobian& right);

    //! Destructor
    virtual ~Electrode_FD_Jacobian();

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
					   bool useDefaultDeltas = true);

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

    //! Calculate the perturbations
    virtual void calc_Perturbations(const std::vector<double>& centerpoint, std::vector<double>& dof_Deltas, 
                                    std::vector<double>& dof_Atol, double base_RelDelta = 1.0E-4);

    //! @copydoc Electrode_Jacobian::print_jacobian
    virtual void print_jacobian(int indentSpaces = 0) const;

private:

    // This helper function handles running a single electrode->integrate call at
    // a given set of dof values.
    /*!
     *  @return              Returns the number of subintegration steps
     */
    int run_electrode_integration(const std::vector<double> & dof_values, double dt);

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

    //! Atols used on the last jacobian delta calculation
    std::vector<double> jac_dof_Atol;
};

} // namespace Cantera

#endif /* ELECTRODE_FD_JACOBIAN_H_ */

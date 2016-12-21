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

    //! Constructor
    /*!
     *  @param[in]             right               Takes a pointer to the underlying electrode object
     */
    Electrode_FD_Jacobian(const Electrode_FD_Jacobian& right);

    //! Assignment Operator
    /*!
     *  @param[in]             right               Object to be copied
     *
     *  @return                                    Returns a reference to the current object
     */
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
     *  @param[out]      centerpoint             On exit, this contains the external variables at t = t_final that are DOFS.
     */
    virtual void default_dofs_fill(std::vector<double>& centerpoint);
  
    //! Calculate a two sided jacobian and store it in the object
    /*!
     *    @param[in]       centerpoint           External varialbes to be used for the calculation
     *    @param[in]       dt                    Delta T
     *    @param[in,out]   dof_Deltas            Input deltas for the dofs, or output dofs.
     *    @param[in]       useDefaultDeltas      Boolean indicating whether deltas are computed or input
     */
    virtual void compute_jacobian(const std::vector<double>& centerpoint, const double dt, double* dof_Deltas = nullptr,
				  bool useDefaultDeltas = true) override;

    //! Calculate a one sided jacobian and store it in the object
    /*!
     *  @param[in]           centerpoint         External varialbes to be used for the calculation
     *  @param[in]           dt                  Delta T
     *  @param[in]           dof_Deltas          Input deltas for the dofs, or output dofs. Defaults to nullptr
     *  @param[in]           useDefaultDeltas    Boolean indicating whether deltas are computed or input. Defaults to true.
     *  @param[in]           baseAlreadyCalculated Boolean indicating whether the base calculation has already been done and storred
     *                                             internally within the object.
     */
    virtual void compute_oneSided_jacobian(const std::vector<double>& centerpoint, const double dt, double* dof_Deltas = nullptr,
					   bool useDefaultDeltas = true, bool baseAlreadyCalculated = false) override;

    //! Store the base calculation if the Jacobian object is being used as a storage vehicle
    //! for an externally evaluated jacobian calculation
    /*!
     *  (THIS IS EMPTY?)
     *
     *  @param[in]           centerpoint         The current vector of values of the external variables to be used for the calculation
     *  @param[in]           dt                  Delta T
     *  @param[in]           speciesSources      Vector of the species sources
     *  @param[in]           enthalpySrc         value of the enthalpy source term
     */
    virtual void store_base_calculation(const std::vector<double>& centerpoint, const double dt,
					const std::vector<double>& speciesSources, double const enthalpySrc);

    //! Add entry to compute
    /*!
     *  @param[in]           entry               DOF_SOURCE_PAIR entry to compute
     */
    virtual void add_entry_to_compute(DOF_SOURCE_PAIR entry) override;

    //! Remove entry to compute
    /*!
     *  @param[in]           entry               DOF_SOURCE_PAIR entry to remove from computation
     */
    virtual void remove_entry_to_compute(DOF_SOURCE_PAIR entry) override;

    //! Set the atol used in the perturbation calculation to a number
    /*!
     *   The default atols used in the perturbation calculation.
     *
     *  @param[in]           dof_Atol            Vector of atols length = number of DOF.
     */
    void set_Atol(const std::vector<double>& dof_Atol);

    //! Calculate default Atol values
    /*!
     *  @param[in]           centerpoint         Default DOF values
     *  @param[out]          dof_Atol            Calculated default values returned to user.
     *                                           The default is zero, in which case the values are not returned to the user
     */
    void calc_dof_Atol(const std::vector<double>& centerpoint, double* const dof_Atol = 0);

    //! Calculate the perturbations that should be used in the calculation of the jacobian
    /*!
     *  There is a lot of heuristic  experience involved with the calculation of perturbation values. There is not much writeup, here
     *  though about that experience.
     *
     *  @param[in]             centerpoint         The current vector of values of the external variables to be used for the calculation
     *  @param[in,out]         dof_Deltas          Input deltas for the dofs, or output dofs.
     *
     *  @param[in]             base_RelDelta       Base relative delta to use in the calculation of the perturbations of the species.
     *                                             The default is 1.0E-4. This is quite a bit larger the sqrt machine precision
     *                                             value usually suggested in textbooks. However, the thought here is that this is a
     *                                             very noisy base calculation. In order to get a good result the perturbation has to be
     *                                             greater than the noise.
     *
     *  @param[in]             dof_Atol            Vector of atols. Defaults to nullptr. Default value of the atol for the perturbation
     *                                             is taken from jac_dof_Atol[].
     */
    virtual void calc_Perturbations(const std::vector<double>& centerpoint, std::vector<double>& dof_Deltas, 
                                    double base_RelDelta = 1.0E-4, double* dof_Atol = nullptr);

    //!  Print jacobians
    /*!
     *  @param[in]           indentSpaces        Indent spacies for printout
     */
    virtual void print_jacobian(int indentSpaces = 0) const override;

private:

    // This helper function handles running a single electrode->integrate() call at a given set of dof values.
    /*!
     *  @param[in]           dof_values          The current vector of values of the external variables to be used for
     *                                           the calculation
     *  @param[in]           dt                  Delta T. Amount of time to integrate 
     *  @param[in]           base                Boolean indicating that this is a base calculation, versus a perturbation
     *                                           calculation. For base calculations, the time step history is kept.
     *                                           Defaults to true.
     *
     *  @return                                  Returns the number of subintegration steps
     */
    int run_electrode_integration(const std::vector<double> & dof_values, double dt, bool base = true);

    //! This helper function sets a specific jacobian entry using a centered difference formula
    //! based on the specified source values and delta in the dof.
    /*!
     *  @param[in]           entry               DOF_SOURCE_PAIR entry to compute
     *  @param[in]           source_values       Vector, length 2, of perturbed plus value and perturbed minus value
     *                                           of the source term
     *  @param[in]           delta               Delta value of the dof
     */
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

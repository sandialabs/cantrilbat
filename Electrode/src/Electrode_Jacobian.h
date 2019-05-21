/**
 *  @file Electrode_Jacobian.h
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

//----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{

//==================================================================================================================================
//! Structure for holding the SOURCE_TYPES and subindex for specifying the pair
struct ST_pair
{
    //! Constructor
    /*!
     *    @param[in]         st                  enum class specifying the SOURCE_TYPES
     *    @param[in]         ksp                 Species index within the lyte thermoPhase class when needed
     */
    ST_pair(SOURCE_TYPES st = SOURCE_TYPES::MAX_SOURCE, size_t ksp = npos) :
        sourceType(st),
        subIndex(ksp)
    {
    };

    //! Equality operator required for find
    /*!
     *  subIndex not pertinent if not SPECIES_SOURCE
     *
     *  @param[in]           alt                 Other ST_pair object to compare to
     *  @return                                  Returns true if (*this == alt);
     */
    bool operator==(ST_pair alt) const {
        if (sourceType != alt.sourceType) return false;
        if (subIndex != alt.subIndex) {
            if (sourceType == SOURCE_TYPES::SPECIES_SOURCE) return false;
        }
        return true;
    }

    //! Not equals operator required for maps
    /*!
     *  subIndex not pertinent if not SPECIES_SOURCE
     *
     *  @param[in]           alt                 Other ST_pair object to compare to
     *  @return                                  Returns true if (*this != alt);
     */
    bool operator!=(ST_pair alt) const {
        return ! (operator==(alt));
    }

    //! Less than operator required for maps
    /*!
     *  subIndex not pertinent if not SPECIES_SOURCE
     *
     *  @param[in]           alt                 Other ST_pair object to compare to
     *  @return                                  Returns true if (*this < alt);
     */
    bool operator<(ST_pair alt) const {
        if (sourceType < alt.sourceType) return true; 
        if (subIndex < alt.subIndex) {
            if (sourceType == alt.sourceType) return true;
        } 
        return false;
    }

    //! enum class specifying the SOURCE_TYPES
    SOURCE_TYPES sourceType;

    //!  Species index within the lyte thermoPhase class when needed
    size_t subIndex;
};

//==================================================================================================================================
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

  //! This pair definition marries an independent variable specified by a DOF enum with a source term specified
  //! by the ST_pair structure, which is the SOURCE_TYPES enum and subIndex pair
  //! The two of them together signifies a Jacobian term (i.e., an entry in a 2D matrix). 
  /*!
   *  The DOFS and SOURCE_TYPES enum is defined in the Electrode.h file, while the ST_pair structure is defined 
   *  at the top of this file
   */
  typedef std::pair<DOFS, ST_pair> DOF_SOURCE_PAIR;

  //! Constructor
  /*!
   *  @param[in]             elect               Takes a pointer to the underlying electrode object
   */
  Electrode_Jacobian(Electrode* elect);

  //! Copy Constructor
  /*!
   *  @param[in]             right               Object to be copied
   */
  Electrode_Jacobian(const Electrode_Jacobian& right);

  //! Assignment Operator
  /*!
   *  @param[in]             right               Object to be copied
   *
   *  @return                                    Returns a reference to the current object
   */
  Electrode_Jacobian& operator=(const Electrode_Jacobian& right);

  //! Destructor
  virtual ~Electrode_Jacobian();

  //! Returns a string given an enum DOFS value
  /*!
   *  The enum DOFS are used to represent degrees of freedom for the independent unknowns, that we are trying to find the jacobian
   *  values for.
   *
   *  @param[in]             dd                  enum DOFS value
   *
   *  @return                                    Returns a string representing the enum
   */
  std::string dofsString(enum DOFS dd) const;

  //! Create a default setup of the Jacobian
  /*!
   *  @param[out]            centerpoint         Value of the dofs
   */ 
  virtual void default_setup(std::vector<double>& centerpoint) = 0;

  //! Compute the Jacobian at the point specified by centerpoint
  /*!
   *  The array centerpoint should contain the value of each dof where
   *  centerpoint[DOFS] = dof_value, where DOFS refers to a value from the DOFS enum.
   *  centerpoint.size() == (MAX_DOF + n_species - 1) (i.e. a value for all possible dofs must be specified even if
   *  Jacobian entries are not being computed for that dof)
   *  Species mole fractions are specified starting at centerpoint[SPECIES] and should be in the order expected by
   *  electrode->setElectrolyteMoleNumbers.
   *
   *  @param[in]             centerpoint         The current vector of values of the external variables to be used for the calculation
   *  @param[in]             dt                  Delta T
   *  @param[in,out]         dof_Deltas          Input deltas for the dofs, or output dofs.
   *  @param[in]             useDefaultDeltas    Boolean indicating whether deltas are computed or input
   */
  virtual void compute_jacobian(const std::vector<double>& centerpoint, const double dt,
                                double* dof_Deltas = nullptr, bool useDefaultDeltas = true) = 0;

  //! Calculate a one sided jacobian and store it in the object
  /*!
   *  @param[in]             centerpoint         External varialbes to be used for the calculation
   *  @param[in]             dt                  Delta T
   *  @param[in]             dof_Deltas          Input deltas for the dofs, or output dofs. Defaults to nullptr
   *  @param[in]             useDefaultDeltas    Boolean indicating whether deltas are computed or input. Defaults to true.
   *  @param[in]             baseAlreadyCalculated Boolean indicating whether the base calculation has already been done and storred
   *                                             internally within the object.
   */
  virtual void compute_oneSided_jacobian(const std::vector<double>& centerpoint, const double dt, double* dof_Deltas = nullptr,
                                         bool useDefaultDeltas = true, bool baseAlreadyCalculated = false) = 0;

  //! Print the jacobian out as a table to stdout
  /*!
   *  @param[in]             indentSpaces        Number of indent spaces
   */
  virtual void print_jacobian(int indentSpaces = 0) const = 0;

  //! Return the partial derivative of the requested source term with respect to the requested dof,
  //! requested source and dof are specifed by dof_source_pair
  /*! 
   *  @param[in]   dof_source_pair              DOF_SOURCE_PAIR pair representing the row and variable
   *
   *  @return                                   Return the jacobian value
   */
  double get_jacobian_value(const DOF_SOURCE_PAIR &dof_source_pair) const
  {
      std::map< DOF_SOURCE_PAIR, double >::const_iterator it = jacobian.find(dof_source_pair);
      if (it == jacobian.end()) {
	  throw Electrode_Error("Electrode_Jacobian::get_jacobian_value", "Jacobian Entry not computed");
      }
      return it->second;
  }

  // These 3 functions enable the user to specify which Jacobian entries need to be calculated.

  //! This function adds a vector of entries into the list of entries that need to be computed
  /*!
   *  Each entry is a <DOF,SOURCES> pair
   *
   *  @param[in]       entries                   Vector of entries that need to be computed
   */
  virtual void add_entries_to_compute(const std::vector<DOF_SOURCE_PAIR> &entries);

  //! This function adds an entry into the list of entries that need to be computed
  /*!
   *  Each entry is a <DOF,SOURCES> pair
   *
   *  @param[in]       entry                     Sing entry that need to be computed
   */
  virtual void add_entry_to_compute(DOF_SOURCE_PAIR entry);

  //! Remove an entry that needs to be computed
  /*!
   *  Each entry is a <DOF,SOURCES> pair
   *
   *  @param[in]       entry                     Sing entry that need to be computed
   */
  virtual void remove_entry_to_compute(DOF_SOURCE_PAIR entry);

  // --------------------------------------------- D A T A ------------------------------------------------------------
protected:

  //! Pointer to the electrode object
  Electrode* electrode;

  //! Index within the electrode object for the start of the electrolyte species
  int electrolytePhaseSpeciesStart;

  // Store the desired Jacobian contributions as a map from [dof, source] -> result
  /*!
   *  The jacobian is d ( source ) / d ( dof )
   *
   *  The map is ok here, since the jacobian is small and all memory is pretty scattered at this point
   */
  std::map< DOF_SOURCE_PAIR, double > jacobian;

  //! Jacobian centerpoint storage for the DOFS
  /*!
   *  This is a vector of the values of the DOFS around which the Jacobian is being calculated.
   *  It has a length equal to the number of degrees of freedom, DOFS, which includes a vector over the electrolyte
   *  species concentrations. The length is taken from the input centerpoint vector from compute_jacobian
   */
  std::vector<double> jac_centerpoint;

  //! Jacobian delta_t
  /*!
   *  This is the global time step value over which we are computing source terms from perhaps multiple local time steps
   *  taken by the integration of the electrode object.
   */
  double jac_dt;

  //! Jacobian t_init_init;
  /*!
   *  Initial time for the global time step. The final time is (jac_t_init_init + jac_dt)
   */
  double jac_t_init_init;

  //! Pointer to the %ThermoPhase object that represents the electrolyte solution phase
  thermo_t_double* tp_solnPhase;

  //! Maximum number of subintegration steps to complete global time integration using Electrode object
  size_t jac_numSubs_Max;

  //! Minimum number of subintegration steps to complete global time integration using Electrode object
  size_t jac_numSubs_Min;

  //! Value of the enhanced enthalpy produced by the Electrode source during the intervale
  /*!
   *  Number of Joules of enthalpy produced during the interval by the electrode object
   *  Units: Joules
   */
  double jac_energySource;

  //!  Value for the electolyte phase source during the interval (Electrode is an extensive variable, so units are straight kmol)
  /*!
   *  Number of moles of electrolyte produced during the interval by the electrode object
   *  Units: kmol
   */
  double jac_electrolytePhaseSource;

  //!  Value for the electron source during the interval (Electrode is an extensive variable, so units are straight kmol)
  /*!
   *  Units: kmol
   */
  double jac_electronSource;

  //!  Value for the electolyte species source during the interval (Electrode is an extensive variable, so units are straight kmol)
  /*!
   *  Number of moles of each species produced during the interval by the electrode object
   *  Units: kmol
   *  Length: number of species in the electrolyte
   */
  std::vector<double> jac_lyteSpeciesSource;


private:
 
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif


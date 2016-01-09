/*
 * Electrode_FD_Jacobian.cpp
 *
 *  Created on: Jun 10, 2013
 *      Author: vebruni
 */

#include "Electrode_FD_Jacobian.h"
#include "Electrode.h"

namespace Cantera {

//===================================================================================================================================
static void drawline(int sp, int ll)
{
  for (int i = 0; i < sp; i++) printf(" ");
  for (int i = 0; i < ll; i++) printf("-");
  printf("\n");
}
//===================================================================================================================================
static void indent(int sp) {
    for (int i = 0; i < sp; i++) printf(" ");
}
//===================================================================================================================================
Electrode_FD_Jacobian::Electrode_FD_Jacobian(Electrode * elect, double baseRelDelta) :
    Electrode_Jacobian(elect),
    base_RelDelta(baseRelDelta)
{
}
//===================================================================================================================================
Electrode_FD_Jacobian::Electrode_FD_Jacobian(const Electrode_FD_Jacobian& right) :
    Electrode_Jacobian(right.electrode),
    base_RelDelta(right.base_RelDelta)
{
    operator=(right);
}
//===================================================================================================================================
Electrode_FD_Jacobian& Electrode_FD_Jacobian::operator=(const Electrode_FD_Jacobian& right)
{
    if (this == &right) {
	return *this;
    }
    Electrode_Jacobian::operator=(right);
    dofs_to_fd = right.dofs_to_fd;
    num_sources_using_dof = right.num_sources_using_dof;
    base_RelDelta = right.base_RelDelta;
    jac_Delta = right.jac_Delta;
    jac_dof_Atol = right.jac_dof_Atol;

    return *this;
}
//===================================================================================================================================
Electrode_FD_Jacobian::~Electrode_FD_Jacobian()
{
}
//===================================================================================================================================
//
//  This creates a full interaction matrix by default
// centerpoint is DOFs in DOFS order (Electrode.h
//
void Electrode_FD_Jacobian::default_setup(std::vector<double>& centerpoint)
{
    // Sources
   
    std::vector<DOFS> dof_vector;
    dof_vector.push_back(SOLID_VOLTAGE);
    dof_vector.push_back(LIQUID_VOLTAGE);
    dof_vector.push_back(TEMPERATURE);
    dof_vector.push_back(PRESSURE);
    size_t ip_solvent = electrode->solnPhaseIndex();
    size_t ip_metal = electrode->metalPhaseIndex();
    ThermoPhase& tpe = electrode->thermo(ip_solvent);
    size_t nsp = tpe.nSpecies();

    //
    // Fill the centerpoint vector
    //
    size_t sz = 4 + nsp;
    centerpoint.resize(4 + nsp);
    double phiMetal = electrode->phaseVoltage(ip_metal);
    centerpoint[SOLID_VOLTAGE] = phiMetal;
    double phiElectrolyte = electrode->phaseVoltage(ip_solvent);
    centerpoint[LIQUID_VOLTAGE] = phiElectrolyte;
    centerpoint[TEMPERATURE] = tpe.temperature();
    centerpoint[PRESSURE] = tpe.pressure();
    tpe.getMoleFractions(& centerpoint[SPECIES]);

     
    for (size_t k = 0; k < nsp; ++k) {
	dof_vector.push_back((DOFS)(SPECIES + k));
    }

    std::vector<SOURCES> src_vector;
    src_vector.push_back(ENTHALPY_SOURCE);
    src_vector.push_back(ELECTROLYTE_PHASE_SOURCE);
    src_vector.push_back(CURRENT_SOURCE);
    for (size_t k = 0; k < nsp; ++k) {
	src_vector.push_back((SOURCES)(SPECIES_SOURCE + k));
    }  
    DOF_SOURCE_PAIR entry;
    for (size_t i = 0; i < src_vector.size(); i++) {
	 for (size_t j = 0; j < dof_vector.size(); j++) {
	     entry = std::make_pair(dof_vector[j], src_vector[i]);
	     add_entry_to_compute(entry);

	 }
    }
    jac_Delta.resize(sz);
    calc_dof_Atol(centerpoint);
}
//===================================================================================================================================
void Electrode_FD_Jacobian::default_dofs_fill(std::vector<double>& centerpoint)
{    
    size_t ip_metal = electrode->metalPhaseIndex();
    size_t ip_solvent = electrode->solnPhaseIndex();
    ThermoPhase& tpe = electrode->thermo(ip_solvent);
    size_t nsp = tpe.nSpecies();
    centerpoint.resize(4 + nsp);
    double phiMetal = electrode->phaseVoltage(ip_metal);
    centerpoint[SOLID_VOLTAGE] = phiMetal;
    double phiElectrolyte = electrode->phaseVoltage(ip_solvent);
    centerpoint[LIQUID_VOLTAGE] = phiElectrolyte;
    centerpoint[TEMPERATURE] = tpe.temperature();
    centerpoint[PRESSURE] = tpe.pressure();
    tpe.getMoleFractions(& centerpoint[SPECIES]);
}
//===================================================================================================================================
void Electrode_FD_Jacobian::compute_jacobian(const std::vector<double> & centerpoint, const double dt,
					     double* dof_Deltas, bool useDefaultDeltas)
{
    // Temporary storage used to do the centered difference calculation:
    /*
     *  speciesSource is vector of length 2 with each entry being a vector<double>
     */
    std::vector< double > speciesSources[2];
    speciesSources[0].resize(electrode->nSpecies());
    speciesSources[1].resize(electrode->nSpecies());
    double energySource[2];
    double electrolytePhaseSource[2];
    double individualSpeciesSource[2];
    int numSubs;

    int electronIndex = electrode->kSpecElectron();

    jac_centerpoint = centerpoint;
    jac_dt = dt;
    jac_t_init_init = electrode->timeInitInit();
    std::vector<double> perturbed_point = centerpoint;
    jac_numSubs_Max = -1;
    jac_numSubs_Min = 10000000;
    // 
    //   Calculate the deltas for the jacobian and then report it back if needed
    //
    if (useDefaultDeltas) {
	calc_Perturbations(centerpoint, jac_Delta, base_RelDelta);
	if (dof_Deltas) {
	    for (size_t i = 0; i < centerpoint.size(); ++i) {
		dof_Deltas[i] = jac_Delta[i];
	    }
	}
    } else {
	for (size_t i = 0; i < centerpoint.size(); ++i) {
	    jac_Delta[i] = dof_Deltas[i];
	}
    }
    //
    // Iterate over the list of dofs that need to be finite differenced with respect to
    // For each dof store the source term values with the dof perturbed +- its centerpoint value.
    // These are then used to set the Jacobian entries using a centered difference formula
    //
    std::list<DOFS>::iterator dof_it = dofs_to_fd.begin();
    for ( ; dof_it != dofs_to_fd.end(); ++dof_it ) {

	// Determine how far to perturb the dof using base_delta as a relative perturbation magnitude
	double dof_delta = base_RelDelta;
	double dof_val = perturbed_point[*dof_it];
	double dof_absval = std::abs(dof_val);
	if (dof_absval > 0.0) {
	    dof_delta *= dof_absval;
	} 
	//dof_delta *= (std::abs(perturbed_point[*dof_it]) > 0.) ? std::abs(perturbed_point[*dof_it]) : 1.0;

	for( int negative=0; negative<2; ++negative) {
	    // Perturb by -1^negative * base_delta and compute
	    // perturbed_point[*dof_it] += dof_delta * std::pow(-1.0, negative);
	    perturbed_point[*dof_it] += jac_Delta[*dof_it] * std::pow(-1.0, negative);
	    run_electrode_integration(perturbed_point, dt);

	    // Store off the source term results in temporary storage
	    energySource[negative] = electrode->getIntegratedSourceTerm(Cantera::ENTHALPY_SOURCE);
	    electrolytePhaseSource[negative] = electrode->getIntegratedSourceTerm(Cantera::ELECTROLYTE_PHASE_SOURCE);
	    numSubs = electrode->integratedSpeciesSourceTerm(&speciesSources[negative][0]);

	    perturbed_point[*dof_it] = centerpoint[*dof_it];
	    jac_numSubs_Max = std::max( jac_numSubs_Max, numSubs);
	    jac_numSubs_Min = std::min( jac_numSubs_Min, numSubs);
	}

	// Use set_jacobian_entry to store the various sensitivities based on a centered difference.

        jac_energySource = 0.5 *( energySource[0] + energySource[1]);
	set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, ENTHALPY_SOURCE), energySource, dof_delta);

        jac_electrolytePhaseSource = 0.5 *( electrolytePhaseSource[0] + electrolytePhaseSource[1]);
	set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, ELECTROLYTE_PHASE_SOURCE), electrolytePhaseSource, dof_delta);

	individualSpeciesSource[0] = (speciesSources[0])[electronIndex];
	individualSpeciesSource[1] = speciesSources[1][electronIndex];

        jac_electronSource =  0.5 *(  individualSpeciesSource[0] + individualSpeciesSource[1]);
	set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, CURRENT_SOURCE), individualSpeciesSource, dof_delta);

	for (int sp = 0; sp < electrode->numSolnPhaseSpecies(); ++sp) {
	    int idx = sp + electrolytePhaseSpeciesStart;
	    individualSpeciesSource[0] = speciesSources[0][idx];
	    individualSpeciesSource[1] = speciesSources[1][idx];
            jac_lyteSpeciesSource[sp] = 0.5 *(  individualSpeciesSource[0] + individualSpeciesSource[1]);
	    set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, (SOURCES)(SPECIES_SOURCE + sp)), individualSpeciesSource, dof_delta);
	}
    }
}
//===================================================================================================================================
void Electrode_FD_Jacobian::compute_oneSided_jacobian(const std::vector<double> & centerpoint, const double dt,
						      double* dof_Deltas, bool useDefaultDeltas, bool baseAlreadyCalculated)
{
    // Temporary storage used to do the one-sided difference calculation:
    /*
     *  speciesSource is vector of length 2 with each entry being a vector<double>
     */
    std::vector< double > speciesSources[2];
    speciesSources[0].resize(electrode->nSpecies());
    speciesSources[1].resize(electrode->nSpecies());
    double energySource[2];
    double electrolytePhaseSource[2];
    double individualSpeciesSource[2];
    int numSubs;
    if (centerpoint.size() > jac_dof_Atol.size()) {
	jac_dof_Atol.resize(centerpoint.size(), 0.0);
    }

    int electronIndex = electrode->kSpecElectron();

    jac_centerpoint = centerpoint;
    jac_dt = dt;
    jac_t_init_init = electrode->timeInitInit();
    std::vector<double> perturbed_point = centerpoint;
   
    // 
    //   Calculate the deltas for the jacobian and then report it back if needed
    //
    if (useDefaultDeltas) {
	calc_Perturbations(centerpoint, jac_Delta, base_RelDelta);
	if (dof_Deltas) {
	    for (size_t i = 0; i < centerpoint.size(); ++i) {
		dof_Deltas[i] = jac_Delta[i];
	    }
	}
    } else {
	for (size_t i = 0; i < centerpoint.size(); ++i) {
	    jac_Delta[i] = dof_Deltas[i];
	}
    }

    //
    //  Calculate the base point
    //
    run_electrode_integration(perturbed_point, dt, true);
    //
    //  Store the results
    //
    energySource[0] = electrode->getIntegratedSourceTerm(Cantera::ENTHALPY_SOURCE);
    electrolytePhaseSource[0] = electrode->getIntegratedSourceTerm(Cantera::ELECTROLYTE_PHASE_SOURCE);
    numSubs = electrode->integratedSpeciesSourceTerm(&speciesSources[0][0]);
 
    jac_numSubs_Max = numSubs;
    jac_numSubs_Min = numSubs;
  
    // Iterate over the list of dofs that need to be finite differenced with respect to
    // For each dof store the source term values with the dof perturbed +- its centerpoint value.
    // These are then used to set the Jacobian entries using a centered difference formula
    std::list<DOFS>::iterator dof_it = dofs_to_fd.begin();
    for ( ; dof_it != dofs_to_fd.end(); ++dof_it ) {
	double dof_delta = jac_Delta[*dof_it];
	for (int negative = 1; negative < 2; ++negative) {
	    // perturb the point
	    perturbed_point[*dof_it] += dof_delta;
	    //
	    //   
	    //
	    run_electrode_integration(perturbed_point, dt, false);
	    //
	    // Store off the source term results in temporary storage
	    //
	    energySource[1] = electrode->getIntegratedSourceTerm(Cantera::ENTHALPY_SOURCE);
	    electrolytePhaseSource[1] = electrode->getIntegratedSourceTerm(Cantera::ELECTROLYTE_PHASE_SOURCE);
	    numSubs = electrode->integratedSpeciesSourceTerm(&speciesSources[1][0]);

	    perturbed_point[*dof_it] = centerpoint[*dof_it];
	    jac_numSubs_Max = std::max( jac_numSubs_Max, numSubs);
	    jac_numSubs_Min = std::min( jac_numSubs_Min, numSubs);
	}

	// Use set_jacobian_entry to store the various sensitivities based on a centered difference.

        jac_energySource = (energySource[0]);
	set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, ENTHALPY_SOURCE), energySource, -0.5 * dof_delta);

        jac_electrolytePhaseSource = (electrolytePhaseSource[0]);
	set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, ELECTROLYTE_PHASE_SOURCE), electrolytePhaseSource, -0.5 * dof_delta);

	individualSpeciesSource[0] = speciesSources[0][electronIndex];
	individualSpeciesSource[1] = speciesSources[1][electronIndex];
        jac_electronSource = (individualSpeciesSource[0]);
	set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, CURRENT_SOURCE), individualSpeciesSource, -0.5 * dof_delta);

	for (int sp = 0; sp < electrode->numSolnPhaseSpecies(); ++sp) {
	    int idx = sp + electrolytePhaseSpeciesStart;
	    individualSpeciesSource[0] = speciesSources[0][idx];
	    individualSpeciesSource[1] = speciesSources[1][idx];
            jac_lyteSpeciesSource[sp] = individualSpeciesSource[0];
	    set_jacobian_entry( DOF_SOURCE_PAIR(*dof_it, (SOURCES)(SPECIES_SOURCE + sp)), individualSpeciesSource, -0.5 * dof_delta);
	}
    }
}
//===================================================================================================================================
void Electrode_FD_Jacobian::print_jacobian(int indentSpaces) const
{
    double val;
    std::vector<bool> electrodeSpeciesSrcInc(electrode->nSpecies(), false);
    std::vector<bool> lyteSpeciesSrcInc(tp_solnPhase->nSpecies(), true);

    indent(indentSpaces);
    int cellNumber = electrode->electrodeCellNumber_;
    int domainNumber = electrode->electrodeDomainNumber_;
    printf(" Electrode jacobian: %d %d, Time_init = %g, Time_final = %g dt = %g",
	   cellNumber, domainNumber, jac_t_init_init,  jac_t_init_init + jac_dt, jac_dt);
    if (  jac_numSubs_Max ==   jac_numSubs_Min) {
	printf(" Equal numSub = %d\n",  jac_numSubs_Max);
    } else {
	printf(" Unequal numSub = %d <= num <= %d\n",  jac_numSubs_Min,  jac_numSubs_Max);
    }

    int cCount = 5 * 17 +  tp_solnPhase->nSpecies() * 17;
    for (size_t i = 0; i < (size_t) electrode->nSpecies(); i++) {
	if (electrodeSpeciesSrcInc[i]) {
	    cCount += 17;
	}
    }
    drawline(indentSpaces, cCount + 5);

    printf("   Entry      | ");
    printf(" DOF_Value   | ");
    printf("  DOF_Delta  |");
    printf("%14.14s|", "CURRENT_SOURCE");
    printf("%14.14s|", "ENTHALPY_SOURCE");
    printf("%14.14s|", "LYTE_PHASE_SRC");
  
    for (size_t i = 0; i < (size_t) electrode->nSpecies(); i++) {
	if (electrodeSpeciesSrcInc[i]) {
	    std::string ss = electrode->speciesName(i);
	    printf("%14.14s|", ss.c_str());
	}
    }
    for (size_t i = 0; i < tp_solnPhase->nSpecies(); i++) {
	if (lyteSpeciesSrcInc[i]) {
	    std::string ss = tp_solnPhase->speciesName(i);
	    printf("%16.16s|", ss.c_str());
	}
    }
    printf("\n");
    drawline(indentSpaces, cCount + 5);
    indent(indentSpaces);
    printf("              |");
    printf("              |");
    printf("              |");
    printf(" % -12.5E |", jac_electronSource);
    printf(" % -12.5E |", jac_energySource);
    printf(" % -12.5E |", jac_electrolytePhaseSource);
    for (size_t i = 0; i < tp_solnPhase->nSpecies(); i++) {
	if (lyteSpeciesSrcInc[i]) {
	    printf("   % -12.5E |", jac_lyteSpeciesSource[i]);
	}
    }
    printf("\n");
    drawline(indentSpaces, cCount + 5);

    std::list<DOFS>::const_iterator dof_it = dofs_to_fd.begin();
    int iDof = 0;
    for ( ; dof_it != dofs_to_fd.end(); ++dof_it ) {
	for (int i = 0; i < indentSpaces; i++) printf(" ");
	std::string ss = dofsString(*dof_it);
	printf("%14.14s|", ss.c_str());
	val = jac_centerpoint[iDof];
	printf(" % -12.5E |", val);
	val = jac_Delta[iDof];
	printf(" % -12.5E |", val);
	val = get_jacobian_value( DOF_SOURCE_PAIR(*dof_it, CURRENT_SOURCE) );
	printf(" % -12.5E |", val);
	val = get_jacobian_value( DOF_SOURCE_PAIR(*dof_it, ENTHALPY_SOURCE) );
	printf(" % -12.5E |", val);
	val = get_jacobian_value( DOF_SOURCE_PAIR(*dof_it, ELECTROLYTE_PHASE_SOURCE) );
	printf(" % -12.5E |", val);
  
	for (size_t i = 0; i < tp_solnPhase->nSpecies(); i++) {
	    if (lyteSpeciesSrcInc[i]) {
		val = get_jacobian_value( DOF_SOURCE_PAIR(*dof_it, (SOURCES)(SPECIES_SOURCE + i)) );
		printf("   % -12.5E |", val);
	    }
	}
	printf("\n");
	iDof++;
    }
    drawline(indentSpaces, cCount + 5);
}
//===================================================================================================================================
void Electrode_FD_Jacobian::add_entry_to_compute(DOF_SOURCE_PAIR entry)
{
    if (jacobian.find(entry) == jacobian.end() ) {
        jacobian[entry] = 0.0;
        if( std::find(dofs_to_fd.begin(), dofs_to_fd.end(), entry.first) == dofs_to_fd.end() ) {
            dofs_to_fd.push_back(entry.first);
            num_sources_using_dof[entry.first] = 1;
        } else {
            ++num_sources_using_dof[entry.first];
        }
    }
}
//===================================================================================================================================
void Electrode_FD_Jacobian::remove_entry_to_compute(DOF_SOURCE_PAIR entry)
{
    std::map< DOF_SOURCE_PAIR, double >::iterator entry_pos = jacobian.find(entry);
    if (entry_pos != jacobian.end()) {
	jacobian.erase(entry_pos);
	--num_sources_using_dof[entry.first];
	if( num_sources_using_dof[entry.first] == 0 )
	{
	    dofs_to_fd.erase( std::find(dofs_to_fd.begin(), dofs_to_fd.end(), entry.first) );
	}
    }
}
//===================================================================================================================================
void Electrode_FD_Jacobian::set_Atol(const std::vector<double>& dof_Atol)
{
    size_t sz = dof_Atol.size();
    jac_dof_Atol = dof_Atol;
    if (jac_Delta.size() < sz) {
	jac_Delta.resize(sz, 0.0);

    }
}
//===================================================================================================================================
void Electrode_FD_Jacobian::calc_dof_Atol(const std::vector<double>& centerpoint, double* const dof_Atol)
{
    //
    // Here we identify the DOF to do finite difference and then associate an atol with it
    // Getting the atol value helps.
    //
    size_t sz = centerpoint.size();
    jac_dof_Atol.resize(sz);
    //
    // Voltages (SOLID_VOLTAGE, LIQUID_VOLTAGE)
    //   This is an absolute quantity related to deltaG. We have a good idea that it must be ~0.01 * millivolt for a good delta
    //
    // double metalPot = centerpoint[SOLID_VOLTAGE);
    jac_dof_Atol[SOLID_VOLTAGE] = 1.0E-8;
    jac_dof_Atol[LIQUID_VOLTAGE] = 1.0E-8;
    //
    // Temperature we have a good idea that it should be a fractional degree ~0.01 
    //
    
    jac_dof_Atol[TEMPERATURE] = 0.00001;
    
    double pres = centerpoint[PRESSURE];
    if (pres > 1.0E-3) {
	jac_dof_Atol[PRESSURE] = 1.0E-5;
    } else {
	jac_dof_Atol[PRESSURE] = pres * 1.0E-3;
    }
    
    //
    // Electrolyte Species Mole Fractions - We don't have diffusion responsibilities here. So we can go really low with the ATOL
    // 
    size_t ip_solvent = electrode->solnPhaseIndex();
    ThermoPhase& tpe = electrode->thermo(ip_solvent);
    size_t nsp = tpe.nSpecies(); 
    for (size_t i = 0; i < nsp; i++) {
	jac_dof_Atol[SPECIES + i] = 1.0E-20;
    }
    //
    // report the result 
    //
    if (dof_Atol) {
	copy(jac_dof_Atol.begin(), jac_dof_Atol.end(), dof_Atol);
    }

}


//===================================================================================================================================
void Electrode_FD_Jacobian::calc_Perturbations(const std::vector<double>& centerpoint, std::vector<double>& dof_Deltas,
					       double base_RelDelta, double* dof_Atol)
 {
     if (dof_Atol == 0) {
	 dof_Atol = DATA_PTR(jac_dof_Atol);
     }
     //
     // Here we identify the DOF to do finite difference and then associate an atol with it
     // Getting the atol value helps.
     //
     //size_t sz = centerpoint.size();
     //
     // Voltages (SOLID_VOLTAGE, LIQUID_VOLTAGE)
     //   This is an absolute quantity related to deltaG. We have a good idea that it must be ~0.01 * millivolt for a good delta
     //
     // double metalPot = centerpoint[SOLID_VOLTAGE);

     dof_Deltas[SOLID_VOLTAGE] = 1.0E-6 + dof_Atol[SOLID_VOLTAGE];
     dof_Deltas[LIQUID_VOLTAGE] = 1.0E-6 + dof_Atol[LIQUID_VOLTAGE];
     //
     // Temperature we have a good idea that it should be a fractional degree ~0.01 
     //
     double temp = centerpoint[TEMPERATURE];
     dof_Deltas[TEMPERATURE] = 0.001 + 0.00001 * temp + dof_Atol[TEMPERATURE];

     //
     // Pressure 
     //
     double pres = centerpoint[PRESSURE];
     dof_Deltas[PRESSURE] = 0.001 + 1.0E-5 * pres + dof_Atol[PRESSURE];

     //
     // Electrolyte Species - We don't have diffusion responsibilities here. So we can really low with the delta 
     // 
     size_t ip_solvent = electrode->solnPhaseIndex();
     ThermoPhase& tpe = electrode->thermo(ip_solvent);
     size_t nsp = tpe.nSpecies(); 
     for (size_t i = 0; i < nsp; i++) {
	 double dof_absval = fabs(centerpoint[SPECIES + i]);
        if (dof_absval > 0.0) {
	    dof_Deltas[SPECIES + i] = dof_absval * base_RelDelta + dof_Atol[SPECIES + i] + 1.0E-21;
	} else {
	    dof_Deltas[SPECIES + i] = dof_Atol[SPECIES + i] + 1.0E-21;
	}
     }
 }

//===================================================================================================================================
int Electrode_FD_Jacobian::run_electrode_integration(const std::vector<double> & dof_values, double dt, bool base)
{
    double  GlobalRtolSrcTerm = 1.0E-3;
    Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS;
    Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR;
    if (!base) {
	subIntegrationType = FVDELTA_TIMEINTEGRATION_SIR;
    }

    electrode->revertToInitialTime(true);
    electrode->setState_TP(dof_values[TEMPERATURE], dof_values[PRESSURE]);
    electrode->setVoltages(dof_values[SOLID_VOLTAGE], dof_values[LIQUID_VOLTAGE]);
    electrode->setFinalStateFromInit();
    electrode->setElectrolyteMoleNumbers(&dof_values[SPECIES], true);

    int numSubs = electrode->integrate(dt, GlobalRtolSrcTerm, fieldInterpolationType, subIntegrationType);
    
    return numSubs;
}
//===================================================================================================================================
void Electrode_FD_Jacobian::set_jacobian_entry(const DOF_SOURCE_PAIR & entry, double source_values[2], double delta)
{
    // Note if the jacobian entry isn't registered, it isn't storred in the mapping.
    if (jacobian.find(entry) != jacobian.end()) {
	jacobian[entry] = (source_values[0] - source_values[1]) / (2*delta);
    }
}
//===================================================================================================================================
} // namespace Cantera

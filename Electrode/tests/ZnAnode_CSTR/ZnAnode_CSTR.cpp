/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode.h"
#include "Electrode_Factory.h"

//=====================================================================================================
int main(int argc, char **argv)
{
  // Set up electrode input file
  Cantera::Electrode_Factory * factory = Cantera::Electrode_Factory::factory();

  ELECTRODE_KEY_INPUT *eki = new ELECTRODE_KEY_INPUT();
  BEInput::BlockEntry *cfa = new BEInput::BlockEntry("command_file");

  // Initial parse of electrode input file to determine what kind of ELECTRODE_KEY_INPUT is needed
  std::string electrodeFile("anode.inp");
  eki->electrode_input(electrodeFile, cfa);
  std::string modelName = eki->electrodeModelName;

  // Create E_K_I of the proper type
  ELECTRODE_KEY_INPUT *eki2 = factory->newElectrodeKeyInputObject(modelName);
  eki2->electrode_input_child(eki->commandFile_, eki->lastBlockEntryPtr_);
  delete eki;
  eki = eki2;
  eki2 = NULL;

  Cantera::Electrode *electrode = factory->newElectrodeObject(modelName);

  // Create the model.  This is all that is required to get the porosity computed inside
  // the electrode object.  There is much more setup required if we want to do any
  // electrochemical calculations with this, but we don't.
  electrode->electrode_model_create(eki);
  electrode->setInitialConditions(eki);
  electrode->setPrintLevel(0);
  electrode->setID(0,0);
  electrode->choiceDeltaTsubcycle_init_ = 0;

  // Set electrode to a unit volume with a known reference porosity.
  electrode->setElectrodeSizeParams(1.0, 1.0, 0.25);
  electrode->setTime(0.);
  // Tell the electrode not to track the number of moles in the electrolyte
  // Setting this causes setElectrolyteMoleNumbers(...) to treat the array it is
  // passed as mole fraction values rather than total mole values
  electrode->turnOffFollowElectrolyteMoles();

  const double cur_T = 294.15;
  const double cur_P = 1.e5;
  const double cur_phi_sol = 0.;
  const double cur_phi_liq = 1.38;
  std::vector<double> moleNums(electrode->numSolnPhaseSpecies(), 0.);
  const ThermoPhase *electrolytePhase = electrode->getPhase("HMW_ZnKOH");
  if( !electrolytePhase )
  {
    std::cout << "Error: HMW_ZnKOH electrolyte phase not found." << std::endl;
    return -1;
  }
  const int k_index = electrolytePhase->speciesIndex("K+");
  const int oh_index = electrolytePhase->speciesIndex("OH-");
  const int zincate_index = electrolytePhase->speciesIndex("Zn(OH)4--");
  const int h2o_index = electrolytePhase->speciesIndex("H2O(L)");
  const double c_oh = 7.e-3;
  const double c_zincate = 5.3e-4;
  const double c_k = c_oh + 2* c_zincate;
  const double c_h2o = 4.63e-2;
  moleNums[oh_index] = c_oh;
  moleNums[zincate_index] = c_zincate;
  moleNums[k_index] = c_k;
  moleNums[h2o_index] = c_h2o;
  double totalMoleNums = 0.0;
  for(unsigned i=0; i < moleNums.size(); ++i)
  {
      totalMoleNums += moleNums[i];
  }
  // Convert the mole numbers to mole fractions since followElectrolyteMoles is off
  for(unsigned i=0; i<moleNums.size(); ++i)
    moleNums[i] /= totalMoleNums;

  double t = 0.;
  const double t_final = 1.e-3;
  const double dt = 1.e-6;
  while( t <= t_final )
  {
    electrode->resetStartingCondition(t);
    electrode->setState_TP(cur_T, cur_P);
    electrode->setVoltages(cur_phi_sol, cur_phi_liq);
    electrode->setElectrolyteMoleNumbers(&moleNums[0],true);

    // Tell Electrode to integrate the next time step
    electrode->integrate(dt);
    electrode->writeCSVData(1);
    t += dt;
  }

  delete cfa;
  delete eki;
  delete electrode;

  Cantera::appdelete();
  return 0;
} 

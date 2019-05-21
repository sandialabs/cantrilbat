/*
 * FDJacobianTests.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: vebruni
 */

#include "gtest/gtest.h"
#include "Electrode.h"
#include "Electrode_FD_Jacobian.h"

namespace Zuzax
{

class MockThermoPhase_lyte : public ThermoPhase
{
public:
  MockThermoPhase_lyte() : 
    ThermoPhase()
  {
     m_kk = 3;
  }
};

//! A mock electrode object to use in unit testing of Electrode_FD_Jacobian
class MockElectrode : public Electrode
{
public:
  MockElectrode() : Electrode()
  {
    m_NumTotSpecies = 4;
    m_NumTotPhases = 1;
    m_PhaseSpeciesStartIndex[0] = 1;
    kElectron_ = 0;
    solnPhase_ = 0;
    metalPhase_ = 1;
    phaseVoltages_.resize(2);
    fake_electrolyte_mole_nums.resize(3);
    m_NumVolPhases = 1;
    VolPhaseList_.push_back( new MockThermoPhase_lyte());
  }

  virtual ~MockElectrode() {}

  virtual void setFinalStateFromInit() override {}
  virtual void updateState() override {}

  virtual void setElectrolyteMoleNumbers(const double* const electrolyteMoleNum, bool setInitial) override
  {
    fake_electrolyte_mole_nums[0] = electrolyteMoleNum[0];
    fake_electrolyte_mole_nums[1] = electrolyteMoleNum[1];
    fake_electrolyte_mole_nums[2] = electrolyteMoleNum[2];
  }

  virtual int integrate(double deltaT, double  GlobalRtolSrcTerm = 1.0E-3,
                        Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS,
                        Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR) override
  {return 0;}

  virtual size_t numSolnPhaseSpecies() const override { return 3; }

  virtual double integratedEnthalpySourceTerm() override
  { return temperature_; }

  virtual void getIntegratedPhaseMoleTransfer(double* const phaseMolesTransfered) override
  {
    phaseMolesTransfered[0] = temperature_ + 2*deltaVoltage_;
    for(int i=0; i<3; ++i)
    {
      phaseMolesTransfered[0] += (i+3) * fake_electrolyte_mole_nums[i];
    }
  }

  virtual size_t integratedSpeciesSourceTerm(double* const spMoleDelta) override
  {
    spMoleDelta[0] = deltaVoltage_;
    spMoleDelta[1] = 2*fake_electrolyte_mole_nums[0];
    spMoleDelta[2] = 3*fake_electrolyte_mole_nums[1];
    spMoleDelta[3] = 4*fake_electrolyte_mole_nums[2];
    return 1;
  }
  
  virtual double getIntegratedSourceType(SOURCE_TYPES sourceType, size_t ksp) override
  {
    double result = 0.0;
    switch( sourceType )
    {
    case SOURCE_TYPES::ELECTROLYTE_PHASE_SOURCE:
      result = temperature_ + 2*deltaVoltage_;
      for(int i=0; i<3; ++i)
      {
        result += (i+3) * fake_electrolyte_mole_nums[i];
      }
      break;
    case SOURCE_TYPES::ENTHALPY_SOURCE:
      result = integratedEnthalpySourceTerm();
      break;
    default:
      result = 0.;
      break;
    }
    return result;
  }

  virtual void revertToInitialTime(bool) override {};

private:
  std::vector<double> fake_electrolyte_mole_nums;
};

//! Unit test fixture for Electrode_FD_Jacobian tests
class FDJacobianTest : public testing::Test
{
public:
  FDJacobianTest() :
    temp_energy_pair(TEMPERATURE, ST_pair(SOURCE_TYPES::ENTHALPY_SOURCE)),
    current_voltage_pair(SOLID_VOLTAGE, ST_pair(SOURCE_TYPES::CURRENT_SOURCE)),
    dt(0.1),
    zero(0.)
  {
    // Create a mock electrode object, an fd_jacobian object that drives it
    // and set up their initial state.
    mock_electrode = new MockElectrode();
    std::vector<Electrode_Jacobian::DOF_SOURCE_PAIR> entries_to_compute;
    entries_to_compute.push_back(temp_energy_pair);
    fd_jacobian = new Electrode_FD_Jacobian(mock_electrode, 1.e-5);
    fd_jacobian->add_entries_to_compute(entries_to_compute);

    for(int idx=0; idx < 3; ++idx)
    {
      species_source_pairs.push_back( Electrode_Jacobian::DOF_SOURCE_PAIR((DOFS)(SPECIES+idx), ST_pair(SOURCE_TYPES::SPECIES_SOURCE, idx)) );
    }

    for(int idx=0; idx < MAX_DOF + 2; ++idx)
    {
      electrolyte_phase_all_dofs.push_back( Electrode_Jacobian::DOF_SOURCE_PAIR((DOFS)(idx), ST_pair(SOURCE_TYPES::ELECTROLYTE_PHASE_SOURCE)));
    }

    point.resize(7);
    std::fill(point.begin(), point.end(), 1.);
  }
  ~FDJacobianTest()
  {
    delete fd_jacobian;
    delete mock_electrode;
  }
protected:
  Electrode_Jacobian::DOF_SOURCE_PAIR temp_energy_pair;
  Electrode_Jacobian::DOF_SOURCE_PAIR current_voltage_pair;
  std::vector<Electrode_Jacobian::DOF_SOURCE_PAIR> species_source_pairs;
  std::vector<Electrode_Jacobian::DOF_SOURCE_PAIR> electrolyte_phase_all_dofs;

  std::vector<double> point;
  double dt;
  double zero;

  Electrode *mock_electrode;
  Electrode_Jacobian *fd_jacobian;
};

//! Test that the expected exception occurs if trying to access a jacobian value that is not being computed.
TEST_F(FDJacobianTest, MissingEntry)
{
  EXPECT_THROW( fd_jacobian->get_jacobian_value(current_voltage_pair), ZuzaxError);
}

//! Test adding an entry to the list of computed jacobian entries.
TEST_F(FDJacobianTest, AddEntry)
{
  fd_jacobian->add_entry_to_compute( current_voltage_pair );
  EXPECT_DOUBLE_EQ( fd_jacobian->get_jacobian_value(current_voltage_pair), zero);
}

//! Test removing an entry from the list of computed entries
TEST_F(FDJacobianTest, RemoveEntry)
{
  ASSERT_DOUBLE_EQ( fd_jacobian->get_jacobian_value(temp_energy_pair), zero);
  fd_jacobian->remove_entry_to_compute( temp_energy_pair );
  EXPECT_THROW( fd_jacobian->get_jacobian_value(temp_energy_pair), ZuzaxError );
}

//! Test computing a single entry of the jacobian
TEST_F(FDJacobianTest, ComputeJacobian)
{
  fd_jacobian->compute_jacobian(point, dt);
  EXPECT_NEAR(1., fd_jacobian->get_jacobian_value(temp_energy_pair), 1.e-9);
}

TEST_F(FDJacobianTest, ComputeJacobianWithDofValueZero)
{
  std::fill(point.begin(), point.end(), 0.);
  fd_jacobian->compute_jacobian(point, dt);
  EXPECT_NEAR(1., fd_jacobian->get_jacobian_value(temp_energy_pair), 1.e-9);
}

//! Test computing a jacobian with multiple entries
TEST_F(FDJacobianTest, ComputeJacobianMultipleEntries)
{
  fd_jacobian->add_entry_to_compute(current_voltage_pair);
  fd_jacobian->compute_jacobian(point, dt);
  EXPECT_NEAR(1., fd_jacobian->get_jacobian_value(temp_energy_pair), 1.e-9);
  EXPECT_NEAR(1., fd_jacobian->get_jacobian_value(current_voltage_pair), 1.e-9);
}

//! Test computing jacobian entries wrt species dofs
TEST_F(FDJacobianTest, SpeciesJacobians)
{
  fd_jacobian->add_entries_to_compute(species_source_pairs);
  fd_jacobian->compute_jacobian(point, dt);
  EXPECT_NEAR(2., fd_jacobian->get_jacobian_value(species_source_pairs[0]), 1.e-9);
  EXPECT_NEAR(3., fd_jacobian->get_jacobian_value(species_source_pairs[1]), 1.e-9);
  EXPECT_NEAR(4., fd_jacobian->get_jacobian_value(species_source_pairs[2]), 1.e-9);
}

//! Test the electrolyte phase source jacobian entries wrt each dof.
TEST_F(FDJacobianTest, ElectrolytePhaseSource)
{
  Electrode_Jacobian::DOF_SOURCE_PAIR electrolyte_phase(SOLID_VOLTAGE, ST_pair(SOURCE_TYPES::ELECTROLYTE_PHASE_SOURCE));
  fd_jacobian->add_entries_to_compute(electrolyte_phase_all_dofs);
  fd_jacobian->compute_jacobian(point, dt);
  // SOLID_VOLTAGE
  EXPECT_NEAR(2., fd_jacobian->get_jacobian_value(electrolyte_phase), 1.e-9);
  electrolyte_phase.first = (DOFS)(electrolyte_phase.first + 1);
  // LIQUID_VOLTAGE
  EXPECT_NEAR(-2., fd_jacobian->get_jacobian_value(electrolyte_phase), 1.e-9);
  electrolyte_phase.first = (DOFS)(electrolyte_phase.first + 1);
  // TEMPERATURE
  EXPECT_NEAR(1., fd_jacobian->get_jacobian_value(electrolyte_phase), 1.e-9);
  electrolyte_phase.first = (DOFS)(electrolyte_phase.first + 1);
  // PRESSURE
  EXPECT_NEAR(0., fd_jacobian->get_jacobian_value(electrolyte_phase), 1.e-9);
  //EXPECT_NEAR(1., fd_jacobian->get_jacobian_value(electrolyte_phase), 1.e-9);
  electrolyte_phase.first = (DOFS)(electrolyte_phase.first + 1);
  // SPECIES 0-2
  EXPECT_NEAR(3., fd_jacobian->get_jacobian_value(electrolyte_phase), 1.e-9);
  electrolyte_phase.first = (DOFS)(electrolyte_phase.first + 1);
  EXPECT_NEAR(4., fd_jacobian->get_jacobian_value(electrolyte_phase), 1.e-9);
  electrolyte_phase.first = (DOFS)(electrolyte_phase.first + 1);
  EXPECT_NEAR(5., fd_jacobian->get_jacobian_value(electrolyte_phase), 1.e-9);
}

} // end of namespace
//----------------------------------------------------------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  printf("Running main() from FDJacobianTests.cpp\n");
  testing::InitGoogleTest(&argc, argv);
  int retnCode = RUN_ALL_TESTS();
  return retnCode;
}


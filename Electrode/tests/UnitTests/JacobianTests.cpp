/*
 * JacobianTests.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: vebruni
 */

#include "gtest/gtest.h"
#include "Electrode.h"
#include "Electrode_FD_Jacobian.h"

namespace Cantera
{

class MockElectrode : public Electrode
{
public:
  MockElectrode() : Electrode()
  {
    m_NumTotSpecies = 1;
    kElectron_ = 0;
    solnPhase_ = 0;
    metalPhase_ = 1;
    phaseVoltages_.resize(2);
  }

  virtual ~MockElectrode() {}

  virtual void setFinalStateFromInit() {}
  virtual void updateState() {}
  virtual void setElectrolyteMoleNumbers(const double* const electrolyteMoleNum, bool setInitial) {}

  virtual int integrate(double deltaT, double  GlobalRtolSrcTerm = 1.0E-3,
                        Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS,
                        Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR)
  {return 0;}

  virtual double energySourceTerm()
  { return temperature_; }

  virtual void getIntegratedPhaseMoleTransfer(doublereal* const phaseMolesTransfered)
  { phaseMolesTransfered[0] = 1.; }

  virtual double integratedSourceTerm(doublereal* const spMoleDelta)
  {
    spMoleDelta[0] = deltaVoltage_;
    return 1.0;
  }

};

class JacobianTest : public testing::Test
{
public:
  JacobianTest() :
    temp_energy_pair(TEMPERATURE, ENTHALPY_SOURCE),
    current_voltage_pair(SOLID_VOLTAGE, CURRENT_SOURCE)
  {
    mock_electrode = new MockElectrode();
    std::vector<Electrode_Jacobian::DOF_SOURCE_PAIR> entries_to_compute;
    entries_to_compute.push_back(temp_energy_pair);
    fd_jacobian = new Electrode_FD_Jacobian(mock_electrode);
    fd_jacobian->add_entries_to_compute(entries_to_compute);
  }
  ~JacobianTest()
  {
    delete fd_jacobian;
    delete mock_electrode;
  }
protected:
  Electrode_Jacobian::DOF_SOURCE_PAIR temp_energy_pair;
  Electrode_Jacobian::DOF_SOURCE_PAIR current_voltage_pair;

  Electrode *mock_electrode;
  Electrode_Jacobian *fd_jacobian;
};

TEST_F(JacobianTest, MissingEntry)
{
  EXPECT_THROW( fd_jacobian->get_jacobian_value(current_voltage_pair), CanteraError);
}

TEST_F(JacobianTest, AddEntry)
{
  fd_jacobian->add_entry_to_compute( current_voltage_pair );
  double zero = 0.;
  EXPECT_DOUBLE_EQ( fd_jacobian->get_jacobian_value(current_voltage_pair), zero);
}

TEST_F(JacobianTest, RemoveEntry)
{
  double zero = 0.;
  ASSERT_DOUBLE_EQ( fd_jacobian->get_jacobian_value(temp_energy_pair), zero);
  fd_jacobian->remove_entry_to_compute( temp_energy_pair );
  EXPECT_THROW( fd_jacobian->get_jacobian_value(temp_energy_pair), CanteraError );
}

TEST_F(JacobianTest, ComputeJacobian)
{
  double dt=0.1;
  std::vector<double> point(5);
  std::fill(point.begin(), point.end(), 1.);
  fd_jacobian->compute_jacobian(point, dt);
  EXPECT_NEAR(1., fd_jacobian->get_jacobian_value(temp_energy_pair), 1e-12);
}

/*
TEST_F(JacobianTest, updateProperties)
{
  double cp_R, h_RT, s_R;

  // Reference values calculated using CHEMKIN II
  // Expect agreement to single-precision tolerance
  set_tpow(298.15);
  poly.updateProperties(&tpow_[0], &cp_R, &h_RT, &s_R);
  EXPECT_NEAR(4.46633496, cp_R, 1e-7);
  EXPECT_NEAR(-158.739244, h_RT, 1e-5);
  EXPECT_NEAR(25.7125777, s_R, 1e-6);

  set_tpow(876.54);
  poly.updateProperties(&tpow_[0], &cp_R, &h_RT, &s_R);
  EXPECT_NEAR(6.33029000, cp_R, 1e-7);
  EXPECT_NEAR(-50.3179924, h_RT, 1e-5);
  EXPECT_NEAR(31.5401226, s_R, 1e-6);
}

TEST_F(JacobianTest, updatePropertiesTemp)
{
  double cp_R1, h_RT1, s_R1;
  double cp_R2, h_RT2, s_R2;
  double T = 481.99;

  set_tpow(T);
  poly.updatePropertiesTemp(T, &cp_R1, &h_RT1, &s_R1);
  poly.updateProperties(&tpow_[0], &cp_R2, &h_RT2, &s_R2);

  EXPECT_DOUBLE_EQ(cp_R1, cp_R2);
  EXPECT_DOUBLE_EQ(h_RT1, h_RT2);
  EXPECT_DOUBLE_EQ(s_R1, s_R2);
}
*/

} // namespace Cantera

int
main(int argc, char** argv)
{
  printf("Running main() from JacobianTests.cpp\n");
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


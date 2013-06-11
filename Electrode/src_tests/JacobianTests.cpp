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
  MockElectrode() : Electrode()
  {}

  virtual int integrate(double deltaT, double  GlobalRtolSrcTerm = 1.0E-3,
                        Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS,
                        Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR)
  {return 0;}
  virtual double energySourceTerm()
  { return temperature_; }
};

class JacobianTest : public testing::Test
{
public:
  JacobianTest() :
    temp_energy_pair(TEMPERATURE, ENTHALPY_SOURCE)
  {
    mock_electrode = new MockElectrode();
    std::vector<Electrode_Jacobian::DOF_SOURCE_PAIR> entries_to_compute;
    fd_jacobian = new Electrode_FD_Jacobian(mock_electrode, entries_to_compute);
  }
  ~JacobianTest()
  {
    delete fd_jacobian;
    delete mock_electrode;
  }
protected:
  Electrode_Jacobian::DOF_SOURCE_PAIR temp_energy_pair;

  Electrode *mock_electrode;
  Electrode_Jacobian *fd_jacobian;
};

TEST_F(JacobianTest, MissingEntry)
{
  EXPECT_THROW( fd_jacobian->get_jacobian_value(temp_energy_pair), CanteraError);
//  EXPECT_EQ(poly.speciesIndex(), (size_t) 0);
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


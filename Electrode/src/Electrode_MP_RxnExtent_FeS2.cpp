/**
 * @file Electrode_MP_RxnExtent.cpp 
 */

#include "Electrode_MP_RxnExtent_FeS2.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
Electrode_MP_RxnExtent_FeS2::Electrode_MP_RxnExtent_FeS2() :
    Electrode_MP_RxnExtent()
{
}
//==================================================================================================================================
Electrode_MP_RxnExtent_FeS2::Electrode_MP_RxnExtent_FeS2(const Electrode_MP_RxnExtent_FeS2& right) :
    Electrode_MP_RxnExtent()
{
    operator=(right);
}
//==================================================================================================================================
Electrode_MP_RxnExtent_FeS2& Electrode_MP_RxnExtent_FeS2::operator=(const Electrode_MP_RxnExtent_FeS2& right)
{
    if (this == &right) {
        return *this;
    }
    Electrode_MP_RxnExtent::operator=(right);
    return *this;
}
//==================================================================================================================================
Electrode_MP_RxnExtent_FeS2::~Electrode_MP_RxnExtent_FeS2()
{
}
//==================================================================================================================================
Electrode_Types_Enum  Electrode_MP_RxnExtent_FeS2::electrodeType() const
{
    return MP_RXNEXTENT_FES2_ET;
}
//==================================================================================================================================
double Electrode_MP_RxnExtent_FeS2::openCircuitVoltageSS_Region_NoCheck(double relativeExtentRxn, int xRegion) const
{
    double voltage, voltage_b, voltage_3, lb;
    switch (xRegion) {
    case 0:
        voltage = 1.889538 + 0.2297e-3*temperature_;
        break;
    case 1:
        voltage = 1.6613762 + 0.4178e-3*temperature_;
        break;
    case 2:
        lb = RegionBoundaries_ExtentRxn_[2];
        voltage_b = 1.6613762 + 0.4178e-3*temperature_;
        voltage_3 = 1.803338 - 0.2355e-3*temperature_;
        voltage = voltage_b + (relativeExtentRxn - lb) / (2.0 - lb) * (voltage_3 - voltage_b);
        break;
    case 3:
        voltage = 1.896438 - 0.3958e-3*temperature_;
        break;
    case -1:
	voltage = -1000;
	break;
    default:
        voltage = 1000.;
        break;
    }
    return voltage;
}
//==================================================================================================================================
} // End of namespace
//----------------------------------------------------------------------------------------------------------------------------------

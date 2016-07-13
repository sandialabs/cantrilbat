/*
 * @file         GFCEO_Electrode.cpp
 */
/*
 * Copywrite 2015 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "GFCEO_Electrode.h"


//-----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//===================================================================================================================================
GFCEO_Electrode::GFCEO_Electrode(Electrode* ee, doublereal atol, int iOwn) :
   ZZCantera::ResidJacEval(atol),
   ee_(ee),
   iOwnObject_(iOwn)
{

}
//===================================================================================================================================
GFCEO_Electrode::~GFCEO_Electrode()
{
    if (iOwnObject_) {
        delete ee_;
	ee_ = 0;
    }
}
//===================================================================================================================================
GFCEO_Electrode::GFCEO_Electrode(const GFCEO_Electrode& right) :
    ZZCantera::ResidJacEval()
{
    operator=(right);
}
//===================================================================================================================================
GFCEO_Electrode& GFCEO_Electrode::operator=(const GFCEO_Electrode& right)
{
    if (this == &right) {
        return *this;
    }
    ZZCantera::ResidJacEval::operator=(right);

    if (right.iOwnObject_) {
        delete ee_;
        ee_ = right.ee_->duplMyselfAsElectrode();
        iOwnObject_ = 1;
    } else {
        ee_ = right.ee_;
        iOwnObject_ = 0;
    }

    return *this;
}
//===================================================================================================================================
Electrode& GFCEO_Electrode::electrode()
{
    return *ee_;
}
//===================================================================================================================================
int GFCEO_Electrode::nEquations() const
{
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::constrain(const int k, const int flag)
{
// COMPLETE
}
//===================================================================================================================================
int GFCEO_Electrode::constraint(const int k) const
{
    // COMPLETE
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::initSizes()
{
// COMPLETE
}
//===================================================================================================================================
void GFCEO_Electrode::setAlgebraic(const int k)
{
// COMPLETE
}
//===================================================================================================================================
int GFCEO_Electrode::evalResidNJ(const doublereal t, const doublereal delta_t,
                            const doublereal* const y,
                            const doublereal* const ydot,
                            doublereal* const resid,
                            const ResidEval_Type_Enum evalType,
                            const int id_x,
                            const doublereal delta_x)
 {
     return 0;
 }
//===================================================================================================================================
int GFCEO_Electrode::eval(const doublereal t, const doublereal* const y, const doublereal* const ydot,
			  doublereal* const resid)
{
    return 0;
}
//===================================================================================================================================
int GFCEO_Electrode::getInitialConditions(const doublereal t0, doublereal* const y, doublereal* const ydot)
{
    return 0;
}
//===================================================================================================================================
doublereal GFCEO_Electrode::filterNewStep(const doublereal t, const doublereal* const ybase,
					  doublereal* const step)
{
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::setAtol(doublereal atol)
{
}
//===================================================================================================================================
int  GFCEO_Electrode::evalTimeTrackingEqns(const doublereal t, const doublereal delta_t, const doublereal* const y,
                                      const doublereal* const ydot)
{
    return 0;
}
//===================================================================================================================================
 bool GFCEO_Electrode::evalStoppingCritera(const doublereal t,
                                     const doublereal delta_t,
                                     const doublereal* const y,
                                     const doublereal* const ydot)
{
    return 0;
}

//===================================================================================================================================
int GFCEO_Electrode::calcDeltaSolnVariables(const doublereal t,
                           const doublereal* const y,
                           const doublereal* const ydot,
                           doublereal* const delta_y,
                           const doublereal* const solnWeights)
{
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::calcSolnScales(const doublereal t, const doublereal* const y,
				     const doublereal* const y_old, doublereal* const yScales)
{
    
}
//===================================================================================================================================
void GFCEO_Electrode::user_out2(const int ifunc, const doublereal t,
				const doublereal delta_t,
				const doublereal* const y,
				const doublereal* const ydot)
{
}
//===================================================================================================================================
void  GFCEO_Electrode::user_out(const int ifunc, const doublereal t,
				const doublereal* y,
				const doublereal* ydot)
{

}
//===================================================================================================================================
int GFCEO_Electrode::matrixConditioning(doublereal* const matrix, const int nrows,doublereal* const rhs)
{
    return 0;
}
//===================================================================================================================================
int GFCEO_Electrode::evalJacobian(const doublereal t, const doublereal delta_t, doublereal cj,
				  const doublereal* const y, const doublereal* const ydot,
				  GeneralMatrix& J, doublereal* const resid)
{
    return 0;
}
//===================================================================================================================================
int GFCEO_Electrode::evalJacobianDP(const doublereal t, const doublereal delta_t, doublereal cj,
				    const doublereal* const y,
				    const doublereal* const ydot,
				    doublereal* const* jacobianColPts,
				    doublereal* const resid)
{
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::writeSolution(int ievent, const double time, const double deltaT,  const int time_step_num,
				    const double* const y, const double* const ydot)
{

}
//===================================================================================================================================
void GFCEO_Electrode::assertOwnership()
{
     if (! iOwnObject_) {
        Electrode* e = ee_->duplMyselfAsElectrode();
        ee_ = e;
        iOwnObject_ = 1;
    } 
}
//===================================================================================================================================
} // End of namespace
//-----------------------------------------------------------------------------------------------------------------------------------

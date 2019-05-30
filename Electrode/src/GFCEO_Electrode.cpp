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

#include "zuzax/numerics/DAE_Solver.h"

//-----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{
//===================================================================================================================================
GFCEO_Electrode::GFCEO_Electrode(Electrode_Integrator* ee, double atol, int iOwn) :
   Zuzax::ResidJacEval(atol),
   ee_(ee),
   iOwnObject_(iOwn),
   integDAE_(nullptr)
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
    Zuzax::ResidJacEval()
{
    operator=(right);
}
//===================================================================================================================================
GFCEO_Electrode& GFCEO_Electrode::operator=(const GFCEO_Electrode& right)
{
    if (this == &right) {
        return *this;
    }
    Zuzax::ResidJacEval::operator=(right);

    if (right.iOwnObject_) {
        delete ee_;
        ee_ = (Electrode_Integrator*) right.ee_->duplMyselfAsElectrode();
        delete integDAE_;
        integDAE_ = right.integDAE_->duplMyselfAs_DAE_Solver(*this);
        iOwnObject_ = 1;
    } else {
        ee_ = right.ee_;
        iOwnObject_ = 0;
        integDAE_ = right.integDAE_;
    }

    return *this;
}
//===================================================================================================================================
Electrode_Integrator& GFCEO_Electrode::electrode()
{
    return *ee_;
}
//===================================================================================================================================
int GFCEO_Electrode::nEquations() const
{
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::set_DAE_Integrator(DAE_Solver* integDAE)
{
    if (integDAE_) {
        if (iOwnObject_) {
            delete integDAE_;
        }
    }
    integDAE_ = integDAE;
    ResidJacEval* r_ptr = integDAE->resid_ptr();
    if (r_ptr != this) {
         throw Electrode_Error("GFCEO_Electrode::set_DAE_Integrator()", "Residual objects are not right");
    }
}
//===================================================================================================================================
void GFCEO_Electrode::setConstraint(const int k, const int flag)
{
// COMPLETE
}
//===================================================================================================================================
int GFCEO_Electrode::isConstraint(const int k) const
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
int GFCEO_Electrode::evalResidNJ(const double t, const double delta_t,
                            const double* const y, const double* const ydot, double* const resid,
                            const ResidEval_Type evalType, const Solve_Type solveType,
                            const int id_x, const double delta_x)
{
     int retn = ee_->GFCEO_evalResidNJ(t, delta_t, y, ydot, resid, evalType, solveType, id_x, delta_x);
     return retn;
}
//===================================================================================================================================
int GFCEO_Electrode::evalResid(const double t, const double* const y, const double* const ydot,
                               double* const resid, const ResidEval_Type evalType, const Solve_Type solveType)
{
    return 0;
}
//===================================================================================================================================
int GFCEO_Electrode::getInitialConditionsWithDot(const double t0, double* const y, double* const ydot)
{
    return 0;
}
//===================================================================================================================================
double GFCEO_Electrode::filterNewStep(const double t, const double* const ybase, double* const step)
{
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::setAtol(double atol)
{
}
//===================================================================================================================================
int  GFCEO_Electrode::evalTimeTrackingEqns(const double t, const double delta_t, const double* const y, const double* const ydot,
                                           const double* const p)
{
    return 0;
}
//===================================================================================================================================
bool GFCEO_Electrode::evalStoppingCritera(const double t, const double delta_t,
                                          const double* const y, const double* const ydot)
{
    return 0;
}

//===================================================================================================================================
int GFCEO_Electrode::calcDeltaSolnVariables(const double t, const double* const y,
                           const double* const ydot,
                           double* const delta_y,
                           const double* const solnWeights)
{
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::calcColumnScales(const double* const y, const double* const y_old, double* const yColScales)
{
    
}
//===================================================================================================================================
void GFCEO_Electrode::user_out2(const int ifunc, const double t,
				const double delta_t, int time_step_num,
				const double* const y,
				const double* const ydot)
{
}
//===================================================================================================================================
int GFCEO_Electrode::matrixConditioning(double* const matrix, const int nrows,double* const rhs)
{
    return 0;
}
//===================================================================================================================================
int GFCEO_Electrode::evalJacobian(const double t, const double delta_t, double cj,
				  const double* const y, const double* const ydot,
				  GeneralMatrix& J, double* const resid, const Solve_Type solveType)
{
    return 0;
}
//===================================================================================================================================
int GFCEO_Electrode::evalJacobianDP(const double t, const double delta_t, double cj,
				    const double* const y,
				    const double* const ydot,
				    double* const* jacobianColPts,
				    double* const resid, const Solve_Type solveType)
{
    return 0;
}
//===================================================================================================================================
void GFCEO_Electrode::writeSolution(int ievent, const bool doTimeDependentResid, const double time, const double deltaT, 
                                    const int time_step_num, const double* const y, const double* const ydot,
                                    const Solve_Type solveType, const double delta_t_np1)
{

}
//===================================================================================================================================
void GFCEO_Electrode::assertOwnership()
{
     if (! iOwnObject_) {
        Electrode_Integrator* e = (Electrode_Integrator*) ee_->duplMyselfAsElectrode();
        ee_ = e;
        iOwnObject_ = 1;
    } 
}
//===================================================================================================================================
} // End of namespace
//-----------------------------------------------------------------------------------------------------------------------------------

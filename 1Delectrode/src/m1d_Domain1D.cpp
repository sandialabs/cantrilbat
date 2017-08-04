/*
 * m1d_BulkDomain1D.cpp
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#include "m1d_Domain1D.h"

#include "m1d_DomainDescription.h"
#include "m1d_exception.h"

//#include "cantera/base/ctml.h"
#include "m1d_globals.h"
#include "m1d_ProblemStatement.h"

#include "mdp_stringUtils.h"


#include <cstdio>
#include <cstdlib>

using namespace std;
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
// Constructor
Domain1D::Domain1D() :
    NumDomainEqns(0),
    coordinateSystemType_(Rectilinear_Coordinates),
    crossSectionalArea_(1.0),
    cylinderLength_(1.0),
    TemperatureReference_(298.15),
    PressureReference_(1.01325e5),
    energyEquationProbType_(0),
    solidMechanicsProbType_(0),
    porosityEquationProbType_(Porosity_EqnType_Status::None),
    residType_Curr_(Zuzax::ResidEval_Type::Base_ResidEval),
    counterResBaseCalcs_(0),
    counterJacBaseCalcs_(0),
    counterJacDeltaCalcs_(0),
    counterResShowSolutionCalcs_(0)
{
    coordinateSystemType_ = PSinput_ptr->coordinateSystemType_;
    crossSectionalArea_ = PSinput_ptr->crossSectionalArea_;
    cylinderLength_ = PSinput_ptr->cylinderLength_;
    TemperatureReference_ = PSinput_ptr->TemperatureReference_;
    PressureReference_ = PSinput_ptr->PressureReference_;
    energyEquationProbType_ = PSinput_ptr->Energy_equation_prob_type_;
#ifdef MECH_MODEL
    solidMechanicsProbType_ = PSinput_ptr->Solid_Mechanics_prob_type_;
#else
    if (PSinput_ptr->Solid_Mechanics_prob_type_ != 0) {
        throw m1d_Error("Domain1D::Domain1D()", "MECH_MODEL define not set, but Solid_Mechanics_prob_type_ nonzero, "
                        + mdpUtil::int2str(PSinput_ptr->Solid_Mechanics_prob_type_));
    }
#endif
}
//==================================================================================================================================
Domain1D::Domain1D(const Domain1D& r) :
    NumDomainEqns(0),
    coordinateSystemType_(Rectilinear_Coordinates),
    crossSectionalArea_(1.0),
    cylinderLength_(1.0),
    TemperatureReference_(298.15),
    PressureReference_(1.01325e5),
    energyEquationProbType_(0),
    solidMechanicsProbType_(0),
    residType_Curr_(Zuzax::ResidEval_Type::Base_ResidEval),
    counterResBaseCalcs_(0),
    counterJacBaseCalcs_(0),
    counterJacDeltaCalcs_(0),
    counterResShowSolutionCalcs_(0)
{
    *this = r;
}
//==================================================================================================================================
Domain1D::~Domain1D()
{
}
//==================================================================================================================================
Domain1D&
Domain1D::operator=(const Domain1D& r)
{
    if (this == &r) {
        return *this;
    }

    NumDomainEqns                   = r.NumDomainEqns;
    coordinateSystemType_           = r.coordinateSystemType_;
    crossSectionalArea_             = r.crossSectionalArea_;
    cylinderLength_                 = r.cylinderLength_;
    residType_Curr_ = r.residType_Curr_;
    counterResBaseCalcs_ = r.counterResBaseCalcs_;
    counterJacBaseCalcs_ = r.counterJacBaseCalcs_;
    counterJacDeltaCalcs_ = r.counterJacDeltaCalcs_;
    counterResShowSolutionCalcs_ = r.counterResShowSolutionCalcs_;
    TemperatureReference_ = r.TemperatureReference_;
    PressureReference_ = r.PressureReference_;
    energyEquationProbType_ = r.energyEquationProbType_;
    solidMechanicsProbType_ = r.solidMechanicsProbType_;
    return *this;
}
//==================================================================================================================================
void Domain1D::setID(const std::string& s)
{
    m_id = s;
}
//==================================================================================================================================
std::string Domain1D::id() const
{
    err("id()");
    return std::string("");
}
//==================================================================================================================================
void Domain1D::domain_prep(LocalNodeIndices* const li_ptr)
{
    err("domain_prep()");
}
//==================================================================================================================================
void Domain1D::residEval(Epetra_Vector& res, const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                         const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr, const double t,
                         const double rdelta_t, const Zuzax::ResidEval_Type residType, const Zuzax::Solve_Type solveType)
{
    residType_Curr_ = residType;
    err("residEval()");
}
//==================================================================================================================================
void 
Domain1D::eval_PostSoln(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                        const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr, const double t,
                        const double rdelta_t)
{
}
//==================================================================================================================================
void 
Domain1D::eval_HeatBalance(const int ifunc, const double t, const double deltaT, const Epetra_Vector* soln_ptr,
                           const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr, 
                           struct globalHeatBalVals& dVals)
{
}
//==================================================================================================================================
void Domain1D::eval_SpeciesElemBalance(const int ifunc, const double t, const double deltaT, const Epetra_Vector* const soln_ptr,
                                       const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr, 
                                       struct globalHeatBalVals& dVals)
{
}
//==================================================================================================================================
void
Domain1D::residEval_PreCalc(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                            const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
                            const double t, const double rdelta_t, const Zuzax::ResidEval_Type residType,
                            const Zuzax::Solve_Type solveType)
{
}
//==================================================================================================================================
void
Domain1D::residEval_PostCalc(Epetra_Vector& res, const bool doTimeDependentResid,
                             const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
                             const Epetra_Vector* const solnOld_ptr, const double t, const double rdelta_t,
                             const Zuzax::ResidEval_Type residType, const Zuzax::Solve_Type solveType)
{
}
//==================================================================================================================================
void
Domain1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                              const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
                              const double t, const double t_old)
{
}
//==================================================================================================================================
void
Domain1D::revertToInitialGlobalTime()
{
}
//==================================================================================================================================
void
Domain1D::saveDomain(ZZCantera::XML_Node& oNode, const Epetra_Vector* const soln_GLALL_ptr,
                     const Epetra_Vector* const solnDot_GLALL_ptr, const double t, bool duplicateOnAllProcs)
{
    err("saveDomain()");
}
//==================================================================================================================================
void
Domain1D::readSimulation(const ZZCantera::XML_Node& simulationNode, Epetra_Vector* const soln_GLALL_ptr,
                         Epetra_Vector* const solnDot_GLALL_ptr)
{
    std::string ida = id();
    ZZCantera::XML_Node* domainNode = simulationNode.findNameID("domain", ida);
    if (!domainNode) {
        throw m1d_Error("Domain1D::readSimulation()", "cant find domain node" + ida);
    }
    readDomain(*domainNode, soln_GLALL_ptr, solnDot_GLALL_ptr, -1.0);
}
//==================================================================================================================================
void
Domain1D::readDomain(const ZZCantera::XML_Node& oNode, Epetra_Vector* const soln_GLALL_ptr,
                     Epetra_Vector* const solnDot_GLALL_ptr, double globalTimeRead)
{
    err("readDomain()");
}
//==================================================================================================================================
void Domain1D::writeSolutionTecplotHeader()
{
    err("writeSolutionTecplotHeader()");
}
//==================================================================================================================================
void Domain1D::writeSolutionTecplot(const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
                                    const double t)
{
    err("writeSolutionTecplot()");
}
//==================================================================================================================================
void
Domain1D::showSolution(const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
                       const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
                       const Epetra_Vector* const solnOld_ptr, const Epetra_Vector_Owned* const residInternal_ptr,
                       const double t, const double rdelta_t, const int indentSpaces, bool duplicateOnAllProcs)
{
    err("showSolution()");
}
//==================================================================================================================================
void
Domain1D::showSolutionVector(std::string& solnVecName, const Epetra_Vector* const solnVector_GlAll_ptr,
                             const Epetra_Vector* const solnVector_ptr, const double t, const double rdelta_t,
                             int indentSpaces, bool duplicateOnAllProcs, FILE* of)
{
    err("showSolutionVector()");
}
//==================================================================================================================================
void
Domain1D::showSolutionIntVector(std::string& solnVecName, const Epetra_IntVector* const solnIntVector_GlAll_ptr,
                                const Epetra_IntVector* const solnIntVector_ptr, const double t, const double rdelta_t,
                                int indentSpaces, bool duplicateOnAllProcs, FILE* of)
{
    err("showSolutionIntVector()");
}
//==================================================================================================================================
int
Domain1D::reportSolutionParam(const std::string& paramID, double* const paramVal) const
{
    paramVal[0] = 0.0;
    return -1;
}
//==================================================================================================================================
int
Domain1D::reportSolutionVector(const std::string& requestID, const int requestType, const Epetra_Vector* soln_ptr,
                               std::vector<double>& vecInfo) const
{
    vecInfo.clear();
    return -1;
}
//==================================================================================================================================
void
Domain1D::setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted* const soln,
                               const Epetra_Vector_Ghosted* const solnDot, const double t, const double delta_t, const double t_old)
{
}
//==================================================================================================================================
void
Domain1D::initialConditions(const bool doTimeDependentResid, Epetra_Vector* const soln,
                            Epetra_Vector* const solnDot, const double t, const double delta_t)
{
}
//==================================================================================================================================
void Domain1D:: calcDeltaSolnVariables(const double t, const Epetra_Vector& soln,
                                       const Epetra_Vector* const solnDot_ptr, Epetra_Vector& deltaSoln,
                                       const Epetra_Vector* const atolVector_ptr,
                                       const Zuzax::Solve_Type solveType,
                                       const  Epetra_Vector* const solnWeights)
{
    printf("Domain1D: function not implemented\n");
    exit(-1);
}
//==================================================================================================================================
void Domain1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
                             Epetra_Vector_Ghosted& atolVector, const Epetra_Vector_Ghosted* const atolV)
{
    printf("Domain1D: function not implemented\n");
    exit(-1);
}
//==================================================================================================================================
void
Domain1D::setAtolDeltaDamping(double atolDefault, double relcoeff,  const Epetra_Vector_Ghosted& soln,
                              Epetra_Vector_Ghosted& atolDeltaDamping, const Epetra_Vector_Ghosted* const atolV)
{
    printf("Domain1D: function not implemented\n");
    exit(-1);
}
//==================================================================================================================================
void
Domain1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted& soln,
                                const Epetra_Vector_Ghosted& solnDot, Epetra_Vector_Ghosted& atolVector_DAEInit,
                                const Epetra_Vector_Ghosted* const atolV)
{
    printf("Domain1D: function not implemented\n");
    exit(-1);
}
//==================================================================================================================================
void
Domain1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, const Epetra_Vector_Ghosted& soln,
                                      const Epetra_Vector_Ghosted& solnDot, Epetra_Vector_Ghosted& atolDeltaDamping,
                                      const Epetra_Vector_Ghosted* const atolV)
{
    printf("Domain1D: function not implemented\n");
    exit(-1);
}
//==================================================================================================================================
void
Domain1D::incrementCounters(const Zuzax::ResidEval_Type residType)
{
    if (residType == Zuzax::ResidEval_Type::Base_ResidEval) {
        counterResBaseCalcs_++;
    } else if (residType == Zuzax::ResidEval_Type::JacBase_ResidEval) {
        counterJacBaseCalcs_++;
    } else if (residType == Zuzax::ResidEval_Type::Base_ShowSolution) {
        counterResShowSolutionCalcs_++;
    } else {
        counterJacDeltaCalcs_++;
    }
}
//==================================================================================================================================
void
Domain1D::err(const char* msg) const
{
    printf("BulkDomain1D: function not implemented: %s\n", msg);
    exit(-1);
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------


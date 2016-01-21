
Work on structuring the Globally Fully Coupled Electrode Object


The Electrode_Integrator uses the multiple inheritance class construct to
make all Electrode objects an inheritor from Cantera::ResidJacEval. 
Therefore that method is not available for the GFCEO problem.

 Electrode_CSTR::evalResidNJ
       Electrode_CSTR::integrateResid()

             unpackNonlinSolnVector(y);
             calcResid(resid, evalType);


The way the integrator is set up, the nonlinear solver refers back to itself
using
        
     pSolve_ = new NonlinearSolver(this);

where pSolve_ is a pointer to a Cantera::NonlinearSolver object.

The new object doesn't necessarily have to be a Cantera::NonlinearSolver
object since we aren't solving the nonlinear problem separately.
However, it probably would save time if it were inherited from
a  Cantera::ResidJacEval object.

SOLUTION 1
------------------------------

One solution is to put the necessary classes as regular virtual function
in Electrode with the prefix GFCEO (globally fully coupled electrode object
    GFCEO_evalresisNJ()
    GFCEO_nEquations()

Then, have an interface class that inherits from Cantera::ResidJacEval 
that takes a pointer to the Electrode object.
Then, we could provide a default interface that attaches the 
new member functions to the Cantera::ResidJacEval member functions.

This may involve the smallest amount of work.




PROBLEM - 

The calcResid() problem is meant to handle discontinuities in functions wrt deltaT by 
breaking up the time step. How is this done in the new problem? Could we add 
significant concepts from Tyler?















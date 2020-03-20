**************************************************************************************************
Work on structuring the Globally Fully Coupled Electrode Object

**************************************************************************************************

The Electrode_Integrator uses the multiple inheritance class construct to
make all Electrode objects an inheritor from Zuzax::ResidJacEval. 
Therefore that method is not available for the GFCEO problem.

 Electrode_CSTR::evalResidNJ
       Electrode_CSTR::integrateResid()

             unpackNonlinSolnVector(y);
             calcResid(resid, evalType);


The way the integrator is set up, the nonlinear solver refers back to itself
using
        
     pSolve_ = new NonlinearSolver(this);

where pSolve_ is a pointer to a Zuzax::NonlinearSolver object.

The new object doesn't necessarily have to be a Zuzax::NonlinearSolver
object since we aren't solving the nonlinear problem separately.
However, it probably would save time if it were inherited from
a  Zuzax::ResidJacEval object.

SOLUTION 1
------------------------------

One solution is to put the necessary classes as regular virtual function
in Electrode with the prefix GFCEO (globally fully coupled electrode object
    GFCEO_evalresisNJ()
    GFCEO_nEquations()

Then, have an interface class that inherits from Zuzax::ResidJacEval 
that takes a pointer to the Electrode object.
Then, we could provide a default interface that attaches the 
new member functions to the Zuzax::ResidJacEval member functions.

This may involve the smallest amount of work.



PROBLEM - 

The calcResid() problem is meant to handle discontinuities in functions wrt deltaT by 
breaking up the time step. How is this done in the new problem? Could we add 
significant concepts from Tyler? If we add Tyler's solution then the solution becomes
smooth. But, we must follow the smoothness if it's below the atol requirements.

This discontinuities problem is big, because the newman equations aren't typically used to
solve this equation system. It may be prudent to make everything continuous.

The equation system to use for multispecies equations is up in the air too.
I could go with a straight mole number equation system. This would be a departure
from the one before, where I went with phase moles / mole fraction of sp >= 1 system.

A straight mole number system has the benefit that all equations go below atol as the phase 
becomes unimportant and small.

We could also see if we can utilize the phase death problem to in a subcycle iteration scheme.


SOLUTION METHOD FOR THE RESULTING TIME EQUATIONS

GFCEO -> To solve these systems, we can fall back on the general DAE solution methods supplied by
         Zuzax using the DAE_Solver class


HANDLING ZEROING PHASES:

Previously I had an extra equation representing the calculate deltaT. I can't maintain this when
coupling everything together, so I need a different solution.

The solution is to use the sundials solve for a root solution using the following equation system:

First Equation:
     dn / dt = sum_j( src_j ) = src^T

     
    n d X_i / dt  = src_i - X_i ( src^T)

    for i = 1, ..., N-1

 
    Note: the i = 0 mole fraction equation is left out, or sum_j X_j = 1 is substituted for that equation.
     (I like the "left out option".

  
This has the advantage that the X_i is exactly correct when n goes through the origin.
This should be stable for n positive and negative. We need to ensure the bounded-ness of X_i for i = 0, ... N-1.

Also from the get-go, we need to ensure that the special species is the species with the largest 
mole fraction. It always occurs that the scheme breaks down when this isn't taken into account.





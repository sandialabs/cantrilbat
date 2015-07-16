
Example Structure of the Electrode XML solution file:


<ctml>
  <electrodeOutput index="1">
    <globalTimeStep index="1">
      <timeIncrement type="global">
        <timeState type="t_init">
          <time> 0.0 </time>
          <electrodeState>



What it All Means
------------------------------

electrodeOutput XML Element

The outer XML element called "electrodeOutput" 
contains all of the data for a single time dependent run of a single
explicit electrode object.  It has an index attribute which is a
numbered entry. If there are more than one electrodeOutput in a file,
then the indexes will be consequuatively numbered.

ElectrodeIdentification XML Element

This element contains the information needed to identify the type of electrode
used and the type of EState file written. 


globalTimeStep  XML Element

Everytime you do a resetStartingCondition(Tinitial) you create a new
global time step and a new globalTimeStep XML element. This occurs for
all global time steps in a global run. The "globalTimeStep" XML element
has an index attribute as well. These will be consequatively numbered starting
with the value 1, usually though restarts may change this value. 
It corresponds to the global time step number in the encompassing calculation. 
The globalTimeStep XML element consists of a series of timeIncrement XML nodes
that contain the subGlobal integration steps. It includes the value that will
be used for the next subGlobal deltaT as well. This is needed to prevent
churning of the calculation.

timeIncrement XML Element

A timeIncrement consists of the results from a single run of the integrate() command
within the Electrode object. Because the data is perhaps discontinuous, the
value of t_init is included in the data. States within the calculation are
labeled as t_init, t_intermediate, and t_final.
Note, if we are formulating a Jacobian and want to dump all of the
data into the solution file, then there will be more than one
timeIncrement XML elements for a single globalTimeStep. This has yet to be
implemented. Also, if we
are trying to find a current and using the root finder, then there may
be more than one timeIncrement XML element for a single
globalTimeStep. Each of these solutions will have different external
solution variables applied to the electrode over the same solution interval.

timeState XML Element

A timeState XML element contains the cellNumber, domain, and type attribute.
The type attribute can be t_init, t_intermediate, and t_final currently.
The timeState XML node contains the time and the electrodeState elements.
The time is the current solution time.
The electrodeState is the current internal state of the electrode. 

electrodeState XML element

The electrodeState is the current internal state of the electrode. It also
includes the values of the external variables.  This
includes the temperature, pressure, and voltages of all phases. The species
numbers are also included. Thus, the external variables are included in the
specification of the the electrodeState, and the voltage of the electrolyte
and the molar composition of the electrode are included in the calculation.

Several history variables are included. The "capacityDischargedtoDate" is
included, which is a historical value for the amount of electrons
created/discharged during the calculation.

The "deltaTsubcycle_init_next" value contains the value to use on the
subintegration time step for the next step.

The electrodeState element varies with the type of electrode.  The
ElectrodeIdentification XML Element contains the information needed to
identify the type of file written.




@TODO

Add external solution values to the XML database.

Add the capability of skipping t_init entries when they are continuous. This
may just involve a pointer to the previous entry which contains the results.








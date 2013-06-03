


<ctml>
  <electrodeOutput index="1">
    <globalTimeStep index="1">
      <timeIncrement type="global">
        <timeState type="t_init">
          <time> 0 </time>
          <ElectrodeState>




What it All Means
------------------------------


The outer XML element called electrodeOutput is a containing that
contains all of the data for a single time dependent run of a single
explicit electrode object.  It has an index attribute which is a
numbered entry. If there are more than one electrodeOutput in a file,
then the indexes will be consequuatively numbered.




Everytime you do a resetStartingCondition(Tinitial) you create a new
global time step and a new globalTimeStep XML element.


A timeIncrement consists of the results from a  single run of the integrate() command
within the Electrode object.

Note, if we are formulating a Jacobian and want to dump all of the
data into the solution file, then there will be more than one
timeIncrement XML elements for a single globalTimeStep. Also, if we
are trying to find a current and using the root finder, then there may
be more than one timeIncrement XML element for a single
globalTimeStep. Each of these solutions will have different external
solution variables applied to the electrode over the same solution interval.






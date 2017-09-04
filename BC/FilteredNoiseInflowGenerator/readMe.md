# Update by Lukas Lebovitz (08/2017)

- This variant of the filtered noise inflow generator uses the proposed method of selective filtering. Selective filtering only filters the virtual grid coordinates which are nearest to actual mesh faces. This can heavily improve performance for grids that are non-uniform (i.e. meshes with refinements). 

- This version fixes a bug that crashes the simulation sometimes when the number of indices on the virtual grid can't be split equally to all processors.

- The `getRandomField()` function now works in parallel and distributes random number generation on all available processors.

- This version fixes a bug that crashes the simulation sometimes when the number of indices on the virtual grid can't be split equally to all processors.

- Furthermore the indices are now distributed to the processors according to their computational load which is depending on the corresponding integral length scales. This only works if the integral length scales are vertically increasing (in z-direction).

- Some minor modifications.

### Instructions:

In order to use this inlet generator with OpenFOAM:

0. Compile this boundary condition with the command wmake.

1. Include this in your controldict:
  libs
  (
      "libfilteredNoiseInflowGenerator.so"
  );

2. Call the inlet generator from your velocity boundary condition for your inlet:
  
```cpp
    <nameOfInletPatch>
    {
        type            filteredNoiseInflowGenerator;
        value           uniform ( 0 0 0 ); // placeholder
        //perturb 1e-6; //optional: used for the interpolation. Change this value if there are artifacts on the inlet patch
        correlationShape gaussian; // other correlationShapes are exp or doubleExp
        //Virtual grid has the size of the smalles inlet patch cell.
        gridFactor 2; //optional: Use the gridFactor to scale the virtual grid size.
    }
```

3. Make sure you include a directory `constant/boundaryData/<nameOfInletPatch>/` which contains the files for `R`, `L` and `ref`. An example case is found at: https://github.com/ETH-BuildingPhysics/ETH-OFTools-2.3.X/tree/master/Tutorials/FilteredNoiseInflowGenerator

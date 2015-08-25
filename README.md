# ETH-OFTools-2.3.X

Chair of Building Physics (ETH Zurich) addons for OpenFOAM-2.3.x.


## Version:


## Authors:
* Marc Immer
* Marcel Vonlanthen


## Installation:
Exectute the following lines in a terminal to install all the addons. Each
addon can be installed separately by executing wmake in its folder.
1. git clone https://github.com/ETH-BuildingPhysics/ETH-OFTools-2.3.X.git
2. cd ETH-OFTools-2.3.X
3. git checkout -b development remotes/origin/development
4. ./Allwmake

## Contents

### Boundary Conditions
* **FilteredNoiseInflowGenerator**     
* **FilteredNoiseInflowGeneratorScalar**     
* **LEMOSinflowGeneratorMod**     
* **nestedBlendedVelocity**   
  Nested boundary condition for the velocity. To use with a nested solver.

### Function Objects
* **CnC**  
* **scalarVelocityProduct**  
* **timeAbort**  

### Python tools

### Sampling libraries
* **foamFileSmall**  
* **hdf5**  

### Solvers
* **nestedImpBlendDsPimpleFoam**  
  One-way, downscale, nested solver with an implicit blending at the nested boundaries. It solves the Navier-Stokes equation with the pimple algorithm. The turbulence model can be laminar, RANS or LES. The solver has nesting cascade capabilities.
* **pisoScalarSourceFoam230**     
  The OpenFOAM pisoFoam solver with passive scalar solver included.


### Utilities
* **cellVolumes**    
* **gradU**    
* **LESdelta**    
* **postAvg**    
* **Q2D**    
* **resLES**  
* **scalarCovariance**  
* **scalarFluxes**  
* **SGSValues**  
* **snGradU**  
* **TKEProd2D**  
* **yPlusMeanLES**      
  Compute yPlus mean for a LES. The fields UMean, nuSgsMean must exist. They can be created with the default OpenFOAM functionObject *fieldAverage*
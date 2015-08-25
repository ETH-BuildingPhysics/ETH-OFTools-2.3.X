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
   Variant of Kleins Filtered Noise turbulent inflow boundary condition
* **FilteredNoiseInflowGeneratorScalar**  
  Variant *FilteredNoiseInflowGenerator* including turbulent scalar generation
* **LEMOSinflowGeneratorMod**     
   Modification of https://github.com/LEMOS-Rostock/LEMOS-2.3.x/tree/master/libLEMOS-2.3.x/boundaryConditions/inflowGenerator to read input files similar to the OpenFOAM BC *TimeVaryingMappedFixedValue*
* **nestedBlendedVelocity**   
  Nested boundary condition for the velocity. To use with a nested solver.

### Function Objects
* **CnC**  
   Job control function object
* **scalarVelocityProduct**  
   Creates the field Us, where s is any scalar field. Must be placed before the OpenFOAM function object *fieldAverage*, to get the averaged UsMean field. Useful to get the fields required for *scalarCovariance*
* **timeAbort**  
   Job control function object

### Python tools

### Sampling libraries
* **foamFileSmall**  
    Similar to the file format *foamFile*, however does not save the points and faces. Useful for runtime sampling of surfaces if mesh does not change.
* **hdf5**  
    New surface file type: hdf5. Saves a surface as a h5 file. Needs the hdf5 library installed

### Solvers
* **nestedImpBlendDsPimpleFoam**  
  One-way, downscale, nested solver with an implicit blending at the nested boundaries. It solves the Navier-Stokes equation with the pimple algorithm. The turbulence model can be laminar, RANS or LES. The solver has nesting cascade capabilities.
* **pisoScalarSourceFoam230**     
  The OpenFOAM pisoFoam solver with passive scalar transport of an arbitrary number of scalars. The list of scalars to be solved is defined in the transportProperties as: scalars ( A B C ); a generic boundary condition can be defined with the file 0/sdefault. Use fvOptions do define sources. Also writes the product field Us, to be used with runtime averaging and the utility *scalarCovariance*


### Utilities
* **cellVolumes**    
  Save the cell volumes into a file *cellVolumes*
* **gradU**    
  saves the fields *gradU* and *gradUMean* of the internal field (not at boundaries)
* **LESdelta**    
  saves the LES filter width sgsModel->delta() in the field *LESdelta*
* **postAvg**    
  Averages scalar fields, saves the averaged field with the postfix *Mean*. E.g *nuSgs* will be saved as *nuSgsMean*
* **Q2D**    
  Similar to *Q*, but sets the off-plane components of the velocity gradient tensor to zero (off-plane is z component). Emulates what a Stereo-PIV system would measure.
* **resLES**  
  Useful for analyzing statistics of LES computations. Needs the fields UMean,UPrime2Mean and kMean. Outputs TKE in *TKEMean*, the ratio of resolved TKE to total TKE in *resLES* and the production of resolved TKE in *TKEMeanProd*
* **scalarCovariance**  
  To compute the scalar covariance, e.g the turbulent temperatur flux <u'T'> is written to *UPrimeTPrimeMean*, calculated by UTMean - (UMean * TMean). Important: the product field UT must be sampled and averaged!
* **scalarFluxes**  
  Used in combination with the ** solver. Not needed with proper runtime averaging and the *scalarCovariance* utility
* **SGSValues**  
  For LES. Saves the production of subgrid TKE in a field *kSgsProd*. Computed using 2.0*nuSgs*magSqr(symm(gradU)));
* **snGradU**  
* **TKEProd2D**  
* **yPlusMeanLES**      
  Compute yPlus mean for a LES. The fields UMean, nuSgsMean must exist. They can be created with the default OpenFOAM functionObject *fieldAverage*. If nuSgsMean does not exist, it can be created with *postAvg*
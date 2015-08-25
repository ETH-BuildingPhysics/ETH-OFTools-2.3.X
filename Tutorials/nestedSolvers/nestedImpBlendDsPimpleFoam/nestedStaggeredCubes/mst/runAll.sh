#!/bin/sh


masterDir=$(pwd)

echo 'Prepare the master case'
echo '======================='
echo ' '
./generateSnappyGrid.sh
decomposePar > output/decomposePar.log 2>&1

echo ' '
echo 'Prepare the nested case'
echo '======================='
echo ' '
cd ../ntd
nestedDir=$(pwd)
./generateSnappyGrid.sh
decomposePar > $nestedDir/output/decomposePar.log 2>&1

echo ' '
echo 'Link mst and ntd'
echo '================'
echo ' '
cd $masterDir
./foamMakeNestedSymLink.py -nd ../ntd -ts 0

echo ' '
echo 'run the case (long!)'
echo '===================='
echo ' '
cd $masterDir
mpirun -np 4 nestedImpBlendDsPimpleFoam -parallel




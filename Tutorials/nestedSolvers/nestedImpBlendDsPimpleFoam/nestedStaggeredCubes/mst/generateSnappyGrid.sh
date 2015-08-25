#!/bin/sh

mkdir -p output

echo 'clean polyMesh directory'
find constant/polyMesh -type f -not -name 'blockMeshDict*' | xargs rm || true


echo 'run surfaceFeatureExtract'
surfaceFeatureExtract > output/surfaceFeatureExtract.log 2>&1
echo 'run blockMesh'
blockMesh > output/blockMesh.log 2>&1
echo 'run snappyHexMesh'
snappyHexMesh -overwrite  > output/snappyHexMesh.log 2>&1
echo 'run checkMesh'
checkMesh  > output/checkMesh.log 2>&1

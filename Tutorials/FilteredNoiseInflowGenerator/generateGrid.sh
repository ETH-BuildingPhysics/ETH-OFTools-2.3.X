#!/bin/sh

mkdir -p output

echo '      create blockMesh'
blockMesh -dict constant/polyMesh/blockMeshDict_withMultiHP > output/blockMesh.out 2>&1

echo '      create refinement gndboxLevel01'
topoSet -dict system/topoSetDict_gndboxLevel01 > output/topoSet_gndBoxLevel01.out 2>&1
refineMesh -dict system/refineMeshDict_gndboxLevel01 -overwrite > output/refineMesh_gndBoxLevel01.out 2>&1

echo '      create refinement gndboxLevel02'
topoSet -dict system/topoSetDict_gndboxLevel02 > output/topoSet_gndBoxLevel02.out 2>&1
refineMesh -dict system/refineMeshDict_gndboxLevel02 -overwrite > output/refineMesh_gndBoxLevel02.out 2>&1

echo '      create refinement gndboxLevel03'
topoSet -dict system/topoSetDict_gndboxLevel03 > output/topoSet_gndBoxLevel03.out 2>&1
refineMesh -dict system/refineMeshDict_gndboxLevel03 -overwrite > output/refineMesh_gndBoxLevel03.out 2>&1

# echo '      create refinement gndboxLevel04'
# topoSet -dict system/topoSetDict_gndboxLevel04 > output/topoSet_gndBoxLevel04.out 2>&1
# refineMesh -dict system/refineMeshDict_gndboxLevel04 -overwrite > output/refineMesh_gndBoxLevel04.out 2>&1

echo '      checkMesh'
checkMesh > output/checkMesh.out 2>&1
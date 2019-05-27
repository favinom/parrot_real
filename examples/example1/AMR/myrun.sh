#!/bin/bash

ol=1;
as=0;
np=6;
mpirun -n ${np} ../../../parrot_real-opt -i refineBlock.i origLevel=${ol} adapSteps=${as}

for (( c=0; c<=as+1; c++ ))
do
rm refinedMesh_${ol}_000$c.xdr
done

for (( c=1; c<=as+1; c++ ))
do
temp=`expr $c - 1`
mv refinedMesh_${ol}_000${c}_mesh.xdr refinedMesh_${ol}_000${temp}_mesh.xdr
done

mpirun -n ${np} ../../../parrot_real-opt -i diffusion.i origLevel=${ol} adapSteps=${as}

mpirun -n ${np} ../../../parrot_real-opt -i advection.i origLevel=${ol} adapSteps=${as}

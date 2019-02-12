This data repository accompanies the paper
"Verification benchmarks for single-phase flow in three-dimensional fractured porous media",
see the corresponding [call for participation](https://arxiv.org/abs/1809.06926).

Four benchmark cases are considered in the following subfolders:

* [single](single):  Refers to Section "4.1 Case 1: single fracture" of the document.

* [regular](regular): Refers to Section "4.2 Case 2: regular fracture network" of the document.

* [small_features](small_features):  Refers to Section "4.3 Case 3: network with small features" of the document.

* [field](field):  Refers to Section "4.4 Case 4: a field case" of the document.

For each benchmark, geometry and result data as well as plotting scripts are provided in corresponding subfolders:

* `geometry`: a geometrical description of the matrix domain and fracture network
in form of a Gmsh input file `gmsh.geo`.

* `results`: usually `csv` files as required by the corresponding case.

* `scripts`: Python and/or Matlab scripts to process the result data.

Explicit instructions on how to contribute are provided in the [call for participation](https://arxiv.org/abs/1809.06926).
# Solving electromechanics of ICC using CHASTE

This code solves electroMechanics for a 2D ICC block:

- electrics solver generates Ca and voltage based on cell model (from CellML)

- contraction model generates active tension based on Ca or voltage

- mechanics solver generates deformation based on the tension and the material law

# How to run:

First install Chaste from:

https://github.com/Chaste/Chaste

Or use the docker version:

https://github.com/Chaste/chaste-docker

Basically, for Docker, you only need to type this in the terminal:

```
sudo docker run --name chaste -it --init --rm -v chaste_data:/home/chaste -v $(pwd)/testoutput:/home/chaste/testoutput chaste/release
```

### Steps:

- clone this repo into the Chaste projects folder.

- copy /fibre/2by3_fibre_s.ortho to the Chaste "testoutput" folder or change the path in line 155 of icc2dsim_Du.hpp.

- if you want to use Nash1998 contraction model (this is the default in the main test, so you **must** do this step), copy /utils/NashContractionModel.cpp (and .hpp) files to /Chaste/heart/src/odes/contractionmodels/ and change the constants if necessary:

```
sudo docker cp /utils/NashContractionModel.cpp chaste:/home/chaste/src/heart/src/odes/contractionmodels/
sudo docker cp /utils/NashContractionModel.hpp chaste:/home/chaste/src/heart/src/odes/contractionmodels/
```
also add **Nash** to line 46 of /Chaste/heart/src/odes/contractionmodels/ContractionModelName.hpp

and add the following to line 44 of Chaste/heart/src/problem/cell_factories/LabelBasedContractionCellFactory.hpp:

```
#include "NashContractionModel.hpp" 
```
then add the following to line 99 of Chaste/heart/src/problem/cell_factories/LabelBasedContractionCellFactory.hpp:

```
case NASH:
            {
                return new NashContractionModel;
                break;
            }

```

- if you want to use an exponential material law, copy /utils/SchmidCostaExponentialLaw2d.cpp (and .hpp) files to /Chaste/continuum_mechanics/src/problem/material_laws/ and change the constants if necessary. This material law is not stable and most of the time the solver does not converge! Consider using Mooney-Rivlin material law (it's the default in the main test code (line 148, 163) at the moment):

```
sudo docker cp /utils/SchmidCostaExponentialLaw2d.cpp chaste:/home/chaste/src/continuum_mechanics/src/problem/material_laws/
sudo docker cp /utils/SchmidCostaExponentialLaw2d.hpp chaste:/home/chaste/src/continuum_mechanics/src/problem/material_laws/
```

- run the following in the terminal (you need to do this only once):

```
  cmake ICC_2D 8 -DCMAKE_BUILD_TYPE:STRING="Release"\
          -DChaste_ERROR_ON_WARNING:BOOL="OFF" \
          -DChaste_UPDATE_PROVENANCE:BOOL="OFF" \
          -H/home/chaste/src \
          -B/home/chaste/lib
```

- go to /home/chaste/lib/projects/ICC_2D/ and run:

```
  make -j8
```

- go to the /test folder and run:

```
  ./icc2dsim_Du 
```

- you can check the results in the testoutput folder and visualise with cmgui
  
  /testoutput/icc2d_Du/deformation/cmgui/LoadSolutions.com
  
- voltage and Ca of the watched node can be checked in /testoutput/icc2d_Du/watched.txt

  it is noted that Ca is not watched by default, so you need to:
  
  uncomment and change line 157 in /Chaste/heart/src/problem/CardioElectroMechanicsProblem.cpp:
  ```
    double Ca = mpElectricsProblem->GetTissue()->GetCardiacCell(mWatchedElectricsNodeIndex)->GetIntracellularCalciumConcentration();
  ```
  add the following (above line 164):
  
  ```
  *mpWatchedLocationFile << Ca << "\n";
  ```
  
  change (remove "\n"):
  ```
  *mpWatchedLocationFile << V << "\n"; 
  ```
  
  to:
  
  ```
  *mpWatchedLocationFile << V << " ";
  ```

- strain and stress of the watched node can be checked in /testoutput/icc2d_Du/deformation/*.strain and *.stress


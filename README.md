## relxill++ - modified and extended version of the relxill model package

Dan Wilkins - wilkins.401@osu.edu
If using these models, please be sure to cite the original Dauser and Garcia papers detailed through the links below!

The master branch of this repository is the unmodified relxill package, forked from upstream. For the different model variants, please checkout one of the branches

- relxill_emis - modifications of fitting accretion disc emissivity profiles, including twice broken power law (relxill3), and an emissivity profile in radial bins for free-fitting (relxillemis)

- relxill_components - provides a switch parameter to isolate the reflection and primary components from a relxill model for visualisation purposes.

To compile each of the model variants, copy all of the files from src/ to build/, then run the XSPEC initpackage command on the lmodel_relxill.dat file.

## relxill - Astrophysics local model for Relativistic Reflection 

Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP

### 1. General
The model description and download of released version can be found on the homepage at 
http://www.sternwarte.uni-erlangen.de/research/relxill/  -- In order to obtain a stable version this
is the best place to get the model from.

### 2. Installation
In order to install and build the model from the GIT distribution, simply run "make model". 
It will be installed in the subdirectory "build/". After the installation, a simple test is executed. If
this test succeed the model should be correctly working. This installation only works if you have a recent
version of Heasoft installed (https://heasarc.gsfc.nasa.gov/lheasoft/).

### 3. Usage
The relxill code is written in C and C++. While its routines can be used directly, its main intention
is to be used as *local model* in X-ray data analysis software such as ISIS, Xspec, or Sherpa. The 
description and meaning of the model parameters can be found in the relxill documentation (https://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/relxill_docu.pdf)
and on the relxill homepage.
If the default installation does not work for you, you might need to check how to install local models for your
data analysis software.

### 4. Contributing

More information on the code and how to contribute can be found at CONTRIBUTE.md.

### 5. Support
For questions or bug reports, please contact thomas.dauser@fau.de and javier@caltech.edu

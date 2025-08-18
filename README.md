# G4DarkBreM

Geant4 Dark Bremmstrahlung from MadGraph

<p align="center">
  <a href="https://www.apache.org/licenses/LICENSE-2.0" alt="Apache 2.0 license">
    <img src="https://img.shields.io/badge/license-Apache%202-blue" />
  </a>
  <a href="https://github.com/LDMX-Software/G4DarkBreM/actions" alt="Actions">
    <img src="https://img.shields.io/github/actions/workflow/status/LDMX-Software/G4DarkBreM/test.yml?branch=main" />
  </a>
  <a href="https://github.com/LDMX-Software/G4DarkBreM/releases" alt="Releases">
    <img src="https://img.shields.io/github/v/release/LDMX-Software/G4DarkBreM" />
  </a>
</p>

The version submitted in 2022 to Computer Physics Communications is [release v1.1.1](https://github.com/LDMX-Software/G4DarkBreM/releases/tag/v1.1.1).
Ongoing development of this package is maintained at on GitHub [LDMX-Software/G4DarkBreM](https://github.com/LDMX-Software/G4DarkBreM).
A seminar was given to provide more detail on this technique - the slides (with speaker notes) are available on 
[the documentation website](https://ldmx-software.github.io/G4DarkBreM/G4DarkBreM_HEP_Seminar.pdf).

## Installation
The only dependencies of G4DarkBreM are Geant4 which has [an extensive installation guide](http://cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/)
and Boost which can be installed [from the website](https://www.boost.org/doc/libs/1_80_0/more/getting_started/unix-variants.html)
or via your package manager (e.g. [on Ubuntu](https://stackoverflow.com/questions/12578499/how-to-install-boost-on-ubuntu)).

As defined by [our CMake infrastructure](CMakeLists.txt), the minimum versions of these dependencies are 1.68 for Boost and 10.2.3 for Geant4.
We use the Boost.Math and Boost.Iostreams subcomponents of Boost if you wish to limit the size of the Boost needed to be installed.

While G4DarkBreM is not explicitly limited to a certain platform, it has only been used on Linux-based operating systems and explicitly uses C++11 standard features.

After installing Geant4, one can build and install G4DarkBreM using tools probably used to install Geant4 (if built from scratch).
```
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=<my-install>
cd build
make install
```

Additionally, G4DarkBreM can be pulled into a more expansive simulation framework as a submodule and included as a subdirectory in CMake 
```
add_subdirectory(G4DarkBreM)
```
This defines the `G4DarkBreM` cmake target which later targets can link to, for example
```
target_link_libraries(MySim PUBLIC G4DarkBreM)
```

## Usage
Getting started can be done by looking at [the configuration](docs/configuration.md) page
and the g4db::example namespace which holds a fully functional (though primitive) Geant4
simulation example.

### Dark Brem Event Libraries
This simulation relies on using reference event libraries generated before running.
The development and usage of this simulation within LDMX uses event libraries generated with
[LDMX-Software/dark-brem-lib-gen](https://github.com/LDMX-Software/dark-brem-lib-gen)
For the paper published in Computer Physics Communications specifically, `v4.3` of this repository was used
(_note:_ the interaction method with this repository in v4.3 was different than the current method).

Follow the [Quick Start](https://github.com/LDMX-Software/dark-brem-lib-gen#quick-start) in that repository
if you wish to use it for generating libraries.

The directories created by runs of `dark-brem-lib-gen` can already be provided as a dark brem event library;
however, if you wish to shrink the size of the files you are carrying around, you can use
`g4db-extract-library` and `gzip` to save disk space.

First, pull out the library into a CSV file.
```
g4db-extract-library <path-to-db-lib>
```
This already saves a factor of ~10 in disk space since the scaling method doesn't use all
of the kinematic information output by MadGraph/MadEvent. Without explicitly providing
an output file, it will simply call the CSV file the same name as the directory of the 
LHE files with the `.csv` extension added. This CSV file can also be loaded by provided
as a dark brem event library.

If you wish to save another factor of ~4 in disk space, you can also compress the CSV.
```
gzip <path-to-db-lib-csv>
```
This `gzip`-compressed CSV file can also be provided as a dark brem event library.

## Validation
Analysis and validation of G4DarkBreM has been studied in another repository 
[tomeichlersmith/ldmx-sim-technique](https://github.com/tomeichlersmith/ldmx-sim-technique).
This repository also stands as an example for integrating G4DarkBreM into a larger simulation and processing framework.

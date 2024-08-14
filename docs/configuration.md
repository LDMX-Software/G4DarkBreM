# Configuration of Simulation

The G4DarkBreM package follows Geant4's structure of separating
particles from processes and specific models of processes from
the process itself. This means configuring the package
requires three steps.
1. Define the dark photon mass with G4APrime::Initialize
2. Define how the model should behave with g4db::G4DarkBreMModel::G4DarkBreMModel
3. Define which model the process should use with G4DarkBremsstrahlung::G4DarkBremsstrahlung

Besides definition of the dark photon mass, all of the parameters
that change how the dark brem physically behaves are defined within
the model. The process contains some functionality (for example,
limiting the simulation to only one dark brem interaction per event);
however, these parameters are second to the model parameters.

It is common to do both of these steps within a `G4VPhysicsConstructor`
since that class can be registered with a Geant4 PhysicsList and handled
later by Geant4. An example G4VPhysicsConstructor is provded by
g4db::example::APrimePhysics.

## Physics Parameters

There are other parameters that tune the functionality of the simulation,
but the parameters listed here are directly related to how the physics
of the simulation is done.

Parameter             | Required | Location              | Description
----------------------|----------|-----------------------|-------------
A' Mass               | Yes      | G4APrime              | Mass of dark photon in MeV
A' tau                | No       | G4APrime              | Proper lifetime of dark photon in seconds, only relevant if decay mode is set to GeantDecay (default -1, meaning no decay)
Decay Mode            | No       | G4APrime              | Choice of whether to allow the A' to decay and the decay distribution. The only decay channel is e+/e- (default NoDecay)
MG/ME Data            | Yes      | g4db::G4DarkBreMModel | Library of dark brem events
Use Muons             | Yes      | g4db::G4DarkBreMModel | Choice to dark brem off muons or electrons
Minimum DB Threshold  | No       | g4db::G4DarkBreMModel | Minimum energy in GeV lepton needs to have to DB (default twice A' mass)
Mixing Strength       | No       | g4db::G4DarkBreMModel | Mixing strength between dark photon and SM photon (default 1.0)
Scaling Method        | No       | g4db::G4DarkBreMModel | Choice of how to scale the outgoing kinematics from the library (default ForwardOnly)
Xsec Method           | No       | g4db::G4DarkBreMModel | Choice of WW approximation variant to calculate the total cross section. Options include Full, Improved, HyperImproved, and Auto (default is Auto where Full is used if the A' to lepton mass ratio is below the max and Improved is used otherwise)
Maximum mass ratio    | No       | g4db::G4DarkBreMModel | If using Auto Xsec Method, the maximum ratio of A' to lepton mass to transition from Full WW to Improved WW (default 50.0)
A' LHE ID             | No       | g4db::G4DarkBreMModel | If MG/ME Data is in the form of LHE files, the PDG ID of the A' inside of those files (default 622)
Scale A'              | No       | g4db::G4DarkBreMModel | Whether to scale the outgoing dark photon energy along with the recoil lepton (default false, meaning the dark photon will conserve total momentum while ignoring the nuclear recoil)
Min A' decay distance | No       | g4db::G4DarkBreMModel | When A' decay mode is set to FlatDecay, this is the minimum lab-frame distance in mm the A' will travel before decaying (default 0)
Max A' decay distance | No       | g4db::G4DarkBreMModel | When A' decay mode is set to FlatDecay, this is the maximum lab-frame distance in mm the A' will travel before decaying (default 1)
Only One per Event    | No       | G4DarkBremsstrahlung  | Limit the dark brem process to only one interaction per event (default false)
Global Bias           | No       | G4DarkBremsstrahlung  | Bias the dark brem process everywhere by this factor (default 1.0)

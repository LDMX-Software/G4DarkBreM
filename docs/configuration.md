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

# Configuration of Simulation

The G4DarkBreM package follows Geant4's structure of separating
particles from processes. This means configuring the package
requires two steps.
1. Initialize the G4APrime object to define its mass
2. Construct the G4DarkBreMModel to define how the model should behave

It is common to do both of these steps within a `G4VPhysicsConstructor`
since that class can be registered with a Geant4 PhysicsList and handled
later by Geant4. An example G4VPhysicsConstructor is provded by
g4db::example::APrimePhysics.

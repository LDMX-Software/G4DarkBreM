#include "G4DarkBreM/G4FractionallyCharged.h"

G4FractionallyCharged* G4FractionallyCharged::theFractionallyCharged = 0;

G4FractionallyCharged* G4FractionallyCharged::FractionallyCharged() {
  if (!theFractionallyCharged) {
    throw std::runtime_error(
        "Attempting to access the FractionallyCharged particle before it has been "
        "initialized.");
  }

  return theFractionallyCharged;
}

void G4FractionallyCharged::Initialize(double mass, int id, float charge) {
  if (theFractionallyCharged)
    throw std::runtime_error(
        "Attempting to initialize the FractionallyCharged particle more than once.");

  /**
   * Here are the properties of the formal Geant4 fcp we define.
   *
   * Property | Value
   * ---|---
   * short name | fcp-
   * mass | **configured**
   * mass width | 0
   * electric charge | **configured**
   * spin | 0
   * parity | 0
   * conjugation | 0
   * isospin | 0
   * isospin3 | 0
   * Gparity | 0
   * long name | FractionallyCharged-
   * lepton number | 0
   * baryon number | 0
   * PDG ID encoding | **configured**
   * is stable (no decay) | true
   * lifetime | stable
   * decay table | since fcp is stable, always nullptr
   */

  // Create the particle (fcp-)
  theFractionallyCharged = new G4FractionallyCharged(
      "fcp-" /* short name */, mass * MeV, 0. /* mass width */,
      -charge /*electric charge */, 0 /* spin */, 0 /* parity */,
      0 /* conjugation */, 0 /* isospine */, 0 /* isospin3 */, 0 /* G parity */,
      "FractionallyCharged-" /* long name */, 0 /* lepton number */, 0 /* baryon number */,
      id, true /* stable? */, -1 /* lifetime */, nullptr /* decay table */
  );

  // Create the antiparticle (fcp+) explicitly
  G4ParticleDefinition* antiFCP = new G4ParticleDefinition(
      "fcp+" /* short name */, mass * MeV, 0. /* mass width */,
      +charge /*electric charge */, 0 /* spin */, 0 /* parity */,
      0 /* conjugation */, 0 /* isospine */, 0 /* isospin3 */, 0 /* G parity */,
      "FractionallyCharged+" /* long name */, 0 /* lepton number */, 0 /* baryon number */,
      -id, true /* stable? */, -1 /* lifetime */, nullptr /* decay table */
  );


  // Link them as particle/antiparticle pair
  theFractionallyCharged->SetAntiPDGEncoding(-id);
  antiFCP->SetAntiPDGEncoding(id);

  return;
}

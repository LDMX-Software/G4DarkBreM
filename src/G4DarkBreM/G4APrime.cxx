#include "G4DarkBreM/G4APrime.h"
#include "G4DarkBreM/G4FractionallyCharged.h"

#include "G4DecayTable.hh"
#include "G4ParticleTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4PhysicalConstants.hh"
#include "globals.hh"

G4APrime* G4APrime::theAPrime = 0;

G4APrime::DecayMode G4APrime::decay_mode_ = G4APrime::DecayMode::NoDecay;

int G4APrime::decay_id_ = 11;

G4APrime* G4APrime::APrime() {
  if (!theAPrime) {
    throw std::runtime_error(
        "Attempting to access the APrime particle before it has been "
        "initialized.");
  }

  return theAPrime;
}

void G4APrime::Initialize(double mass, int id, double tau,
                          G4APrime::DecayMode decay_mode, int decay_id) {
  if (theAPrime)
    throw std::runtime_error(
        "Attempting to initialize the APrime particle more than once.");

  if (decay_mode == G4APrime::DecayMode::GeantDecay && tau < 0.0)
    throw std::runtime_error(
        "Invalid configuration: DecayMode set to GeantDecay but tau is "
        "negative.");

  G4APrime::decay_mode_ = decay_mode;
  G4APrime::decay_id_ = decay_id;

  /**
   * Here are the properties of the formal Geant4 dark photon we define.
   *
   * Property | Value
   * ---|---
   * short name | A^1
   * mass | **configured**
   * mass width | 0
   * electric charge | 0
   * spin | 0
   * parity | 0
   * conjugation | 0
   * isospin | 0
   * isospin3 | 0
   * Gparity | 0
   * long name | APrime
   * lepton number | 0
   * baryon number | 0
   * PDG ID encoding | **configured**
   * is stable (no decay) | depends on DecayMode
   * lifetime | depends on DecayMode
   * decay table | depends on DecayMode
   */

  theAPrime = new G4APrime(
      "A^1" /* short name */, mass * MeV, 0. /* mass width */,
      0. /*electric charge */, 0 /* spin */, 0 /* parity */,
      0 /* conjugation */, 0 /* isospine */, 0 /* isospin3 */, 0 /* G parity */,
      "APrime" /* long name */, 0 /* lepton number */, 0 /* baryon number */,
      id, true /* stable? */, -1 /* lifetime */, nullptr /* decay table */
  );

  if (decay_mode != G4APrime::DecayMode::NoDecay) {
    G4DecayTable* table = new G4DecayTable();
    G4VDecayChannel* mode = nullptr;

    if (decay_id == 11) {
        mode = new G4PhaseSpaceDecayChannel("A^1", 1.0, 2, "e-", "e+");
    } else if (decay_id == 17) {
        mode = new G4PhaseSpaceDecayChannel("A^1", 1.0, 2, "fcp-", "fcp+");
    } else {
        throw std::runtime_error(
            "Invalid configuration: Decay ID not supported for APrime decay.");
    }
    table->Insert(mode);

    theAPrime->SetPDGStable(false);
    theAPrime->SetPDGLifeTime(tau * second);
    if (decay_mode == G4APrime::DecayMode::FlatDecay)
      theAPrime->SetPDGLifeTime(0.0);  // decay configured in G4DarkBreMModel
    theAPrime->SetDecayTable(table);
  }

  return;
} // end of Initialize

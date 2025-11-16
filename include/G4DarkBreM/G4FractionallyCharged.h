/**
 * @file G4FractionallyCharged.h
 * @brief Class creating the milli-charged particle in Geant.
 * @author Tamas Almos Vami (UCSB)
 */

#ifndef SIMCORE_DARKBREM_G4FRACTIONALLYCHARGED_H_
#define SIMCORE_DARKBREM_G4FRACTIONALLYCHARGED_H_

// Geant
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayTable.hh"
#include "G4ParticleTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4PhysicalConstants.hh"
#include "globals.hh"

class G4String;
class G4DecayTable;

/**
 * Formal class representing the milli-charged particle in Geant4.
 *
 * This class follows the standard prototype for all G4ParticleDefinitions
 * where a static private member holds onto the single instance of this
 * particle definition to be shared by everyone. In this case,
 * G4FractionallyCharged::Initialize needs to be called before any calls to
 * G4FractionallyCharged::G4FractionallyCharged so that the simulation has a defined FCM mass.
 *
 * The G4FractionallyCharged::Initialize call should be done in the
 * ConstructParticle function of a physics constructor.
 */
class G4FractionallyCharged : public G4ParticleDefinition {
 public:
  /**
   * @enum DecayMode
   *
   * How to handle FractionallyCharged decays.
   */
  enum class DecayMode {
    /**
     * No decay/stable -- this is the default.
     */
    NoDecay = 1,
  };

  /**
   * Accessor for FractionallyCharged definition
   *
   * @throws std::runtime_error if the FractionallyCharged is not initialized yet.
   *
   * @see Initialize for configuring and constructing the FractionallyCharged
   * at the start of a run.
   */
  static G4FractionallyCharged* FractionallyCharged();

  /**
   * Initialize the FractionallyCharged particle with the passed configuration
   *
   * @throws std::runtime_error if the FractionallyCharged has already been initialized
   *
   * @param[in] mass The mass of the FractionallyCharged in MeV
   * @param[in] id The PDG ID number to use for the FractionallyCharged particle
   * The default value for the PDG ID is set to 17 as a 4th generation lepton 
   * whose charge is a free parameter which we can set to be a small fraction of e.
   */
  static void Initialize(double mass, int id = 17, float charge = 0.1);

  /// Get the G4FractionallyCharged::DecayMode that was provided to G4FractionallyCharged::Initialize
  static DecayMode getDecayMode() { return DecayMode::NoDecay; }

 private:
  /** Reference to single particle definition of A' */
  static G4FractionallyCharged* theFractionallyCharged;

  /**
   * Constructor
   *
   * Passes all parameters to the base class constructor
   * to register this particle definition with Geant4.
   */
  G4FractionallyCharged(const G4String& name, G4double mass, G4double width, G4double charge,
           G4int iSpin, G4int iParity, G4int iConjugation, G4int iIsospin,
           G4int iIsospin3, G4int gParity, const G4String& pType, G4int lepton,
           G4int baryon, G4int encoding, G4bool stable, G4double lifetime,
           G4DecayTable* decaytable)
      : G4ParticleDefinition(name, mass, width, charge, iSpin, iParity,
                             iConjugation, iIsospin, iIsospin3, gParity, pType,
                             lepton, baryon, encoding, stable, lifetime,
                             decaytable) {}
  /**
   * Destructor
   *
   * Does nothing on purpose.
   */
  virtual ~G4FractionallyCharged() = default;
};

#endif  // SIMCORE_DARKBREM_G4FRACTIONALLYCHARGED_H_

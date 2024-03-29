/**
 * @file PrototypeModel.h
 * @brief Class providing a prototype model for dark brem
 * @author Michael Revering, University of Minnesota
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef G4DARKBREM_PROTOTYPEMODEL_H
#define G4DARKBREM_PROTOTYPEMODEL_H

#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "G4Track.hh"

/**
 * G4DarkBreM internal namespace
 */
namespace g4db {

/**
 * Abstract class representing a model for dark brem.
 *
 * The model is what actually determines two important things:
 *  1. How the cross section is calculated
 *  2. What the particle change is when the process happens
 *
 * This class is the base class that shows what is necessary
 * for the model to function properly.
 */
class PrototypeModel {
 public:
  /**
   * Constructor
   *
   * Configures the model based on the passed parameters
   *
   * Names the logger after the name for this model.
   */
  PrototypeModel(bool muons) : muons_{muons} {}

  /// Destructor, nothing on purpose
  virtual ~PrototypeModel() = default;

  /**
   * Define the verbose level
   *
   * This is called after the model is passed into the process
   * so that the model and the process have the same verbose level.
   * Only call this if you wish the model and the process to not have
   * the same verbose level.
   *
   * @param[in] l verbose leve
   */
  void SetVerboseLevel(int l) { verbose_level_ = l; }

  /**
   * Get the verbose level of the model
   *
   * Use this in model functions to decide if a printout should
   * be done or not.
   *
   * @return integer level (higher means more detail)
   */
  int GetVerboseLevel() const { return verbose_level_; }

  /**
   * Check if we are bremming off muons
   *
   * If false, assume we are bremming off electrons.
   *
   * @return true if dark bremming of muons
   */
  bool DarkBremOffMuons() const { return muons_; }

  /**
   * Print the configuration of this model
   *
   * Helpful for debugging and keeping the process compliant
   * with the other Geant4 processes.
   */
  virtual void PrintInfo() const = 0;

  /**
   * Calculate the cross section given the input parameters
   *
   * @see G4DarkBremmstrahlung::GetMeanFreePath
   * @param[in] leptonKE current lepton kinetic energy
   * @param[in] atomicA atomic-mass number for the element the lepton is in
   * @param[in] atomicZ atomic-number for the element the lepton is in
   * @returns cross section with units incorporated as a G4double
   */
  virtual G4double ComputeCrossSectionPerAtom(G4double leptonKE,
                                              G4double atomicA,
                                              G4double atomicZ) = 0;

  /**
   * Generate the change in the particle now that we can assume the interaction
   * is occuring
   *
   * @note The input particleChange has already been cleared and then
   * initialized, so there is no need for the model to do those steps.
   *
   * @see G4DarkBremmstrahlung::PostStepDoIt
   * @param[in,out] particleChange particle change class that stores information
   * @param[in] track current track that needs the change
   * @param[in] step current step of the track
   * @param[in] element the element selected to dark brem off of
   */
  virtual void GenerateChange(G4ParticleChange& particleChange,
                              const G4Track& track, const G4Step& step,
                              const G4Element& element) = 0;

 protected:
  /// whether muons (true) or electrons (false) are dark bremming
  bool muons_;
  /// verbose level for this model
  int verbose_level_{0};
};  // PrototypeModel

}  // namespace g4db

#endif

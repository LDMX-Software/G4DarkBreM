/**
 * @file G4DarkBremsstrahlung.cxx
 * @brief Class providing the Dark Bremsstrahlung process class.
 * @author Michael Revering, University of Minnesota
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "G4DarkBreM/G4DarkBremsstrahlung.h"

#include "G4Electron.hh"      //for electron definition
#include "G4MuonMinus.hh"     //for muon definition
#include "G4MuonPlus.hh"      //for muon definition
#include "G4EventManager.hh"  //for EventID number
#include "G4ProcessTable.hh"  //for deactivating dark brem process
#include "G4ProcessType.hh"   //for type of process
#include "G4RunManager.hh"    //for VerboseLevel

#include "G4DarkBreM/G4APrime.h"

const std::string G4DarkBremsstrahlung::PROCESS_NAME = "DarkBrem";

G4DarkBremsstrahlung::G4DarkBremsstrahlung(
    std::shared_ptr<g4db::PrototypeModel> the_model,
    bool only_one_per_event, double global_bias, 
    bool cache_xsec, int verbose_level)
    : G4VDiscreteProcess(G4DarkBremsstrahlung::PROCESS_NAME,
                         fElectromagnetic),
      only_one_per_event_{only_one_per_event},
      global_bias_{global_bias}, cache_xsec_{cache_xsec}, model_{the_model} {
  // we need to pretend to be an EM process so 
  // the biasing framework recognizes us
  SetProcessSubType(63);  // needs to be different from the other Em Subtypes

  SetVerboseLevel(verbose_level);
  model_->SetVerboseLevel(verbose_level);

  /*
   * In G4 speak, a "discrete" process is one that only happens at the end of
   * steps. we want the DB to be discrete because it is not a "slow braking"
   * like ionization, the lepton suddenly has the interaction and loses a
   * lot of its energy.
   *
   * The first argument to this function is the process we are adding.
   *      The process manager handles cleaning up the processes,
   *      so we just give it a new pointer.
   * The second argument is the "ordering" index.
   *      This index determines when the process is called w.r.t. the other
   * processes that could be called at the end of the step. Not providing the
   * second argument means that the ordering index is given a default value of
   * 1000 which seems to be safely above all the internal/default processes.
   */
  G4ParticleDefinition* particle_def{G4Electron::ElectronDefinition()};
  if (model_->DarkBremOffMuons()) {
    particle_def = G4MuonMinus::Definition();
  }
  if (GetVerboseLevel() > 0) {
    G4cout << "[ G4DarkBremsstrahlung ] : Connecting dark brem to " 
      << particle_def->GetParticleName() << " "
      << particle_def->GetPDGEncoding() << G4endl;
  }
  G4int ret = particle_def->GetProcessManager()->AddDiscreteProcess(this);
  if (ret < 0) {
    throw std::runtime_error("Particle process manager returned a non-zero status "
        + std::to_string(ret) + " when attempting to register dark brem to it.");
  } else if (GetVerboseLevel() > 0) {
    G4cout
      << "[ G4DarkBremsstrahlung ] : successfully put dark brem in index " 
      << ret << " of process table." << G4endl;
  }
  /**
   * have our custom dark brem process go first in any process ordering
   */
  particle_def->GetProcessManager()->SetProcessOrderingToFirst(this,
      G4ProcessVectorDoItIndex::idxAll);
  if (GetVerboseLevel() > 0) {
    G4cout << "[ G4DarkBremsstrahlung ] : set dark brem process ordering to first" << G4endl;
  }

  if (cache_xsec_) {
    element_xsec_cache_ = g4db::ElementXsecCache(model_);
  }
}

G4bool G4DarkBremsstrahlung::IsApplicable(const G4ParticleDefinition& p) {
  if (model_->DarkBremOffMuons()) return &p == G4MuonMinus::Definition() or &p == G4MuonPlus::Definition();
  else return &p == G4Electron::Definition();
}

void G4DarkBremsstrahlung::PrintInfo() {
  G4cout 
    << " Muons              : " << model_->DarkBremOffMuons() << "\n"
    << " Only One Per Event : " << only_one_per_event_ << "\n"
    << " Global Bias        : " << global_bias_ << "\n"
    << " Cache Xsec         : " << cache_xsec_
    << G4endl;
  model_->PrintInfo();
}

G4VParticleChange* G4DarkBremsstrahlung::PostStepDoIt(const G4Track& track,
                                                       const G4Step& step) {
  // Debugging Purposes: Check if track we get is the configured lepton
  if (not IsApplicable(*track.GetParticleDefinition()))
    throw std::runtime_error("Dark brem process received a track that isn't applicable."); 

  /*
   * Geant4 has decided that it is our time to interact,
   * so we are going to change the particle
   */
  if (GetVerboseLevel() > 2) G4cout << "A dark brem occurred!" << G4endl;

  if (only_one_per_event_) {
    /**
     * Deactivate the process after one dark brem if we restrict ourselves to 
     * only one per event. If this is in the stepping action instead, more than 
     * one brem can occur within each step. Reactivated in RunManager::TerminateOneEvent 
     *
     * Both biased and unbiased process could be in the run (but not at the same time),
     * so we turn off both while silencing the warnings from the process table.
    std::cout << "Deactivating dark brem process" << std::endl;
     */
    std::vector<G4String> db_process_name_options = {
        "biasWrapper(" + PROCESS_NAME + ")", PROCESS_NAME};
    G4ProcessManager* pman = track.GetDefinition()->GetProcessManager();
    for (std::size_t i_proc{0}; i_proc < pman->GetProcessList()->size(); i_proc++) {
      G4VProcess* p{(*(pman->GetProcessList()))[i_proc]};
      if (p->GetProcessName().contains(PROCESS_NAME)) {
        pman->SetProcessActivation(p, false);
        break;
      }
    }
  }

  if (GetVerboseLevel() > 2) G4cout << "Initializing track" << G4endl;
  aParticleChange.Initialize(track);

  if (GetVerboseLevel() > 2) G4cout << "Calling model's GenerateChange" << G4endl;
  model_->GenerateChange(aParticleChange, track, step);

  /*
   * Parent class has some internal counters that need to be reset,
   * so we call it before returning. It will return our shared
   * protected member variable aParticleChange that we have been modifying
   */
  if (GetVerboseLevel() > 2) G4cout << "Calling parent's PostStepDoIt" << G4endl;
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}

G4double G4DarkBremsstrahlung::GetMeanFreePath(const G4Track& track, G4double,
                                                G4ForceCondition*) {
  // won't happen if it isn't applicable
  if (not IsApplicable(*track.GetParticleDefinition())) return DBL_MAX;

  G4double energy = track.GetDynamicParticle()->GetKineticEnergy();
  G4double SIGMA = 0;
  G4Material* materialWeAreIn = track.GetMaterial();
  const G4ElementVector* theElementVector = materialWeAreIn->GetElementVector();
  const G4double* NbOfAtomsPerVolume = materialWeAreIn->GetVecNbOfAtomsPerVolume();
  
  for (size_t i = 0; i < materialWeAreIn->GetNumberOfElements(); i++) {
    G4double AtomicZ = (*theElementVector)[i]->GetZ();
    G4double AtomicA = (*theElementVector)[i]->GetA() / (g / mole);
  
    G4double element_xsec;
  
    if (cache_xsec_)
      element_xsec = element_xsec_cache_.get(energy, AtomicA, AtomicZ);
    else
      element_xsec =
          model_->ComputeCrossSectionPerAtom(energy, AtomicA, AtomicZ);
  
    SIGMA += NbOfAtomsPerVolume[i] * element_xsec;
  }
  SIGMA *= global_bias_;
  if (GetVerboseLevel() > 3) {
    G4cout << "G4DBrem : sigma = " << SIGMA 
      << " initIntLenLeft = " << theInitialNumberOfInteractionLength
      << " nIntLenLeft = " << theNumberOfInteractionLengthLeft << G4endl;
  }
  return SIGMA > DBL_MIN ? 1. / SIGMA : DBL_MAX;
}

#include "G4DarkBreM/ElementXsecCache.h"

namespace g4db {

G4double ElementXsecCache::get(G4double energy, G4double A, G4double Z) {
  key_t key = computeKey(energy, A, Z);
  if (the_cache_.find(key) == the_cache_.end()) {
    if (model_.get() == nullptr) {
      throw std::runtime_error(
          "ElementXsecCache not given a model to calculate cross "
          "sections with.");
    }
    the_cache_[key] = model_->ComputeCrossSectionPerAtom(energy, A, Z);
  }
  return the_cache_.at(key);
}

void ElementXsecCache::stream(std::istream& i) {
  // Buffer for each line
  std::string buffer{};
  // Buffer for ',' characters. Char type so that the stream only extracts one
  // character.
  char comma{};
  // Skip the header
  std::getline(i, buffer);
  double small{1e-3};
  while (std::getline(i, buffer)) {
    // Skip empty lines (avoid parsing as 0s)
    if (buffer == "") {
      continue;
    }
    std::stringstream ss{buffer};
    ss >> std::setprecision(std::numeric_limits<double>::digits10 + 1);
    double E{};
    double A{};
    double Z{};
    double xsec{};
    // Manually parsing CSV sure is fun and not bugprone at all
    ss >> A >> charbuffer;
    ss >> Z >> charbuffer;
    ss >> E >> charbuffer;
    ss >> xsec;
    key_t key{computeKey(E, A, Z)};
    the_cache_[key] = xsec * CLHEP::picobarn;
  }
}
void ElementXsecCache::stream(std::ostream& o) const {
  o << "A [au],Z [protons],Energy [MeV],Xsec [pb]\n"
    << std::setprecision(std::numeric_limits<double>::digits10 +
                         1);  // maximum precision
  for (auto const& cache_entry : the_cache_) {
    const key_t& key = cache_entry.first;
    const double& xsec = cache_entry.second;
    key_t E = key % MAX_E;
    key_t A = ((key - E) / MAX_E) % MAX_A;
    key_t Z = ((key - E) / MAX_E - A) / MAX_A;
    o << A << "," << Z << "," << E << "," << xsec / CLHEP::picobarn << "\n";
  }
  o << std::endl;
}

ElementXsecCache::key_t ElementXsecCache::computeKey(G4double energy,
                                                     G4double A,
                                                     G4double Z) const {
  /**
   * @note This implicitly converts the input double to a unsigned long int.
   * Since the internal units for energy are MeV, this means the cache is
   * binned at a 1 MeV level.
   * The atomic inputs A and Z also undergo the implicit conversion; however,
   * this is less of a worry since different elements are separated by at
   * least one whole unit in A and Z.
   */
  key_t energyKey = energy;
  key_t AKey = A;
  key_t ZKey = Z;
  return (ZKey * MAX_A + AKey) * MAX_E + energyKey;
}

}  // namespace g4db

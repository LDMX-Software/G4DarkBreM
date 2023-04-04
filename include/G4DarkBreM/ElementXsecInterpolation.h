#ifndef G4DarkBREM_ELEMENTXSECINTERPOLATION_H
#define G4DarkBREM_ELEMENTXSECINTERPOLATION_H

#include <memory>

#include "G4DarkBreM/PrototypeModel.h"

namespace g4db {

/**
 * Interpolation of cross sections by element
 *
 * The interpolation method is a modified quadratic interpolation
 * where the main modification is - since we have a way to calculate
 * new sample points - we expand the set of sample points if a
 * cross section is requested for a energy and/or element outside
 * of our current sample range.
 */
class ElementXsecInterpolation {
 public:
  /**
   * Default constructor
   *
   * Does nothing interesting, but no model for calculating cross section has
   * been set. If interpolation is attempted with a default-constructed
   * interpolator, an exception will be thrown.
   */
  ElementXsecInterpolation() = default;

  /**
   * Constructor with a model to calculate the cross section.
   */
  ElementXsecInterpolation(std::shared_ptr<PrototypeModel> model)
      : model_{model} {}

  /**
   * Get the value of the cross section for the input variables.
   *
   * For the interpolation to be successful, we impose the kinematic
   * limit immediately - i.e. if energy is less than twice the A'
   * mass, the cross section is zero.
   *
   * @throws std::runtime_error if no model is available for calculating cross sections
   * @param[in] energy Energy of incident lepton [MeV]
   * @param[in] A atomic mass of element [atomic mass units]
   * @param[in] Z atomic number of element [num protons]
   * @returns cross section corresponding to the input parameters (including
   * units Geant4 style)
   */
  G4double get(G4double energy, G4double A, G4double Z);

 private:
  /**
   * A sample set is two parallel vectors limited to operate for 
   * our interpolation goal
   */
  class SampleSet {
   public:
    /**
     * Initialize the sample set with the two input vectors
     */
    SampleSet(const std::vector<double>& init_x,
              const std::vector<double>& inti_y);

    /**
     * Get the minimum x currently in sample set
     */
    const double& min() const;

    /**
     * Get the maximum x currently in sample set
     */
    const double& max() const;

    /**
     * prepend the sample set with the input point
     */
    void prepend(double x, double y);

    /**
     * append the sample set with the input point
     */
    void append(double x, double y);

    /**
     * Calculate an interpolation for the point x
     */
    double interpolate(double x) const;

   private:
    /// x points
    std::vector<double> x_;
    /// y points
    std::vector<double> y_;
  };

 private:
  /// each element has its own set of sample points
  std::map<int, SampleSet> the_samples_;

  /// shared pointer to the model for calculating cross sections
  std::shared_ptr<PrototypeModel> model_;

};  // ElementXsecInterpolation

}

#endif

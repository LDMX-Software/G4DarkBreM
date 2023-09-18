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
   * @throws std::runtime_error if no model is available for calculating cross
   * sections
   * @param[in] energy Energy of incident lepton [MeV]
   * @param[in] A atomic mass of element [atomic mass units]
   * @param[in] Z atomic number of element [num protons]
   * @returns cross section corresponding to the input parameters (including
   * units Geant4 style)
   */
  G4double get(G4double energy, G4double A, G4double Z);

  /**
   * Stream the table of sample points into the output stream.
   *
   * @param[in,out] o ostream to write to
   */
  void stream(std::ostream& o) const;

  /**
   * Overload the streaming operator for ease
   *
   * @param[in] o ostream to write to
   * @param[in] c interpolation to write out
   * @returns modified ostream
   */
  friend std::ostream& operator<<(std::ostream& o,
                                  const ElementXsecInterpolation c) {
    c.stream(o);
    return o;
  }

 private:
  /**
   * A sample set is two parallel vectors limited to operate for
   * our interpolation goal
   */
  class SampleSet {
   public:
    /**
     * Initialize the sample set with the two input vectors
     *
     * @note We assume that the input vectors are ordered properly
     * i.e. init_x are ordered by their values and init_y are paired
     * with their x.
     *
     * @param[in] init_x x values ordered by x
     * @param[in] init_y corresponding y values
     */
    SampleSet(const std::vector<double>& init_x,
              const std::vector<double>& init_y);

    /**
     * Get the minimum x currently in sample set
     * @return minimum x in sample set
     */
    const double& min() const;

    /**
     * Get the maximum x currently in sample set
     * @return maximum x in sample set
     */
    const double& max() const;

    /**
     * prepend the sample set with the input point
     *
     * We assume that the input x is less than the current
     * minimum x. In debug build mode, this assumption is
     * enforced by assert.
     *
     * @param[in] x new lowest x sample point
     * @param[in] y corresponding y value
     */
    void prepend(double x, double y);

    /**
     * append the sample set with the input point
     *
     * We assume that the input x is greater than the current
     * maximum x. In debug build mode, this assumption is
     * enforced by assert.
     *
     * @param[in] x new highest x sample point
     * @param[in] y corresponding y value
     */
    void append(double x, double y);

    /**
     * Calculate an interpolation for the point x
     *
     * The interpolation is a quadratic interpolation over
     * three points whose range encloses the requested sample x.
     * The first such three points that satisfies this criteria
     * is used, so this usually selects the three points where
     * two points are below x and one above.
     *
     * @note We assume that the input x is within the range
     * of the sample set. In debug build mode, this assumption
     * is enforced with assert.
     *
     * @return interpolated value y for the input x
     */
    double interpolate(double x) const;

    /**
     * Dump our sample points to the input stream.
     *
     * We also want to label each row with the Z of
     * the nucleus the sample is for, so that is a parameter,
     * moreover, we want to add the kinematic_min to the x
     * value we used to recover the actual energy of
     * the particle we have the sampled cross section for.
     *
     * @param[in] o output stream to write to
     * @param[in] z atomic z to write as first entry in CSV
     * @param[in] kinematic_min minimum energy to add to x value
     */
    void dump(std::ostream& o, int z, double kinematic_min) const;

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

}  // namespace g4db

#endif

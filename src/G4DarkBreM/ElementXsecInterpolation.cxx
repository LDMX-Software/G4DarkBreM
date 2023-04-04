#include "G4DarkBreM/ElementXsecInterpolation.h"
#include "G4DarkBreM/G4APrime.h"

namespace g4db {

G4double ElementXsecInterpolation::get(G4double energy, G4double A, G4double Z) {
  static const double kinematic_min = 2*G4APrime::APrime()->GetPDGMass();
  if (energy <= kinematic_min) return 0;
  double x = energy - kinematic_min;
  /**
   * @note Implicit conversion to integer from the double form of Z.
   * This may not function properly if users have a custom element that
   * is an averaged Z, but since elements close in Z behave similarly, 
   * this may not matter either.
   */
  int Z_key = Z;

  auto sample_set = the_samples_.find(Z_key);
  if (sample_set == the_samples_.end()) {
    if (model_.get() == nullptr) {
      throw std::runtime_error(
                      "ElementXsecInterpolation not given a model to calculate cross "
                      "sections with.");
    }
    /**
     * When there is no entry for an element yet, we
     * initialize the samples by calculating the cross
     * section for the requested energy and a sample above
     * and below.
     */
    double xsec = model_->ComputeCrossSectionPerAtom(energy, A, Z);
    std::vector<double> init_x{x/2, x, 2*x};
    std::vector<double> init_y{ 
      model_->ComputeCrossSectionPerAtom(x/2+kinematic_min, A, Z),
      xsec,
      model_->ComputeCrossSectionPerAtom(2*x+kinematic_min, A, Z)
    };
    the_samples_.emplace(std::piecewise_construct,
        std::forward_as_tuple(Z_key),
        std::forward_as_tuple(init_x, init_y));
    return xsec;
  } else if (x < sample_set->second.min()) {
    /**
     * When the requested energy is below the lowest sample point,
     * we determine two more samples. One at the requested energy
     * and one below it.
     */
    double xsec = model_->ComputeCrossSectionPerAtom(energy, A, Z);
    sample_set->second.prepend(x, xsec);
    sample_set->second.prepend(
        x/2, 
        model_->ComputeCrossSectionPerAtom(x/2+kinematic_min, A, Z)
        );
    return xsec;
  } else if (x > sample_set->second.max()) {
    /**
     * When the requested energy is above the highest sample point,
     * we determine two more samples. One at the requested energy
     * and one above it.
     */
    double xsec = model_->ComputeCrossSectionPerAtom(energy, A, Z);
    sample_set->second.append(x, xsec);
    sample_set->second.append(
        2*x, 
        model_->ComputeCrossSectionPerAtom(2*x+kinematic_min, A, Z)
        );
    return xsec;
  } else {
    /**
     * When the requested energy is within the range set by our
     * samples, we do a quadratic interpolation using the three
     * points that enclose the requested energy.
     *
     * We do the interpolation with the _first_ set of three points
     * that encloses the requested energy, so this prefers the point
     * set that has two samples below and one above.
     */
    return sample_set->second.interpolate(x);
  }
}

ElementXsecInterpolation::SampleSet::SampleSet(const std::vector<double>& x, 
    const std::vector<double>& y) : x_{x}, y_{y} {
  assert(x_.size() == y_.size());
}

const double& ElementXsecInterpolation::SampleSet::min() const {
  return x_.front();
}

const double& ElementXsecInterpolation::SampleSet::max() const {
  return x_.back();
}

void ElementXsecInterpolation::SampleSet::prepend(double x, double y) {
  assert(x < x_.front());
  x_.insert(x_.begin(), x);
  y_.insert(y_.begin(), y);
}

void ElementXsecInterpolation::SampleSet::append(double x, double y) {
  assert(x > x_.back());
  x_.push_back(x);
  y_.push_back(y);
}

double ElementXsecInterpolation::SampleSet::interpolate(double x) const {
  assert((x >= x_.front() && x <= x_.back()));
  if (x == x_.front()) {
    return y_.front();
  }
  if (x == x_.back()) {
    return y_.back();
  }
  std::size_t i_0{0};
  for (; i_0 < x_.size()-2; ++i_0) {
    if (x > x_[i_0] and x < x_[i_0+2]) {
      break;
    } else if (x == x_[i_0]) {
      return y_[i_0];
    }
  }
  double x_0{x_[i_0  ]}, y_0{y_[i_0  ]},
         x_1{x_[i_0+1]}, y_1{y_[i_0+1]},
         x_2{x_[i_0+2]}, y_2{y_[i_0+2]};
  return (
       (x - x_1)*(x - x_2)/(x_0 - x_2)/(x_0 - x_1)*y_0
      +(x - x_0)*(x - x_2)/(x_1 - x_2)/(x_1 - x_0)*y_1
      +(x - x_0)*(x - x_1)/(x_2 - x_0)/(x_2 - x_1)*y_2
      );
}

}

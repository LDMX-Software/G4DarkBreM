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
     * section in 10% energy steps from the kinematic minimum
     * up to the requested energy.
     */
    std::vector<double> init_x{0.}, init_y{0.};
    double curr_x{0.1*kinematic_min};
    while (curr_x+kinematic_min < energy) {
      init_x.push_back(curr_x);
      init_y.push_back(model_->ComputeCrossSectionPerAtom(curr_x+kinematic_min, A, Z));
      curr_x *= 1.1;
    }
    init_x.push_back(energy-kinematic_min);
    double xsec = model_->ComputeCrossSectionPerAtom(energy, A, Z);
    init_y.push_back(xsec);
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
        x/1.1, 
        model_->ComputeCrossSectionPerAtom(x/1.1+kinematic_min, A, Z)
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
        1.1*x, 
        model_->ComputeCrossSectionPerAtom(1.1*x+kinematic_min, A, Z)
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

void ElementXsecInterpolation::stream(std::ostream& o) const {
  static const double kinematic_min = 2*G4APrime::APrime()->GetPDGMass();
  o << "Z [protons],Energy [MeV],Xsec [pb]\n"
    << std::setprecision(std::numeric_limits<double>::digits10 +
                         1);  // maximum precision
  for (const auto& sample : the_samples_) {
    sample.second.dump(o, sample.first, kinematic_min);
  }
  o << std::endl;
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

void ElementXsecInterpolation::SampleSet::dump(std::ostream& o, int z, double kinematic_min) const {
  for (std::size_t i_sample{0}; i_sample < x_.size(); ++i_sample) {
    o << z << "," << kinematic_min+x_[i_sample] << "," << y_[i_sample] / CLHEP::picobarn << "\n";
  }
}

}

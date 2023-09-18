/**
 * @file xsec_calc.cxx
 * definition of g4db-xsec-calc executable
 */

#include <unistd.h>

#include <fstream>
#include <iostream>

#include "G4DarkBreM/ElementXsecInterpolation.h"
#include "G4DarkBreM/G4APrime.h"
#include "G4DarkBreM/G4DarkBreMModel.h"

/**
 * print out how to use g4db-xsec-calc
 */
void usage() {
  std::cout
      << "USAGE:\n"
         "  g4db-xsec-calc [options]\n"
         "\n"
         "Calculate dark brem cross sections and write them out to a CSV "
         "table\n"
         "\n"
         "OPTIONS\n"
         "  -h,--help    : produce this help and exit\n"
         "  -o,--output  : output file to write scaled events to\n"
         "  -M,--ap-mass : mass of dark photon in GeV\n"
         "  --muons      : pass to set lepton to muons (otherwise electrons)\n"
         "  --energy     : python-like arange for input energies in GeV (stop, "
         "start stop, start stop step)\n"
         "                 default start is 0 and default step is 0.1 GeV\n"
         "  --target     : define target material with two parameters (atomic "
         "units): Z A\n"
         "  --method     : method to calculate xsec, one of 'fullww', 'hiww', "
         "or 'iww'\n"
         "  --interpolate : run the expanding interpolation instead of the "
         "full xsec\n"
         "                  if an argument is provided to this option, then "
         "the table of sampling points\n"
         "                  used for the interpolation will be dumped into a "
         "file named after that arg.\n"
      << std::flush;
}

/**
 * The names for the different xsec_methods that are acceptable on the command
 * line
 */
static const std::map<std::string, g4db::G4DarkBreMModel::XsecMethod>
    xsec_methods = {{"fullww", g4db::G4DarkBreMModel::XsecMethod::Full},
                    {"hiww", g4db::G4DarkBreMModel::XsecMethod::HyperImproved},
                    {"iww", g4db::G4DarkBreMModel::XsecMethod::Improved}};

/**
 * definition of g4db-xsec-calc
 *
 * We use the cross section caching table used within the G4DarkBremsstrahlung
 * process.
 */
int main(int argc, char* argv[]) try {
  std::string output_filename{"xsec.csv"};
  double ap_mass{0.1};
  double min_energy{0.};
  double max_energy{4.};
  double energy_step{0.1};
  double target_Z{74.};
  double target_A{183.84};
  bool muons{false};
  bool interpolate{false};
  std::string interpo_samples{};
  std::string method{};
  for (int i_arg{1}; i_arg < argc; ++i_arg) {
    std::string arg{argv[i_arg]};
    if (arg == "-h" or arg == "--help") {
      usage();
      return 0;
    } else if (arg == "--muons") {
      muons = true;
    } else if (arg == "--interpolate") {
      interpolate = true;
      if (i_arg + 1 < argc and argv[i_arg + 1][0] != '-') {
        interpo_samples = argv[++i_arg];
      }
    } else if (arg == "--method") {
      if (i_arg + 1 >= argc or argv[i_arg + 1][0] == '-') {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      method = argv[++i_arg];
      if (xsec_methods.find(method) == xsec_methods.end()) {
        std::cerr << method << " is not a recognized xsec method" << std::endl;
        return 1;
      }
    } else if (arg == "-o" or arg == "--output") {
      if (i_arg + 1 >= argc or argv[i_arg + 1][0] == '-') {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      output_filename = argv[++i_arg];
    } else if (arg == "-M" or arg == "--ap-mass") {
      if (i_arg + 1 >= argc or argv[i_arg + 1][0] == '-') {
        std::cerr << arg << " requires an argument after it" << std::endl;
        return 1;
      }
      ap_mass = std::stod(argv[++i_arg]);
    } else if (arg == "--energy") {
      std::vector<std::string> args;
      while (i_arg + 1 < argc and argv[i_arg + 1][0] != '-') {
        args.push_back(argv[++i_arg]);
      }
      if (args.size() == 0) {
        std::cerr << arg << " requires arguments after it" << std::endl;
        return 1;
      } else if (args.size() == 1) {
        max_energy = std::stod(args[0]);
      } else if (args.size() == 2) {
        min_energy = std::stod(args[0]);
        max_energy = std::stod(args[1]);
      } else if (args.size() == 3) {
        min_energy = std::stod(args[0]);
        max_energy = std::stod(args[1]);
        energy_step = std::stod(args[2]);
      }
    } else if (arg == "--target") {
      std::vector<std::string> args;
      while (i_arg + 1 < argc and argv[i_arg + 1][0] != '-') {
        args.push_back(argv[++i_arg]);
      }
      if (args.size() != 2) {
        std::cerr << arg << " requires two arguments: Z A" << std::endl;
        return 1;
      }
      target_Z = std::stod(args[0]);
      target_A = std::stod(args[1]);
    } else {
      std::cout << arg << " is an unrecognized option" << std::endl;
      return 1;
    }
  }

  std::ofstream table_file(output_filename);
  if (!table_file.is_open()) {
    std::cerr << "File '" << output_filename << "' was not able to be opened."
              << std::endl;
    return 2;
  }

  // start at max and work our way down
  //    this mimics the actual progress of a simulation slightly better
  max_energy *= GeV;
  energy_step *= GeV;
  min_energy *= GeV;

  if (method.empty()) {
    if (muons)
      method = "fullww";
    else
      method = "hiww";
  }

  std::cout << "Parameter         : Value\n"
            << "Mass A' [MeV]     : " << ap_mass * GeV << "\n"
            << "Min Energy [MeV]  : " << min_energy << "\n"
            << "Max Energy [MeV]  : " << max_energy << "\n"
            << "Energy Step [MeV] : " << energy_step << "\n"
            << "Lepton            : " << (muons ? "Muons" : "Electrons") << "\n"
            << "Xsec Method       : " << method << "\n"
            << "Target A [amu]    : " << target_A << "\n"
            << "Target Z [amu]    : " << target_Z << "\n"
            << std::flush;

  // the process accesses the A' mass from the G4 particle
  G4APrime::Initialize(ap_mass * GeV);
  auto model = std::make_shared<g4db::G4DarkBreMModel>(
      "LIBRARY NOT NEEDED", muons,
      0.0,  // threshold for non-zero xsec
      1.0,  // epsilon
      g4db::G4DarkBreMModel::ScalingMethod::Undefined,  // scaling method
      xsec_methods.at(method),  // xsec calculation method
      50.0,                     // maximum R for XsecMethod::Auto
      622,                      // ID of dark photon in event library
      false                     // load event library
  );

  // wrap the model in an interpolation object
  //   wrapping is almost no cost, cross section calculations
  //   are only made if the `get` method is called
  g4db::ElementXsecInterpolation interpolation(model);

  table_file << "A [au],Z [protons],Energy [MeV],Xsec [pb]\n";

  G4double current_energy = max_energy;
  int energy_width = (max_energy - min_energy);
  int bar_width = 80;
  int pos = 0;
  bool is_redirected = (isatty(STDOUT_FILENO) == 0);
  while (current_energy > min_energy - energy_step and current_energy > 0) {
    double xsec = interpolate
                      ? interpolation.get(current_energy, target_A, target_Z)
                      : model->ComputeCrossSectionPerAtom(current_energy,
                                                          target_A, target_Z);
    table_file << target_A << "," << target_Z << "," << current_energy << ","
               << xsec / CLHEP::picobarn << "\n";
    current_energy -= energy_step;
    if (not is_redirected) {
      int old_pos{pos};
      pos = bar_width * (max_energy - current_energy) / energy_width;
      if (pos != old_pos) {
        std::cout << "[";
        for (int i{0}; i < bar_width; ++i) {
          if (i < pos)
            std::cout << "=";
          else if (i == pos)
            std::cout << ">";
          else
            std::cout << " ";
        }
        std::cout << "] "
                  << int((max_energy - current_energy) / energy_width * 100.0)
                  << " %\r";
        std::cout.flush();
      }
    }
  }
  if (not is_redirected) std::cout << std::endl;

  table_file.flush();
  table_file.close();

  if (not interpo_samples.empty()) {
    std::ofstream sample_file{interpo_samples};
    if (not sample_file.is_open()) {
      throw std::runtime_error(
          "Unable to open file '" + interpo_samples +
          "' to store table of interpolation sample points.");
    }
    sample_file << interpolation;
    sample_file.flush();
    sample_file.close();
  }

  return 0;
} catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << std::endl;
  return 127;
}

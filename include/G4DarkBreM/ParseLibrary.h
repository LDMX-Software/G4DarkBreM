/**
 * @file ParseLibrary.h
 * Declaration of library parsing function
 */

#ifndef G4DARKBREM_PARSELIBRARY_H
#define G4DARKBREM_PARSELIBRARY_H

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "CLHEP/Vector/LorentzVector.h"

namespace g4db {

/**
 * Data frame to store necessary information from LHE files
 */
struct OutgoingKinematics {
  /// 4-momentum of lepton in center of momentum frame for electron-A'
  /// system
  CLHEP::HepLorentzVector lepton;
  /// 4-vector pointing to center of momentum frame
  CLHEP::HepLorentzVector centerMomentum;
  /// energy of lepton before brem (used as key in mad graph data map)
  double E;
};  // OutgoingKinematics

/**
 * parse the input library and return the in-memory kinematics library
 *
 * ## Library Specification
 * A "library" includes a single file or a directory holding multiple files.
 * No subdirectories are inspected, so the directory must be "flat".
 * Each file in a library can either be a CSV file or an LHE file,
 * both of which can be compressed by gzip if desired.
 *
 * ### CSV
 * The CSV file is expected to have a **single** header line which
 * names the columns. These column names have no requirements
 * (besides the existence of this line).
 *
 * The CSV is required to have 10 columns on all non-empty lines of the file.
 * The 10 columns of the CSV all are in MeV and _in order_ are
 * 1. The target Z
 * 2. The incident lepton energy
 * 3. The total energy of the recoil
 * 4. The x-component of the recoil momentum
 * 5. The y-component of the recoil momentum
 * 6. The z-component of the recoil momentum
 * 7. The total energy of the A'
 * 8. The x-component of the A' momentum
 * 9. The y-component of the A' momentum
 * 10. The z-component of the A' momentum
 *
 * ### LHE
 * The LHE files must have dark brem events in it where "dark brem event"
 * in this context is defined below.
 * ```
 *   lepton_id -1 <skip> <skip> <skip> <skip> px py pz E m
 *   <skip-line>
 *   <skip-line>
 *   lepton_id 1 <skip> <skip> <skip> <skip> px py pz E m
 *   <skip-line>
 *   <skip-line>
 *   aprime_id 1 <skip> <skip> <skip> <skip> px py pz E m
 * ```
 *
 * We also check each line if it contains a 'Znuc' string. Looking
 * for the target Z value from this line matching the following format.
 * ```
 *     <num>   <Z> # Znuc <other comments>
 * ```
 * This is supposed to match a line of the param_card dumped into the header
 * of the LHE. In other words, this parser should function properly if your
 * MG/ME model has a configurable parameter in param_card.dat for the target
 * Z and that parameter is called 'Znuc'.
 *
 * This matches a subcomponent of the LHE scheme written by MadGraph/MadEvent
 * (hence the reason this is the "lhe" parser); however, a lot of information
 * is skipped and additional assumptions are made in order to increase the
 * parsing speed.
 *
 * The `lepton_id` is allowed to be _either_ 11 or 13 _everywhere_. No
 * consistency checking is done.
 *
 * The `E` from the first line is used as the incident lepton energy.
 * The four-momentum from the middle line is the recoil lepton's four momentum,
 * and the four-momentum from the last line is used in conjuction with the
 * recoil four-momentum to calculate the center of momentum vector.
 *
 * @param[in] path path to library to parse
 * @param[in] aprime_lhe_id the ID number for the A' (aka dark photon) in the
 * library being parsed
 * @param[in,out] lib map of target Z and incident energy keys to set of
 * outgoing kinematics
 */
void parseLibrary(
    const std::string& path, int aprime_lhe_id,
    std::map<int, std::map<double, std::vector<OutgoingKinematics>>>& lib);

/**
 * Dump the input library to the input output stream
 *
 * @note This is only helpful for our extraction executable
 * and for testing that the parsing is performing correctly.
 *
 * @param[in,out] o output stream to write CSV to
 * @param[in] lib library to write out
 */
void dumpLibrary(
    std::ostream& o,
    const std::map<int, std::map<double, std::vector<OutgoingKinematics>>>&
        lib);

}  // namespace g4db

#endif

/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "oops/util/stringFunctions.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "eckit/config/Configuration.h"

namespace util {
namespace stringfunctions {

// -----------------------------------------------------------------------------
void swapNameMember(const eckit::Configuration & conf, std::string & filename) {
  if (conf.has("member")) {
    const int mymember = conf.getInt("member");
    std::ostringstream mm;
    mm << std::setw(3) << std::setfill('0') << mymember;
    // Construct the output file name
    std::string str_member = "%{member}%";
    std::size_t member_index = filename.find(str_member);
    if (member_index != std::string::npos) {
      filename = filename.replace(member_index, str_member.length(), mm.str());
    }
  }
}
// -----------------------------------------------------------------------------

}  // namespace stringfunctions
}  // namespace util

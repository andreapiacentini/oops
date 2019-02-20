/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_SETUP_H_
#define TEST_SETUP_H_

#include "eckit/runtime/Main.h"

namespace testing {

//----------------------------------------------------------------------------------------------------------------------

template<bool useBoost>
struct Setup {
  Setup() {
  }
};

#ifndef OOPS_TEST_NO_BOOST
template<>
struct Setup<true> {
  Setup() {
    eckit::Main::initialise(boost::unit_test::framework::master_test_suite().argc,
                                boost::unit_test::framework::master_test_suite().argv);
  }
};
#endif

//----------------------------------------------------------------------------------------------------------------------

}  // namespace testing

#endif  // TEST_SETUP_H_

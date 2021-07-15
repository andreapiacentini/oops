/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ErrorCovarianceQG.h"

#include <cmath>

#include "eckit/config/Configuration.h"
#include "model/FieldsQG.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
//AQ #include "model/QgFortran.h"
#include "model/StateQG.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
ErrorCovarianceQG::ErrorCovarianceQG(const GeometryQG & resol, const oops::Variables & vars,
                                     const eckit::Configuration & conf,
                                     const StateQG &, const StateQG &) {
  oops::Log::trace() << "ErrorCovarianceQG created" << std::endl;
}
// -----------------------------------------------------------------------------
ErrorCovarianceQG::~ErrorCovarianceQG() {
  oops::Log::trace() << "ErrorCovarianceQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceQG::multiply(const IncrementQG & dxin, IncrementQG & dxout) const {
}
// -----------------------------------------------------------------------------
void ErrorCovarianceQG::inverseMultiply(const IncrementQG & dxin, IncrementQG & dxout) const {
  oops::IdentityMatrix<IncrementQG> Id;
  dxout.zero();
  GMRESR(dxout, dxin, *this, Id, 20, 1.0e-5);
}
// -----------------------------------------------------------------------------
void ErrorCovarianceQG::randomize(IncrementQG & dx) const {
}
// -----------------------------------------------------------------------------
void ErrorCovarianceQG::print(std::ostream & os) const {
  os << "ErrorCovarianceQG::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace qg

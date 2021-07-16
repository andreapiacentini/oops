/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ChangeVarTLADQG.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/StateQG.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVarTLADQG::ChangeVarTLADQG(const StateQG &, const StateQG &,
                                 const GeometryQG & resol, const eckit::Configuration & conf) {}
// -----------------------------------------------------------------------------
ChangeVarTLADQG::~ChangeVarTLADQG() {}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::multiply(const IncrementQG & dxa, IncrementQG & dxm) const {
  //AQ qg_change_var_tl_f90(dxa.fields().toFortran(), dxm.fields().toFortran());
  oops::Log::debug() << "ChangeVarTLADQG::multiply" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::multiplyInverse(const IncrementQG & dxm, IncrementQG & dxa) const {
  //AQ qg_change_var_tl_f90(dxm.fields().toFortran(), dxa.fields().toFortran());
  oops::Log::debug() << "ChangeVarTLADQG::multiplyInverse" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::multiplyAD(const IncrementQG & dxm, IncrementQG & dxa) const {
  //AQ qg_change_var_ad_f90(dxm.fields().toFortran(), dxa.fields().toFortran());
  oops::Log::debug() << "ChangeVarTLADQG::multiplyAD" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::multiplyInverseAD(const IncrementQG & dxa, IncrementQG & dxm) const {
  //AQ qg_change_var_ad_f90(dxa.fields().toFortran(), dxm.fields().toFortran());
  oops::Log::debug() << "ChangeVarTLADQG::multiplyInverseAD" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::print(std::ostream & os) const {
  os << "QG linear change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg


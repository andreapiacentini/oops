/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
//AQ #include "model/instantiateQgChangeVarFactory.h"
#include "model/QgTraits.h"
#include "oops/runs/EnsVariance.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  //AQ qg::instantiateQgChangeVarFactory();
  oops::EnsVariance<qg::QgTraits> var;
  return run.execute(var);
}

/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSFILTERS_H_
#define OOPS_BASE_OBSFILTERS_H_

#include <memory>
#include <set>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Printable.h"

namespace oops {

/// Holds observation filters (usually QC) for one observation type

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsFilters : public util::Printable,
                   private boost::noncopyable {
  typedef GeoVaLs<OBS>            GeoVaLs_;
  typedef ObsDiagnostics<OBS>     ObsDiags_;
  typedef ObsSpace<OBS>           ObsSpace_;
  typedef ObsVector<OBS>          ObsVector_;
  typedef std::shared_ptr<ObsFilterBase<OBS> >  ObsFilterPtr_;
  typedef std::shared_ptr<ObsDataVector<OBS, int> >  ObsDataPtr_;

 public:
  /// Initialize all filters for \p obspace, from parameters, using
  /// \p qcflags and \p obserr (observation error variances)
  /// \p iteration argument indicates outer loop iteration in the variational
  /// assimilation
  ObsFilters(const ObsSpace_ &, const std::vector<ObsFilterParametersWrapper<OBS>> &,
             ObsDataPtr_ qcflags, ObsVector_ & obserr, const int iteration = 0);

  void preProcess() const;
  void priorFilter(const GeoVaLs_ &) const;
  void postFilter(const ObsVector_ &, const ObsDiags_ &) const;

  Variables requiredVars() const {return geovars_;}
  Variables requiredHdiagnostics() const {return diagvars_;}

 private:
  void print(std::ostream &) const override;

  std::vector<ObsFilterPtr_> filters_;
  Variables geovars_;
  Variables diagvars_;
  ObsDataPtr_ qcflags_;
  ObsVector_ & obserr_;
  std::shared_ptr<ObsDataVector<OBS, float> > obserrtmp_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsFilters<OBS>::ObsFilters(const ObsSpace_ & os,
                            const std::vector<ObsFilterParametersWrapper<OBS>> & filtersParams,
                            ObsDataPtr_ qcflags, ObsVector_ & obserr, const int iteration)
  : filters_(), geovars_(), diagvars_(), qcflags_(qcflags), obserr_(obserr),
    obserrtmp_(new ObsDataVector<OBS, float>(obserr)) {
  Log::trace() << "ObsFilters::ObsFilters starting:\n";
  for (const ObsFilterParametersWrapper<OBS> &filterParams : filtersParams)
    Log::trace() << "  " << filterParams << std::endl;

// Prepare QC handling and statistics if any filters are present
  if (filtersParams.size() > 0) {
    eckit::LocalConfiguration preconf;
    preconf.set("filter", "QCmanager");
    filters_.push_back(FilterFactory<OBS>::create(os, preconf, qcflags_, obserrtmp_));
  }

// Create the filters, only at 0-th iteration, or at iterations specified in "apply at iterations"
  for (const ObsFilterParametersWrapper<OBS> &filterParams : filtersParams) {
    // Only create filters for the 0-th iteration by default
    bool apply = (iteration == 0);
    // If "apply at iterations" is set, check if this is the right iteration
    if (filterParams.applyAtIterations.value() != boost::none) {
      std::set<int> iters = parseIntSet(*filterParams.applyAtIterations.value());
      apply = contains(iters, iteration);
    }
    if (apply) {
      ObsFilterPtr_ tmp(FilterFactory<OBS>::create(os, filterParams.filterParameters,
                                                   qcflags_, obserrtmp_));
      geovars_ += tmp->requiredVars();
      diagvars_ += tmp->requiredHdiagnostics();
      filters_.push_back(tmp);
    }
  }

  Log::trace() << "ObsFilters::ObsFilters done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::preProcess() const {
  for (const auto & filter : filters_) {
    filter->preProcess();
  }
  obserrtmp_->mask(*qcflags_);
  obserr_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::priorFilter(const GeoVaLs_ & gv) const {
  for (const auto & filter : filters_) {
    filter->priorFilter(gv);
  }
  obserrtmp_->mask(*qcflags_);
  obserr_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::postFilter(const ObsVector_ & hofx, const ObsDiags_ & diags) const {
  for (const auto & filter : filters_) {
    filter->postFilter(hofx, diags);
  }
  obserrtmp_->mask(*qcflags_);
  obserr_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsFilters<OBS>::print(std::ostream & os) const {
  os << "ObsFilters: " << filters_.size() << " elements:" << std::endl;
  for (const auto & filter : filters_) {
    os << *filter << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTERS_H_
